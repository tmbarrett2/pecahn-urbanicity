# ============================================================
# rebuild_mirror_and_run_sites.ps1
#
# Unattended weekend run (2026-11):
#   1. preflight  - Docker, disk, package install, site list + Earth Engine check
#   2. 3a         - buffers for reference points + study sites (02_generate_local_buffers.R)
#   3. 3c-3d      - clip planet, filter local layers, merge with EXISTING global layers
#   4. 3e         - wipe old Overpass DB, init new one (~29 h last time), restart, verify
#   5. 06         - compute urbanicity metrics for the study sites (06_compute_study_sites.R)
#
# Same pinned planet (2026-06-08) and global layers as the June build, so the
# reference points see identical OSM data; the only change is that the local
# 1 km layers now also cover the study sites.
#
# Launch from an ordinary PowerShell window (NOT inside an IDE that may close):
#   powershell -ExecutionPolicy Bypass -File "C:\Users\tmbar\OneDrive\Documents\GitHub\pecahn-urbanicity\reference_data\osm_backend\rebuild_mirror_and_run_sites.ps1"
#
# Resume after a failure by skipping finished steps, e.g.:
#   ... -File rebuild_mirror_and_run_sites.ps1 -StartAt init     (steps: preflight, buffers, extract, init, sites)
#
# Log: C:\osm\logs\weekend_<timestamp>.log
# ============================================================

param(
	[ValidateSet("preflight", "buffers", "extract", "init", "wait", "sites")]
	[string]$StartAt = "preflight"
)

# Continue (not Stop): in PS 5.1, Stop + native stderr (docker/osmium progress) throws spuriously.
# Native failures are caught via exit codes in Check(); cmdlets that matter use -ErrorAction Stop.
$ErrorActionPreference = "Continue"

$REPO    = "C:\Users\tmbar\OneDrive\Documents\GitHub\pecahn-urbanicity"
$OSM     = "C:\osm"
$DATA    = "C:\osm\data"
$DBDIR   = "C:\osm\overpass_db"
$VOL     = "C:/osm/data:/data"
$IMG     = "wiktorn/overpass-api"
$RSCRIPT = "C:\Program Files\R\R-4.4.3\bin\Rscript.exe"
$API     = "http://localhost:12345/api/interpreter?data=%5Bout%3Ajson%5D%3Bnode%281%29%3Bout%3B"
$CONTAINER = "overpass_local"

New-Item -ItemType Directory -Force "$OSM\logs" | Out-Null
$LOG = "$OSM\logs\weekend_{0}.log" -f (Get-Date -Format "yyyyMMdd_HHmmss")
Start-Transcript -Path $LOG | Out-Null

# Keep the PC awake while THIS process runs (released automatically on exit).
Add-Type -Namespace Win32 -Name Power -MemberDefinition '[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);'
[Win32.Power]::SetThreadExecutionState([uint32]"0x80000001") | Out-Null   # ES_CONTINUOUS | ES_SYSTEM_REQUIRED

$steps = @("preflight", "buffers", "extract", "init", "wait", "sites")
$startIdx = [array]::IndexOf($steps, $StartAt)
function Run-Step($name) { return ([array]::IndexOf($steps, $name) -ge $startIdx) }

function Say($msg, $color = "Cyan") { Write-Host ("[{0}] {1}" -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"), $msg) -ForegroundColor $color }
function Fail($msg) { Say "FAILED: $msg" "Red"; Say "Log: $LOG" "Red"; Stop-Transcript | Out-Null; exit 1 }
# Run a native exe so its stdout+stderr land in the transcript as plain text.
function X { $exe, $rest = $args; & $exe @rest 2>&1 | ForEach-Object { "$_" } | Out-Host }
function Check($what) { if ($LASTEXITCODE -ne 0) { Fail "$what (exit $LASTEXITCODE)" } }

function Test-Api {
	try { $r = Invoke-WebRequest -UseBasicParsing -TimeoutSec 30 $API; return ($r.Content -match '"elements"') }
	catch { return $false }
}

# ------------------------------------------------------------
if (Run-Step "preflight") {
	Say "=== PREFLIGHT ==="

	# Docker Desktop up?
	docker info *> $null
	if ($LASTEXITCODE -ne 0) {
		Say "Docker not running - starting Docker Desktop..."
		Start-Process "C:\Program Files\Docker\Docker\Docker Desktop.exe"
		$ok = $false
		for ($i = 0; $i -lt 36; $i++) { Start-Sleep 10; docker info *> $null; if ($LASTEXITCODE -eq 0) { $ok = $true; break } }
		if (-not $ok) { Fail "Docker Desktop did not come up within 6 minutes" }
	}
	Say ("Docker OK. Image digest: " + (docker image inspect $IMG --format "{{index .RepoDigests 0}}"))

	# Inputs that are reused from the June build
	foreach ($f in @("$DATA\planet-latest.osm.pbf", "$DATA\global_layers.osm.pbf")) {
		if (-not (Test-Path $f)) { Fail "missing $f" }
	}

	# Disk: need ~160 GB peak during init AFTER the old DB (~73 GB) is removed
	$freeGB = [math]::Round((Get-PSDrive C).Free / 1GB)
	Say "Free on C: $freeGB GB (plus ~73 GB freed when old DB is wiped)"
	if ($freeGB -lt 100) { Fail "less than 100 GB free on C:" }

	# Install the package from the repo so the run uses the current fixes
	Say "Installing pecahnurbanicity from repo..."
	X $RSCRIPT -e "devtools::install('$($REPO -replace '\\','/')/pecahnurbanicity', upgrade = FALSE, quiet = TRUE)"
	Check "package install"

	# Site list parses + Earth Engine authenticates (Overpass is expected to be down here)
	X $RSCRIPT "$REPO\reference_data\06_compute_study_sites.R" --check
	Check "06 --check (site list / Earth Engine)"

	Say "Preflight OK - you can walk away now." "Green"
}

# ------------------------------------------------------------
if (Run-Step "buffers") {
	Say "=== 3a: buffers (reference points + study sites) ==="
	Push-Location "$REPO\reference_data"
	# Retry: R/sf has crashed intermittently here (0xC0000005) yet succeeds on rerun.
	for ($try = 1; $try -le 3; $try++) {
		X $RSCRIPT "02_generate_local_buffers.R"
		$code = $LASTEXITCODE
		if ($code -eq 0) { break }
		Say "  buffers attempt $try failed (exit $code); retrying..." "Yellow"
	}
	Pop-Location
	if ($code -ne 0) { Fail "02_generate_local_buffers.R (exit $code)" }
}

# ------------------------------------------------------------
if (Run-Step "extract") {
	Say "=== 3c: extract (clip planet to buffers) ==="
	X docker run --rm --entrypoint osmium -v $VOL $IMG extract -s smart -p /data/local_buffers.geojson /data/planet-latest.osm.pbf -o /data/local_clip.osm.pbf --overwrite
	Check "osmium extract"

	Say "=== 3c: filter local heavy layers ==="
	X docker run --rm --entrypoint osmium -v $VOL $IMG tags-filter /data/local_clip.osm.pbf nwr/building nwr/shop nwr/public_transport nwr/railway=station,halt nwr/communication=mobile_phone -o /data/local_layers.osm.pbf --overwrite
	Check "osmium tags-filter"

	Say "=== 3d: merge with existing global layers ==="
	X docker run --rm --entrypoint osmium -v $VOL $IMG merge /data/global_layers.osm.pbf /data/local_layers.osm.pbf -o /data/planet-pecahn.osm.pbf --overwrite
	Check "osmium merge"

	Remove-Item "$DATA\local_clip.osm.pbf" -Force   # intermediate; frees disk for init
	Say ("Merged extract: {0:N1} GB" -f ((Get-Item "$DATA\planet-pecahn.osm.pbf").Length / 1GB))
}

# ------------------------------------------------------------
if (Run-Step "init") {
	Say "=== 3e: wipe old Overpass DB and rebuild ==="
	X docker compose -f "$OSM\docker-compose.yml" down
	if (-not (Test-Path "$DATA\planet-pecahn.osm.pbf")) { Fail "new $DATA\planet-pecahn.osm.pbf not found (did 'extract' run?)" }

	Get-ChildItem $DBDIR -Force | Remove-Item -Recurse -Force -ErrorAction Stop
	Move-Item "$DATA\planet-pecahn.osm.pbf" "$DBDIR\planet-pecahn.osm.pbf" -ErrorAction Stop

	X docker compose -f "$OSM\docker-compose.yml" up -d
	Check "docker compose up (init)"
	Say "Overpass init started (took ~29 h in June)."
}

# ------------------------------------------------------------
# "wait" is safe to resume into while an init is already running: it never wipes
# anything, only waits for init_done, restarts the container to serve, and verifies.
if (Run-Step "wait") {
	Say "=== 3e: waiting for Overpass init (polling every 10 min) ==="
	if (Test-Path "$DBDIR\init_done") {
		Say "init_done already present."
	}

	# Init mode exits the container when the DB is built; init_done marks success.
	$deadline = (Get-Date).AddHours(72)
	while (-not (Test-Path "$DBDIR\init_done") -or ((docker inspect -f "{{.State.Running}}" $CONTAINER 2>$null) -eq "true" -and -not (Test-Api))) {
		Start-Sleep 600
		$running = (docker inspect -f "{{.State.Running}}" $CONTAINER 2>$null) -eq "true"
		$done = Test-Path "$DBDIR\init_done"
		if (-not $running) {
			if ($done) { Say "Init finished; restarting container to serve." "Green"; break }
			X docker logs --tail 50 $CONTAINER
			Fail "Overpass container stopped before init_done appeared (see log tail above)"
		}
		if ($done -and (Test-Api)) { Say "Init finished and API already serving." "Green"; break }
		if ((Get-Date) -gt $deadline) { Fail "Overpass init exceeded 72 h" }
		Say ("  still initialising... db size {0:N1} GB" -f ((Get-ChildItem "$DBDIR\db" -Recurse -File -ErrorAction SilentlyContinue | Measure-Object Length -Sum).Sum / 1GB)) "Gray"
	}

	X docker compose -f "$OSM\docker-compose.yml" up -d
	Check "docker compose up (serve)"
	$ok = $false
	for ($i = 0; $i -lt 60; $i++) { if (Test-Api) { $ok = $true; break }; Start-Sleep 30 }
	if (-not $ok) { Fail "Overpass API not answering 30 min after restart" }
	Say "Overpass mirror live at http://localhost:12345/api/interpreter" "Green"
}

# ------------------------------------------------------------
if (Run-Step "sites") {
	Say "=== 06: study-site urbanicity metrics ==="
	if (-not (Test-Api)) { Fail "Overpass API not responding" }
	X $RSCRIPT "$REPO\reference_data\06_compute_study_sites.R"
	Check "06_compute_study_sites.R"
}

Say "ALL DONE." "Green"
Say "Results: C:\Users\tmbar\OneDrive - Duke University\Active Projects\pecahn\study_sites_results\study_sites_urbanicity.csv" "Green"
Stop-Transcript | Out-Null
