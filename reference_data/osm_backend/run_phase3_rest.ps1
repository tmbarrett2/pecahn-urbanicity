# ============================================================
# run_phase3_rest.ps1
#
# Runs Phase 3c-3d of the PECAHN Overpass extract sequentially,
# stopping immediately if any step fails. Start it AFTER 3b
# (global_layers.osm.pbf) has finished.
#
# Launch from C:\osm :
#   powershell -ExecutionPolicy Bypass -File .\run_phase3_rest.ps1
# ============================================================

$VOL = "C:/osm/data:/data"
$IMG = "wiktorn/overpass-api"

Write-Host "=== 3c: extract (clip planet to point buffers) ===" -ForegroundColor Cyan
docker run --rm --entrypoint osmium -v $VOL $IMG extract -s smart -p /data/local_buffers.geojson /data/planet-latest.osm.pbf -o /data/local_clip.osm.pbf
if ($LASTEXITCODE -ne 0) { Write-Host "extract FAILED (exit $LASTEXITCODE); stopping." -ForegroundColor Red; exit 1 }

Write-Host "=== 3c: filter local heavy layers ===" -ForegroundColor Cyan
docker run --rm --entrypoint osmium -v $VOL $IMG tags-filter /data/local_clip.osm.pbf nwr/building nwr/shop nwr/public_transport nwr/railway=station,halt nwr/communication=mobile_phone -o /data/local_layers.osm.pbf
if ($LASTEXITCODE -ne 0) { Write-Host "filter FAILED (exit $LASTEXITCODE); stopping." -ForegroundColor Red; exit 1 }

Write-Host "=== 3d: merge into planet-pecahn.osm.pbf ===" -ForegroundColor Cyan
docker run --rm --entrypoint osmium -v $VOL $IMG merge /data/global_layers.osm.pbf /data/local_layers.osm.pbf -o /data/planet-pecahn.osm.pbf
if ($LASTEXITCODE -ne 0) { Write-Host "merge FAILED (exit $LASTEXITCODE); stopping." -ForegroundColor Red; exit 1 }

Write-Host "Phase 3c-3d complete." -ForegroundColor Green
Get-Item C:\osm\data\planet-pecahn.osm.pbf | Select-Object Name, @{N='GB';E={[math]::Round($_.Length/1GB,1)}}
