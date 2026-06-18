# Self-hosted Overpass backend — build notes

The urbanicity run (`04_compute_urbanicity_reference.Rmd`) issues thousands of progressive OSM feature queries.
To avoid public-API rate limits **and** to pin every query to one reproducible OSM
snapshot, the run is served by a **local Overpass instance** built from a
tag-filtered + buffer-clipped extract of the OSM planet.

Served at: `http://localhost:12345/api/interpreter`
(set in `04_compute_urbanicity_reference.Rmd` via `osmdata::set_overpass_url(...)`).

---

## Pinned inputs & tooling

| Item | Value |
|---|---|
| OSM planet snapshot | **2026-06-08** (downloaded as `planet-260608.osm.pbf`, renamed `planet-latest.osm.pbf`) |
| Planet MD5 | `eb419b0f06089ef97e4b49f937ebf4a8` |
| `osmium-tool` | 1.15.0 (libosmium 2.19.0) |
| Docker image | `wiktorn/overpass-api@sha256:04d4e189ba98aa42619414b8dab2599e852d3e7add9b84001b6de1f22319c01c` |
| Docker | Docker Desktop (Docker Engine 29.5.3) |

> The `docker-compose.yml` in this folder references the image as
> `wiktorn/overpass-api` (i.e. the floating `:latest` tag) — that is the artifact
> as it was run. For a byte-reproducible rebuild, pin it to the digest above:
> `image: wiktorn/overpass-api@sha256:04d4e189ba98aa42619414b8dab2599e852d3e7add9b84001b6de1f22319c01c`.

All `osmium` steps run inside the same image via `--entrypoint osmium`, so no
separate `osmium` install is required. Working directory is `C:\osm`; the
`data/` subfolder is bind-mounted at `/data`.

---

## Split-layer strategy

A whole-planet Overpass instance is unnecessary and enormous. The extract keeps:

- **Global layers** (kept worldwide so progressive distance searches never clip
  at a buffer edge): all `highway` ways + `amenity=hospital,school,bank,atm`.
- **Local 1 km layers** (heavy, only needed within the analysis boxes, so clipped
  to the 3 km point buffers first): `building`, `shop`, `public_transport`,
  `railway=station,halt`, `communication=mobile_phone`.

The two are merged into a single `planet-pecahn.osm.pbf` that backs the instance.

---

## Build steps

### Phase 3a — buffers (R, see `../02_generate_local_buffers.R`)
Buffers the 1,400 reference points 3 km in World Mollweide, dissolves, and writes
`local_buffers.geojson` (WGS84) for `osmium extract -p`.

### Phase 3b — global layers (kept worldwide)
```sh
docker run --rm --entrypoint osmium -v C:/osm/data:/data wiktorn/overpass-api \
  tags-filter /data/planet-latest.osm.pbf \
  nwr/highway nwr/amenity=hospital,school,bank,atm \
  -o /data/global_layers.osm.pbf
```

### Phase 3c — clip planet to the point buffers, then filter local layers
```sh
# clip (smart strategy keeps referenced nodes/ways/relations complete)
docker run --rm --entrypoint osmium -v C:/osm/data:/data wiktorn/overpass-api \
  extract -s smart -p /data/local_buffers.geojson \
  /data/planet-latest.osm.pbf -o /data/local_clip.osm.pbf

# filter the heavy local layers from the clip
docker run --rm --entrypoint osmium -v C:/osm/data:/data wiktorn/overpass-api \
  tags-filter /data/local_clip.osm.pbf \
  nwr/building nwr/shop nwr/public_transport nwr/railway=station,halt \
  nwr/communication=mobile_phone \
  -o /data/local_layers.osm.pbf
```

### Phase 3d — merge global + local
```sh
docker run --rm --entrypoint osmium -v C:/osm/data:/data wiktorn/overpass-api \
  merge /data/global_layers.osm.pbf /data/local_layers.osm.pbf \
  -o /data/planet-pecahn.osm.pbf

# optional sanity check on the merged extent
docker run --rm --entrypoint osmium -v C:/osm/data:/data wiktorn/overpass-api \
  fileinfo -e /data/planet-pecahn.osm.pbf
```

> Phases 3c–3d are scripted in `run_phase3_rest.ps1` (run after 3b finishes):
> `powershell -ExecutionPolicy Bypass -File .\run_phase3_rest.ps1`

### Phase 3e — build & serve the Overpass instance
```sh
# the compose file expects the merged PBF inside the DB volume
move C:\osm\data\planet-pecahn.osm.pbf C:\osm\overpass_db\planet-pecahn.osm.pbf

# init mode builds the Overpass DB from the PBF on first start
docker compose up -d
```

`docker-compose.yml` runs in `OVERPASS_MODE: init` with a PBF→bz2 preprocess
(`OVERPASS_PLANET_PREPROCESS`), because the `wiktorn/overpass-api` init path
expects a bzip2-compressed planet. `OVERPASS_META: "no"` (metadata not needed)
and `OVERPASS_COMPRESSION: gz` keep the built DB compact. Port `12345:80`
exposes the API at the URL above.

### Verify
```sh
curl "http://localhost:12345/api/interpreter?data=[out:json];node(1);out;"
```
A well-formed (even if empty) JSON response confirms the instance is live before
launching the urbanicity run.
