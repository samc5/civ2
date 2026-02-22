#!/usr/bin/env python3
"""
inject_assets.py — Inject coloured asset-marker points into LAZ tiles.

For each asset (from assets_with_tiles.geojson) that falls within a
configurable radius of the tile centre, a small vertical pillar of
brightly-coloured points is injected into the point cloud so the assets
are clearly visible when the output file is opened in Potree Desktop,
CloudCompare, or any other LAS/LAZ viewer.

Usage
-----
    # Download tile #0 from the CDN and inject nearby assets
    python inject_assets.py --tile 0

    # Process a locally-downloaded LAZ file
    python inject_assets.py myfile.laz

    # Process all 47 tiles (downloads each one — ~14 GB total!)
    python inject_assets.py --all

    # Custom radius (default 200 m)
    python inject_assets.py --tile 0 --radius 300

Output goes to ./output_laz/<original_name>_assets.laz

Requirements
------------
    pip install laspy[lazrs] pyproj numpy requests
"""

import argparse
import json
import math
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    import laspy
except ImportError:
    sys.exit("Missing dependency — run:  pip install laspy[lazrs]")

try:
    from pyproj import CRS, Transformer
except ImportError:
    sys.exit("Missing dependency — run:  pip install pyproj")

try:
    import requests
except ImportError:
    requests = None  # only needed for --tile / --all downloads


# ═══════════════════════════════════════════════════════════════════════════
# PATHS
# ═══════════════════════════════════════════════════════════════════════════
SCRIPT_DIR = Path(__file__).resolve().parent
TILES_JSON = SCRIPT_DIR / "pointcloudjp.json"
ASSETS_FILE = SCRIPT_DIR / "assets_with_tiles.geojson"
OUTPUT_DIR = SCRIPT_DIR / "output_laz"

# ═══════════════════════════════════════════════════════════════════════════
# TUNING
# ═══════════════════════════════════════════════════════════════════════════
DEFAULT_RADIUS_M = 200       # show assets within this distance of tile centre
PILLAR_HEIGHT_M = 4.0        # height of the marker pillars
PILLAR_NPTS = 40             # number of points per pillar (vertical)
DISC_RADIUS_M = 0.5          # radius of the base disc
DISC_NPTS = 12               # points around the disc
LINE_SPACING_M = 0.5         # spacing of markers along LineString assets
LINE_PILLAR_NPTS = 5         # mini-pillar points per line position
MARKER_CLASSIFICATION = 31   # classification value (max for 5-bit format ≤5)
Z_GRID_CELL_M = 2.0          # cell size for ground-elevation grid
Z_GRID_MAX_SAMPLE = 2_000_000  # max points to sample for Z grid
CHUNK_SIZE = 1_000_000       # points per streaming chunk (keeps RAM ~100 MB)

# ═══════════════════════════════════════════════════════════════════════════
# COLOUR MAP  —  asset_type → (R, G, B) in 0–255
# ═══════════════════════════════════════════════════════════════════════════
COLOR_MAP = {
    "VALVE_COVER":            (231,  76,  60),
    "MANHOLE_COVER":          ( 52, 152, 219),
    "CATCH_BASIN":            ( 41, 128, 185),
    "RAMP":                   ( 46, 204, 113),
    "DRIVEWAY":               (243, 156,  18),
    "CURB":                   (230, 126,  34),
    "SIDEWALK":               (149, 165, 166),
    "PLANTER":                ( 39, 174,  96),
    "TRASH_BIN":              (127, 140, 141),
    "BUS_STOP":               (241, 196,  15),
    "CABINET":                (211,  84,   0),
    "HYDRANT":                (192,  57,  43),
    "CCTV":                   (142,  68, 173),
    "TREE":                   ( 30, 132,  73),
    "BUS_SHELTER":            (233,  30,  99),
    "BIKE_RACK":              (  0, 188, 212),
    "BOX":                    (121,  85,  72),
    "BENCH":                  (160,  82,  45),
    "UTILITY_POLE":           (158, 158, 158),
    "SIGNAL_POLE":            (255, 152,   0),
    "FENCE":                  (188, 143, 143),
    "RETAINING_WALL":         ( 96, 125, 139),
    "PEDESTRIAN_PUSH_BUTTON": (  0, 172, 193),
    "TRAFFIC_SIGNAL":         (244,  67,  54),
    "FLASHER":                (255,  64, 129),
    "LUMINARIES":             (255, 213,  79),
    "JERSEY_BARRIER":         (189, 189, 189),
}
DEFAULT_COLOR = (255, 0, 255)  # magenta for unknown types


# ═══════════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════════

def haversine_m(lat1, lon1, lat2, lon2):
    """Great-circle distance in metres."""
    R = 6_371_000
    to_rad = math.pi / 180
    dlat = (lat2 - lat1) * to_rad
    dlon = (lon2 - lon1) * to_rad
    a = (math.sin(dlat / 2) ** 2
         + math.cos(lat1 * to_rad) * math.cos(lat2 * to_rad)
         * math.sin(dlon / 2) ** 2)
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def detect_crs(las):
    """Return a pyproj CRS for the LAZ file.

    Tries VLR parsing first, then falls back to a coordinate-range heuristic
    appropriate for the Brookline / Jamaica Plain, MA area.
    Works with both a LasData object and a LasHeader object.
    """
    header = las.header if hasattr(las, 'header') else las
    return detect_crs_from_header(header)


def detect_crs_from_header(header):
    """Return a pyproj CRS from a LasHeader."""
    # Attempt 1: parse CRS from file VLRs
    try:
        crs = header.parse_crs()
        if crs is not None:
            epsg = crs.to_epsg()
            print(f"  CRS from VLR: {crs.name} (EPSG:{epsg})")
            return crs
    except Exception:
        pass

    # Attempt 2: coordinate-range heuristic
    x_mid = (header.x_min + header.x_max) / 2
    y_mid = (header.y_min + header.y_max) / 2

    if 600_000 < x_mid < 950_000:
        # MA State Plane NAD83 US feet (EPSG:2249)
        epsg = 2249
    elif 200_000 < x_mid < 280_000:
        # MA State Plane NAD83 metres (EPSG:26986)
        epsg = 26986
    elif 280_000 < x_mid < 420_000:
        # UTM Zone 19N (EPSG:32619)
        epsg = 32619
    else:
        # Safe default for this dataset
        epsg = 32619

    print(f"  CRS heuristic: EPSG:{epsg}  (x_mid={x_mid:.0f}, y_mid={y_mid:.0f})")
    return CRS.from_epsg(epsg)


def build_z_grid(x, y, z, cell=Z_GRID_CELL_M):
    """Build a coarse 2-D median-Z grid for ground-elevation lookup.

    Accepts plain numpy arrays (already subsampled by caller).
    """
    n = len(x)
    x_min, y_min = float(x.min()), float(y.min())
    ix = ((x - x_min) / cell).astype(np.int32)
    iy = ((y - y_min) / cell).astype(np.int32)

    cells = defaultdict(list)
    for i in range(n):
        cells[(int(ix[i]), int(iy[i]))].append(float(z[i]))

    grid = {}
    for k, vs in cells.items():
        grid[k] = float(np.median(vs))

    return grid, x_min, y_min, cell


def lookup_z(grid, gx_min, gy_min, gcell, px, py, fallback):
    """Look up ground Z from the grid, searching nearby cells if needed."""
    ix = int((px - gx_min) / gcell)
    iy = int((py - gy_min) / gcell)

    # Spiral outward up to 5 cells
    for r in range(6):
        for dx in range(-r, r + 1):
            for dy in range(-r, r + 1):
                if max(abs(dx), abs(dy)) != r:
                    continue
                z = grid.get((ix + dx, iy + dy))
                if z is not None:
                    return z
    return fallback


# ═══════════════════════════════════════════════════════════════════════════
# MARKER GENERATORS
# ═══════════════════════════════════════════════════════════════════════════

def make_pillar(cx, cy, cz, rgb8):
    """Vertical pillar + base disc of coloured points."""
    r, g, b = rgb8
    pts = []

    # Vertical column
    for i in range(PILLAR_NPTS):
        z = cz + (PILLAR_HEIGHT_M * i / max(1, PILLAR_NPTS - 1))
        pts.append((cx, cy, z, r, g, b))

    # Base disc
    for i in range(DISC_NPTS):
        angle = 2 * math.pi * i / DISC_NPTS
        dx = DISC_RADIUS_M * math.cos(angle)
        dy = DISC_RADIUS_M * math.sin(angle)
        pts.append((cx + dx, cy + dy, cz, r, g, b))

    return pts


def make_line_markers(coords_proj, cz, rgb8):
    """Points along a projected LineString with small pillars."""
    r, g, b = rgb8
    pts = []

    for k in range(len(coords_proj) - 1):
        x1, y1 = coords_proj[k]
        x2, y2 = coords_proj[k + 1]
        seg_len = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        if seg_len < 0.01:
            continue
        n_pts = max(2, int(seg_len / LINE_SPACING_M))
        for i in range(n_pts):
            t = i / max(1, n_pts - 1)
            x = x1 + t * (x2 - x1)
            y = y1 + t * (y2 - y1)
            # Small pillar at each sample
            for j in range(LINE_PILLAR_NPTS):
                z = cz + (PILLAR_HEIGHT_M * j / max(1, LINE_PILLAR_NPTS - 1))
                pts.append((x, y, z, r, g, b))

    return pts


# ═══════════════════════════════════════════════════════════════════════════
# TILE DOWNLOAD
# ═══════════════════════════════════════════════════════════════════════════

def download_tile(url, dest_path):
    """Stream-download a LAZ tile with progress."""
    if requests is None:
        sys.exit("Install requests for CDN downloads:  pip install requests")
    resp = requests.get(url, stream=True)
    resp.raise_for_status()
    total = int(resp.headers.get("content-length", 0))
    dl = 0
    with open(dest_path, "wb") as f:
        for chunk in resp.iter_content(chunk_size=131072):
            f.write(chunk)
            dl += len(chunk)
            if total:
                pct = dl * 100 // total
                mb = dl // 1_000_000
                tot_mb = total // 1_000_000
                print(f"\r  Downloading… {pct}%  ({mb}/{tot_mb} MB)", end="", flush=True)
    print()


# ═══════════════════════════════════════════════════════════════════════════
# MAIN PROCESSING
# ═══════════════════════════════════════════════════════════════════════════

def process_tile(laz_path, tile_lat, tile_lon, assets_geojson, output_path,
                 radius_m=DEFAULT_RADIUS_M):
    """Read LAZ in streaming chunks, inject markers, write output.

    Never holds the entire point cloud in memory — peak RAM is ~1 chunk
    (~1 M points ≈ 40–60 MB) plus the small Z-grid and marker array.
    """
    print(f"\n{'=' * 60}")
    print(f"Processing: {laz_path}")
    t0 = time.time()

    # ── 1. Read header only ──────────────────────────────────────────
    print("  Reading header…")
    with laspy.open(str(laz_path)) as reader:
        hdr = reader.header
        n_orig = hdr.point_count
        fmt_id = hdr.point_format.id
        has_color = "red" in hdr.point_format.dimension_names
        print(f"  {n_orig:,} points, format {fmt_id}, "
              f"bbox x=[{hdr.x_min:.1f}, {hdr.x_max:.1f}]")
        if not has_color:
            print("  ⚠  Point format has no RGB — markers will appear but without colour.")

    # ── 2. Detect CRS & build transformer ────────────────────────────
    # Use an unopened read of just the header for CRS detection
    with laspy.open(str(laz_path)) as reader:
        crs_pc = detect_crs_from_header(reader.header)
    crs_wgs = CRS.from_epsg(4326)
    transformer = Transformer.from_crs(crs_wgs, crs_pc, always_xy=True)

    # ── 3. Filter nearby assets ──────────────────────────────────────
    nearby = []
    for f in assets_geojson["features"]:
        geom = f["geometry"]
        if geom["type"] == "Point":
            flon, flat = geom["coordinates"][:2]
        elif geom["type"] == "LineString":
            coords = geom["coordinates"]
            mid = coords[len(coords) // 2]
            flon, flat = mid[:2]
        else:
            continue
        if haversine_m(tile_lat, tile_lon, flat, flon) <= radius_m:
            nearby.append(f)

    if not nearby:
        print(f"  No assets within {radius_m} m — copying file unchanged.")
        import shutil
        shutil.copy2(str(laz_path), str(output_path))
        print(f"  → {output_path.name}")
        return

    type_counts = defaultdict(int)
    for f in nearby:
        type_counts[f["properties"].get("asset_type", "UNKNOWN")] += 1
    print(f"  {len(nearby)} assets within {radius_m} m:")
    for t, c in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"      {t}: {c}")

    # ── 4. Streaming pass 1: build Z grid from subsampled points ─────
    print("  Building elevation grid (streaming)…")
    sample_xs, sample_ys, sample_zs = [], [], []
    sampled = 0
    subsample_step = max(1, n_orig // Z_GRID_MAX_SAMPLE)

    with laspy.open(str(laz_path)) as reader:
        point_idx = 0
        for chunk in reader.chunk_iterator(CHUNK_SIZE):
            n_chunk = len(chunk)
            # Compute which indices in this chunk to sample
            # We want global indices 0, step, 2*step, ... that fall in [point_idx, point_idx+n_chunk)
            first_in_chunk = ((point_idx + subsample_step - 1) // subsample_step) * subsample_step
            local_indices = np.arange(first_in_chunk - point_idx, n_chunk, subsample_step)
            if len(local_indices) > 0:
                x_c = np.asarray(chunk.x)
                y_c = np.asarray(chunk.y)
                z_c = np.asarray(chunk.z)
                sample_xs.append(x_c[local_indices])
                sample_ys.append(y_c[local_indices])
                sample_zs.append(z_c[local_indices])
                sampled += len(local_indices)
            point_idx += n_chunk
            print(f"\r    Scanned {point_idx:,} / {n_orig:,} pts, sampled {sampled:,}",
                  end="", flush=True)

    print()
    sx = np.concatenate(sample_xs); del sample_xs
    sy = np.concatenate(sample_ys); del sample_ys
    sz = np.concatenate(sample_zs); del sample_zs
    z_grid, gx_min, gy_min, gcell = build_z_grid(sx, sy, sz)
    z_fallback = float(np.median(sz[::max(1, len(sz) // 10_000)]))
    del sx, sy, sz
    print(f"  Z grid: {len(z_grid):,} cells")

    # ── 5. Generate marker points ────────────────────────────────────
    print("  Generating markers…")
    all_pts = []

    for feat in nearby:
        geom = feat["geometry"]
        atype = feat["properties"].get("asset_type", "UNKNOWN")
        rgb = COLOR_MAP.get(atype, DEFAULT_COLOR)

        if geom["type"] == "Point":
            lon, lat = geom["coordinates"][:2]
            px, py = transformer.transform(lon, lat)
            pz = lookup_z(z_grid, gx_min, gy_min, gcell, px, py, z_fallback)
            all_pts.extend(make_pillar(px, py, pz, rgb))

        elif geom["type"] == "LineString":
            coords_wgs = geom["coordinates"]
            coords_proj = [transformer.transform(c[0], c[1]) for c in coords_wgs]
            mid = coords_proj[len(coords_proj) // 2]
            pz = lookup_z(z_grid, gx_min, gy_min, gcell, mid[0], mid[1], z_fallback)
            all_pts.extend(make_line_markers(coords_proj, pz, rgb))

    del z_grid  # free memory before the write pass

    n_mark = len(all_pts)
    if n_mark == 0:
        print("  No marker points generated — copying unchanged.")
        import shutil
        shutil.copy2(str(laz_path), str(output_path))
        return

    print(f"  {n_mark:,} marker points for {len(nearby)} assets")

    # ── 6. Streaming pass 2: copy original + append markers ──────────
    print(f"  Writing {output_path.name} (streaming)…")

    with laspy.open(str(laz_path)) as reader:
        out_header = laspy.LasHeader(
            point_format=reader.header.point_format,
            version=reader.header.version,
        )
        out_header.scales = reader.header.scales
        out_header.offsets = reader.header.offsets
        for vlr in reader.header.vlrs:
            out_header.vlrs.append(vlr)

        with laspy.open(str(output_path), mode="w", header=out_header) as writer:
            # Stream-copy original points
            written = 0
            for chunk in reader.chunk_iterator(CHUNK_SIZE):
                writer.write_points(chunk)
                written += len(chunk)
                print(f"\r    Copied {written:,} / {n_orig:,} pts",
                      end="", flush=True)
            print()

            # Build and write marker points
            print("  Appending markers…")
            # Need a sample raw array to get the dtype
            sample_raw = np.zeros(1, dtype=chunk.array.dtype)
            marker_raw = np.zeros(n_mark, dtype=sample_raw.dtype)

            arr = np.array(all_pts, dtype=np.float64)
            del all_pts
            marker_x, marker_y, marker_z = arr[:, 0], arr[:, 1], arr[:, 2]
            marker_r = (arr[:, 3] * 256).astype(np.uint16)
            marker_g = (arr[:, 4] * 256).astype(np.uint16)
            marker_b = (arr[:, 5] * 256).astype(np.uint16)

            scales = np.array(reader.header.scales)
            offsets = np.array(reader.header.offsets)
            marker_raw["X"] = np.round((marker_x - offsets[0]) / scales[0]).astype(np.int32)
            marker_raw["Y"] = np.round((marker_y - offsets[1]) / scales[1]).astype(np.int32)
            marker_raw["Z"] = np.round((marker_z - offsets[2]) / scales[2]).astype(np.int32)

            if has_color:
                marker_raw["red"] = marker_r
                marker_raw["green"] = marker_g
                marker_raw["blue"] = marker_b

            if "raw_classification" in marker_raw.dtype.names:
                marker_raw["raw_classification"] = np.uint8(MARKER_CLASSIFICATION)
            elif "classification" in marker_raw.dtype.names:
                marker_raw["classification"] = np.uint8(MARKER_CLASSIFICATION)

            if "intensity" in marker_raw.dtype.names:
                marker_raw["intensity"] = np.uint16(65535)

            rec = laspy.PackedPointRecord(marker_raw, reader.header.point_format)
            writer.write_points(rec)

    elapsed = time.time() - t0
    total_pts = n_orig + n_mark
    out_size_mb = output_path.stat().st_size / 1_000_000
    print(f"  ✓ {output_path.name}  ({total_pts:,} pts, +{n_mark:,} markers, "
          f"{out_size_mb:.0f} MB)  in {elapsed:.1f} s")


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Inject coloured asset markers into LAZ point clouds.")
    parser.add_argument(
        "file", nargs="?", default=None,
        help="Path to a local .laz file to process.")
    parser.add_argument(
        "--tile", type=int, default=None,
        help="Tile index (0-based) from pointcloudjp.json to download & process.")
    parser.add_argument(
        "--all", action="store_true",
        help="Process all tiles (downloads each from CDN).")
    parser.add_argument(
        "--radius", type=float, default=DEFAULT_RADIUS_M,
        help=f"Asset search radius in metres (default: {DEFAULT_RADIUS_M}).")
    parser.add_argument(
        "--output-dir", type=str, default=None,
        help=f"Output directory (default: {OUTPUT_DIR}).")
    args = parser.parse_args()

    out_dir = Path(args.output_dir) if args.output_dir else OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load assets
    if not ASSETS_FILE.exists():
        sys.exit(f"Assets file not found: {ASSETS_FILE}\n"
                 f"Run spatial_join.py first to generate it.")
    with open(ASSETS_FILE) as f:
        assets = json.load(f)
    print(f"Loaded {len(assets['features'])} assets from {ASSETS_FILE.name}")

    # Load tile index
    if not TILES_JSON.exists():
        sys.exit(f"Tile index not found: {TILES_JSON}")
    with open(TILES_JSON) as f:
        tiles = json.load(f)
    all_tiles = tiles["features"]
    print(f"Loaded {len(all_tiles)} tiles from {TILES_JSON.name}")

    # Determine which tiles to process
    if args.file:
        # Local file mode — match filename to tile metadata for lat/lon
        laz_path = Path(args.file)
        if not laz_path.exists():
            sys.exit(f"File not found: {laz_path}")

        fname = laz_path.name
        tile_props = None
        for t in all_tiles:
            if t["properties"]["filename"] == fname:
                tile_props = t["properties"]
                break

        if tile_props is None:
            # Can't match by name — read LAZ header only to estimate lat/lon
            print(f"  Filename '{fname}' not in tile index — detecting location from LAZ header…")
            with laspy.open(str(laz_path)) as tmp:
                crs_pc = detect_crs_from_header(tmp.header)
                crs_wgs = CRS.from_epsg(4326)
                rev = Transformer.from_crs(crs_pc, crs_wgs, always_xy=True)
                cx = (tmp.header.x_min + tmp.header.x_max) / 2
                cy = (tmp.header.y_min + tmp.header.y_max) / 2
                lon, lat = rev.transform(cx, cy)
            tile_props = {"lat": lat, "lon": lon, "filename": fname}
            print(f"  Estimated centre: lat={lat:.6f}, lon={lon:.6f}")

        out_name = laz_path.stem + "_assets.laz"
        process_tile(laz_path, tile_props["lat"], tile_props["lon"],
                     assets, out_dir / out_name, args.radius)

    elif args.tile is not None:
        idx = args.tile
        if idx < 0 or idx >= len(all_tiles):
            sys.exit(f"Tile index {idx} out of range (0–{len(all_tiles)-1})")
        tile = all_tiles[idx]
        p = tile["properties"]
        url = p["baseUrl"] + p["lasPath"]

        # Download to a temp location
        dl_dir = out_dir / "downloads"
        dl_dir.mkdir(exist_ok=True)
        dl_path = dl_dir / p["filename"]
        if dl_path.exists():
            print(f"  Using cached download: {dl_path}")
        else:
            print(f"  Downloading tile {idx}: {p['filename']} …")
            download_tile(url, dl_path)

        out_name = dl_path.stem + "_assets.laz"
        process_tile(dl_path, p["lat"], p["lon"],
                     assets, out_dir / out_name, args.radius)

    elif args.all:
        print(f"\nProcessing all {len(all_tiles)} tiles…")
        dl_dir = out_dir / "downloads"
        dl_dir.mkdir(exist_ok=True)

        for idx, tile in enumerate(all_tiles):
            p = tile["properties"]
            url = p["baseUrl"] + p["lasPath"]
            dl_path = dl_dir / p["filename"]

            if dl_path.exists():
                print(f"\n[{idx+1}/{len(all_tiles)}] Cached: {p['filename']}")
            else:
                print(f"\n[{idx+1}/{len(all_tiles)}] Downloading: {p['filename']}")
                download_tile(url, dl_path)

            out_name = dl_path.stem + "_assets.laz"
            process_tile(dl_path, p["lat"], p["lon"],
                         assets, out_dir / out_name, args.radius)

    else:
        parser.print_help()
        print("\nExamples:")
        print("  python inject_assets.py --tile 0        # download & process tile 0")
        print("  python inject_assets.py myfile.laz      # process a local LAZ file")
        print("  python inject_assets.py --all           # process all 47 tiles")
        sys.exit(1)

    print("\n✓ Done — output in", out_dir)


if __name__ == "__main__":
    main()
