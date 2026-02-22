#!/usr/bin/env python3
"""
compute_chokepoints.py

Calculates sidewalk chokepoints caused by obstacles between paired
curb + sidewalk line features.

Algorithm
---------
1. Pair each CURB LineString with the nearest SIDEWALK LineString 
   (centroid proximity < 15 m).
2. For each pair, compute the "corridor width" (perpendicular distance
   between curb and sidewalk at sampled stations along the curb).
3. Find every point obstacle (hydrant, tree, utility pole, …) that
   lies inside the corridor (within the curb–sidewalk strip).
4. For each obstacle, measure the *effective clear width* — the
   shortest perpendicular distance from the obstacle to the nearer
   edge (curb or sidewalk).  That is the "remaining passable width"
   for a pedestrian on the narrower side.
5. The chokepoint of the block is the obstacle with the smallest
   effective clear width.

Outputs
-------
  chokepoints.geojson   — A FeatureCollection with:
    • LineString features representing the measurement line from the
      obstacle to its nearest point on the closer edge.
    • Properties: corridor_width_m, clear_width_m, obstacle_type,
      reduction_pct, curb_fid, sidewalk_fid, etc.

Usage
-----
    python compute_chokepoints.py
"""

import json, math, sys, os

# --------------------------------------------------------------------------- #
#  Configuration
# --------------------------------------------------------------------------- #
INPUT_ASSETS = "aboveGroundAssets.geojson"
OUTPUT_FILE  = "chokepoints.geojson"

MAX_PAIR_DIST_M = 15        # max centroid distance to pair curb ↔ sidewalk
CORRIDOR_BUFFER_M = 10      # obstacle must be within this of BOTH lines
CORRIDOR_SAMPLE_STEP = 5    # sample corridor width every N vertices on curb
CHOKEPOINT_THRESHOLD_FT = 4 # only report obstacles where clear width < this (feet)
M_TO_FT = 3.28084

# Obstacle types that can cause chokepoints (point geometry only)
OBSTACLE_TYPES = {
    "HYDRANT", "TREE", "UTILITY_POLE", "SIGNAL_POLE", "TRASH_BIN",
    "BIKE_RACK", "BENCH", "BUS_SHELTER", "BUS_STOP", "CABINET", "BOX",
    "PLANTER", "PEDESTRIAN_PUSH_BUTTON", "FLASHER", "CATCH_BASIN",
    "MANHOLE_COVER", "VALVE_COVER",
}

# --------------------------------------------------------------------------- #
#  Geometry helpers
# --------------------------------------------------------------------------- #
def haversine_m(lat1, lon1, lat2, lon2):
    """Distance in metres between two WGS-84 points."""
    R = 6_371_000
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2
         + math.cos(math.radians(lat1))
         * math.cos(math.radians(lat2))
         * math.sin(dlon / 2) ** 2)
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def point_to_line(plat, plon, line_coords):
    """
    Nearest point on a LineString (list of [lon, lat]) to a query point.
    Returns (distance_m, nearest_lat, nearest_lon).
    """
    best_d = float("inf")
    best_lat = best_lon = None
    for i in range(len(line_coords) - 1):
        ax, ay = line_coords[i]       # lon, lat
        bx, by = line_coords[i + 1]
        dx, dy = bx - ax, by - ay
        if dx == 0 and dy == 0:
            t = 0.0
        else:
            t = max(0.0, min(1.0, ((plon - ax) * dx + (plat - ay) * dy)
                                   / (dx * dx + dy * dy)))
        nx, ny = ax + t * dx, ay + t * dy
        d = haversine_m(plat, plon, ny, nx)
        if d < best_d:
            best_d, best_lat, best_lon = d, ny, nx
    return best_d, best_lat, best_lon


def line_centroid(coords):
    """Return (lat, lon) centroid of a coordinate list [[lon,lat], …]."""
    lat = sum(c[1] for c in coords) / len(coords)
    lon = sum(c[0] for c in coords) / len(coords)
    return lat, lon


def sample_corridor_width(curb_coords, sw_coords, step=CORRIDOR_SAMPLE_STEP):
    """
    Sample perpendicular corridor width by projecting vertices of
    the curb onto the sidewalk and averaging the distances.
    Returns (mean_width_m, max_width_m).
    """
    widths = []
    for i in range(0, len(curb_coords), max(1, step)):
        clat, clon = curb_coords[i][1], curb_coords[i][0]
        d, _, _ = point_to_line(clat, clon, sw_coords)
        widths.append(d)
    if not widths:
        return 0.0, 0.0
    return sum(widths) / len(widths), max(widths)


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #
def main():
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), INPUT_ASSETS)
    if not os.path.isfile(src):
        sys.exit(f"ERROR: {src} not found.  Run spatial_join.py first.")

    with open(src) as f:
        gj = json.load(f)

    # ── Partition features ------------------------------------------------- #
    curbs, sidewalks, obstacles = [], [], []
    for feat in gj["features"]:
        at = feat["properties"].get("asset_type") or ""
        geom = feat["geometry"]
        if not geom:
            continue
        if at == "CURB" and geom["type"] == "LineString":
            curbs.append(feat)
        elif at == "SIDEWALK" and geom["type"] == "LineString":
            sidewalks.append(feat)
        elif at in OBSTACLE_TYPES and geom["type"] == "Point":
            obstacles.append(feat)

    print(f"Loaded {len(curbs)} curbs, {len(sidewalks)} sidewalks, "
          f"{len(obstacles)} point obstacles")

    # ── Pair curbs ↔ sidewalks -------------------------------------------- #
    pairs = []  # (curb_feat, sidewalk_feat, centroid_dist)
    for c in curbs:
        cc = c["geometry"]["coordinates"]
        clat, clon = line_centroid(cc)
        best = None
        for s in sidewalks:
            sc = s["geometry"]["coordinates"]
            slat, slon = line_centroid(sc)
            d = haversine_m(clat, clon, slat, slon)
            if d < MAX_PAIR_DIST_M and (best is None or d < best[0]):
                best = (d, s)
        if best:
            pairs.append((c, best[1], best[0]))

    print(f"Paired {len(pairs)} curb–sidewalk corridors")

    # ── For each pair, find obstacles and compute chokepoints -------------- #
    chokepoint_features = []

    for curb, sw, _ in pairs:
        cc = curb["geometry"]["coordinates"]
        sc = sw["geometry"]["coordinates"]
        c_fid = curb["properties"].get("feature_id", "?")
        s_fid = sw["properties"].get("feature_id", "?")

        mean_w, max_w = sample_corridor_width(cc, sc)
        if mean_w < 0.3:
            continue  # degenerate pair

        # Find obstacles inside this corridor
        corridor_obs = []
        for obs in obstacles:
            olat = obs["geometry"]["coordinates"][1]
            olon = obs["geometry"]["coordinates"][0]

            dc, c_nlat, c_nlon = point_to_line(olat, olon, cc)
            ds, s_nlat, s_nlon = point_to_line(olat, olon, sc)

            # Obstacle must be reasonably close to BOTH edges
            if dc > CORRIDOR_BUFFER_M or ds > CORRIDOR_BUFFER_M:
                continue
            # Obstacle should be between the two lines (sum ≈ corridor width)
            # Allow some tolerance (1.5× mean corridor width)
            if dc + ds > mean_w * 2.5:
                continue

            # The "clear width" is the remaining passable corridor width
            # at the obstacle's station.  The curb and sidewalk lines in the
            # data bound the "furniture zone", not the full walkable area
            # (the sidewalk extends further toward buildings).  Pedestrians
            # can route around the obstacle on either side, so the effective
            # passable width = dc + ds (full corridor cross-section minus
            # the obstacle's near-zero footprint).
            clear_width = dc + ds

            # Only keep obstacles that actually create a chokepoint
            # (clear width under threshold)
            if clear_width * M_TO_FT >= CHOKEPOINT_THRESHOLD_FT:
                continue

            # Draw measurement line across the full corridor (curb → sidewalk)
            # to show the passable width

            corridor_obs.append({
                "obs": obs,
                "olat": olat, "olon": olon,
                "dc": dc, "ds": ds,
                "clear_width": clear_width,
                # edge projections for the measurement line
                "curb_lat": c_nlat, "curb_lon": c_nlon,
                "sw_lat": s_nlat, "sw_lon": s_nlon,
            })

        if not corridor_obs:
            continue

        # Sort by clear width ascending (tightest first)
        corridor_obs.sort(key=lambda x: x["clear_width"])

        # Emit all obstacles as features (tightest = the chokepoint)
        for rank, co in enumerate(corridor_obs):
            otype = (co["obs"]["properties"].get("asset_type")
                     or co["obs"]["properties"].get("type") or "UNKNOWN")
            is_chokepoint = (rank == 0)
            clear_ft = round(co["clear_width"] * M_TO_FT, 2)
            corridor_ft = round(mean_w * M_TO_FT, 2)
            reduction = ((mean_w - co["clear_width"]) / mean_w * 100
                         if mean_w > 0 else 0)

            # Feature 1: measurement line spanning curb → sidewalk
            # (the full passable corridor width at this obstacle)
            chokepoint_features.append({
                "type": "Feature",
                "properties": {
                    "kind": "chokepoint_line",
                    "is_chokepoint": is_chokepoint,
                    "obstacle_type": otype,
                    "clear_width_m": round(co["clear_width"], 2),
                    "clear_width_ft": clear_ft,
                    "corridor_width_m": round(mean_w, 2),
                    "corridor_width_ft": corridor_ft,
                    "max_corridor_width_m": round(max_w, 2),
                    "reduction_pct": round(reduction, 1),
                    "curb_fid": c_fid,
                    "sidewalk_fid": s_fid,
                    "rank": rank + 1,
                    "image_url": co["obs"]["properties"].get("image_url"),
                    "condition": co["obs"]["properties"].get("condition"),
                },
                "geometry": {
                    "type": "LineString",
                    "coordinates": [
                        [co["curb_lon"], co["curb_lat"]],
                        [co["sw_lon"], co["sw_lat"]],
                    ],
                },
            })

            # Feature 2: obstacle point marker
            chokepoint_features.append({
                "type": "Feature",
                "properties": {
                    "kind": "obstacle_point",
                    "is_chokepoint": is_chokepoint,
                    "obstacle_type": otype,
                    "clear_width_m": round(co["clear_width"], 2),
                    "clear_width_ft": clear_ft,
                    "corridor_width_m": round(mean_w, 2),
                    "corridor_width_ft": corridor_ft,
                    "reduction_pct": round(reduction, 1),
                    "rank": rank + 1,
                    "curb_fid": c_fid,
                    "sidewalk_fid": s_fid,
                    "image_url": co["obs"]["properties"].get("image_url"),
                    "condition": co["obs"]["properties"].get("condition"),
                },
                "geometry": {
                    "type": "Point",
                    "coordinates": [co["olon"], co["olat"]],
                },
            })

    # ── Write output ------------------------------------------------------- #
    out = {
        "type": "FeatureCollection",
        "features": chokepoint_features,
    }
    dst = os.path.join(os.path.dirname(os.path.abspath(__file__)), OUTPUT_FILE)
    with open(dst, "w") as f:
        json.dump(out, f)

    n_corridors = len(set(
        (ft["properties"]["curb_fid"], ft["properties"]["sidewalk_fid"])
        for ft in chokepoint_features
    ))
    n_chokepoints = sum(
        1 for ft in chokepoint_features
        if ft["properties"].get("kind") == "obstacle_point"
           and ft["properties"].get("is_chokepoint")
    )
    print(f"Wrote {len(chokepoint_features)} features to {OUTPUT_FILE}")
    print(f"  {n_corridors} corridors with obstacles")
    print(f"  {n_chokepoints} chokepoints (tightest per corridor)")


if __name__ == "__main__":
    main()
