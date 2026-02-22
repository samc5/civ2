#!/usr/bin/env python3
"""
Compute road centerlines from paired curb geometries.

Reads curb LineStrings from aboveGroundAssets.geojson, chains them into two
continuous polylines (one per road side), then averages them to produce a
centerline.  The centerline is split at each traffic intersection to produce
one LineString per inter-intersection segment.

Output: road_centerlines.geojson
  - One Feature per segment between consecutive intersections
  - Properties carry the intersection names and indices so the HTML can
    apply the traffic gradient.
"""

import json, math, csv
from typing import List, Tuple

# ── Config ───────────────────────────────────────────────────────────────────
ASSETS_FILE = "aboveGroundAssets.geojson"
TRAFFIC_CSV = "combined_intersection_15min_counts.csv"
OUTPUT_FILE = "road_centerlines.geojson"

CHAIN_GAP_M = 15          # max gap (m) when linking chain fragments
RESAMPLE_STEP_M = 5.0     # resample interval for averaging two sides
SNAP_DIST_M = 80          # max dist to snap an intersection to a chain

# ── Helpers ──────────────────────────────────────────────────────────────────
Pt = Tuple[float, float]  # (lat, lon)

def haversine(lat1, lon1, lat2, lon2) -> float:
    """Distance in metres between two WGS-84 points."""
    R = 6_371_000
    p = math.pi / 180
    a = (math.sin((lat2 - lat1) * p / 2) ** 2
         + math.cos(lat1 * p) * math.cos(lat2 * p)
         * math.sin((lon2 - lon1) * p / 2) ** 2)
    return 2 * R * math.asin(min(1, math.sqrt(a)))


def chain_length(pts: List[Pt]) -> float:
    """Total arc-length of a polyline in metres."""
    return sum(haversine(*pts[i], *pts[i + 1]) for i in range(len(pts) - 1))


def resample(pts: List[Pt], step_m: float) -> List[Pt]:
    """Resample a polyline at uniform arc-length intervals."""
    if len(pts) < 2:
        return list(pts)
    total = chain_length(pts)
    if total < 1e-6:
        return [pts[0]]
    n_out = max(2, int(round(total / step_m)) + 1)
    out: List[Pt] = []
    seg_idx = 0
    seg_start = 0.0  # arc-length at start of current segment
    seg_len = haversine(*pts[0], *pts[1])
    for i in range(n_out):
        target = total * i / (n_out - 1)
        # Advance to the segment that contains 'target'
        while seg_idx < len(pts) - 2 and target > seg_start + seg_len + 1e-9:
            seg_start += seg_len
            seg_idx += 1
            seg_len = haversine(*pts[seg_idx], *pts[seg_idx + 1])
        # Interpolate within segment
        if seg_len < 1e-9:
            out.append(pts[seg_idx])
        else:
            t = (target - seg_start) / seg_len
            t = max(0.0, min(1.0, t))
            lat = pts[seg_idx][0] + t * (pts[seg_idx + 1][0] - pts[seg_idx][0])
            lon = pts[seg_idx][1] + t * (pts[seg_idx + 1][1] - pts[seg_idx][1])
            out.append((lat, lon))
    return out


def nearest_index(pts: List[Pt], lat: float, lon: float) -> int:
    """Index of the point in *pts* nearest to (lat, lon)."""
    best_i = 0
    best_d = haversine(pts[0][0], pts[0][1], lat, lon)
    for i in range(1, len(pts)):
        d = haversine(pts[i][0], pts[i][1], lat, lon)
        if d < best_d:
            best_d = d
            best_i = i
    return best_i


def midline(chain_a: List[Pt], chain_b: List[Pt],
            start_lat: float, end_lat: float) -> List[Pt]:
    """
    Compute the midline between two chains within a latitude band.

    1. Extract sub-chains within [start_lat, end_lat] (with small buffer).
    2. Resample both to the same length.
    3. Average corresponding points.
    """
    buf = 0.00005  # ~5 m latitude buffer
    lo = min(start_lat, end_lat) - buf
    hi = max(start_lat, end_lat) + buf

    def extract(chain):
        sub = [(lat, lon) for lat, lon in chain if lo <= lat <= hi]
        if len(sub) < 2:
            return sub
        # Ensure south→north ordering
        if sub[-1][0] < sub[0][0]:
            sub.reverse()
        return sub

    sub_a = extract(chain_a)
    sub_b = extract(chain_b)

    if len(sub_a) < 2 or len(sub_b) < 2:
        # Fallback: straight line
        return []

    # Resample both at RESAMPLE_STEP_M so they have similar density
    ra = resample(sub_a, RESAMPLE_STEP_M)
    rb = resample(sub_b, RESAMPLE_STEP_M)

    # Make both the same length (resample shorter to match longer count)
    n = min(len(ra), len(rb))
    if len(ra) > len(rb):
        rb = resample(sub_b, chain_length(sub_b) / (len(ra) - 1)) if len(ra) > 1 else rb
        n = min(len(ra), len(rb))
    elif len(rb) > len(ra):
        ra = resample(sub_a, chain_length(sub_a) / (len(rb) - 1)) if len(rb) > 1 else ra
        n = min(len(ra), len(rb))
    n = min(len(ra), len(rb))

    mid = []
    for i in range(n):
        mid.append(((ra[i][0] + rb[i][0]) / 2,
                     (ra[i][1] + rb[i][1]) / 2))
    return mid


# ── Load data ────────────────────────────────────────────────────────────────
with open(ASSETS_FILE) as f:
    assets = json.load(f)

curbs_raw = [feat for feat in assets["features"]
             if feat["properties"].get("asset_type") == "CURB"]

# Convert each curb to list of (lat, lon)
curbs = []
for feat in curbs_raw:
    coords = feat["geometry"]["coordinates"]
    pts = [(c[1], c[0]) for c in coords]  # lon,lat → lat,lon
    closed = (len(pts) > 2
              and abs(pts[0][0] - pts[-1][0]) < 1e-6
              and abs(pts[0][1] - pts[-1][1]) < 1e-6)
    curbs.append({"pts": pts, "closed": closed})

print(f"Loaded {len(curbs)} curbs ({sum(1 for c in curbs if c['closed'])} closed)")

# ── Chain curbs by endpoint proximity ────────────────────────────────────────
open_idx = [i for i, c in enumerate(curbs) if not c["closed"]]

# Build adjacency: end[i] → start[j]
next_map = {}  # i → j
for i in open_idx:
    end = curbs[i]["pts"][-1]
    best_j, best_d = -1, 1e18
    for j in open_idx:
        if j == i:
            continue
        start = curbs[j]["pts"][0]
        d = haversine(*end, *start)
        if d < best_d:
            best_d = d
            best_j = j
    if best_d < CHAIN_GAP_M:
        next_map[i] = best_j

has_pred = set(next_map.values())
chain_starts = [i for i in open_idx if i not in has_pred]

chains = []
for s in chain_starts:
    chain = [s]
    cur = s
    visited = {s}
    while cur in next_map:
        nxt = next_map[cur]
        if nxt in visited:
            break
        chain.append(nxt)
        visited.add(nxt)
        cur = nxt
    # Concatenate all points
    pts = []
    for ci in chain:
        pts.extend(curbs[ci]["pts"])
    chains.append({"indices": chain, "pts": pts})

# Sort chains by total point count (the two big ones are the road sides)
chains.sort(key=lambda c: len(c["pts"]), reverse=True)

print(f"Built {len(chains)} chains:")
for i, ch in enumerate(chains):
    print(f"  Chain {i}: {len(ch['indices'])} curbs, {len(ch['pts'])} pts, "
          f"lat {ch['pts'][0][0]:.5f}→{ch['pts'][-1][0]:.5f}")

# The longest chain is one road side (Chain A).
# All other non-closed curbs form the other side (Chain B).
chain_a_set = set(chains[0]["indices"])
chain_a_pts = chains[0]["pts"]

# Ensure Chain A runs south → north
if chain_a_pts[-1][0] < chain_a_pts[0][0]:
    chain_a_pts.reverse()

# Gather all remaining open curb points (not in chain A) for Chain B
chain_b_pts: List[Pt] = []
for ch in chains[1:]:
    pts = ch["pts"]
    # Ensure each fragment runs south→north
    if len(pts) >= 2 and pts[-1][0] < pts[0][0]:
        pts = list(reversed(pts))
    chain_b_pts.extend(pts)

# Sort chain B by latitude so it runs south → north monotonically
chain_b_pts.sort(key=lambda p: p[0])

print(f"\nChain A: {len(chain_a_pts)} pts, lat {chain_a_pts[0][0]:.5f}→{chain_a_pts[-1][0]:.5f}")
print(f"Chain B: {len(chain_b_pts)} pts, lat {chain_b_pts[0][0]:.5f}→{chain_b_pts[-1][0]:.5f}")

# ── Load intersections ───────────────────────────────────────────────────────
rush_cols = ["7:00 AM", "7:15 AM", "7:30 AM", "7:45 AM",
             "5:00 PM", "5:15 PM", "5:30 PM", "5:45 PM"]

imap = {}
with open(TRAFFIC_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        name = row["Intersection name"]
        lat = float(row["latitude of intersection"])
        lon = float(row["longitude of intersection"])
        if name not in imap:
            imap[name] = {"lat": lat, "lon": lon, "total": 0}
        for col in rush_cols:
            v = int(row.get(col, 0) or 0)
            imap[name]["total"] += v

intersections = sorted(
    [{"name": n, **d} for n, d in imap.items()],
    key=lambda x: x["lat"]
)
print(f"\nIntersections ({len(intersections)}):")
for ix in intersections:
    print(f"  {ix['lat']:.6f}, {ix['lon']:.6f}  {ix['name']}  (traffic={ix['total']})")

# ── Build centerline segments between consecutive intersections ──────────────
features = []

for idx in range(len(intersections) - 1):
    A = intersections[idx]
    B = intersections[idx + 1]

    mid_pts = midline(chain_a_pts, chain_b_pts, A["lat"], B["lat"])

    if len(mid_pts) < 2:
        # Fallback: straight line
        print(f"  Segment {idx}: {A['name']} → {B['name']} — FALLBACK (straight line)")
        mid_pts = [(A["lat"], A["lon"]), (B["lat"], B["lon"])]
    else:
        # Snap first/last point to intersection coords
        mid_pts[0] = (A["lat"], A["lon"])
        mid_pts[-1] = (B["lat"], B["lon"])
        print(f"  Segment {idx}: {A['name']} → {B['name']} — {len(mid_pts)} centerline pts")

    # Build GeoJSON Feature
    features.append({
        "type": "Feature",
        "properties": {
            "segment_index": idx,
            "from_name": A["name"],
            "to_name": B["name"],
            "from_lat": A["lat"],
            "from_lon": A["lon"],
            "to_lat": B["lat"],
            "to_lon": B["lon"],
            "from_traffic": A["total"],
            "to_traffic": B["total"],
        },
        "geometry": {
            "type": "LineString",
            "coordinates": [[lon, lat] for lat, lon in mid_pts],
        },
    })

geojson = {
    "type": "FeatureCollection",
    "features": features,
}

with open(OUTPUT_FILE, "w") as f:
    json.dump(geojson, f)

total_pts = sum(len(feat["geometry"]["coordinates"]) for feat in features)
print(f"\nWrote {OUTPUT_FILE}: {len(features)} segments, {total_pts} total points")
