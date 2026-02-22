"""
spatial_join.py — Join JP assets to point cloud tile polygons.

Reads:  abovegroundjp.geojson   (1208 asset features, WGS84)
        pointcloudjp.json        (47 tile polygon features, WGS84)

Writes: assets_with_tiles.geojson
        Adds properties to each asset feature:
          tile_file      — LAZ filename  (e.g. global_xyz_rgb_icgu_1_0_2000.laz)
          las_url        — full CDN URL  (baseUrl + lasPath)
          street_dataset — dataset name from tile props

Strategy:
  Points     — point-in-polygon test; fallback to nearest tile centroid
  LineStrings — midpoint coordinate used for the join
"""

import json
import math
from shapely.geometry import shape, Point


def midpoint(coords):
    """Return midpoint of a LineString coordinate list."""
    n = len(coords)
    if n == 1:
        return coords[0]
    # Use fractional mid-index interpolation
    half = (n - 1) / 2.0
    lo = int(half)
    hi = lo + 1 if lo + 1 < n else lo
    frac = half - lo
    lon = coords[lo][0] + frac * (coords[hi][0] - coords[lo][0])
    lat = coords[lo][1] + frac * (coords[hi][1] - coords[lo][1])
    return [lon, lat]


def get_test_point(feature):
    """Return [lon, lat] representative point for a feature."""
    geom = feature["geometry"]
    gtype = geom["type"]
    if gtype == "Point":
        return geom["coordinates"][:2]
    elif gtype == "LineString":
        return midpoint(geom["coordinates"])
    elif gtype == "MultiPoint":
        return geom["coordinates"][0][:2]
    elif gtype == "MultiLineString":
        return midpoint(geom["coordinates"][0])
    else:
        # Fallback: try to extract something usable
        coords = geom.get("coordinates", [])
        if coords:
            c = coords[0]
            while isinstance(c[0], list):
                c = c[0]
            return c[:2]
    return None


def nearest_tile(lon, lat, tile_shapes):
    """Return tile with centroid nearest to (lon, lat)."""
    min_dist = float("inf")
    nearest = None
    for tile in tile_shapes:
        cx, cy = tile["centroid"]
        dist = math.sqrt((lon - cx) ** 2 + (lat - cy) ** 2)
        if dist < min_dist:
            min_dist = dist
            nearest = tile
    return nearest


def main():
    # --- Load data ---
    with open("abovegroundjp.geojson", encoding="utf-8") as f:
        assets_fc = json.load(f)
    with open("pointcloudjp.json", encoding="utf-8") as f:
        tiles_fc = json.load(f)

    # --- Build tile index ---
    tile_shapes = []
    for tf in tiles_fc["features"]:
        props = tf["properties"]
        tile_shapes.append(
            {
                "shape": shape(tf["geometry"]),
                "props": props,
                "centroid": (props["lon"], props["lat"]),
            }
        )

    # --- Join ---
    tagged = 0
    untagged = 0
    fallback = 0
    output_features = []

    for feature in assets_fc["features"]:
        test_pt = get_test_point(feature)

        new_feature = {
            "type": "Feature",
            "geometry": feature["geometry"],
            "properties": dict(feature.get("properties") or {}),
        }

        if test_pt is None:
            untagged += 1
            output_features.append(new_feature)
            continue

        lon, lat = test_pt[0], test_pt[1]
        pt = Point(lon, lat)

        # Primary: containment test
        matched_tile = None
        for tile in tile_shapes:
            if tile["shape"].contains(pt):
                matched_tile = tile
                break

        # Fallback: nearest centroid
        if matched_tile is None:
            matched_tile = nearest_tile(lon, lat, tile_shapes)
            fallback += 1

        if matched_tile is not None:
            p = matched_tile["props"]
            new_feature["properties"]["tile_file"] = p["filename"]
            new_feature["properties"]["las_url"] = p["baseUrl"] + p["lasPath"]
            new_feature["properties"]["street_dataset"] = p["dataset"]
            tagged += 1
        else:
            untagged += 1

        output_features.append(new_feature)

    # --- Write output ---
    output = {
        "type": "FeatureCollection",
        "features": output_features,
    }
    with open("assets_with_tiles.geojson", "w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False)

    print(f"Tagged {tagged} / Untagged {untagged} (fallback-nearest: {fallback})")
    print(f"Wrote assets_with_tiles.geojson with {len(output_features)} features")


if __name__ == "__main__":
    main()
