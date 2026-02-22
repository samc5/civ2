"""
viewer_app.py — Flask server for the Potree point cloud + asset viewer.

Usage:
    python viewer_app.py

Then open http://localhost:5001/viewer.html in your browser.

Prerequisites:
    - pip install flask
    - A local Potree 1.8 installation. Set POTREE_PATH in viewer.html
      to point at it (default: "../potree-1.8").
    - Run spatial_join.py first to generate assets_with_tiles.geojson.

The app serves all project files (GeoJSON, JSON, HTML) from the working
directory so the viewer can fetch them.
"""

from flask import Flask, send_from_directory
import os
import sys

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

app = Flask(__name__, static_folder=PROJECT_DIR)


@app.route("/")
def index():
    """Redirect root to the viewer."""
    return send_from_directory(PROJECT_DIR, "viewer.html")


@app.route("/<path:filename>")
def serve_file(filename):
    """Serve any file from the project directory (and parent for Potree)."""
    # Allow traversal up one level so that ../potree-1.8/... works
    full = os.path.normpath(os.path.join(PROJECT_DIR, filename))
    directory = os.path.dirname(full)
    basename = os.path.basename(full)
    return send_from_directory(directory, basename)


if __name__ == "__main__":
    # Quick sanity checks
    required = ["pointcloudjp.json", "viewer.html"]
    optional = ["assets_with_tiles.geojson"]

    for f in required:
        if not os.path.exists(os.path.join(PROJECT_DIR, f)):
            print(f"WARNING: {f} not found in {PROJECT_DIR}")

    for f in optional:
        path = os.path.join(PROJECT_DIR, f)
        if not os.path.exists(path):
            print(f"NOTE: {f} not found — run 'python spatial_join.py' first to generate it.")

    port = 5001
    print(f"\n  Serving viewer at http://localhost:{port}/viewer.html")
    print(f"  Project dir: {PROJECT_DIR}\n")
    app.run(debug=True, port=port)
