"""
viewer_app.py — Flask server for the Potree point cloud + asset viewer.

Usage:
    python viewer_app.py

Then open http://localhost:5001/viewer.html in your browser.

Prerequisites:
    - pip install flask
    - A Potree installation at ../potree (sibling to civ2).
      The POTREE_PATH in viewer.html should be "../potree".
    - Run spatial_join.py first to generate assets_with_tiles.geojson.

The app serves civ2 project files and the sibling potree/ directory so
the viewer can load all its dependencies.
"""

from flask import Flask, send_from_directory, abort
import os

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(PROJECT_DIR)  # one level up (Projects/)

app = Flask(__name__)


@app.route("/")
def index():
    """Serve the viewer as the landing page."""
    return send_from_directory(PROJECT_DIR, "viewer.html")


@app.route("/<path:filename>")
def serve_file(filename):
    """
    Serve files from the project dir.
    Because POTREE_PATH = "../potree", the browser resolves requests to
    /potree/... (relative to the site root).  We handle that by serving
    from the parent directory.
    """
    # First, try inside the project directory (civ2/)
    local = os.path.normpath(os.path.join(PROJECT_DIR, filename))
    if local.startswith(PROJECT_DIR) and os.path.isfile(local):
        return send_from_directory(os.path.dirname(local), os.path.basename(local))

    # Then, try from the parent directory (for ../potree/... → /potree/...)
    parent = os.path.normpath(os.path.join(PARENT_DIR, filename))
    if parent.startswith(PARENT_DIR) and os.path.isfile(parent):
        return send_from_directory(os.path.dirname(parent), os.path.basename(parent))

    abort(404)


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

    potree_dir = os.path.join(PARENT_DIR, "potree")
    if os.path.isdir(potree_dir):
        potree_js = os.path.join(potree_dir, "build", "potree", "potree.js")
        if os.path.isfile(potree_js):
            print(f"  Potree found at {potree_dir} ✓")
        else:
            print(f"  WARNING: Potree dir exists but build/potree/potree.js is missing — run the Potree build")
    else:
        print(f"  WARNING: Potree not found at {potree_dir}")

    port = 5001
    print(f"\n  Serving viewer at http://localhost:{port}/viewer.html")
    print(f"  Project dir: {PROJECT_DIR}")
    print(f"  Parent dir:  {PARENT_DIR}\n")
    app.run(debug=True, port=port)
