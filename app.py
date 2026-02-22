from flask import Flask, send_from_directory
import os

app = Flask(__name__, static_folder=os.getcwd())

@app.route("/")
def index():
    return send_from_directory(app.static_folder, "rush_hour_distress_map.html")

@app.route("/pointcloud")
def pointcloud():
    return send_from_directory(app.static_folder, "potree_viewer.html")

@app.route("/<path:filename>")
def static_files(filename):
    return send_from_directory(app.static_folder, filename)

if __name__ == "__main__":
    app.run(debug=True, port=5000)
