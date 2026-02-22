# civ2

### How to use Inject assets

# Download tile #0 from CDN and inject nearby assets
python inject_assets.py --tile 0

# Process a local LAZ file you already downloaded
python inject_assets.py path/to/scan.laz

# Custom search radius (default 200m)
python inject_assets.py --tile 0 --radius 300

# Process all 47 tiles (warning: ~14 GB of downloads)
python inject_assets.py --all