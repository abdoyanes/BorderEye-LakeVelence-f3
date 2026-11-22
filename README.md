# Border Eye – UAV-LiDAR analysis of the Lake Velence reed wetland (flight f3)

R code for the Border Eye course project (University of Sopron). It analyzes a UAV-LiDAR survey of
Lake Velence, Hungary, flown on 15 October 2026 with a Quantum Systems Trinity Pro and a Qube 640 LiDAR.

## Scripts

| Script | What it does | Run time |
| --- | --- | --- |
| `01_trajectory_analysis.R` | Reads the 200 Hz flight trajectory and derives survey lines, altitude, speed, roll/pitch/heading, navigation accuracy, and wind | ~1 min |
| `02_pointcloud_light.R` | Processes the whole point cloud in one pass (300 m tiles, in memory): noise removal, ground classification (CSF), and 1 m rasters of density, DTM, DSM, and vegetation height | ~30 min |

## Data (not included)

The data are too large for GitHub. Put the f3 files in one folder:

- `Velence_261015_f3.laz` – point cloud (81.8 million points)
- `Velence_261015_f3_trajectory.txt` – flight trajectory

Both files use EOV coordinates (EPSG:23700), and heights are in EOMA 1980.

## How to run

1. Install R and the packages: `lidR`, `RCSF`, `terra`, `sf`, `data.table`, `ggplot2`, `patchwork`, `future`.
2. Set your data folder at the top of each script:
   ```r
   BASE_DIR <- "G:/PhD Courses/7th Semester/Border_Eye"
   ```
3. Run `01_trajectory_analysis.R`, then `02_pointcloud_light.R`.

The results are written to `outputs/trajectory/` and `outputs/pointcloud/` inside the data folder: figures (PNG), tables (CSV), and rasters (GeoTIFF).

## Main results

- 14 survey lines, 81.7 m apart, flown about 118 m above the water at 17.7 m/s.
- Water level 121.35 m; reed covers 59.9 ha with a median height of 1.96 m above the water.
- Open water returns laser points only close to the flight lines, which leaves stripe-shaped gaps in the data.

## Note

The script `Start.R` declares the trajectory with `crs = 4237` (degrees). Use `crs = 23700` (EOV metres); otherwise, the bounding-box crop returns wrong points.

## Author

Abdelbagi Yanes Fadlalmwlla Adam—PhD candidate, University of Sopron
