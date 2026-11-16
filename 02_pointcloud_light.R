# =============================================================================
# Border Eye – Step 2 (light): full-area point-cloud products in ONE pass
# Every chunk is read once (thinned, 7 attributes only) and processed in memory:
# noise -> ground (CSF) -> DTM, DSM, normalised heights -> 1 m rasters.
# No intermediate LAZ files are written, so it is much faster and lighter than
# the old full pipeline. Starts with a coverage check of the LAZ files.
# Author : Abdelbagi Yanes Fadlalmwlla Adam
# Note   : save as UTF-8 (° ² symbols in labels)
# =============================================================================
pkgs <- c("lidR", "RCSF", "terra", "sf", "data.table", "ggplot2", "future")
miss <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(miss)) install.packages(miss)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

# ---- 0. Settings ------------------------------------------------------------
args     <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args)) args[1] else if (exists("BASE_DIR")) BASE_DIR else
  "G:/PhD Courses/7th Semester/Border_Eye"               # <-- set your data folder
LAS_DIR     <- base_dir   # folder with the f3 LAZ files (searched recursively)
LAS_PATTERN <- ""         # e.g. "f3" to use only files whose name contains f3
EPSG        <- 23700      # EOV; used only if the LAZ header has no CRS
RES         <- 1          # m, cell size of all rasters
TARGET_D    <- 25         # pts/m² kept after random thinning (speed)
CHUNK       <- 300        # m, chunk size
BUFFER      <- 15         # m, chunk buffer (avoids edge artefacts)
WORKERS     <- max(1, min(4, parallel::detectCores() - 1))   # lower it if RAM is short

out_dir <- file.path(base_dir, "outputs", "pointcloud")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fmt  <- function(x, d = 2) formatC(x, format = "f", digits = d, big.mark = ",")
step <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
t0   <- Sys.time()

# PNG sized to the map's aspect ratio (no large white margins)
png_map <- function(f, r, w = 8, panels = 1) {
  e <- terra::ext(r); asp <- (e[4] - e[3]) / (e[2] - e[1])
  png(file.path(out_dir, f), w, min(max((w / panels) * asp * 0.85 + 0.9, 3), 11), "in", res = 300)
}
# Clamp to the 2–98 % range so outliers take the end colours instead of being blank
clp <- function(r, p = c(0.02, 0.98)) {
  q <- quantile(terra::values(r, na.rm = TRUE), p); terra::clamp(r, q[1], q[2], values = TRUE)
}

# ---- 1. Files and coverage check --------------------------------------------
files <- list.files(LAS_DIR, "\\.(las|laz)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
files <- files[!grepl("/outputs/", normalizePath(files, winslash = "/", mustWork = FALSE))]
if (nzchar(LAS_PATTERN)) files <- files[grepl(LAS_PATTERN, basename(files))]
stopifnot("No .las/.laz files found - check LAS_DIR / LAS_PATTERN" = length(files) > 0)

ctg <- readLAScatalog(files)
if (is.na(st_crs(ctg))) st_crs(ctg) <- EPSG
hdr <- as.data.table(ctg@data)[, .(file = basename(filename), points = Number.of.point.records,
                                   xmin = Min.X, xmax = Max.X, ymin = Min.Y, ymax = Max.Y)]
print(hdr)
fwrite(hdr, file.path(out_dir, "P0_files.csv"))
area_m2  <- as.numeric(st_area(ctg))
dens_hdr <- npoints(ctg) / area_m2
step("Catalogue:", length(files), "file(s),", fmt(npoints(ctg), 0), "points,",
     fmt(area_m2 / 1e4, 1), "ha,", fmt(dens_hdr, 1), "pts/m²")

# Compare with the trajectory-based survey area (from 04_survey_boundary.R)
survey_shp <- file.path(base_dir, "outputs", "boundary", "survey_area_trajectory.shp")
lines_shp  <- file.path(base_dir, "outputs", "boundary", "survey_lines.shp")
survey <- if (file.exists(survey_shp)) st_read(survey_shp, quiet = TRUE) else NULL
tiles  <- st_as_sf(ctg)
if (!is.null(survey)) {
  cov_pct <- 100 * as.numeric(st_area(st_intersection(st_union(tiles), st_union(survey)))) /
             as.numeric(st_area(survey))
  step("LAZ files cover", fmt(cov_pct, 1), "% of the trajectory-based survey area")
  if (cov_pct < 80) warning("The LAZ files cover only part of the survey area - check LAS_DIR / LAS_PATTERN")
}
png(file.path(out_dir, "P1_coverage.png"), 6, 8, "in", res = 200)
plot(st_geometry(if (!is.null(survey)) st_union(st_union(tiles), st_union(survey)) else tiles),
     border = NA, axes = TRUE, main = "Coverage: LAZ tiles (grey) vs survey area (blue)")
plot(st_geometry(tiles), col = "grey85", border = "grey40", add = TRUE)
if (!is.null(survey)) plot(st_geometry(survey), border = "#2c7fb8", lwd = 2, add = TRUE)
if (file.exists(lines_shp)) plot(st_geometry(st_read(lines_shp, quiet = TRUE)), col = "black", add = TRUE)
invisible(dev.off())

# ---- 2. One-pass processing -------------------------------------------------
frac <- min(1, TARGET_D / dens_hdr)
opt_select(ctg)       <- "xyzirnc"           # coordinates, intensity, returns, class
opt_filter(ctg)       <- paste("-drop_withheld", if (frac < 1) paste("-keep_random_fraction", round(frac, 4)))
opt_chunk_size(ctg)   <- CHUNK
opt_chunk_buffer(ctg) <- BUFFER
opt_progress(ctg)     <- TRUE
plan(multisession, workers = WORKERS)
set_lidr_threads(1)
step("Processing in", CHUNK, "m chunks with", WORKERS, "worker(s); thinning fraction", round(frac, 3))

process_chunk <- function(las, res, frac) {
  if (npoints(las) < 500) return(NULL)
  ca   <- res^2 * frac                                   # converts counts to true pts/m²
  tmpl <- pixel_metrics(las, ~list(density = length(Z)), res = res)
  fit  <- function(x) terra::crop(terra::extend(x, tmpl), tmpl)
  las  <- classify_noise(las, ivf(res = 3, n = 10))
  noise <- fit(pixel_metrics(las, ~list(noise = sum(Classification == 18L)), res = res))
  las  <- filter_poi(las, Classification != 18L & Classification != 7L)
  las  <- classify_ground(las, csf(sloop_smooth = FALSE, class_threshold = 0.3,
                                   cloth_resolution = 1, rigidness = 3L), last_returns = FALSE)
  n_ground <- sum(las$Classification == 2L)
  dsm  <- fit(rasterize_canopy(las, res = res, algorithm = p2r(0.2)))
  names(dsm) <- "dsm"
  gcnt <- fit(pixel_metrics(las, ~list(ground = sum(Classification == 2L)), res = res))
  if (n_ground < 10) {                                   # e.g. a chunk of open water only
    na  <- tmpl * NA
    out <- c(tmpl / ca, noise / ca, gcnt / ca, dsm, setNames(na, "dtm"),
             setNames(na, "h95"), setNames(na, "hmax"), setNames(na, "below"))
    return(out)
  }
  dtm  <- rasterize_terrain(las, res = res, algorithm = tin())
  nlas <- normalize_height(las, dtm)
  hm   <- fit(pixel_metrics(nlas, ~list(h95 = quantile(Z, 0.95), hmax = max(Z),
                                        below = sum(Z < -0.3)), res = res))
  dtm  <- fit(dtm); names(dtm) <- "dtm"
  c(tmpl / ca, noise / ca, gcnt / ca, dsm, dtm, hm$h95, hm$hmax, hm$below / ca)
}

R <- catalog_map(ctg, process_chunk, res = RES, frac = frac)
names(R) <- c("density", "noise", "ground", "dsm", "dtm", "h95", "hmax", "below")
terra::writeRaster(R, file.path(out_dir, "pointcloud_layers_1m.tif"), overwrite = TRUE)
step("Rasters done:", terra::ncell(R), "cells")

# ---- 3. Derived layers ------------------------------------------------------
dens  <- R$density
med_d <- terra::global(dens, median, na.rm = TRUE)[[1]]
foot  <- if (!is.null(survey)) terra::rasterize(terra::vect(survey), dens) else !is.na(dens)
foot  <- terra::ifel(foot == 1, 1, NA)
d5    <- terra::resample(terra::aggregate(terra::ifel(is.na(dens), 0, dens), 5, mean, na.rm = TRUE), dens, "near")
drop  <- terra::mask(terra::ifel(d5 < 0.1 * med_d, 1, 0), foot)
slope <- terra::terrain(R$dtm, "slope", unit = "degrees")
hs    <- terra::shade(terra::terrain(R$dtm, "slope", unit = "radians"),
                      terra::terrain(R$dtm, "aspect", unit = "radians"), 45, 315)
h95   <- terra::clamp(R$h95, lower = 0, values = TRUE)
hcls  <- terra::classify(h95, rbind(c(-Inf, 0.5, 1), c(0.5, 2, 2), c(2, 4, 3), c(4, Inf, 4)))
hcls  <- terra::ifel(terra::resample(drop, hcls, "near") == 1, 0, hcls)
levels(hcls) <- data.frame(id = 0:4, class = c("No returns (water)", "< 0.5 m", "0.5–2 m", "2–4 m", "> 4 m"))
for (n in c("drop", "slope", "h95", "hcls"))
  terra::writeRaster(get(n), file.path(out_dir, paste0(n, "_1m.tif")), overwrite = TRUE)

# ---- 4. Figures -------------------------------------------------------------
png_map("P2_density_dropout.png", dens, 11, 2); par(mfrow = c(1, 2))
terra::plot(clp(dens), col = hcl.colors(50, "viridis"), main = "Point density (pts/m²)")
terra::plot(drop, col = c("grey90", "#2166ac"), legend = FALSE, main = "Return dropout (blue = no returns)")
invisible(dev.off())

png_map("P3_DTM_DSM.png", R$dtm, 11, 2); par(mfrow = c(1, 2))
terra::plot(hs, col = gray.colors(50, 0.3, 1), legend = FALSE, main = "DTM (m a.s.l.) on hillshade")
terra::plot(clp(R$dtm), col = hcl.colors(50, "Terrain2", alpha = 0.55), add = TRUE)
terra::plot(clp(R$dsm), col = hcl.colors(50, "Terrain2"), main = "DSM (m a.s.l.)")
invisible(dev.off())

png_map("P4_slope.png", slope)
terra::plot(clp(slope, c(0, 0.98)), col = hcl.colors(50, "YlOrRd", rev = TRUE), main = "Terrain slope (°)")
invisible(dev.off())

png_map("P5_reed_height.png", h95, 12, 2); par(mfrow = c(1, 2))
terra::plot(clp(h95, c(0, 0.99)), col = hcl.colors(50, "viridis"), main = "Vegetation height, 95th percentile (m)")
terra::plot(hcls, col = c("#2166ac", "#f7f7f7", "#d9f0a3", "#78c679", "#006837"), main = "Height classes",
            mar = c(3.1, 3.1, 2.1, 9), plg = list(cex = 0.8))
invisible(dev.off())

veg <- terra::values(h95, na.rm = TRUE); veg <- veg[veg >= 0.5]
ggsave(file.path(out_dir, "P6_height_histogram.png"),
       ggplot(data.frame(h = veg), aes(h)) + geom_histogram(binwidth = 0.1, fill = "#1b7837", colour = NA) +
         geom_vline(xintercept = median(veg), linetype = 2) +
         annotate("text", x = median(veg), y = Inf, vjust = 1.5, hjust = -0.05,
                  label = paste("median", fmt(median(veg)), "m")) +
         labs(x = "Vegetation height (m), cells ≥ 0.5 m", y = "1 m cells",
              title = "Vegetation (reed) height distribution") + theme_bw(),
       width = 6, height = 3.8, dpi = 300)

# Cross-sections: three 4 m wide W–E strips at 25 %, 50 % and 75 % of the extent
e  <- terra::ext(dens)
ys <- e[3] + c(0.25, 0.5, 0.75) * (e[4] - e[3])
sec <- rbindlist(lapply(seq_along(ys), function(i) {
  s <- clip_rectangle(ctg, e[1], ys[i] - 2, e[2], ys[i] + 2)
  if (is.empty(s) || npoints(s) < 100) return(NULL)
  s <- classify_noise(s, ivf(res = 3, n = 10)); s <- filter_poi(s, Classification != 18L)
  s <- classify_ground(s, csf(sloop_smooth = FALSE, class_threshold = 0.3, cloth_resolution = 1, rigidness = 3L))
  data.table(x = s$X, z = s$Z, cls = ifelse(s$Classification == 2L, "Ground", "Non-ground"),
             prof = sprintf("Profile %s (EOV Y = %s)", LETTERS[i], fmt(ys[i], 0)))
}))
if (nrow(sec)) ggsave(file.path(out_dir, "P7_cross_sections.png"),
  ggplot(sec, aes(x, z, colour = cls)) + geom_point(size = 0.1) + facet_wrap(~prof, ncol = 1) +
    scale_colour_manual(values = c("Ground" = "#8c510a", "Non-ground" = "#1b7837"), name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(x = "EOV X (m), west → east", y = "Z (m a.s.l.)", title = "W–E cross-sections (4 m wide)") +
    theme_bw(), width = 10, height = 7, dpi = 300)

# ---- 5. Summary table -------------------------------------------------------
gsum <- function(r, f = "sum") terra::global(r, f, na.rm = TRUE)[[1]]
v    <- function(r) terra::values(r, na.rm = TRUE)
pts  <- gsum(R$density) * RES^2
cls_area <- as.data.table(terra::freq(hcls))[, .(class = value, ha = count * RES^2 / 1e4)]
cls_area[, pct := 100 * ha / sum(ha)]
summ <- data.table(
  Metric = c("Files / points", "Catalogue area (ha)", "Point density, mean / median (pts/m²)",
             "Return dropout in survey area (%)", "Noise points (%)", "Ground points (%)",
             "Points > 0.3 m below DTM (%)", "DTM median / range (m a.s.l.)",
             "DSM range (m a.s.l.)", "Vegetated cells ≥ 0.5 m (%)",
             "Vegetation height median / mean ± SD (m)", "Vegetation height P95 / max (m)",
             "Thinning fraction / run time (min)"),
  Value = c(paste(length(files), "/", fmt(npoints(ctg), 0)), fmt(area_m2 / 1e4, 1),
            paste(fmt(gsum(dens, "mean"), 1), "/", fmt(med_d, 1)),
            fmt(100 * gsum(drop, "mean"), 1),
            fmt(100 * gsum(R$noise) * RES^2 / pts, 3), fmt(100 * gsum(R$ground) * RES^2 / pts, 1),
            fmt(100 * gsum(R$below) * RES^2 / pts, 2),
            paste(fmt(median(v(R$dtm))), "/", fmt(min(v(R$dtm))), "–", fmt(max(v(R$dtm)))),
            paste(fmt(min(v(R$dsm))), "–", fmt(max(v(R$dsm)))),
            fmt(100 * length(veg) / length(v(h95)), 1),
            paste0(fmt(median(veg)), " / ", fmt(mean(veg)), " ± ", fmt(sd(veg))),
            paste0(fmt(quantile(veg, 0.95)), " / ", fmt(max(v(R$hmax)))),
            paste(round(frac, 3), "/", fmt(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))))
fwrite(summ, file.path(out_dir, "PT1_summary.csv"))
fwrite(cls_area, file.path(out_dir, "PT2_height_classes.csv"))
# Surface level for 04_survey_boundary.R and 05_scan_geometry.R (replaces the assumed 104 m)
fwrite(data.table(surface_z = median(v(R$dtm))), file.path(out_dir, "surface_level.csv"))
print(summ); print(cls_area)
plan(sequential)
step("Finished. Outputs in", out_dir)
