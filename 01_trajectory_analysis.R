# =============================================================================
# Border Eye – Part A: Flight trajectory (SBET-type) analysis
# Data   : Velence_261015_f3_trajectory.txt  (Lake Velence, 15 Oct 2026, flight 3)
# CRS    : EOV / HD72 (EPSG:23700); Z as delivered (m)
# Author : Abdelbagi Yanes Fadlalmwlla Adam
# Note   : save this file as UTF-8 (° ± symbols in labels)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

# ---- 0. Paths ---------------------------------------------------------------
base_dir <- "G:/PhD Courses/7th Semester/Border_Eye"
if (!dir.exists(base_dir)) base_dir <- getwd()           # fallback for testing
traj_file <- file.path(base_dir, "Velence_261015_f3_trajectory.txt")
out_dir   <- file.path(base_dir, "outputs", "trajectory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_bw(base_size = 11) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", size = 11)))

# ---- 1. Read ----------------------------------------------------------------
tr <- fread(traj_file)
setnames(tr, c("Time","X","Y","Z","Roll","Pitch","Heading",
               "sdX","sdY","sdZ","sdRoll","sdPitch","sdHeading"))
wrap180 <- function(a) ((a + 180) %% 360) - 180

# GPS seconds-of-week -> clock time (date taken from file name, GPS-UTC = 18 s)
flight_date <- as.Date("2026-10-15")
week_start  <- flight_date - as.integer(format(flight_date, "%w"))  # Sunday
t0_utc <- as.POSIXct(paste(week_start, "00:00:00"), tz = "UTC") + tr$Time[1] - 18
t1_utc <- t0_utc + diff(range(tr$Time))

# ---- 2. Kinematics ----------------------------------------------------------
dt <- diff(tr$Time)
tr[, `:=`(
  t_rel  = Time - Time[1],
  v_h    = c(NA, sqrt(diff(X)^2 + diff(Y)^2) / dt),       # m/s
  v_z    = c(NA, diff(Z) / dt),                           # m/s
  course = c(NA, atan2(diff(X), diff(Y)) * 180 / pi),     # deg from grid north
  sd_h   = sqrt(sdX^2 + sdY^2)                            # horizontal 1-sigma
)]

# ---- 3. Flight-line segmentation (on 1 Hz subsample) ------------------------
# A survey line = sustained flight (>= 60 s) with low bank angle along the
# dominant survey axis; transit legs and turns are excluded.
s <- tr[seq(1, .N, by = 200)]
k <- 10                                                   # 10-s smoothing window
s[, crs := atan2(shift(X, -k/2) - shift(X, k/2), shift(Y, -k/2) - shift(Y, k/2)) * 180 / pi]
dom <- atan2(median(sin(2 * s$crs * pi/180), na.rm = TRUE),
             median(cos(2 * s$crs * pi/180), na.rm = TRUE)) / 2 * 180 / pi
ax_dev <- function(a) abs(wrap180(2 * (a - dom))) / 2     # axial deviation (0-90)
s[, straight := !is.na(crs) & ax_dev(crs) < 15 & abs(Roll) < 10]
s[, run := rleid(straight)]
runs <- s[straight == TRUE, .(t_start = min(Time), t_end = max(Time),
                              dur = max(Time) - min(Time)), by = run][dur >= 60]

tr[, line_id := NA_integer_]
for (i in seq_len(nrow(runs)))
  tr[Time >= runs$t_start[i] & Time <= runs$t_end[i], line_id := runs$run[i]]

lines <- tr[!is.na(line_id), .(
  t0_sow   = min(Time), t1_sow = max(Time),                # GPS seconds of week
  t_start  = min(t_rel), duration_s = max(Time) - min(Time),
  length_m = sum(sqrt(diff(X)^2 + diff(Y)^2)),
  course   = atan2(last(X) - first(X), last(Y) - first(Y)) * 180 / pi,
  heading  = atan2(mean(sin(Heading*pi/180)), mean(cos(Heading*pi/180)))*180/pi,
  v_mean   = mean(v_h, na.rm = TRUE),
  z_mean   = mean(Z), z_sd = sd(Z),
  roll_sd  = sd(Roll), pitch_mean = mean(Pitch),
  xm = mean(X), ym = mean(Y),
  x0 = first(X), y0 = first(Y), x1 = last(X), y1 = last(Y)), by = line_id]
lines <- lines[length_m >= 500]                           # keep survey lines
# Heading convention check: if the IMU heading is ~180 deg from the direction of
# travel, the IMU x-axis points aft; correct before computing crab angle.
yaw_offset <- if (median(abs(wrap180(lines$heading - lines$course))) > 90) 180 else 0
lines[, crab := wrap180(heading - course - yaw_offset)]
setorder(lines, t_start); lines[, line_no := .I]
tr[!line_id %in% lines$line_id, line_id := NA]

# Wind estimate from the wind triangle: Vg = Va*(sin psi, cos psi) + W
# (linear least squares, on-line samples at 1 Hz; assumes constant airspeed)
w <- tr[!is.na(line_id)][seq(1, .N, by = 200)]
w[, `:=`(psi = (Heading - yaw_offset) * pi / 180,
         vgx = c(NA, diff(X)) / c(NA, diff(Time)), vgy = c(NA, diff(Y)) / c(NA, diff(Time)))]
w <- w[-1][is.finite(vgx)]
A <- rbind(cbind(sin(w$psi), 1, 0), cbind(cos(w$psi), 0, 1))
b <- c(w$vgx, w$vgy)
wind_fit <- qr.solve(A, b)                                # Va, Wx, Wy
Va <- wind_fit[1]; Wx <- wind_fit[2]; Wy <- wind_fit[3]
wind_speed <- sqrt(Wx^2 + Wy^2)
wind_from  <- (atan2(-Wx, -Wy) * 180 / pi) %% 360         # meteorological direction

# Line spacing: perpendicular offsets of parallel lines along dominant axis
ax  <- atan2(mean(sin(2*lines$course*pi/180)), mean(cos(2*lines$course*pi/180)))/2
ux  <- sin(ax); uy <- cos(ax)                             # along-track unit vector
lines[, perp := xm * uy - ym * ux]                        # cross-track coordinate
perp_sorted <- sort(lines$perp)
spacing <- diff(perp_sorted)
spacing <- spacing[spacing > 5]                           # drop re-flown lines

# ---- 4. Summary tables ------------------------------------------------------
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d, big.mark = ",")
on_line <- !is.na(tr$line_id)
summ <- data.table(
  Parameter = c("Records", "Sampling rate", "Start (UTC)", "End (UTC)", "Duration",
                "2-D path length", "Extent X (EOV)", "Extent Y (EOV)",
                "Take-off elevation", "Survey altitude (on lines), mean ± SD",
                "Ground speed (on lines), mean ± SD", "Survey lines detected",
                "Mean line length", "Dominant line azimuth",
                "Median line spacing", "Time on survey lines",
                "Heading–course offset (IMU axis)", "Mean |crab angle|",
                "Estimated airspeed", "Estimated wind (from / speed)"),
  Value = c(fmt(nrow(tr), 0),
            paste0(round(1 / median(dt)), " Hz"),
            format(t0_utc, "%Y-%m-%d %H:%M:%S"), format(t1_utc, "%H:%M:%S"),
            paste0(fmt(diff(range(tr$Time)) / 60, 1), " min"),
            paste0(fmt(sum(sqrt(diff(tr$X)^2 + diff(tr$Y)^2)) / 1000, 1), " km"),
            paste(fmt(min(tr$X), 0), "–", fmt(max(tr$X), 0), "m"),
            paste(fmt(min(tr$Y), 0), "–", fmt(max(tr$Y), 0), "m"),
            paste0(fmt(tr$Z[1], 1), " m"),
            paste0(fmt(mean(tr$Z[on_line]), 1), " ± ", fmt(sd(tr$Z[on_line]), 1), " m"),
            paste0(fmt(mean(tr$v_h[on_line], na.rm = TRUE), 1), " ± ",
                   fmt(sd(tr$v_h[on_line], na.rm = TRUE), 1), " m/s"),
            nrow(lines),
            paste0(fmt(mean(lines$length_m), 0), " m"),
            paste0(fmt((ax * 180 / pi) %% 180, 1), "° / ",
                   fmt((ax * 180 / pi) %% 180 + 180, 1), "°"),
            paste0(fmt(median(spacing), 1), " m"),
            paste0(fmt(100 * mean(on_line), 1), " %"),
            paste0(yaw_offset, "°"),
            paste0(fmt(mean(abs(lines$crab)), 1), "°"),
            paste0(fmt(Va, 1), " m/s"),
            paste0(fmt(wind_from, 0), "° / ", fmt(wind_speed, 1), " m/s")))

acc <- tr[, .(Variable = c("sdX (m)","sdY (m)","sdZ (m)","sdRoll (°)","sdPitch (°)","sdHeading (°)"),
              Mean = sapply(.SD, mean), Median = sapply(.SD, median),
              P95 = sapply(.SD, quantile, 0.95), Max = sapply(.SD, max)),
          .SDcols = c("sdX","sdY","sdZ","sdRoll","sdPitch","sdHeading")]

att <- tr[on_line, .(Variable = c("Roll (°)","Pitch (°)","Vertical speed (m/s)","Ground speed (m/s)"),
                     Mean = c(mean(Roll), mean(Pitch), mean(v_z, na.rm=TRUE), mean(v_h, na.rm=TRUE)),
                     SD   = c(sd(Roll), sd(Pitch), sd(v_z, na.rm=TRUE), sd(v_h, na.rm=TRUE)),
                     Min  = c(min(Roll), min(Pitch), min(v_z, na.rm=TRUE), min(v_h, na.rm=TRUE)),
                     Max  = c(max(Roll), max(Pitch), max(v_z, na.rm=TRUE), max(v_h, na.rm=TRUE)))]

fwrite(summ,  file.path(out_dir, "T1_trajectory_summary.csv"))
fwrite(acc,   file.path(out_dir, "T2_trajectory_accuracy.csv"))
fwrite(att,   file.path(out_dir, "T3_attitude_on_lines.csv"))
fwrite(lines[, .(line_no, t_start, duration_s, length_m, course, heading, crab,
                 v_mean, z_mean, z_sd, roll_sd)],
       file.path(out_dir, "T4_flight_lines.csv"))
saveRDS(list(flight_date = flight_date, wind = c(Va = Va, Wx = Wx, Wy = Wy), yaw_offset = yaw_offset, summary = summ, accuracy = acc, attitude = att, lines = lines,
             spacing = spacing), file.path(out_dir, "trajectory_results.rds"))
print(summ); print(acc); print(att)
cat("\nCrab by direction:\n"); print(lines[, .(crab = mean(crab), v = mean(v_mean)), by = .(dir = ifelse(cos(course*pi/180) > 0, "northbound", "southbound"))])
cat("Line spacing (m):", paste(fmt(spacing, 1), collapse = ", "), "\n")

# ---- 5. Figures -------------------------------------------------------------
sub <- tr[seq(1, .N, by = 20)]                            # 10 Hz for plotting

# Fig A: plan view of the trajectory coloured by altitude, survey lines numbered
fA <- ggplot(sub, aes(X, Y)) +
  geom_path(aes(colour = Z), linewidth = 0.5) +
  geom_segment(data = lines, aes(x = x0, y = y0, xend = x1, yend = y1),
               colour = "black", linewidth = 0.25, alpha = 0.5) +
  geom_label(data = lines, aes(xm, ym, label = line_no), size = 2.3,
             label.padding = unit(0.08, "lines")) +
  geom_point(data = tr[1], aes(X, Y), shape = 24, fill = "red", size = 3) +
  scale_colour_viridis_c(name = "Z (m)") + coord_equal() +
  labs(x = "EOV X / Easting (m)", y = "EOV Y / Northing (m)",
       title = "Flight trajectory and detected survey lines")
ggsave(file.path(out_dir, "Fig_T1_trajectory_map.png"), fA, width = 6.5, height = 7.5, dpi = 300)

# Fig B: time series
ts <- sub[, .(t_min = t_rel / 60, Z, v_h, Roll, Pitch, Heading, on = !is.na(line_id))]
p1 <- ggplot(ts, aes(t_min, Z)) + geom_line(linewidth = 0.3) +
  geom_point(data = ts[on == TRUE], colour = "#2c7fb8", size = 0.05) +
  labs(x = NULL, y = "Z (m)")
p2 <- ggplot(ts, aes(t_min, v_h)) + geom_line(linewidth = 0.3) + labs(x = NULL, y = "Speed (m/s)")
p3 <- ggplot(melt(ts[, .(t_min, Roll, Pitch)], id.vars = "t_min"),
             aes(t_min, value, colour = variable)) + geom_line(linewidth = 0.3) +
  scale_colour_manual(values = c("#d95f02", "#1b9e77"), name = NULL) +
  labs(x = NULL, y = "Angle (°)") + theme(legend.position = "right")
p4 <- ggplot(ts, aes(t_min, Heading)) + geom_point(size = 0.05) +
  labs(x = "Time since start (min)", y = "Heading (°)")
fB <- (p1 / p2 / p3 / p4) + plot_annotation(title = "Platform kinematics and attitude (blue = on survey line)")
ggsave(file.path(out_dir, "Fig_T2_kinematics.png"), fB, width = 7.5, height = 8, dpi = 300)

# Fig C: navigation accuracy (1-s means)
s1 <- tr[, lapply(.SD, mean), by = .(t_min = floor(t_rel) / 60),
         .SDcols = c("sdX","sdY","sdZ","sdRoll","sdPitch","sdHeading")]
c1 <- ggplot(melt(s1[, .(t_min, sdX, sdY, sdZ)], id.vars = "t_min"),
             aes(t_min, value * 1000, colour = variable)) + geom_line(linewidth = 0.35) +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3"), name = NULL) +
  labs(x = NULL, y = "Position σ (mm)")
c2 <- ggplot(melt(s1[, .(t_min, sdRoll, sdPitch, sdHeading)], id.vars = "t_min"),
             aes(t_min, value, colour = variable)) + geom_line(linewidth = 0.35) +
  scale_colour_manual(values = c("#e7298a", "#66a61e", "#e6ab02"), name = NULL) +
  labs(x = "Time since start (min)", y = "Attitude σ (°)")
fC <- (c1 / c2) + plot_annotation(title = "Estimated navigation uncertainty (1σ, 1-s means)")
ggsave(file.path(out_dir, "Fig_T3_accuracy.png"), fC, width = 7.5, height = 5, dpi = 300)

# Fig D: line geometry – course rose, crab angle and ground speed by direction
lines[, direction := ifelse(cos(course * pi / 180) > 0, "Northbound", "Southbound")]
fD1 <- ggplot(lines, aes(course %% 360)) +
  geom_histogram(binwidth = 10, boundary = 0, fill = "#2c7fb8", colour = "white") +
  coord_polar(start = 0) + scale_x_continuous(limits = c(0, 360), breaks = seq(0, 315, 45)) +
  labs(x = NULL, y = "Lines", title = "(a) Line course (°)")
fD2 <- ggplot(lines, aes(direction, crab, fill = direction)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(width = 0.5, alpha = 0.6) +
  geom_jitter(width = 0.08, size = 1.5) + scale_fill_manual(values = c("#1b9e77", "#d95f02")) +
  labs(x = NULL, y = "Crab angle (°)", title = "(b) Crab angle") + theme(legend.position = "none")
fD3 <- ggplot(lines, aes(direction, v_mean, fill = direction)) +
  geom_boxplot(width = 0.5, alpha = 0.6) + geom_jitter(width = 0.08, size = 1.5) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02")) +
  labs(x = NULL, y = "Ground speed (m/s)", title = "(c) Ground speed") +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "Fig_T4_line_geometry.png"), fD1 + fD2 + fD3,
       width = 9, height = 3.6, dpi = 300)

cat("\nOutputs written to:", out_dir, "\n")
