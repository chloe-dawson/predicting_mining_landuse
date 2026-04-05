library(sf)
library(dplyr)
library(readxl)
library(writexl)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(terra, warn.conflicts = FALSE)

# ============================================================
# PART 1: DATA LOADING & SETUP
# ============================================================

excel_path <- "your-file-path/S&P_nickel_17022025.xlsx"
points_data <- read_excel(excel_path) %>%
  filter(STAGE_AGG == "Exploration") %>%
  filter(DEPOSIT_TYPE %in% c("Laterite", "Magmatic Sulphide")) %>%
  select(LONGITUDE, LATITUDE, CONTAINED_R_AND_R_PCT_TONNE, DEPOSIT_TYPE) %>%
  mutate(CONTAINED_R_AND_R_PCT_TONNE = as.numeric(CONTAINED_R_AND_R_PCT_TONNE))

crs_moll <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"
pts_orig <- st_as_sf(points_data,
                     coords = c("LONGITUDE", "LATITUDE"),
                     crs = 4326) %>%
  st_transform(crs = crs_moll)

land_union <- ne_countries(scale = "large", returnclass = "sf") %>%
  st_transform(crs = crs_moll) %>%
  st_union()

# Pre-clip land to areas near mining points for faster lookups
cat("Pre-clipping land to mining regions...\n")
all_pts_buf <- st_buffer(st_union(st_geometry(pts_orig)), dist = 100000)  # 100km around all points
land_clipped <- st_intersection(land_union, all_pts_buf)
cat("Done. Land clipped from global to mining regions only.\n")

# ============================================================
# PART 2: DEPOSIT-TYPE PARAMETERS
# ============================================================

ltf_params <- list(
  "Magmatic Sulphide" = list(median = 30, min = 4, max = 398, log_sd = 1.39),
  "Laterite"          = list(median = 50, min = 7, max = 229, log_sd = 0.99)
)

mine_size_caps <- list(
  "Magmatic Sulphide" = list(min_area = 0.40e6, max_area = 79.45e6),
  "Laterite"          = list(min_area = 2.73e6, max_area = 70.96e6)
)

# ============================================================
# PART 3: FLOOD-FILL FUNCTION
# ============================================================

fill_mine_area <- function(shifted_point, need_cells, land_rast) {
  seed <- cellFromXY(land_rast, sf::st_coordinates(shifted_point))
  vals <- values(land_rast)
  ncells_total <- ncell(land_rast)
  if (is.na(seed) || seed < 1 || seed > ncells_total || vals[seed] == 0) return(integer(0))

  region <- integer(0)
  queue  <- seed

  while (length(region) < need_cells && length(queue) > 0) {
    cell <- queue[1]; queue <- queue[-1]
    if (cell %in% region) next
    if (is.na(cell) || cell < 1 || cell > ncells_total) next
    xy_check <- xyFromCell(land_rast, cell)
    if (any(is.na(xy_check))) next
    region <- c(region, cell)
    nbs   <- terra::adjacent(land_rast, cells = cell,
                             directions = 8, pairs = FALSE)
    nbs   <- nbs[!is.na(nbs)]
    keeps <- nbs[nbs >= 1 & nbs <= ncells_total & !is.na(vals[nbs]) & vals[nbs] == 1]
    newq  <- setdiff(keeps, c(region, queue))
    queue <- c(queue, newq)
  }
  head(region, need_cells)
}

# ============================================================
# HELPER: safe geometry checks
# ============================================================

safe_intersection <- function(x, y) {
  tryCatch(st_intersection(x, y), error = function(e) NULL)
}

safe_is_usable <- function(geom) {
  if (is.null(geom)) return(FALSE)
  if (length(geom) == 0) return(FALSE)
  tryCatch({ !all(st_is_empty(geom)) }, error = function(e) FALSE)
}

safe_area <- function(geom) {
  tryCatch({ as.numeric(st_area(st_union(geom))) }, error = function(e) 0)
}

# ============================================================
# PART 4: SIMULATION FUNCTION (ONE ITERATION)
# ============================================================

simulate_union_area <- function() {
  n            <- nrow(pts_orig)
  all_geoms    <- vector("list", n)
  flagged_min  <- 0
  snap_count   <- 0
  flood_count  <- 0
  max_shift    <- 15000

  for (i in seq_len(n)) {
    orig_pt   <- st_geometry(pts_orig[i, ])
    dep_type  <- pts_orig$DEPOSIT_TYPE[i]
    tonnage   <- pts_orig$CONTAINED_R_AND_R_PCT_TONNE[i]
    params    <- ltf_params[[dep_type]]
    caps      <- mine_size_caps[[dep_type]]

    # --- A) SHIFT POINT WITHIN 15 KM ---
    shifted <- NULL
    for (attempt in 1:50) {
      dx <- runif(1, -max_shift, max_shift)
      dy <- runif(1, -max_shift, max_shift)
      candidate <- orig_pt + c(dx, dy)
      st_crs(candidate) <- crs_moll
      hit <- tryCatch(
        as.numeric(st_intersects(candidate, land_clipped, sparse = FALSE)),
        error = function(e) 0
      )
      if (!is.na(hit) && hit > 0) {
        shifted <- candidate
        break
      }
    }
    if (is.null(shifted)) {
      snap_count <- snap_count + 1
      shifted <- orig_pt
    }

    # --- B) DRAW LAND TRANSFORMATION FACTOR ---
    ltf_log_mean <- log(params$median)
    ltf_log_sd   <- params$log_sd
    ltf_raw      <- rlnorm(1, meanlog = ltf_log_mean, sdlog = ltf_log_sd)
    ltf          <- max(params$min, min(params$max, ltf_raw))

    # --- C) CALCULATE TARGET AREA ---
    target_area <- ltf * tonnage
    target_area <- max(caps$min_area, min(caps$max_area, target_area))

    # --- D) CREATE LAND-ONLY BUFFER ---
    buffer_radius <- sqrt(target_area / pi)
    raw_buffer    <- st_buffer(shifted, dist = buffer_radius)
    land_buffer   <- safe_intersection(raw_buffer, land_clipped)

    lb_area <- 0
    if (safe_is_usable(land_buffer)) {
      lb_area <- safe_area(land_buffer)
    }

    if (lb_area >= caps$min_area) {
      all_geoms[[i]] <- land_buffer

    } else {
      # FLOOD-FILL: grow square raster cells outward on land
      flood_count <- flood_count + 1
      cell_size   <- 100
      need_cells  <- ceiling(target_area / (cell_size^2))
      bbox_expand <- buffer_radius * 2
      bbox_buf    <- st_buffer(shifted, dist = bbox_expand)
      land_clip   <- safe_intersection(bbox_buf, land_clipped)

      if (safe_is_usable(land_clip)) {
        r <- tryCatch({
          rr <- rast(ext(st_bbox(land_clip)), res = cell_size, crs = crs_moll)
          rasterize(vect(land_clip), rr, field = 1, background = 0)
        }, error = function(e) NULL)

        filled <- if (!is.null(r)) fill_mine_area(shifted, need_cells, r) else integer(0)

        if (length(filled) > 0) {
          # Convert raster cells to square polygons (NOT circular buffers)
          filled_rast <- r
          values(filled_rast) <- NA
          filled_rast[filled] <- 1
          filled_poly <- as.polygons(filled_rast, dissolve = TRUE)
          filled_sf   <- st_as_sf(filled_poly)
          st_crs(filled_sf) <- crs_moll
          all_geoms[[i]] <- st_geometry(filled_sf)
        } else {
          flagged_min <- flagged_min + 1
          fallback_radius <- sqrt(caps$min_area / pi)
          fallback_buf    <- st_buffer(shifted, dist = fallback_radius)
          fallback_land   <- safe_intersection(fallback_buf, land_clipped)
          if (safe_is_usable(fallback_land)) {
            all_geoms[[i]] <- fallback_land
          } else {
            all_geoms[[i]] <- fallback_buf
          }
        }
      } else {
        flagged_min <- flagged_min + 1
        fallback_radius <- sqrt(caps$min_area / pi)
        all_geoms[[i]] <- st_buffer(shifted, dist = fallback_radius)
      }
    }
  }

  valid    <- all_geoms[!sapply(all_geoms, is.null)]
  combined <- do.call(c, lapply(valid, st_geometry))
  st_crs(combined) <- crs_moll
  ua <- as.numeric(st_area(st_union(combined)))

  list(union_area = ua, fallback_count = flagged_min,
       snap_count = snap_count, flood_count = flood_count)
}

# ============================================================
# PART 5: RUN PILOT SIMULATIONS
# ============================================================

pilot_n <- 1000
cat("=== CONFIDENCE INTERVAL HALF-WIDTH ANALYSIS ===\n")
cat("Running", pilot_n, "pilot simulations...\n")
cat("Number of exploration points:", nrow(pts_orig), "\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

start_time <- Sys.time()

results <- vector("list", pilot_n)
for (s in seq_len(pilot_n)) {
  results[[s]] <- simulate_union_area()
  if (s %% 50 == 0) {
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    cat(sprintf("  Simulation %d/%d complete (%.1f min elapsed)\n",
                s, pilot_n, as.numeric(elapsed)))
  }
}

end_time <- Sys.time()
cat("Total runtime:", round(difftime(end_time, start_time, units = "hours"), 2), "hours\n\n")

areas      <- sapply(results, `[[`, "union_area")
fallbacks  <- sapply(results, `[[`, "fallback_count")
snaps      <- sapply(results, `[[`, "snap_count")
floods     <- sapply(results, `[[`, "flood_count")
areas_km2  <- areas / 1e6

# ============================================================
# PART 6: CONFIDENCE INTERVAL HALF-WIDTH CALCULATION
# ============================================================

meanA <- mean(areas)
sdA   <- sd(areas)
cvA   <- (sdA / meanA) * 100

cat("=== PILOT SIMULATION RESULTS ===\n")
cat(sprintf("Mean total footprint: %.2f km2\n", meanA / 1e6))
cat(sprintf("Standard deviation:   %.2f km2\n", sdA / 1e6))
cat(sprintf("Coefficient of variation: %.2f%%\n", cvA))
cat(sprintf("Range: %.2f - %.2f km2\n", min(areas_km2), max(areas_km2)))

cat("\n=== SPATIAL ISSUE TRACKING ===\n")
cat(sprintf("Points using original location (all 50 shifts in water):\n"))
cat(sprintf("  Mean per sim: %.2f / %d points (%.2f%%)\n",
            mean(snaps), nrow(pts_orig), mean(snaps)/nrow(pts_orig)*100))
cat(sprintf("  Range: %d - %d per simulation\n", min(snaps), max(snaps)))
cat(sprintf("Buffers needing flood-fill (coastal/water overlap):\n"))
cat(sprintf("  Mean per sim: %.2f / %d points (%.2f%%)\n",
            mean(floods), nrow(pts_orig), mean(floods)/nrow(pts_orig)*100))
cat(sprintf("  Range: %d - %d per simulation\n", min(floods), max(floods)))
cat(sprintf("Fallback buffers (flood-fill also failed):\n"))
cat(sprintf("  Mean per sim: %.2f / %d points (%.2f%%)\n",
            mean(fallbacks), nrow(pts_orig), mean(fallbacks)/nrow(pts_orig)*100))
cat(sprintf("  Range: %d - %d per simulation\n", min(fallbacks), max(fallbacks)))

z <- 1.96
precision_levels <- c(0.01, 0.02, 0.05, 0.10)

ci_calcs <- data.frame(
  precision_percent = precision_levels * 100,
  margin_of_error_km2 = (meanA * precision_levels) / 1e6,
  required_n = NA_integer_
)

for (j in seq_along(precision_levels)) {
  E <- meanA * precision_levels[j]
  n_required <- ceiling((z * sdA / E)^2)
  ci_calcs$required_n[j] <- n_required
}

cat("\n=== REQUIRED SIMULATIONS (95% CI Half-Width) ===\n")
for (j in seq_len(nrow(ci_calcs))) {
  cat(sprintf("  +/-%d%% of mean (+/-%.2f km2): %d simulations\n",
              ci_calcs$precision_percent[j],
              ci_calcs$margin_of_error_km2[j],
              ci_calcs$required_n[j]))
}
cat("\nRecommendation: Use the 5% precision level for your final Monte Carlo.\n")
cat("This means the 95% CI around your mean estimate will be within\n")
cat(sprintf("+/-%.2f km2 of the true mean.\n",
            ci_calcs$margin_of_error_km2[ci_calcs$precision_percent == 5]))

# ============================================================
# PART 7: SAVE RESULTS
# ============================================================

summary_df <- data.frame(
  metric = c("Pilot simulations", "Mean area (km2)", "SD (km2)",
             "CV (%)", "Min area (km2)", "Max area (km2)",
             "Avg original-location per sim", "Original-location rate (%)",
             "Avg flood-fills per sim", "Flood-fill rate (%)",
             "Avg fallbacks per sim", "Fallback rate (%)",
             "Confidence level", "Recommended precision",
             "Recommended n (5%)"),
  value  = c(pilot_n, round(meanA/1e6, 2), round(sdA/1e6, 2),
             round(cvA, 2), round(min(areas_km2), 2), round(max(areas_km2), 2),
             round(mean(snaps), 2),
             round(mean(snaps)/nrow(pts_orig)*100, 2),
             round(mean(floods), 2),
             round(mean(floods)/nrow(pts_orig)*100, 2),
             round(mean(fallbacks), 2),
             round(mean(fallbacks)/nrow(pts_orig)*100, 2),
             "95%", "5%",
             ci_calcs$required_n[ci_calcs$precision_percent == 5])
)

sim_results_df <- data.frame(
  simulation  = 1:pilot_n,
  area_m2     = areas,
  area_km2    = areas_km2,
  snaps       = snaps,
  flood_fills = floods,
  fallbacks   = fallbacks
)

write_xlsx(
  list(summary = summary_df,
       ci_calculations = ci_calcs,
       simulation_results = sim_results_df),
  path = "your-file-path/ci_halfwidth_results.xlsx"
)

cat("\nResults saved to ci_halfwidth_results.xlsx\n")
