library(sf)
library(dplyr)
library(readxl)
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

# Mollweide for the simulation (equal-area, accurate global distances)
crs_moll <- "+proj=moll +lon_0=0 +datum=WGS84 +units=m +no_defs"
# World Behrmann for final output (to match habitat maps)
crs_behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

pts_orig <- st_as_sf(points_data,
                     coords = c("LONGITUDE", "LATITUDE"),
                     crs = 4326) %>%
  st_transform(crs = crs_moll)

land_union <- ne_countries(scale = "large", returnclass = "sf") %>%
  st_transform(crs = crs_moll) %>%
  st_union()

cat("Pre-clipping land to mining regions...\n")
all_pts_buf <- st_buffer(st_union(st_geometry(pts_orig)), dist = 100000)
land_clipped <- st_intersection(land_union, all_pts_buf)
cat("Done.\n")

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
# PART 4: SIMULATION FUNCTION (RETURNS GEOMETRIES)
# ============================================================

simulate_one <- function() {
  n            <- nrow(pts_orig)
  all_geoms    <- vector("list", n)
  point_info   <- data.frame(
    point_id     = 1:n,
    deposit_type = pts_orig$DEPOSIT_TYPE,
    tonnage      = pts_orig$CONTAINED_R_AND_R_PCT_TONNE,
    ltf          = NA_real_,
    target_km2   = NA_real_,
    actual_km2   = NA_real_,
    buffer_type  = NA_character_,
    used_orig    = FALSE
  )
  max_shift <- 15000
  
  for (i in seq_len(n)) {
    orig_pt   <- st_geometry(pts_orig[i, ])
    dep_type  <- pts_orig$DEPOSIT_TYPE[i]
    tonnage   <- pts_orig$CONTAINED_R_AND_R_PCT_TONNE[i]
    params    <- ltf_params[[dep_type]]
    caps      <- mine_size_caps[[dep_type]]
    
    # --- A) SHIFT POINT ---
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
      point_info$used_orig[i] <- TRUE
      shifted <- orig_pt
    }
    
    # --- B) DRAW LTF ---
    ltf_raw <- rlnorm(1, meanlog = log(params$median), sdlog = params$log_sd)
    ltf     <- max(params$min, min(params$max, ltf_raw))
    point_info$ltf[i] <- ltf
    
    # --- C) TARGET AREA ---
    target_area <- ltf * tonnage
    target_area <- max(caps$min_area, min(caps$max_area, target_area))
    point_info$target_km2[i] <- target_area / 1e6
    
    # --- D) BUFFER CONSTRUCTION ---
    buffer_radius <- sqrt(target_area / pi)
    raw_buffer    <- st_buffer(shifted, dist = buffer_radius)
    land_buffer   <- safe_intersection(raw_buffer, land_clipped)
    
    lb_area <- if (safe_is_usable(land_buffer)) safe_area(land_buffer) else 0
    
    if (lb_area >= caps$min_area) {
      all_geoms[[i]] <- land_buffer
      point_info$buffer_type[i] <- "Circular"
    } else {
      cell_size   <- 100
      need_cells  <- ceiling(target_area / (cell_size^2))
      bbox_buf    <- st_buffer(shifted, dist = buffer_radius * 2)
      land_clip   <- safe_intersection(bbox_buf, land_clipped)
      
      if (safe_is_usable(land_clip)) {
        r <- tryCatch({
          rr <- rast(ext(st_bbox(land_clip)), res = cell_size, crs = crs_moll)
          rasterize(vect(land_clip), rr, field = 1, background = 0)
        }, error = function(e) NULL)
        
        filled <- if (!is.null(r)) fill_mine_area(shifted, need_cells, r) else integer(0)
        
        if (length(filled) > 0) {
          filled_rast <- r
          values(filled_rast) <- NA
          filled_rast[filled] <- 1
          filled_poly <- as.polygons(filled_rast, dissolve = TRUE)
          filled_sf   <- st_as_sf(filled_poly)
          st_crs(filled_sf) <- crs_moll
          all_geoms[[i]] <- st_geometry(filled_sf)
          point_info$buffer_type[i] <- "FloodFill"
        } else {
          fallback_radius <- sqrt(caps$min_area / pi)
          fallback_buf    <- st_buffer(shifted, dist = fallback_radius)
          fallback_land   <- safe_intersection(fallback_buf, land_clipped)
          all_geoms[[i]]  <- if (safe_is_usable(fallback_land)) fallback_land else fallback_buf
          point_info$buffer_type[i] <- "Fallback"
        }
      } else {
        fallback_radius <- sqrt(caps$min_area / pi)
        all_geoms[[i]]  <- st_buffer(shifted, dist = fallback_radius)
        point_info$buffer_type[i] <- "Fallback"
      }
    }
    
    point_info$actual_km2[i] <- safe_area(all_geoms[[i]]) / 1e6
  }
  
  # Build sf object for this simulation
  sim_sf <- do.call(rbind, lapply(seq_len(n), function(i) {
    st_sf(
      point_id     = point_info$point_id[i],
      deposit_type = point_info$deposit_type[i],
      tonnage      = point_info$tonnage[i],
      ltf          = point_info$ltf[i],
      target_km2   = point_info$target_km2[i],
      actual_km2   = point_info$actual_km2[i],
      buffer_type  = point_info$buffer_type[i],
      used_orig    = point_info$used_orig[i],
      geometry     = st_geometry(all_geoms[[i]])
    )
  }))
  
  st_crs(sim_sf) <- crs_moll
  sim_sf
}

# ============================================================
# PART 5: RUN 69 SIMULATIONS, SAVE EACH AS GEOPACKAGE
# ============================================================

n_sims  <- 69
out_dir <- "your-file-path/simulations"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("Running", n_sims, "simulations and saving to:", out_dir, "\n")
start_time <- Sys.time()

summary_list <- vector("list", n_sims)

for (s in seq_len(n_sims)) {
  sim_sf <- simulate_one()
  
  # Reproject to World Behrmann to match habitat maps
  sim_behrmann <- st_transform(sim_sf, crs = crs_behrmann)
  
  # Save as GeoPackage (preferred for complex sf objects)
  out_file <- file.path(out_dir, sprintf("sim_%03d.gpkg", s))
  st_write(sim_behrmann, out_file, delete_dsn = TRUE, quiet = TRUE)
  
  # Summary row
  summary_list[[s]] <- data.frame(
    simulation        = s,
    total_area_km2    = sum(sim_sf$actual_km2),
    mean_area_km2     = mean(sim_sf$actual_km2),
    n_circular        = sum(sim_sf$buffer_type == "Circular"),
    n_floodfill       = sum(sim_sf$buffer_type == "FloodFill"),
    n_fallback        = sum(sim_sf$buffer_type == "Fallback"),
    n_used_orig       = sum(sim_sf$used_orig)
  )
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat(sprintf("  Sim %d/%d saved | %.1f min elapsed\n", s, n_sims, elapsed))
}

end_time <- Sys.time()
cat("\nTotal runtime:",
    round(difftime(end_time, start_time, units = "hours"), 2), "hours\n")

# ============================================================
# PART 6: SAVE SUMMARY
# ============================================================

summary_df <- do.call(rbind, summary_list)
write.csv(summary_df,
          file.path(out_dir, "simulation_summary.csv"),
          row.names = FALSE)

cat("\n=== SUMMARY ACROSS 69 SIMULATIONS ===\n")
cat(sprintf("Mean total footprint: %.2f km2 (SD: %.2f)\n",
            mean(summary_df$total_area_km2), sd(summary_df$total_area_km2)))
cat(sprintf("Range: %.2f - %.2f km2\n",
            min(summary_df$total_area_km2), max(summary_df$total_area_km2)))

n_pts <- 373
cat("\n=== SPATIAL ISSUE TRACKING ===\n")
cat(sprintf("Points using original location (all 50 shifts in water):\n"))
cat(sprintf("  Mean per sim: %.2f / %d points (%.2f%%)\n",
            mean(summary_df$n_used_orig), n_pts,
            mean(summary_df$n_used_orig)/n_pts*100))
cat(sprintf("  Range: %d - %d per simulation\n",
            min(summary_df$n_used_orig), max(summary_df$n_used_orig)))
cat(sprintf("Buffers needing flood-fill (coastal/water overlap):\n"))
cat(sprintf("  Mean per sim: %.2f / %d points (%.2f%%)\n",
            mean(summary_df$n_floodfill), n_pts,
            mean(summary_df$n_floodfill)/n_pts*100))
cat(sprintf("  Range: %d - %d per simulation\n",
            min(summary_df$n_floodfill), max(summary_df$n_floodfill)))
cat(sprintf("Fallback buffers (flood-fill also failed):\n"))
cat(sprintf("  Mean per sim: %.2f / %d points (%.2f%%)\n",
            mean(summary_df$n_fallback), n_pts,
            mean(summary_df$n_fallback)/n_pts*100))
cat(sprintf("  Range: %d - %d per simulation\n",
            min(summary_df$n_fallback), max(summary_df$n_fallback)))

cat("\nAll simulations saved to:", out_dir, "\n")
cat("Each file is a GeoPackage (.gpkg) in World Behrmann projection.\n")
cat("File naming: sim_001.gpkg through sim_069.gpkg\n")
