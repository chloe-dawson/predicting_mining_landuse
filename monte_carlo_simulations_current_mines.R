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

excel_path     <- "your-file-path/data/S&P_nickel_17022025.xlsx"
shapefile_path <- "your-file-path/data/current_mines/current_ni_mines.shp"

# CRS: World Behrmann throughout (equal-area, matches habitat / AOH rasters).
crs_behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# --- Excel data: Operating + Expansion mines ---
points_data <- read_excel(excel_path) %>%
  filter(DEV_STAGE %in% c("Expansion", "Operating")) %>%
  filter(DEPOSIT_TYPE %in% c("Laterite", "Magmatic Sulphide")) %>%
  select(PROP_ID, LONGITUDE, LATITUDE, CONTAINED_R_AND_R_PCT_TONNE, DEPOSIT_TYPE) %>%
  mutate(CONTAINED_R_AND_R_PCT_TONNE = as.numeric(CONTAINED_R_AND_R_PCT_TONNE))

pts_orig <- st_as_sf(points_data,
                     coords = c("LONGITUDE", "LATITUDE"),
                     crs = 4326) %>%
  st_transform(crs = crs_behrmann)

# --- Current mines shapefile (existing footprints) ---
cat("Loading current mines shapefile...\n")
current_mines_in <- st_read(shapefile_path, quiet = TRUE)
if (is.na(st_crs(current_mines_in))) {
  stop("Shapefile has no CRS defined (missing or broken .prj). ",
       "Set it explicitly with st_crs(current_mines_in) <- <EPSG> before transforming.")
}
current_mines <- st_transform(current_mines_in, crs = crs_behrmann)
current_mines$Area_km2 <- as.numeric(current_mines$Area_km2)

# Sanity check: PROP_ID should be unique. If not, keep first occurrence and warn.
dup_ids <- current_mines$PROP_ID[duplicated(current_mines$PROP_ID)]
if (length(dup_ids) > 0) {
  warning("Duplicate PROP_IDs found in shapefile: ",
          paste(unique(dup_ids), collapse = ", "),
          ". Only the first occurrence of each will be used.")
  current_mines <- current_mines[!duplicated(current_mines$PROP_ID), ]
}

# Lookup: PROP_ID (as character) -> row index in current_mines
shapefile_index <- setNames(seq_len(nrow(current_mines)),
                            as.character(current_mines$PROP_ID))

# ============================================================
# PART 1b: PER-MINE AUDIT (which method each mine will use)
# ============================================================

mine_audit <- data.frame(
  PROP_ID      = pts_orig$PROP_ID,
  deposit_type = pts_orig$DEPOSIT_TYPE,
  tonnage      = pts_orig$CONTAINED_R_AND_R_PCT_TONNE,
  longitude    = points_data$LONGITUDE,
  latitude     = points_data$LATITUDE,
  has_polygon  = as.character(pts_orig$PROP_ID) %in% names(shapefile_index),
  stringsAsFactors = FALSE
)
mine_audit$existing_km2 <- ifelse(
  mine_audit$has_polygon,
  current_mines$Area_km2[shapefile_index[as.character(mine_audit$PROP_ID)]],
  NA_real_
)
mine_audit$method <- ifelse(
  mine_audit$has_polygon,
  "Shapefile centroid + expansion subtraction",
  "Random location shift (no existing polygon)"
)

# Shapefile mines that don't appear in Excel (won't be simulated, flagged for audit)
shp_only <- setdiff(as.character(current_mines$PROP_ID),
                    as.character(pts_orig$PROP_ID))
shp_only_df <- if (length(shp_only) > 0) {
  data.frame(
    PROP_ID = shp_only,
    note    = "In shapefile but no matching Operating/Expansion record in Excel",
    stringsAsFactors = FALSE
  )
} else {
  data.frame(PROP_ID = character(0), note = character(0))
}

cat(sprintf("Excel mines (Operating + Expansion):              %d\n", nrow(pts_orig)))
cat(sprintf("Unique PROP_IDs in shapefile:                     %d\n", nrow(current_mines)))
cat(sprintf("  Excel mines WITH a shapefile polygon:           %d\n",
            sum(mine_audit$has_polygon)))
cat(sprintf("  Excel mines WITHOUT a polygon (fallback):       %d\n",
            sum(!mine_audit$has_polygon)))
cat(sprintf("  Shapefile mines NOT in Excel (no sim run):      %d\n",
            length(shp_only)))

# --- Land for clipping ---
land_union <- ne_countries(scale = "large", returnclass = "sf") %>%
  st_transform(crs = crs_behrmann) %>%
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
# PART 3: FLOOD-FILL FUNCTION (UNCHANGED)
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
# HELPERS
# ============================================================

safe_intersection <- function(x, y) {
  tryCatch(st_intersection(x, y), error = function(e) NULL)
}

safe_difference <- function(x, y) {
  tryCatch(st_difference(x, y), error = function(e) x)
}

safe_is_usable <- function(geom) {
  if (is.null(geom)) return(FALSE)
  if (length(geom) == 0) return(FALSE)
  tryCatch({ !all(st_is_empty(geom)) }, error = function(e) FALSE)
}

safe_area <- function(geom) {
  if (!safe_is_usable(geom)) return(0)
  tryCatch({ as.numeric(st_area(st_union(geom))) }, error = function(e) 0)
}

# Clean possible GEOMETRYCOLLECTION output from st_difference (keep polygons only)
clean_polys <- function(g) {
  if (!safe_is_usable(g)) return(g)
  tryCatch({
    types <- st_geometry_type(g)
    if (any(types == "GEOMETRYCOLLECTION")) {
      g <- st_collection_extract(g, "POLYGON")
    }
    g <- st_make_valid(g)
    g
  }, error = function(e) g)
}

# Secondary rescue for polygon mines whose centroid-based circle (even at
# target-area radius) sits entirely inside the existing polygon - typical
# for elongated / concave / coastal shapes. Buffers the EXISTING polygon
# outward by a distance estimated from the area-from-perimeter identity
estimate_buffer_distance <- function(geom, expansion_area_m2) {
  perim <- tryCatch({
    as.numeric(st_length(st_cast(st_geometry(geom), "MULTILINESTRING")))
  }, error = function(e) NA_real_)
  if (length(perim) == 0 || !is.finite(perim) || perim == 0) {
    # Square-root approximation fallback if perimeter can't be computed
    existing_area <- safe_area(geom)
    perim <- 4 * sqrt(max(existing_area, 1))
  }
  d <- (-perim + sqrt(perim^2 + 4 * pi * expansion_area_m2)) / (2 * pi)
  if (!is.finite(d) || d <= 0) d <- sqrt(expansion_area_m2 / pi)
  d
}

buffer_existing_outward <- function(existing_geom, target_area) {
  existing_area <- safe_area(existing_geom)
  if (target_area <= existing_area) return(NULL)
  d <- estimate_buffer_distance(existing_geom, target_area - existing_area)
  buf <- tryCatch(st_buffer(existing_geom, dist = d), error = function(e) NULL)
  if (!safe_is_usable(buf)) return(NULL)
  buf_land <- safe_intersection(buf, land_clipped)
  # If land-clip produces nothing usable, accept the unclipped buffer (may
  # extend into water) as last resort so we at least emit *some* expansion.
  if (!safe_is_usable(buf_land)) buf_land <- buf
  expansion <- safe_difference(buf_land, existing_geom)
  clean_polys(expansion)
}

# Build a footprint geometry from a seed point and target area.
# Returns list(geom, type) where type is "Circular" / "FloodFill" / "Fallback".
build_footprint <- function(seed_pt, target_area, caps) {
  buffer_radius <- sqrt(target_area / pi)
  raw_buffer    <- st_buffer(seed_pt, dist = buffer_radius)
  land_buffer   <- safe_intersection(raw_buffer, land_clipped)
  lb_area       <- if (safe_is_usable(land_buffer)) safe_area(land_buffer) else 0
  
  if (lb_area >= caps$min_area) {
    return(list(geom = land_buffer, type = "Circular"))
  }
  
  # Try flood-fill on a 100 m raster
  cell_size  <- 100
  need_cells <- ceiling(target_area / (cell_size^2))
  bbox_buf   <- st_buffer(seed_pt, dist = buffer_radius * 2)
  land_clip  <- safe_intersection(bbox_buf, land_clipped)
  
  if (safe_is_usable(land_clip)) {
    r <- tryCatch({
      rr <- rast(ext(st_bbox(land_clip)), res = cell_size, crs = crs_behrmann)
      rasterize(vect(land_clip), rr, field = 1, background = 0)
    }, error = function(e) NULL)
    
    filled <- if (!is.null(r)) fill_mine_area(seed_pt, need_cells, r) else integer(0)
    
    if (length(filled) > 0) {
      filled_rast <- r
      values(filled_rast) <- NA
      filled_rast[filled] <- 1
      filled_poly <- as.polygons(filled_rast, dissolve = TRUE)
      filled_sf   <- st_as_sf(filled_poly)
      st_crs(filled_sf) <- crs_behrmann
      return(list(geom = st_geometry(filled_sf), type = "FloodFill"))
    }
  }
  
  # Fallback: sized to target_area (not caps$min_area) so the circle has a
  # chance of extending beyond the existing polygon on polygon-mines. Still
  # clipped to land; raw buffer used only if land clip fails.
  fallback_radius <- sqrt(target_area / pi)
  fallback_buf    <- st_buffer(seed_pt, dist = fallback_radius)
  fallback_land   <- safe_intersection(fallback_buf, land_clipped)
  geom <- if (safe_is_usable(fallback_land)) fallback_land else fallback_buf
  list(geom = geom, type = "Fallback")
}

# ============================================================
# PART 4: SIMULATION FUNCTION
# ============================================================

simulate_one <- function() {
  n  <- nrow(pts_orig)
  all_geoms  <- vector("list", n)
  point_info <- data.frame(
    point_id      = 1:n,
    PROP_ID       = pts_orig$PROP_ID,
    deposit_type  = pts_orig$DEPOSIT_TYPE,
    tonnage       = pts_orig$CONTAINED_R_AND_R_PCT_TONNE,
    has_polygon   = FALSE,
    existing_km2  = NA_real_,
    ltf           = NA_real_,
    target_km2    = NA_real_,
    new_total_km2 = NA_real_,
    expansion_km2 = NA_real_,
    buffer_type   = NA_character_,
    zero_reason   = NA_character_,
    used_orig     = FALSE
  )
  max_shift <- 15000
  
  for (i in seq_len(n)) {
    orig_pt   <- st_geometry(pts_orig[i, ])
    dep_type  <- pts_orig$DEPOSIT_TYPE[i]
    tonnage   <- pts_orig$CONTAINED_R_AND_R_PCT_TONNE[i]
    prop_id   <- as.character(pts_orig$PROP_ID[i])
    params    <- ltf_params[[dep_type]]
    caps      <- mine_size_caps[[dep_type]]
    
    # --- A) Decide seed point: shapefile centroid (no shift) or random shift ---
    sf_idx        <- shapefile_index[prop_id]
    has_polygon   <- !is.na(sf_idx)
    existing_geom <- NULL
    existing_area <- 0  # in m^2
    
    if (has_polygon) {
      existing_geom <- st_geometry(current_mines)[sf_idx]
      existing_area <- current_mines$Area_km2[sf_idx] * 1e6   # km^2 -> m^2
      point_info$has_polygon[i]  <- TRUE
      point_info$existing_km2[i] <- existing_area / 1e6
      
      # Centroid of dissolved polygon as the fixed seed (no random shift)
      seed_pt <- st_centroid(existing_geom)
      st_crs(seed_pt) <- crs_behrmann
    } else {
      # ---- ORIGINAL RANDOM-SHIFT LOGIC (unchanged) ----
      shifted <- NULL
      for (attempt in 1:50) {
        dx <- runif(1, -max_shift, max_shift)
        dy <- runif(1, -max_shift, max_shift)
        candidate <- orig_pt + c(dx, dy)
        st_crs(candidate) <- crs_behrmann
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
      seed_pt <- shifted
    }
    
    # --- B) Random LTF (always varies across simulations) ---
    ltf_raw <- rlnorm(1, meanlog = log(params$median), sdlog = params$log_sd)
    ltf     <- max(params$min, min(params$max, ltf_raw))
    point_info$ltf[i] <- ltf
    
    # --- C) Target area, with min raised to existing footprint when applicable ---
    target_area_raw <- ltf * tonnage
    min_for_this_mine <- if (has_polygon) {
      max(existing_area, caps$min_area)
    } else {
      caps$min_area
    }
    # Bounded between min_for_this_mine and the deposit-type max cap.
    # If min_for_this_mine > caps$max_area (existing already exceeds the cap),
    # max() wins -> target = existing_area, i.e. no expansion this sim.
    target_area <- max(min_for_this_mine, min(caps$max_area, target_area_raw))
    point_info$target_km2[i] <- target_area / 1e6
    
    # --- D) Build the new footprint (full new mine, including existing area) ---
    fp        <- build_footprint(seed_pt, target_area, caps)
    new_geom  <- fp$geom
    point_info$buffer_type[i]   <- fp$type
    point_info$new_total_km2[i] <- safe_area(new_geom) / 1e6
    
    # --- E) For shapefile mines: subtract existing polygon -> expansion only ---
    if (has_polygon && safe_is_usable(new_geom)) {
      expansion <- safe_difference(new_geom, existing_geom)
      expansion <- clean_polys(expansion)
      
      # Secondary rescue: if the primary expansion is empty (or effectively
      # zero) but the model wanted growth (target > existing), buffer the
      # existing polygon outward. Catches pathological cases where the
      # centroid circle sits entirely inside an elongated / concave polygon.
      exp_area_m2 <- safe_area(expansion)
      if (target_area > existing_area && exp_area_m2 < 1000) {   # < 0.001 km^2
        alt <- buffer_existing_outward(existing_geom, target_area)
        if (safe_is_usable(alt) && safe_area(alt) > exp_area_m2) {
          expansion <- alt
          point_info$buffer_type[i] <- "OutwardBuffer"
        }
      }
      all_geoms[[i]] <- expansion
    } else {
      all_geoms[[i]] <- new_geom
    }
    
    point_info$expansion_km2[i] <- safe_area(all_geoms[[i]]) / 1e6
    
    # --- F) Classify zero / near-zero expansion for audit ---
    if (point_info$expansion_km2[i] < 0.001) {
      # Tolerance of 1 m^2 to account for floating-point equality
      if (target_area <= existing_area + 1) {
        point_info$zero_reason[i] <- "by_design"      # ltf draw or max-cap made expansion unnecessary
      } else if (has_polygon) {
        point_info$zero_reason[i] <- "geometry_failure"   # wanted expansion, code couldn't render
      }
      # (Non-polygon mines with zero expansion would be a different bug path;
      # leave zero_reason = NA and let buffer_type surface it.)
    }
  }
  
  # Build sf object - drop rows whose geometry is empty (no expansion this sim)
  sim_rows <- lapply(seq_len(n), function(i) {
    g <- all_geoms[[i]]
    if (!safe_is_usable(g)) return(NULL)
    st_sf(
      point_id      = point_info$point_id[i],
      PROP_ID       = point_info$PROP_ID[i],
      deposit_type  = point_info$deposit_type[i],
      tonnage       = point_info$tonnage[i],
      has_polygon   = point_info$has_polygon[i],
      existing_km2  = point_info$existing_km2[i],
      ltf           = point_info$ltf[i],
      target_km2    = point_info$target_km2[i],
      new_total_km2 = point_info$new_total_km2[i],
      expansion_km2 = point_info$expansion_km2[i],
      buffer_type   = point_info$buffer_type[i],
      used_orig     = point_info$used_orig[i],
      geometry      = st_geometry(g)
    )
  })
  sim_rows <- sim_rows[!vapply(sim_rows, is.null, logical(1))]
  sim_sf   <- do.call(rbind, sim_rows)
  st_crs(sim_sf) <- crs_behrmann
  
  attr(sim_sf, "point_info") <- point_info
  sim_sf
}

# ============================================================
# PART 5: RUN 125 SIMULATIONS, SAVE EACH AS GEOPACKAGE
# ============================================================

n_sims  <- 125
out_dir <- "your-file-path/outputs/simulations_current"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("Running", n_sims, "simulations and saving to:", out_dir, "\n")
start_time <- Sys.time()

summary_list <- vector("list", n_sims)

for (s in seq_len(n_sims)) {
  sim_sf <- simulate_one()
  pinfo  <- attr(sim_sf, "point_info")
  
  out_file <- file.path(out_dir, sprintf("sim_%03d.gpkg", s))
  st_write(sim_sf, out_file, delete_dsn = TRUE, quiet = TRUE)
  
  summary_list[[s]] <- data.frame(
    simulation             = s,
    n_with_polygon         = sum(pinfo$has_polygon),
    n_without_polygon      = sum(!pinfo$has_polygon),
    total_expansion_km2    = sum(pinfo$expansion_km2, na.rm = TRUE),
    mean_expansion_km2     = mean(pinfo$expansion_km2, na.rm = TRUE),
    total_new_km2          = sum(pinfo$new_total_km2, na.rm = TRUE),
    n_circular             = sum(pinfo$buffer_type == "Circular",      na.rm = TRUE),
    n_floodfill            = sum(pinfo$buffer_type == "FloodFill",     na.rm = TRUE),
    n_fallback             = sum(pinfo$buffer_type == "Fallback",      na.rm = TRUE),
    n_outward_buffer       = sum(pinfo$buffer_type == "OutwardBuffer", na.rm = TRUE),
    n_zero_by_design       = sum(pinfo$zero_reason == "by_design",        na.rm = TRUE),
    n_zero_geom_failure    = sum(pinfo$zero_reason == "geometry_failure", na.rm = TRUE),
    n_used_orig            = sum(pinfo$used_orig)
  )
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat(sprintf("  Sim %d/%d saved | %.1f min elapsed\n", s, n_sims, elapsed))
}

end_time <- Sys.time()
cat("\nTotal runtime:",
    round(difftime(end_time, start_time, units = "hours"), 2), "hours\n")

# ============================================================
# PART 6: SUMMARY
# ============================================================

summary_df <- do.call(rbind, summary_list)
n_pts      <- nrow(pts_orig)

overall_summary <- data.frame(
  metric = c(
    "Number of simulations",
    "Mean total EXPANSION footprint (km2)",
    "SD total EXPANSION footprint (km2)",
    "Min total EXPANSION footprint (km2)",
    "Max total EXPANSION footprint (km2)",
    "Mean per-mine expansion area (km2)",
    "Mean total NEW footprint, incl. existing (km2)",
    "Mines with shapefile polygon (count)",
    "Mines without polygon -> random location (count)",
    "Mean original-location per sim (no-polygon subset)",
    "Mean flood-fills per sim",
    "Mean fallbacks per sim",
    "Mean OutwardBuffer rescues per sim",
    "Mean zero-expansion BY DESIGN per sim",
    "Mean zero-expansion GEOMETRY FAILURE per sim"
  ),
  value = c(
    nrow(summary_df),
    round(mean(summary_df$total_expansion_km2), 2),
    round(sd(summary_df$total_expansion_km2),   2),
    round(min(summary_df$total_expansion_km2),  2),
    round(max(summary_df$total_expansion_km2),  2),
    round(mean(summary_df$mean_expansion_km2),  2),
    round(mean(summary_df$total_new_km2),       2),
    summary_df$n_with_polygon[1],
    summary_df$n_without_polygon[1],
    round(mean(summary_df$n_used_orig),          2),
    round(mean(summary_df$n_floodfill),          2),
    round(mean(summary_df$n_fallback),           2),
    round(mean(summary_df$n_outward_buffer),     2),
    round(mean(summary_df$n_zero_by_design),     2),
    round(mean(summary_df$n_zero_geom_failure),  2)
  )
)

write_xlsx(
  list(overall_summary      = overall_summary,
       mine_audit           = mine_audit,
       shapefile_only_mines = shp_only_df,
       per_simulation       = summary_df),
  path = file.path(out_dir, "simulation_summary.xlsx")
)

cat("\n=== SUMMARY ACROSS", n_sims, "SIMULATIONS ===\n")
cat(sprintf("Mean total EXPANSION footprint: %.2f km2 (SD: %.2f)\n",
            mean(summary_df$total_expansion_km2),
            sd(summary_df$total_expansion_km2)))
cat(sprintf("Range: %.2f - %.2f km2\n",
            min(summary_df$total_expansion_km2),
            max(summary_df$total_expansion_km2)))
cat(sprintf("Mines with existing polygon (fixed location):    %d / %d\n",
            summary_df$n_with_polygon[1], n_pts))
cat(sprintf("Mines without polygon (random-location fallback): %d / %d\n",
            summary_df$n_without_polygon[1], n_pts))
cat("\nMethod breakdown (full per-mine list saved to mine_audit sheet):\n")
cat(sprintf("  Shapefile centroid + expansion subtraction: %d mines\n",
            sum(mine_audit$has_polygon)))
cat(sprintf("  Random location shift (no polygon):         %d mines\n",
            sum(!mine_audit$has_polygon)))
cat(sprintf("\nPer-sim averages:\n"))
cat(sprintf("  OutwardBuffer rescues (2ndary strategy)   : %.2f\n",
            mean(summary_df$n_outward_buffer)))
cat(sprintf("  Zero expansion BY DESIGN (no growth drawn): %.2f\n",
            mean(summary_df$n_zero_by_design)))
cat(sprintf("  Zero expansion GEOMETRY FAILURE           : %.2f  <- investigate if > 0\n",
            mean(summary_df$n_zero_geom_failure)))
if (nrow(shp_only_df) > 0) {
  cat(sprintf("\n  Shapefile-only mines (in .shp but not Excel): %d (see shapefile_only_mines sheet)\n",
              nrow(shp_only_df)))
}
cat("\nAll simulations saved to:", out_dir, "\n")
cat("Each file is a GeoPackage (.gpkg) in World Behrmann projection.\n")
cat("File naming: sim_001.gpkg through sim_125.gpkg\n")
