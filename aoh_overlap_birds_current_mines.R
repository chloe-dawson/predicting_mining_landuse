# =============================================================================
# Per-mine mining x AOH overlap  +  Per-species AOH change — BIRDS
#                                                      (current / expansion)
#
# Merged pipeline. For each sim_XXX.gpkg, computes overlap area (km^2)
# between EACH individual mine and every selected species' AOH, then in the
# same pass records the original AOH area per species while the mosaicked
# raster is still loaded. At the end, aggregates per (species, simulation)
# and writes two CSVs whose contents are identical to running the original
# per-mine overlap script followed by the per-species aoh-change script.
#
# Birds have multiple tifs per species (e.g. [id]_1.tif, [id]_2.tif, ...
# for seasonal / range variants) spread across six BOTW subfolders. All
# tifs for a species are mosaicked into one layer (max where they overlap)
# before any per-mine overlap or original-AOH computation, so a cell shared
# between seasonal/range tifs is counted once per species.
#
# Mines are NOT unioned, so each mine's standalone habitat impact is
# recorded — required for scenario toggling. If two mines physically overlap
# and both touch a habitat cell, both get credit (intentional double-counting
# at the per-mine level).
# =============================================================================

library(sf)
library(terra)
library(readxl)
library(readr)
library(dplyr)
library(tools)

# Send terra temp files to H: drive instead of C:
terra_tmp <- "H:/scenario_analysis/terra_tmp"
dir.create(terra_tmp, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = terra_tmp)

# ---- 1. File paths -----------------------------------------------------------

sim_dir         <- "H:/scenario_analysis/outputs/simulations_current"
species_path    <- "H:/scenario_analysis/data/species_selection/birds_species_selection.xlsx"
output_path     <- "H:/scenario_analysis/outputs/birds_per_mine_overlap_current.csv"
aoh_change_path <- "H:/scenario_analysis/outputs/birds_aoh_change_current.csv"

aoh_dirs <- c(
  "H:/scenario_analysis/data/aoh/birds/BOTW_2023_1-part-1",
  "H:/scenario_analysis/data/aoh/birds/BOTW_2023_1-part-2",
  "H:/scenario_analysis/data/aoh/birds/BOTW_2023_1-part-3",
  "H:/scenario_analysis/data/aoh/birds/BOTW_2023_1-part-4",
  "H:/scenario_analysis/data/aoh/birds/BOTW_2023_1-part-5",
  "H:/scenario_analysis/data/aoh/birds/BOTW_2023_1-part-6"
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

# ---- 2. Load selected species list ------------------------------------------

cat("Loading selected species list...\n")
species_df <- read_excel(species_path)
species_df$species_id <- as.character(species_df$species_id)

# ---- 3. Build AOH lookup: species_id -> list of raster paths ----------------

cat("Scanning AOH raster folders...\n")

all_tifs <- unlist(lapply(aoh_dirs, function(d) {
  if (!dir.exists(d)) return(character(0))
  list.files(d, pattern = "\\.tif$", full.names = TRUE,
             recursive = FALSE, ignore.case = TRUE)
}))

cat(sprintf("  Total raster files found: %d\n", length(all_tifs)))

basenames   <- file_path_sans_ext(basename(all_tifs))
species_ids <- sub("_[0-9]+$", "", basenames)

tif_df <- data.frame(species_id = species_ids, raster_path = all_tifs,
                     stringsAsFactors = FALSE)

aoh_lookup <- tif_df %>%
  group_by(species_id) %>%
  summarise(raster_paths = list(raster_path), n_files = n(), .groups = "drop")

species_to_check <- species_df %>%
  inner_join(aoh_lookup, by = "species_id")

cat(sprintf("  Selected species          : %d\n", nrow(species_df)))
cat(sprintf("  With AOH raster available : %d\n", nrow(species_to_check)))

# ---- 4. Load simulations as INDIVIDUAL mines (not unioned) ------------------
# Each sim is loaded as the full sf object so we keep per-mine attributes
# (point_id, PROP_ID) for the per-mine attribution.

sim_files <- sort(list.files(sim_dir, pattern = "^sim_\\d+\\.gpkg$", full.names = TRUE))
if (length(sim_files) == 0) stop("No sim_XXX.gpkg files found in ", sim_dir)
cat(sprintf("Found %d simulation files.\n", length(sim_files)))

cat("Loading simulation geometries (per-mine, not unioned)...\n")
sims <- lapply(sim_files, function(f) {
  s <- st_read(f, quiet = TRUE)
  list(
    id = as.integer(sub("^sim_(\\d+)\\.gpkg$", "\\1", basename(f))),
    sf = s
  )
})

# ---- 5. Per-species, per-mine overlap  +  original AOH capture --------------
# For each species:
#   1. Mosaic all tifs into one layer (fun = "max") so overlapping seasonal
#      cells are not double-counted within the species.
#   2. Compute ORIGINAL AOH from the mosaicked raster (one global() call).
#   3. For each sim, intersect the mine geometry with the species' bbox to
#      get only the candidate mines.
#   4. For each candidate mine, run crop/mask/global on its OWN geometry
#      and record one row (sim, species, point_id, PROP_ID, overlap_km2).
# Species with no overlap against any mine in a sim contribute one zero
# row per sim (point_id and PROP_ID set to NA), so every species has at
# least one row per simulation for downstream summary stats.

cat("Computing per-mine overlap per species x simulation...\n")

results <- vector("list", nrow(species_to_check))
n_sp <- nrow(species_to_check)

# Per-species total AOH area, captured while the mosaicked raster is in memory.
original_aoh_km2 <- setNames(rep(NA_real_, n_sp), species_to_check$species_id)

for (k in seq_len(n_sp)) {
  sid    <- species_to_check$species_id[k]
  rpaths <- species_to_check$raster_paths[[k]]
  
  if (k %% 50 == 0 || k == 1) {
    cat(sprintf("  Species %d / %d  (species_id: %s, %d tifs)\n",
                k, n_sp, sid, length(rpaths)))
  }
  
  sp_rows <- list()
  
  # Outer try only protects raster loading + mosaic. If any tif is
  # missing/corrupt, the species gets no rows at all (and a warning prints).
  raster_loaded <- tryCatch({
    # Mosaic all tifs for this species into one layer (max where they overlap).
    rasters <- lapply(rpaths, rast)
    r <- if (length(rasters) == 1) rasters[[1]] else
      do.call(terra::mosaic, c(rasters, list(fun = "max")))
    cell_km2 <- prod(res(r)) / 1e6
    r_bbox  <- st_as_sfc(st_bbox(r))
    
    # Original AOH area — one global() pass on the MOSAICKED raster, so
    # cells shared between seasonal tifs are counted once. 
    n_habitat_total <- global(r > 0, "sum", na.rm = TRUE)[1, 1]
    original_aoh_km2[sid] <- if (is.na(n_habitat_total)) 0 else n_habitat_total * cell_km2
    
    TRUE
  }, error = function(e) {
    cat(sprintf("  WARNING: could not load/mosaic rasters for species %s -- %s\n",
                sid, conditionMessage(e)))
    FALSE
  })
  
  if (raster_loaded) {
    for (sm in sims) {
      sim_rows <- tryCatch({
        mines_proj <- st_transform(sm$sf, crs(r))
        
        # Spatial filter: only mines whose geometry intersects the raster bbox
        hits <- suppressWarnings(st_intersects(mines_proj, r_bbox, sparse = FALSE)[, 1])
        
        if (!any(hits)) {
          # No candidate mines this sim -- record one zero row so the species has a complete per-sim record for downstream stats
          list(data.frame(
            simulation  = sm$id,
            species_id  = sid,
            point_id    = NA_integer_,
            PROP_ID     = NA_integer_,
            overlap_km2 = 0,
            n_tifs      = length(rpaths),
            failed      = FALSE,
            stringsAsFactors = FALSE
          ))
        } else {
          candidate_mines <- mines_proj[hits, ]
          out_rows <- list()
          sim_had_overlap <- FALSE
          
          for (m in seq_len(nrow(candidate_mines))) {
            mine    <- candidate_mines[m, ]
            pid     <- mine$point_id
            prop_id <- mine$PROP_ID
            
            ov_km2 <- tryCatch({
              v  <- vect(mine)
              rc <- crop(r, v, snap = "out")
              rm <- mask(rc, v)
              n_habitat <- global(rm > 0, "sum", na.rm = TRUE)[1, 1]
              if (is.na(n_habitat)) 0 else n_habitat * cell_km2
            }, error = function(e) 0)
            
            if (ov_km2 > 0) {
              sim_had_overlap <- TRUE
              out_rows[[length(out_rows) + 1]] <- data.frame(
                simulation  = sm$id,
                species_id  = sid,
                point_id    = pid,
                PROP_ID     = prop_id,
                overlap_km2 = ov_km2,
                n_tifs      = length(rpaths),
                failed      = FALSE,
                stringsAsFactors = FALSE
              )
            }
          }
          
          # Bbox intersected mines but none actually contained habitat cells -- still record a zero row
          if (!sim_had_overlap) {
            out_rows[[length(out_rows) + 1]] <- data.frame(
              simulation  = sm$id,
              species_id  = sid,
              point_id    = NA_integer_,
              PROP_ID     = NA_integer_,
              overlap_km2 = 0,
              n_tifs      = length(rpaths),
              failed      = FALSE,
              stringsAsFactors = FALSE
            )
          }
          out_rows
        }
      }, error = function(e) {
        # Sim-level failure -- record one flagged zero row so the species still has an entry for this sim, and continue with the next sim
        cat(sprintf("  WARNING: species %s sim %d failed -- %s\n",
                    sid, sm$id, conditionMessage(e)))
        list(data.frame(
          simulation  = sm$id,
          species_id  = sid,
          point_id    = NA_integer_,
          PROP_ID     = NA_integer_,
          overlap_km2 = 0,
          n_tifs      = length(rpaths),
          failed      = TRUE,
          stringsAsFactors = FALSE
        ))
      })
      
      sp_rows <- c(sp_rows, sim_rows)
    }
  }
  
  if (length(sp_rows) > 0) {
    results[[k]] <- bind_rows(sp_rows)
  }
}

# ---- 6. Write per-mine overlap CSV (Script 1 output) -------------------------

out <- bind_rows(results)
write.csv(out, output_path, row.names = FALSE)
cat(sprintf("\nDone -- wrote %d rows to %s\n", nrow(out), output_path))
cat(sprintf("Unique species with non-zero overlap: %d\n",
            length(unique(out$species_id[out$overlap_km2 > 0]))))
cat(sprintf("Unique mines responsible for habitat loss: %d\n",
            length(unique(out$PROP_ID[!is.na(out$PROP_ID)]))))

# ---- 7. Aggregate into AOH change CSV (Script 2 output) ----------------------

cat("\nAggregating overlap per species x simulation...\n")

species_with_overlap <- out %>%
  filter(overlap_km2 > 0) %>%
  pull(species_id) %>%
  unique()

agg <- out %>%
  filter(species_id %in% species_with_overlap) %>%
  group_by(simulation, species_id) %>%
  summarise(total_overlap_km2 = sum(overlap_km2, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(original_aoh_km2 = original_aoh_km2[species_id],
         final_aoh_km2    = original_aoh_km2 - total_overlap_km2)

write_csv(agg, aoh_change_path)
cat(sprintf("Saved %d rows to: %s\n", nrow(agg), aoh_change_path))
cat(sprintf("Species with no AOH raster found: %d\n",
            sum(is.na(original_aoh_km2[species_with_overlap]))))