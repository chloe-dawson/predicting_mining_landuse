# =============================================================================
# Per-mine mining x AOH overlap  +  Per-species AOH change — REPTILES
#                                                      (current / expansion)
#
# Merged pipeline. For each sim_XXX.gpkg, computes overlap area (km^2)
# between EACH individual mine and every selected species' AOH, then in the
# same pass records the original AOH area per species while the raster is
# still loaded. At the end, aggregates per (species, simulation) and writes
# two CSVs whose contents are identical to running the original per-mine
# overlap script followed by the per-species aoh-change script.
#
# Mines are NOT unioned, so each mine's standalone habitat impact is
# recorded — required for scenario toggling. If two mines physically overlap
# and both touch a habitat cell, both get credit (intentional double-counting
# at the per-mine level).
#
# Updated rasters (modified_[id].tif) take priority, falling back to
# originals ([id]_1.tif) across the reptiles folder.
# =============================================================================

library(sf)
library(terra)
library(readxl)
library(readr)
library(dplyr)
library(tools)

# Send terra temp files to H: drive instead of C:
terra_tmp <- "your-file-path/terra_tmp"
dir.create(terra_tmp, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = terra_tmp)

# ---- 1. File paths -----------------------------------------------------------

sim_dir       <- "your-file-path/outputs/simulations_current"
species_path  <- "your-file-path/data/species_selection/reptiles_species_selection.xlsx"
output_path   <- "your-file-path/outputs/reptiles_per_mine_overlap_current.csv"
aoh_change_path <- "your-file-path/outputs/reptiles_aoh_change_current.csv"

updated_dir   <- "your-file-path/data/aoh/updated_reptiles"
original_dirs <- "your-file-path/data/aoh/reptiles"

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

# ---- 2. Load selected species list ------------------------------------------

cat("Loading selected species list...\n")
species_df <- read_excel(species_path)
species_df$species_id <- as.character(species_df$species_id)

# ---- 3. Build AOH lookup: updated priority ----------------------------------

build_lookup <- function(folder, pattern_type = c("updated", "original")) {
  pattern_type <- match.arg(pattern_type)
  if (!dir.exists(folder)) return(data.frame(species_id = character(),
                                             raster_path = character(),
                                             stringsAsFactors = FALSE))
  tifs <- list.files(folder, pattern = "\\.tif$", full.names = TRUE,
                     recursive = FALSE, ignore.case = TRUE)
  basenames <- file_path_sans_ext(basename(tifs))
  ids <- if (pattern_type == "updated") sub("^modified_", "", basenames) else sub("_1$", "", basenames)
  data.frame(species_id = ids, raster_path = tifs, stringsAsFactors = FALSE)
}

cat("Scanning AOH raster folders...\n")
updated_lkp  <- build_lookup(updated_dir, "updated")
original_lkp <- bind_rows(lapply(original_dirs, build_lookup, pattern_type = "original"))

aoh_lookup <- bind_rows(
  updated_lkp  %>% mutate(source = "updated"),
  original_lkp %>% mutate(source = "original")
) %>%
  distinct(species_id, .keep_all = TRUE)

species_to_check <- species_df %>%
  inner_join(aoh_lookup, by = "species_id")

cat(sprintf("  Selected species          : %d\n", nrow(species_df)))
cat(sprintf("  With AOH raster available : %d\n", nrow(species_to_check)))
cat(sprintf("    updated : %d\n", sum(species_to_check$source == "updated")))
cat(sprintf("    original: %d\n", sum(species_to_check$source == "original")))

# ---- 4. Load simulations as INDIVIDUAL mines (not unioned) ------------------

sim_files <- sort(list.files(sim_dir, pattern = "^sim_\\d+\\.gpkg$", full.names = TRUE))
if (length(sim_files) == 0) stop("No sim_XXX.gpkg files found in ", sim_dir)
cat(sprintf("Found %d simulation files.\n", length(sim_files)))

cat("Loading simulation geometries (per-mine, not unioned)...\n")
sims <- lapply(sim_files, function(f) {
  s <- st_read(f, quiet = TRUE)
  list(
    id = as.integer(sub("^sim_(\\d+)\\.gpkg$", "\\1", basename(f))),
    sf = s   # keep full sf with point_id + PROP_ID
  )
})

# ---- 5. Per-species, per-mine overlap  +  original AOH capture --------------
# For each species:
#   1. Load raster, compute cell area, get raster bbox.
#   2. Compute ORIGINAL AOH in the same raster load (one global() call).
#   3. For each sim, intersect mine geometry with species' raster bbox, then
#      for each candidate mine run crop/mask/global on its OWN geometry and
#      record one row (sim, species, point_id, PROP_ID, overlap_km2).
# Species with no overlap against any mine in a sim contribute one zero
# row per sim (point_id and PROP_ID set to NA).

cat("Computing per-mine overlap per species x simulation...\n")

results <- vector("list", nrow(species_to_check))
n_sp <- nrow(species_to_check)

# Per-species total AOH area, captured while the raster is already in memory.
original_aoh_km2 <- setNames(rep(NA_real_, n_sp), species_to_check$species_id)

for (k in seq_len(n_sp)) {
  sid   <- species_to_check$species_id[k]
  rpath <- species_to_check$raster_path[k]
  src   <- species_to_check$source[k]
  
  if (k %% 50 == 0 || k == 1) {
    cat(sprintf("  Species %d / %d  (species_id: %s)\n", k, n_sp, sid))
  }
  
  sp_rows <- list()
  
  # Outer try only protects raster loading itself. If the raster file is
  # missing/corrupt, the species gets no rows at all (and a warning prints).
  raster_loaded <- tryCatch({
    r <- rast(rpath)
    cell_km2 <- prod(res(r)) / 1e6
    r_bbox  <- st_as_sfc(st_bbox(r))   # in raster CRS
    
    # Original AOH area — one global() pass while the raster is loaded.
    n_habitat_total <- global(r > 0, "sum", na.rm = TRUE)[1, 1]
    original_aoh_km2[sid] <- if (is.na(n_habitat_total)) 0 else n_habitat_total * cell_km2
    
    TRUE
  }, error = function(e) {
    cat(sprintf("  WARNING: could not load raster for species %s -- %s\n",
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
            source      = src,
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
                source      = src,
                failed      = FALSE,
                stringsAsFactors = FALSE
              )
            }
          }
          
          # Bbox intersected mines but none actually contained habitat cells --
          # still record a zero row.
          if (!sim_had_overlap) {
            out_rows[[length(out_rows) + 1]] <- data.frame(
              simulation  = sm$id,
              species_id  = sid,
              point_id    = NA_integer_,
              PROP_ID     = NA_integer_,
              overlap_km2 = 0,
              source      = src,
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
          source      = src,
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
