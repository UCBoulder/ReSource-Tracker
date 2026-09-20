## ============================================================================
## dar_comparison_and_loocv.R
##
## ONE script, three methods, two questions. Section 1 (data prep) is your
## actual working data-prep code, unmodified except for one addition (HUC8
## detection, needed for the regional-exponent method) - it produces
## `base_flow_order` exactly as your existing pipeline does.
##
##   METHODS COMPARED:
##     (a) baseline   - current app.R method: QE_MA 50% similarity gate, no scaling
##     (b) dar_state  - drainage-area-ratio gate + ONE statewide fitted exponent
##     (c) dar_region - drainage-area-ratio gate + per-HUC8 fitted exponent
##                       (falls back to the statewide exponent where a HUC8 has
##                       fewer than `min_gages_per_region` gages)
##
##   QUESTION 1 - "how different are the estimates?" (Part A)
##     Full-network run with ALL real gages in place. Compares network
##     coverage, released flow, and %WWef across the three methods.
##
##   QUESTION 2 - "which one is actually more accurate?" (Part B)
##     Leave-one-gage-out cross-validation (LOOCV). Hides each gage in turn,
##     re-propagates from the rest, and checks how close the reconstructed
##     flow at that segment comes to the gage's real reading.
##
##   Standalone script - does not read, source, or modify app.R.
##   Delivered as code for you to run; nothing here executed against your
##   real data on this end.
## ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

out_dir <- "dar_comparison_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

data_dir <- "data-raw"
paths <- list(
  matched_dmr   = file.path(data_dir, "matchedDMRwwtpupdated6.csv"),
  flowlines     = file.path(data_dir, "Flowlines_Colorado_Filtered.gpkg"),
  monthly_flows = file.path(data_dir, "CO_Seasonal_Flow_Summary_Feb2026.csv"),
  gages_sf      = file.path(data_dir, "processed_gages_sf_combinedV4.gpkg")
)
for (p in paths) if (!file.exists(p)) stop("Missing required data file: ", p, call. = FALSE)

mon_codes <- sprintf("%02d", 1:12)
compute_pct <- function(ww, denom) {
  ifelse(is.na(denom) | is.na(ww) | denom <= 0, NA_real_, pmin(100, pmax(0, 100 * ww / denom)))
}

HIGHLIGHT_MONTH      <- "06"   # month used for the headline flow/%WWef/LOOCV comparison
SIMILARITY_THRESHOLD <- 0.50   # current app.R gate (QE_MA-based); also used as the DA-ratio gate
min_gages_per_region <- 15

# ==============================================================================
# 1. DATA PREP  (your existing pipeline, unchanged, + HUC8 detection added)
# ==============================================================================
message("Loading and preparing network data...")

matched_dmr <- read_csv(paths$matched_dmr,
                        col_types = cols(.default = col_guess(), geometry = col_skip(), `...1` = col_skip()),
                        show_col_types = FALSE)

dmr_monthly_wide <- matched_dmr %>%
  mutate(matched_COMID = as.character(matched_COMID), mon_code = sprintf("%02d", as.integer(mon_code))) %>%
  filter(mon_code %in% mon_codes) %>%
  group_by(Facility, matched_COMID, mon_code) %>%
  summarise(mean_dmr_cfs_fac = mean(mean_dmr_cfs, na.rm = TRUE),
            sd_dmr_cfs_fac   = mean(sd_dmr_cfs,   na.rm = TRUE), .groups = "drop") %>%
  mutate(var_dmr_cfs_fac = sd_dmr_cfs_fac^2) %>%
  group_by(matched_COMID, mon_code) %>%
  summarise(WW_mean = sum(mean_dmr_cfs_fac, na.rm = TRUE), WW_var = sum(var_dmr_cfs_fac, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = mon_code, values_from = c(WW_mean, WW_var), names_sep = "_")

flowlines_raw <- st_read(paths$flowlines, quiet = TRUE)

# --- Identify the drainage-area column -------------------------------------
da_candidates <- c("TotDASqKM", "TOTDASQKM", "DrainageArea", "DrainageAreaSqKM",
                   "DA_SqKM", "DASqKM", "AreaSqKM", "CumDrainAreaSqKm", "totdasqkm")
da_field <- da_candidates[da_candidates %in% names(flowlines_raw)][1]

if (is.na(da_field)) {
  numeric_cols <- names(flowlines_raw)[sapply(st_drop_geometry(flowlines_raw), is.numeric)]
  stop(
    "Could not find a drainage-area column automatically. None of these matched: ",
    paste(da_candidates, collapse = ", "),
    "\n\nNumeric columns actually present in your flowlines file:\n  ",
    paste(numeric_cols, collapse = ", "),
    "\n\nSet `da_field <- \"<the correct column name>\"` near the top of this",
    "\nscript (right after this check) and rerun.", call. = FALSE
  )
}
message("Using '", da_field, "' as the drainage-area field.")

# --- Identify the HUC8 field (needed for the per-region exponent method) ----
# NHDPlus flowlines almost universally carry REACHCODE, whose first 8 digits
# are the HUC8 by USGS convention; used only if no direct HUC8 column exists.
huc_candidates <- c("HUC8", "huc8", "REACHCODE", "ReachCode", "reachcode")
huc_field <- huc_candidates[huc_candidates %in% names(flowlines_raw)][1]

if (is.na(huc_field)) {
  stop(
    "Could not find a HUC8 or REACHCODE column for regional grouping. Columns present:\n  ",
    paste(names(flowlines_raw), collapse = ", "),
    "\n\nSet `huc_field <- \"<the correct column name>\"` near the top of this",
    "\nscript and rerun.", call. = FALSE
  )
}
derive_huc8 <- function(x) {
  x <- as.character(x)
  if (all(nchar(na.omit(x)) == 8, na.rm = TRUE)) return(x)
  substr(x, 1, 8)
}
message("Using '", huc_field, "' (derived to HUC8) for the regional-exponent method.")

flowlines <- flowlines_raw %>%
  mutate(COMID = as.character(COMID), Hydroseq = as.character(Hydroseq), DnHydroseq = as.character(DnHydroseq),
         DA = as.numeric(.data[[da_field]]),
         HUC8 = derive_huc8(.data[[huc_field]])) %>%
  left_join(dmr_monthly_wide, by = c("COMID" = "matched_COMID"))

flow_map <- flowlines %>%
  st_drop_geometry() %>%
  select(Hydroseq, DnHydroseq, starts_with("WW_mean_"), starts_with("WW_var_")) %>%
  mutate(across(starts_with("WW_"), as.numeric)) %>%
  group_by(Hydroseq) %>%
  summarise(DnHydroseq = first(DnHydroseq),
            across(starts_with("WW_mean_"), ~sum(.x, na.rm = TRUE)),
            across(starts_with("WW_var_"),  ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(desc(as.numeric(Hydroseq)))

idx <- setNames(seq_len(nrow(flow_map)), flow_map$Hydroseq)
dn_idx <- idx[flow_map$DnHydroseq]
for (mm in mon_codes) {
  flow_map[[paste0("WW_cum_mean_", mm)]] <- coalesce(flow_map[[paste0("WW_mean_", mm)]], 0)
  flow_map[[paste0("WW_cum_var_",  mm)]] <- coalesce(flow_map[[paste0("WW_var_",  mm)]],  0)
}
for (i in seq_len(nrow(flow_map))) {
  j <- dn_idx[i]
  if (!is.na(j)) {
    for (mm in mon_codes) {
      flow_map[[paste0("WW_cum_mean_", mm)]][j] <- flow_map[[paste0("WW_cum_mean_", mm)]][j] + flow_map[[paste0("WW_cum_mean_", mm)]][i]
      flow_map[[paste0("WW_cum_var_",  mm)]][j] <- flow_map[[paste0("WW_cum_var_",  mm)]][j] + flow_map[[paste0("WW_cum_var_",  mm)]][i]
    }
  }
}
ww_local_wide <- flow_map %>% select(Hydroseq, starts_with("WW_mean_"), starts_with("WW_var_"))
ww_cum_wide   <- flow_map %>% select(Hydroseq, starts_with("WW_cum_mean_"))

monthly_flows <- read_csv(paths$monthly_flows, show_col_types = FALSE)
gage_monthlies <- monthly_flows %>%
  select(site_no, starts_with("mean_cfs_"), starts_with("sd_cfs_"), start_year, end_year, n_obs) %>%
  rename_with(~ gsub("mean_cfs_", "Qgage_", .x), starts_with("mean_cfs_")) %>%
  rename_with(~ gsub("sd_cfs_", "Qgage_sd_", .x), starts_with("sd_cfs_"))
for (mm in mon_codes) {
  sd_name <- paste0("Qgage_sd_", mm); var_name <- paste0("Qgage_var_", mm)
  if (sd_name %in% names(gage_monthlies)) gage_monthlies[[var_name]] <- gage_monthlies[[sd_name]]^2
}

gages_sf <- st_read(paths$gages_sf, layer = "gages_sf", quiet = TRUE)
gages_sf$best_COMID <- as.character(gages_sf$best_COMID)
gages_per_comid <- gages_sf %>%
  st_drop_geometry() %>% select(site_no, best_COMID) %>%
  left_join(gage_monthlies, by = "site_no") %>%
  group_by(best_COMID) %>% arrange(desc(n_obs), desc(end_year), .by_group = TRUE) %>% slice(1) %>% ungroup()

flowlines <- flowlines %>% left_join(gages_per_comid %>% rename(COMID = best_COMID), by = "COMID")

flow_order <- flowlines %>%
  st_drop_geometry() %>%
  select(Hydroseq, DnHydroseq, site_no, QE_MA, DA, HUC8, starts_with("Qgage_")) %>%
  arrange(as.numeric(Hydroseq)) %>%
  mutate(low_flow = if_all(all_of(paste0("Qgage_", mon_codes)), ~ . < 0.5),
         site_no = if_else(low_flow, NA_character_, site_no))
for (mm in mon_codes) {
  mcol <- paste0("Qgage_", mm); sdcol <- paste0("Qgage_sd_", mm); vcol <- paste0("Qgage_var_", mm)
  if (mcol  %in% names(flow_order)) flow_order[[mcol]]  <- ifelse(flow_order$low_flow, NA_real_, flow_order[[mcol]])
  if (sdcol %in% names(flow_order)) flow_order[[sdcol]] <- ifelse(flow_order$low_flow, NA_real_, flow_order[[sdcol]])
  if (vcol  %in% names(flow_order)) flow_order[[vcol]]  <- ifelse(flow_order$low_flow, NA_real_, flow_order[[vcol]])
}
flow_order <- flow_order %>% select(-low_flow)

forced_hseq <- c("550025960","550020503","390005386","390025451","390005505")
same_seg_ww_adj_seed <- c("550022475","550015536","550062530","390045364","550025946",
                          "550019184","760015713","390005515","550018746")

flow_order <- flow_order %>%
  mutate(force_gage = Hydroseq %in% forced_hseq,
         gage_hseq  = if_else(!is.na(site_no), Hydroseq, NA_character_))
dn_index <- match(flow_order$DnHydroseq, flow_order$Hydroseq)
flow_order$gage_used <- flow_order$site_no
for (mm in mon_codes) {
  mcol <- paste0("Qgage_", mm); vcol <- paste0("Qgage_var_", mm)
  flow_order[[paste0(mcol, "_prop")]] <- as.numeric(flow_order[[mcol]])
  flow_order[[paste0(vcol, "_prop")]] <- as.numeric(flow_order[[vcol]])
}
flow_order$is_confluence <- flow_order$Hydroseq %in% flow_order$DnHydroseq[duplicated(flow_order$DnHydroseq)]
flow_order <- flow_order %>%
  left_join(ww_local_wide, by = "Hydroseq") %>%
  mutate(across(starts_with("WW_mean_"), ~coalesce(as.numeric(.x), 0)),
         across(starts_with("WW_var_"),  ~coalesce(as.numeric(.x), 0))) %>%
  ## WW_cum_mean_MM = wastewater accumulated from every upstream WWTP (built
  ## once via flow_map above) - this, not the LOCAL WW_mean_MM, is what your
  ## app.R actually divides by when computing %WW. WW_mean_MM stays in the
  ## table only because propagation needs it (self-adjustment at gaged seeds,
  ## and the local addition at each newly-filled segment).
  left_join(ww_cum_wide, by = "Hydroseq") %>%
  mutate(across(starts_with("WW_cum_mean_"), ~coalesce(as.numeric(.x), 0)))
for (mm in mon_codes) {
  mcol <- paste0("Qgage_", mm, "_prop"); ww_col <- paste0("WW_mean_", mm)
  seed_idx <- which(flow_order$Hydroseq %in% same_seg_ww_adj_seed)
  for (i in seed_idx) flow_order[[mcol]][i] <- flow_order[[mcol]][i] + coalesce(as.numeric(flow_order[[ww_col]][i]), 0)
}
base_flow_order <- flow_order
message("Data prep complete: ", nrow(base_flow_order), " network segments.")

## Note: this reproduces your BFS-fill data prep exactly through
## `base_flow_order`. Confluence handling (below, Section 2b) is ported
## faithfully from your app.R code, verbatim in logic - see notes there.

# ==============================================================================
# 2. SHARED PROPAGATION ENGINE (matrix-based, preallocated queue)
# ==============================================================================
# gate_mode: "qema" -> baseline app.R gate on QE_MA, no scaling
#            "dar"  -> gate on DA, scale flow by (DA_j/DA_i)^b_vec[j]
# Reproduces the same allow_propagate / should_add_ww logic as your BFS above:
#   - a forced hydroseq (force_gage) always propagates regardless of the gate
#   - WW is added at j only when j starts out ungaged (seed rows are already
#     marked visited, so the loop only ever reaches j through propagation when
#     j wasn't already a gage) - the same-segment WW self-adjustment for
#     `same_seg_ww_adj_seed` is already baked into base_flow_order's `_prop`
#     seed values above, so it does not need to be repeated here.

n_seg     <- nrow(base_flow_order)
QE_MA_vec <- base_flow_order$QE_MA
DA_vec    <- base_flow_order$DA
force_vec <- base_flow_order$force_gage

## ---- confluence-handling support structures (used by Section 2b) ----------
Hydroseq_vec      <- base_flow_order$Hydroseq
DnHydroseq_vec    <- base_flow_order$DnHydroseq
site_no_vec       <- base_flow_order$site_no
is_confluence_vec <- base_flow_order$is_confluence
hseq_to_idx       <- setNames(seq_len(n_seg), Hydroseq_vec)
confluence_rows   <- which(is_confluence_vec)

single_parent_bfs <- function(seed_rows, mon, gate_mode = c("qema", "dar"),
                              threshold = SIMILARITY_THRESHOLD, b_vec = NULL) {
  gate_mode <- match.arg(gate_mode)
  
  Qgage_col <- paste0("Qgage_", mon, "_prop")
  Vgage_col <- paste0("Qgage_var_", mon, "_prop")
  WWm_col   <- paste0("WW_mean_", mon)
  WWv_col   <- paste0("WW_var_", mon)
  
  Qm <- rep(NA_real_, n_seg)
  Vm <- rep(NA_real_, n_seg)
  Qm[seed_rows] <- base_flow_order[[Qgage_col]][seed_rows]
  Vm[seed_rows] <- base_flow_order[[Vgage_col]][seed_rows]
  
  visited <- rep(FALSE, n_seg)
  visited[seed_rows] <- TRUE
  is_seed <- visited
  
  buf <- integer(n_seg)
  n0  <- length(seed_rows)
  buf[seq_len(n0)] <- seed_rows
  head_p <- 1L; tail_p <- n0
  
  WWm <- base_flow_order[[WWm_col]]
  WWv <- base_flow_order[[WWv_col]]
  
  while (head_p <= tail_p) {
    i <- buf[head_p]; head_p <- head_p + 1L
    j <- dn_index[i]
    if (is.na(j) || visited[j]) next
    
    if (gate_mode == "qema") {
      qi <- QE_MA_vec[i]; qj <- QE_MA_vec[j]
      gate_ok <- is.finite(qi) && is.finite(qj) && abs(qi - qj) / max(qj, 1e-6) < threshold
      scale_factor <- 1
    } else {
      ai <- DA_vec[i]; aj <- DA_vec[j]
      gate_ok <- is.finite(ai) && is.finite(aj) && ai > 0 && aj > 0 &&
        abs(ai - aj) / max(aj, 1e-6) < threshold
      scale_factor <- if (is.finite(ai) && is.finite(aj) && ai > 0) (aj / ai) ^ b_vec[j] else 1
    }
    allow_propagate <- gate_ok || isTRUE(force_vec[j])
    if (!allow_propagate) next
    if (!gate_ok && isTRUE(force_vec[j])) scale_factor <- 1  # forced hop: pass through unscaled
    
    ww_local <- if (!is_seed[j]) WWm[j] else 0
    wv_local <- if (!is_seed[j]) WWv[j] else 0
    
    Qm[j] <- Qm[i] * scale_factor + ww_local
    Vm[j] <- Vm[i] * (scale_factor ^ 2) + wv_local
    visited[j] <- TRUE
    
    tail_p <- tail_p + 1L
    buf[tail_p] <- j
  }
  
  list(Q = Qm, V = Vm)
}

# ==============================================================================
# 2b. CONFLUENCE HANDLING - ported verbatim (in logic) from app.R
# ==============================================================================
# Your app.R does this as a SEPARATE, UNGATED pass after the main BFS:
#   for each confluence Hydroseq with >=2 upstream branches that already carry
#   a value, sum those branches' mean+variance, then walk downstream from the
#   confluence - NO similarity/DA gate check at all - overwriting each
#   segment's mean+variance with the running sum (+ that segment's own local
#   WW, added unconditionally since by construction every segment on this walk
#   is ungaged) until hitting a segment with a real direct gage, where the
#   walk stops without overwriting it.
#
# This is intentionally NOT gated and NOT scaled - it mirrors your code
# exactly, including that a confluence's downstream walk can overwrite values
# the main BFS above already assigned via a different upstream path.
#
# Applied identically after EITHER gate_mode's main BFS - the branches
# reaching the confluence already carry whatever scaling that method applied
# to get there; the merge/summation step itself does not add further scaling,
# same as your app.R (which only ever has one propagation method).

apply_confluences_fast <- function(Qm, Vm, mon, confluence_rows_use = confluence_rows) {
  WWm <- base_flow_order[[paste0("WW_mean_", mon)]]
  WWv <- base_flow_order[[paste0("WW_var_",  mon)]]
  
  for (crow in confluence_rows_use) {
    hseq_c <- Hydroseq_vec[crow]
    upstream_idxs <- which(DnHydroseq_vec == hseq_c & !is.na(Qm))
    if (length(upstream_idxs) < 2) next
    
    summed_mean <- sum(Qm[upstream_idxs], na.rm = TRUE)
    summed_var  <- sum(Vm[upstream_idxs], na.rm = TRUE)
    
    current_hseq <- hseq_c
    repeat {
      idx <- unname(hseq_to_idx[current_hseq])   # NA if current_hseq not found (network edge)
      if (is.na(idx) || !is.na(site_no_vec[idx])) break
      
      ww_local <- coalesce(WWm[idx], 0)
      wv_local <- coalesce(WWv[idx], 0)
      summed_mean <- summed_mean + ww_local
      summed_var  <- summed_var  + wv_local
      
      Qm[idx] <- summed_mean
      Vm[idx] <- summed_var
      
      current_hseq <- DnHydroseq_vec[idx]
      if (is.na(current_hseq)) break
    }
  }
  
  list(Q = Qm, V = Vm)
}

## propagate_old() = your CURRENT app.R pipeline, unchanged: gated
## single-parent BFS, THEN the ungated confluence merge/overwrite pass.
## Kept exactly as-is so you can compare it against the fix below.
## `confluence_rows_use` lets a caller (LOOCV) restrict which confluence nodes
## get checked - see the ancestor-restriction note in Section 2c below.
propagate_old <- function(seed_rows, mon, gate_mode = c("qema", "dar"),
                          threshold = SIMILARITY_THRESHOLD, b_vec = NULL,
                          confluence_rows_use = confluence_rows) {
  main <- single_parent_bfs(seed_rows, mon, gate_mode, threshold, b_vec)
  apply_confluences_fast(main$Q, main$V, mon, confluence_rows_use)
}

# ==============================================================================
# 2c. NEW: confluence handling using the SAME gate as ordinary propagation
# ==============================================================================
# Fix discussed: a confluence stops being a special ungated case. Instead:
#   - Sum upstream branches' mean+variance at the confluence node itself
#     (mass balance - still unconditional, no gate needed for a true merge
#     of two already-validated tributary estimates).
#   - Everywhere else - single-predecessor nodes, and continuing downstream
#     FROM a confluence - use the exact same gated single-hop rule as the
#     rest of the network (QE_MA or DA-ratio, whichever gate_mode is active).
#
# This is implemented as ONE topological pass instead of two separate passes
# (BFS + patch): NHDPlus Hydroseq strictly decreases downstream (checked
# earlier in your data prep), so processing rows in descending-Hydroseq order
# visits every node only after all of its upstream neighbors are already
# finalized. That makes confluence merging and single-hop gating fall out of
# the same loop, with no separate walk needed - and it naturally stops the
# unconditional propagation problem, since anything downstream of a merge is
# gated exactly like anywhere else.

## Precompute once (independent of seed_rows / month / method):
##   pred_list[[j]] = row indices of every direct upstream neighbor of j
##   order_idx      = row processing order, most-upstream-first
valid_edges <- which(!is.na(dn_index))
pred_split  <- split(valid_edges, dn_index[valid_edges])
pred_list   <- vector("list", n_seg)
pred_list[as.integer(names(pred_split))] <- pred_split
pred_list[sapply(pred_list, is.null)] <- list(integer(0))

order_idx <- order(as.numeric(Hydroseq_vec), decreasing = TRUE)

## ---- ancestor restriction (performance) -----------------------------------
## LOOCV only needs the reconstructed value at ONE segment (the held-out
## gage) per call - everything that isn't upstream of it can't possibly
## affect that value. get_ancestors() walks pred_list backward from a target
## row to find its full upstream watershed once; Part B then restricts the
## topological order (propagate_new) and confluence check list (propagate_old)
## to that watershed instead of scanning all ~66,000 segments on every one of
## the ~28,000 LOOCV propagation calls. This changes nothing about the math -
## it just skips segments that were always going to be no-ops for that call.
get_ancestors <- function(target_row) {
  visited <- logical(n_seg)
  visited[target_row] <- TRUE
  queue <- target_row
  head_p <- 1L
  while (head_p <= length(queue)) {
    node <- queue[head_p]; head_p <- head_p + 1L
    preds <- pred_list[[node]]
    new_preds <- preds[!visited[preds]]
    if (length(new_preds) > 0) {
      visited[new_preds] <- TRUE
      queue <- c(queue, new_preds)
    }
  }
  which(visited)
}

propagate_new <- function(seed_rows, mon, gate_mode = c("qema", "dar"),
                          threshold = SIMILARITY_THRESHOLD, b_vec = NULL,
                          order_idx_use = order_idx) {
  gate_mode <- match.arg(gate_mode)
  
  Qgage_col <- paste0("Qgage_", mon, "_prop")
  Vgage_col <- paste0("Qgage_var_", mon, "_prop")
  WWm <- base_flow_order[[paste0("WW_mean_", mon)]]
  WWv <- base_flow_order[[paste0("WW_var_",  mon)]]
  
  Qm <- rep(NA_real_, n_seg)
  Vm <- rep(NA_real_, n_seg)
  Qm[seed_rows] <- base_flow_order[[Qgage_col]][seed_rows]
  Vm[seed_rows] <- base_flow_order[[Vgage_col]][seed_rows]
  
  is_seed <- rep(FALSE, n_seg)
  is_seed[seed_rows] <- TRUE
  
  for (j in order_idx_use) {
    if (is_seed[j]) next   # real (non-held-out) gage - keep its own value, never overwrite
    
    preds     <- pred_list[[j]]
    if (length(preds) == 0) next
    preds_val <- preds[is.finite(Qm[preds])]
    n_val     <- length(preds_val)
    if (n_val == 0) next
    
    ww_local <- coalesce(WWm[j], 0)
    wv_local <- coalesce(WWv[j], 0)
    
    if (n_val == 1) {
      ## ordinary single-hop propagation - same gate as everywhere else
      i <- preds_val
      if (gate_mode == "qema") {
        qi <- QE_MA_vec[i]; qj <- QE_MA_vec[j]
        gate_ok <- is.finite(qi) && is.finite(qj) && abs(qi - qj) / max(qj, 1e-6) < threshold
        scale_factor <- 1
      } else {
        ai <- DA_vec[i]; aj <- DA_vec[j]
        gate_ok <- is.finite(ai) && is.finite(aj) && ai > 0 && aj > 0 &&
          abs(ai - aj) / max(aj, 1e-6) < threshold
        scale_factor <- if (is.finite(ai) && is.finite(aj) && ai > 0) (aj / ai) ^ b_vec[j] else 1
      }
      allow_propagate <- gate_ok || isTRUE(force_vec[j])
      if (!allow_propagate) next
      if (!gate_ok && isTRUE(force_vec[j])) scale_factor <- 1
      
      Qm[j] <- Qm[i] * scale_factor + ww_local
      Vm[j] <- Vm[i] * (scale_factor ^ 2) + wv_local
      
    } else {
      ## true confluence: >=2 branches carrying a value - sum unconditionally
      ## (mass balance; no gate on the merge itself), then j is available to
      ## be gated normally by whatever comes after it, same as any other node
      Qm[j] <- sum(Qm[preds_val]) + ww_local
      Vm[j] <- sum(Vm[preds_val]) + wv_local
    }
  }
  
  list(Q = Qm, V = Vm)
}

## Both pipelines share this signature - callers pick which to run:
##   propagate_old(seed_rows, mon, gate_mode, threshold, b_vec)  - current app.R
##   propagate_new(seed_rows, mon, gate_mode, threshold, b_vec)  - proposed fix

# ==============================================================================
# 3. EXPONENT FITTING (statewide + per-HUC8), refittable excluding one gage
# ==============================================================================

qgage_raw_cols <- paste0("Qgage_", mon_codes)   # raw gage values, NOT _prop, for fitting
base_flow_order$ann_mean_flow <- rowMeans(as.matrix(base_flow_order[, qgage_raw_cols]), na.rm = TRUE)
base_flow_order$ann_mean_flow[is.nan(base_flow_order$ann_mean_flow)] <- NA_real_

gage_rows <- which(!is.na(base_flow_order$site_no) & is.finite(base_flow_order$ann_mean_flow))
cat(sprintf("Gaged segments usable for exponent fit + LOOCV: %d\n", length(gage_rows)))

fit_b_lookup <- function(exclude_row = NA_integer_) {
  fit_rows <- gage_rows
  if (!is.na(exclude_row)) fit_rows <- setdiff(fit_rows, exclude_row)
  
  fit_data <- base_flow_order[fit_rows, ] %>%
    filter(is.finite(ann_mean_flow), is.finite(DA), ann_mean_flow > 0, DA > 0)
  
  statewide_fit <- lm(log(ann_mean_flow) ~ log(DA), data = fit_data)
  b_statewide   <- unname(coef(statewide_fit)[2])
  r2_statewide  <- summary(statewide_fit)$r.squared
  
  region_b <- fit_data %>%
    group_by(HUC8) %>%
    filter(n() >= min_gages_per_region) %>%
    group_modify(~ {
      m <- tryCatch(lm(log(ann_mean_flow) ~ log(DA), data = .x), error = function(e) NULL)
      tibble(
        n_gages  = nrow(.x),
        b_region = if (is.null(m)) NA_real_ else unname(coef(m)[2]),
        r2       = if (is.null(m)) NA_real_ else summary(m)$r.squared
      )
    }) %>%
    ungroup()
  
  b_lookup <- setNames(region_b$b_region, region_b$HUC8)
  b_vec <- rep(b_statewide, n_seg)
  has_region <- base_flow_order$HUC8 %in% names(b_lookup)
  b_vec[has_region] <- b_lookup[base_flow_order$HUC8[has_region]]
  b_vec[!is.finite(b_vec)] <- b_statewide
  
  list(b_vec = b_vec, b_statewide = b_statewide, r2_statewide = r2_statewide,
       region_table = region_b, n_gages_fit = nrow(fit_data))
}

full_fit <- fit_b_lookup(exclude_row = NA_integer_)
cat(sprintf("Statewide exponent (full data, n=%d): b = %.3f, R^2 = %.3f\n",
            full_fit$n_gages_fit, full_fit$b_statewide, full_fit$r2_statewide))
write_csv(full_fit$region_table, file.path(out_dir, "0_region_exponent_fits.csv"))

# ==============================================================================
# PART A - full-network comparison (all real gages in place, no hold-out)
# ==============================================================================
# Now run for BOTH confluence pipelines (old = current app.R, new = gated fix)
# crossed with all 3 methods, so you can see each effect separately:
#   - method (baseline / dar_statewide / dar_regional) changes the DENOMINATOR
#     assumption everywhere
#   - confluence pipeline (old / new) changes only how confluences behave
# 6 combinations total.

## variant definitions: one row per (method) - confluence pipeline is looped separately
method_variants <- list(
  list(label = "baseline",      gate_mode = "qema", b_vec = NULL),
  list(label = "dar_statewide", gate_mode = "dar",  b_vec = rep(full_fit$b_statewide, n_seg)),
  list(label = "dar_regional",  gate_mode = "dar",  b_vec = full_fit$b_vec)
)
confluence_variants <- list(old = propagate_old, new = propagate_new)

run_full_network <- function(propagate_fn, gate_mode, b_vec, label, confl_label) {
  per_month <- map_dfr(mon_codes, function(mm) {
    run <- propagate_fn(gage_rows, mm, gate_mode = gate_mode, b_vec = b_vec)
    tibble(
      month           = mm,
      n_gaged_segs    = sum(is.finite(base_flow_order[[paste0("Qgage_", mm)]][gage_rows])),
      n_propagated    = sum(is.finite(run$Q)) - length(gage_rows),
      pct_network_gaged_or_filled = round(100 * sum(is.finite(run$Q)) / n_seg, 2),
      released_flow_median = median(run$Q, na.rm = TRUE),
      released_flow_mean   = mean(run$Q, na.rm = TRUE)
    )
  })
  per_month$method            <- label
  per_month$confluence_method <- confl_label
  per_month
}

coverage_comparison <- map_dfr(names(confluence_variants), function(confl_label) {
  pfn <- confluence_variants[[confl_label]]
  map_dfr(method_variants, function(v) {
    run_full_network(pfn, v$gate_mode, v$b_vec, v$label, confl_label)
  })
}) %>%
  select(method, confluence_method, everything())
write_csv(coverage_comparison, file.path(out_dir, "1_network_coverage_by_method.csv"))

compute_pctWWef <- function(run, mon) {
  ## WW_cum_mean_MM = wastewater accumulated from every upstream WWTP (matches
  ## your app.R's actual %WW denominator logic) - NOT the local-only WW_mean_MM.
  ww_col <- paste0("WW_cum_mean_", mon)
  compute_pct(base_flow_order[[ww_col]], run$Q)
}

pctWWef_summary <- map_dfr(names(confluence_variants), function(confl_label) {
  pfn <- confluence_variants[[confl_label]]
  map_dfr(method_variants, function(v) {
    run <- pfn(gage_rows, HIGHLIGHT_MONTH, v$gate_mode, b_vec = v$b_vec)
    tibble(method = v$label, confluence_method = confl_label,
           pctWWef = compute_pctWWef(run, HIGHLIGHT_MONTH))
  })
}) %>%
  filter(is.finite(pctWWef)) %>%
  group_by(method, confluence_method) %>%
  summarize(
    n      = n(),
    median = round(median(pctWWef), 3),
    mean   = round(mean(pctWWef), 3),
    p90    = round(quantile(pctWWef, 0.90), 3),
    .groups = "drop"
  )
write_csv(pctWWef_summary, file.path(out_dir, sprintf("2_pctWWef_summary_month_%s.csv", HIGHLIGHT_MONTH)))

cat("\n=== PART A: full-network comparison (method x confluence pipeline) ===\n")
print(coverage_comparison, n = 100)
cat(sprintf("\n%%WWef summary (month %s):\n", HIGHLIGHT_MONTH))
print(pctWWef_summary)

# ==============================================================================
# PART B - leave-one-gage-out cross-validation (accuracy check), ALL 12 MONTHS
# ==============================================================================
# Run across every month, not just HIGHLIGHT_MONTH, since month-to-month flow
# can differ a lot and a single-month accuracy number could be an outlier.
# Now ALSO crossed with both confluence pipelines (old vs new), so each of the
# 3 methods gets tested both ways: 6 prediction columns per gage x month.
#
# Cost-saving note: the DAR exponent is fit on ANNUAL mean flow (ann_mean_flow),
# which doesn't depend on month or confluence pipeline - so for a given
# held-out gage, the "exclude this gage" exponent fit is done once here and
# reused for every month AND both confluence pipelines.
#
# PERFORMANCE: LOOCV only needs the reconstructed value at the held-out
# segment, not the whole 66,000-segment network, on every call. For each
# held-out gage, we precompute its upstream ancestor set ONCE (before the
# month loop) and restrict every propagate_old()/propagate_new() call for
# that gage to just that watershed - the topological loop in propagate_new()
# and the confluence-check loop in propagate_old() then only touch the
# handful of segments that could possibly matter for this gage, instead of
# scanning the entire state 28,000+ times. This changes nothing about the
# math, only what gets skipped as a guaranteed no-op.
#
# If it's still too slow after this, the next levers are commenting out a
# method in `method_variants` above, or shrinking mon_codes to a few
# representative months.

safe_pct_err <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs) & obs != 0
  abs(pred[ok] - obs[ok]) / abs(obs[ok]) * 100
}

loocv_results <- vector("list", length(gage_rows) * length(mon_codes))
r <- 0L

for (k in seq_along(gage_rows)) {
  hold_row  <- gage_rows[k]
  
  ## restrict this gage's work to its own upstream watershed (see note above)
  ancestors_k        <- get_ancestors(hold_row)
  seed_rows          <- intersect(setdiff(gage_rows, hold_row), ancestors_k)
  order_idx_k        <- order_idx[order_idx %in% ancestors_k]
  confluence_rows_k  <- confluence_rows[confluence_rows %in% ancestors_k]
  
  ## exponent fit excluding this gage - same for every month, so done once here
  fit_k <- fit_b_lookup(exclude_row = hold_row)
  b_state_vec_k  <- rep(fit_k$b_statewide, n_seg)
  b_region_vec_k <- fit_k$b_vec
  
  method_variants_k <- list(
    list(label = "baseline",      gate_mode = "qema", b_vec = NULL),
    list(label = "dar_state",     gate_mode = "dar",  b_vec = b_state_vec_k),
    list(label = "dar_region",    gate_mode = "dar",  b_vec = b_region_vec_k)
  )
  
  for (mm in mon_codes) {
    obs_val <- base_flow_order[[paste0("Qgage_", mm)]][hold_row]
    
    row_vals <- list(Hydroseq = base_flow_order$Hydroseq[hold_row],
                     site_no  = base_flow_order$site_no[hold_row],
                     month    = mm,
                     observed = obs_val)
    
    for (v in method_variants_k) {
      row_vals[[paste0("pred_", v$label, "_old")]] <-
        propagate_old(seed_rows, mm, v$gate_mode, b_vec = v$b_vec,
                      confluence_rows_use = confluence_rows_k)$Q[hold_row]
      row_vals[[paste0("pred_", v$label, "_new")]] <-
        propagate_new(seed_rows, mm, v$gate_mode, b_vec = v$b_vec,
                      order_idx_use = order_idx_k)$Q[hold_row]
    }
    
    r <- r + 1L
    loocv_results[[r]] <- as_tibble(row_vals)
  }
  
  if (k %% 25 == 0) cat(sprintf("  LOOCV progress: %d / %d gages (all 12 months, both confluence pipelines)\n", k, length(gage_rows)))
}

loocv_raw <- bind_rows(loocv_results)
write_csv(loocv_raw, file.path(out_dir, "3_loocv_raw_predictions_all_months.csv"))

pred_cols <- grep("^pred_", names(loocv_raw), value = TRUE)

summarize_method <- function(df, pred_col) {
  pred <- df[[pred_col]]
  obs  <- df$observed
  reconstructed <- is.finite(pred)
  pe <- safe_pct_err(pred, obs)
  err_cfs <- (pred - obs)[reconstructed & is.finite(obs)]
  
  tibble(
    method               = sub("^pred_(.*)_(old|new)$", "\\1", pred_col),
    confluence_method    = sub("^pred_.*_(old|new)$", "\\1", pred_col),
    n_gage_months_total  = nrow(df),
    n_reconstructed      = sum(reconstructed),
    pct_reconstructed    = round(100 * sum(reconstructed) / nrow(df), 1),
    median_abs_pct_error = round(median(pe, na.rm = TRUE), 1),
    rmse_cfs             = round(sqrt(mean(err_cfs^2, na.rm = TRUE)), 2),
    median_bias_cfs      = round(median(err_cfs, na.rm = TRUE), 2),
    pearson_r = {
      ok <- reconstructed & is.finite(obs)
      if (sum(ok) < 2 || sd(pred[ok]) == 0 || sd(obs[ok]) == 0) NA_real_
      else round(suppressWarnings(cor(pred[ok], obs[ok])), 3)
    }
  )
}

## Overall summary, pooling all 12 months together - one row per method x confluence pipeline
loocv_summary_overall <- map_dfr(pred_cols, ~ summarize_method(loocv_raw, .x))
write_csv(loocv_summary_overall, file.path(out_dir, "4_loocv_summary_overall_all_months.csv"))

## Per-month breakdown, so you can see whether accuracy (or the old-vs-new gap) varies seasonally
loocv_summary_by_month <- map_dfr(mon_codes, function(mm) {
  df_mm <- loocv_raw %>% filter(month == mm)
  map_dfr(pred_cols, ~ summarize_method(df_mm, .x)) %>% mutate(month = mm, .before = 1)
})
write_csv(loocv_summary_by_month, file.path(out_dir, "5_loocv_summary_by_month.csv"))

plot_data <- loocv_raw %>%
  pivot_longer(cols = starts_with("pred_"), names_to = "variant", values_to = "predicted") %>%
  filter(is.finite(predicted), is.finite(observed))

p <- ggplot(plot_data, aes(x = observed, y = predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ variant, ncol = 3) +
  labs(
    title = "Leave-one-gage-out accuracy check - all 12 months pooled",
    subtitle = "Rows/cols: method x confluence pipeline (_old = current app.R, _new = gated fix)",
    x = "Observed flow (cfs, log scale)",
    y = "Reconstructed flow (cfs, log scale)"
  ) +
  theme_bw()
ggsave(file.path(out_dir, "6_loocv_observed_vs_predicted_all_months.png"), p, width = 13, height = 8, dpi = 150)

p_month <- ggplot(loocv_summary_by_month,
                  aes(x = month, y = median_abs_pct_error, color = method,
                      linetype = confluence_method, group = interaction(method, confluence_method))) +
  geom_line() + geom_point() +
  labs(
    title = "Median absolute % error by month, method, and confluence pipeline",
    x = "Month", y = "Median absolute % error",
    linetype = "Confluence handling"
  ) +
  theme_bw()
ggsave(file.path(out_dir, "7_loocv_error_by_month.png"), p_month, width = 9, height = 5.5, dpi = 150)

cat("\n=== PART B: LOOCV accuracy check (all 12 months pooled, method x confluence pipeline) ===\n")
print(loocv_summary_overall, n = 20)
cat("\nPer-month breakdown (first 6 rows = January, all method x confluence combos):\n")
print(loocv_summary_by_month, n = 6)

# ---- Final summary -----------------------------------------------------------

cat(sprintf("\nAll outputs written to ./%s/:\n", out_dir))
cat("  0_region_exponent_fits.csv                 - per-HUC8 sample size, b, R^2\n")
cat("  1_network_coverage_by_method.csv           - segments gaged/filled, by month, method x confluence pipeline\n")
cat(sprintf("  2_pctWWef_summary_month_%s.csv          - %%WWef distribution, method x confluence pipeline, month %s\n",
            HIGHLIGHT_MONTH, HIGHLIGHT_MONTH))
cat("  3_loocv_raw_predictions_all_months.csv     - every held-out gage x month: observed vs. 6 predictions\n")
cat("  4_loocv_summary_overall_all_months.csv     - reconstruction rate, error, RMSE, r - pooled, per combo\n")
cat("  5_loocv_summary_by_month.csv               - the same metrics, broken out by month\n")
cat("  6_loocv_observed_vs_predicted_all_months.png - scatter, one panel per method x confluence combo\n")
cat("  7_loocv_error_by_month.png                  - median abs %% error by month; line color = method,\n")
cat("                                                 line type = confluence pipeline (old vs new)\n\n")
cat("Reading this for the reviewer response:\n")
cat("  - Part A tells you HOW MUCH method and confluence handling each change coverage/flow/%WWef.\n")
cat("  - Part B tells you WHICH COMBINATION is most accurate via held-out accuracy, across the full year.\n")
cat("  - Compare _old vs _new within the SAME method first (e.g. dar_state_old vs dar_state_new) to\n")
cat("    isolate the confluence fix's effect on its own, separate from the method choice.\n")
cat("  - Expect _new to reconstruct SLIGHTLY FEWER gage-months than _old (it can now legitimately stop\n")
cat("    at a failed gate downstream of a confluence) - that is the intended, more conservative behavior,\n")
cat("    not a regression. What matters is whether accuracy (median_abs_pct_error, rmse_cfs, bias) improves\n")
cat("    on the gage-months it does reconstruct.\n")

# ==============================================================================
# PART C - SENSITIVITY ANALYSIS FOR THE SELECTED METHOD
#          (DAR, statewide exponent, gated/new confluence handling)
# ==============================================================================
# This is the piece that directly answers the reviewers' original request:
# quantify how much released flows and %WWef move under plausible alternative
# values of the two key assumptions baked into the FINAL adopted method -
#   (C1) the DA-ratio similarity threshold - the direct replacement for the
#        old 50% QE_MA threshold the reviewers asked you to justify
#   (C2) the fitted scaling exponent b - tested at its point estimate AND its
#        95% CI bounds, so uncertainty in the exponent itself is represented,
#        not just uncertainty in the threshold
#
# This is a "how much do the estimates move" sensitivity check (like Part A),
# not an accuracy check (that's Part B) - it uses the full gaged network with
# ALL real gages in place, no LOOCV hold-out.

SELECTED_GATE_MODE <- "dar"
b_statewide_vec    <- rep(full_fit$b_statewide, n_seg)

## ---- C1. Threshold sensitivity ---------------------------------------------
threshold_grid   <- c(0.10, 0.25, 0.50, 0.75, 1.00, 1.50, Inf)
threshold_labels <- c("10%", "25%", "50% (adopted)", "75%", "100%", "150%", "No gate (Inf)")

threshold_coverage_flow <- map2_dfr(threshold_grid, threshold_labels, function(th, lbl) {
  map_dfr(mon_codes, function(mm) {
    run <- propagate_new(gage_rows, mm, gate_mode = SELECTED_GATE_MODE,
                         threshold = th, b_vec = b_statewide_vec)
    tibble(
      threshold_label       = lbl,
      threshold_value       = th,
      month                 = mm,
      pct_network_filled    = round(100 * sum(is.finite(run$Q)) / n_seg, 2),
      released_flow_median  = median(run$Q, na.rm = TRUE),
      released_flow_mean    = mean(run$Q, na.rm = TRUE)
    )
  })
})
write_csv(threshold_coverage_flow, file.path(out_dir, "8_sensitivity_threshold_coverage_flow.csv"))

threshold_pctWWef <- map2_dfr(threshold_grid, threshold_labels, function(th, lbl) {
  run <- propagate_new(gage_rows, HIGHLIGHT_MONTH, gate_mode = SELECTED_GATE_MODE,
                       threshold = th, b_vec = b_statewide_vec)
  tibble(threshold_label = lbl, threshold_value = th,
         pctWWef = compute_pctWWef(run, HIGHLIGHT_MONTH))
}) %>%
  filter(is.finite(pctWWef)) %>%
  group_by(threshold_label, threshold_value) %>%
  summarize(n = n(), median = round(median(pctWWef), 3), mean = round(mean(pctWWef), 3),
            p90 = round(quantile(pctWWef, 0.90), 3), .groups = "drop") %>%
  arrange(threshold_value)
write_csv(threshold_pctWWef, file.path(out_dir, sprintf("9_sensitivity_threshold_pctWWef_month_%s.csv", HIGHLIGHT_MONTH)))

p_threshold <- ggplot(
  threshold_pctWWef %>% mutate(threshold_label = factor(threshold_label, levels = threshold_labels)),
  aes(x = threshold_label, y = median)
) +
  geom_col(fill = "#4B6CB7") +
  geom_text(aes(label = median), vjust = -0.4, size = 3.5) +
  labs(
    title = "Sensitivity of median %WWef to the DA-ratio similarity threshold",
    subtitle = sprintf("DAR statewide exponent, gated confluence handling, month %s", HIGHLIGHT_MONTH),
    x = "Similarity threshold", y = "Median %WWef"
  ) +
  theme_bw()
ggsave(file.path(out_dir, "11_sensitivity_threshold_plot.png"), p_threshold, width = 8, height = 5, dpi = 150)

## ---- C2. Exponent sensitivity (point estimate + 95% CI bounds) -------------
full_fit_data <- base_flow_order[gage_rows, ] %>%
  filter(is.finite(ann_mean_flow), is.finite(DA), ann_mean_flow > 0, DA > 0)
full_lm <- lm(log(ann_mean_flow) ~ log(DA), data = full_fit_data)
b_point <- unname(coef(full_lm)[2])
b_ci    <- unname(confint(full_lm)["log(DA)", ])   # c(lower, upper)

exponent_grid   <- c(b_ci[1], b_point, b_ci[2])
exponent_labels <- c(
  sprintf("Lower 95%% CI (b=%.3f)", b_ci[1]),
  sprintf("Fitted (b=%.3f, adopted)", b_point),
  sprintf("Upper 95%% CI (b=%.3f)", b_ci[2])
)

exponent_coverage_flow <- map2_dfr(exponent_grid, exponent_labels, function(b, lbl) {
  map_dfr(mon_codes, function(mm) {
    run <- propagate_new(gage_rows, mm, gate_mode = SELECTED_GATE_MODE,
                         threshold = SIMILARITY_THRESHOLD, b_vec = rep(b, n_seg))
    tibble(
      exponent_label        = lbl,
      exponent_value        = b,
      month                 = mm,
      pct_network_filled    = round(100 * sum(is.finite(run$Q)) / n_seg, 2),
      released_flow_median  = median(run$Q, na.rm = TRUE),
      released_flow_mean    = mean(run$Q, na.rm = TRUE)
    )
  })
})
write_csv(exponent_coverage_flow, file.path(out_dir, "12_sensitivity_exponent_coverage_flow.csv"))

exponent_pctWWef <- map2_dfr(exponent_grid, exponent_labels, function(b, lbl) {
  run <- propagate_new(gage_rows, HIGHLIGHT_MONTH, gate_mode = SELECTED_GATE_MODE,
                       threshold = SIMILARITY_THRESHOLD, b_vec = rep(b, n_seg))
  tibble(exponent_label = lbl, exponent_value = b,
         pctWWef = compute_pctWWef(run, HIGHLIGHT_MONTH))
}) %>%
  filter(is.finite(pctWWef)) %>%
  group_by(exponent_label, exponent_value) %>%
  summarize(n = n(), median = round(median(pctWWef), 3), mean = round(mean(pctWWef), 3),
            p90 = round(quantile(pctWWef, 0.90), 3), .groups = "drop") %>%
  arrange(exponent_value)
write_csv(exponent_pctWWef, file.path(out_dir, sprintf("13_sensitivity_exponent_pctWWef_month_%s.csv", HIGHLIGHT_MONTH)))

p_exponent <- ggplot(
  exponent_pctWWef %>% mutate(exponent_label = factor(exponent_label, levels = exponent_labels)),
  aes(x = exponent_label, y = median)
) +
  geom_col(fill = "#54B887") +
  geom_text(aes(label = median), vjust = -0.4, size = 3.5) +
  labs(
    title = "Sensitivity of median %WWef to the fitted scaling exponent (95% CI)",
    subtitle = sprintf("DAR statewide, 50%% threshold, gated confluence handling, month %s", HIGHLIGHT_MONTH),
    x = "Exponent", y = "Median %WWef"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(out_dir, "14_sensitivity_exponent_plot.png"), p_exponent, width = 8, height = 5, dpi = 150)

cat("\n=== PART C: sensitivity analysis for the adopted method ===\n")
cat(sprintf("Fitted exponent: b = %.3f (95%% CI %.3f - %.3f)\n", b_point, b_ci[1], b_ci[2]))
cat("\nThreshold sensitivity (%WWef, month ", HIGHLIGHT_MONTH, "):\n", sep = "")
print(threshold_pctWWef)
cat("\nExponent sensitivity (%WWef, month ", HIGHLIGHT_MONTH, "):\n", sep = "")
print(exponent_pctWWef)

cat(sprintf("\nPart C outputs written to ./%s/:\n", out_dir))
cat("  8_sensitivity_threshold_coverage_flow.csv - coverage/flow by month, across threshold values\n")
cat(sprintf("  9_sensitivity_threshold_pctWWef_month_%s.csv - %%WWef distribution, across threshold values\n", HIGHLIGHT_MONTH))
cat("  11_sensitivity_threshold_plot.png          - median %WWef vs. threshold, bar chart\n")
cat("  12_sensitivity_exponent_coverage_flow.csv  - coverage/flow by month, across exponent values (95% CI)\n")
cat(sprintf("  13_sensitivity_exponent_pctWWef_month_%s.csv - %%WWef distribution, across exponent values\n", HIGHLIGHT_MONTH))
cat("  14_sensitivity_exponent_plot.png           - median %WWef vs. exponent, bar chart\n\n")
cat("Reading this for the reviewer response:\n")
cat("  - This is the direct answer to 'was a sensitivity analysis conducted on the threshold,\n")
cat("    and how much do alternative values change the released flows and %WWef estimates' -\n")
cat("    now scoped to the actual adopted method (DAR statewide + gated confluence), not the\n")
cat("    old QE_MA-based method.\n")
cat("  - If %WWef stays roughly flat across most of the threshold grid and only shifts at the\n")
cat("    extremes (very tight or no gate at all), that supports 50% as a reasonable, non-arbitrary\n")
cat("    choice sitting in a stable part of the sensitivity curve.\n")
cat("  - The exponent panel shows whether the uncertainty already reported in your fit\n")
cat("    (the 95% CI) translates into a meaningfully different %WWef estimate - if the three\n")
cat("    bars are close together, that's evidence the conclusions are not fragile to exponent\n")
cat("    uncertainty; if they diverge a lot, that uncertainty needs to be stated as a limitation.\n")

# ==============================================================================
# PART D - IS THE DA-RATIO EXPONENT SEASONAL?
# ==============================================================================
# Everything above fits ONE exponent on ANNUAL MEAN flow (ann_mean_flow) and
# applies it to every month. That's an assumption, not a tested fact - DAR
# exponents can plausibly differ between baseflow-dominated months and
# snowmelt/storm-driven months. This section tests that assumption directly:
#
#   (D1) Refit b separately for each month, using that month's own raw gage
#        flow (Qgage_MM) instead of the annual mean - same regression form,
#        same gage_rows, just done 12 times instead of once.
#   (D2) Ask whether the log(DA) slope actually differs by month in a
#        statistically defensible way, via a nested-model F-test on an
#        interaction term (log(Q) ~ log(DA) * month vs. log(Q) ~ log(DA) +
#        month) on the pooled monthly data - this is the real answer to
#        "is there seasonal change in the exponent," not just eyeballing
#        12 point estimates that will never be numerically identical anyway.
#
# This does NOT feed back into Parts A-C or change the adopted method - it's
# a standalone check of whether a fixed annual exponent is defensible, for
# the reviewer response and/or as a documented limitation.

cat("\n=== PART D: is the DA-ratio exponent seasonal? ===\n")

## ---- D1. Per-month exponent fits -------------------------------------------
fit_b_month <- function(mon) {
  qcol <- paste0("Qgage_", mon)
  fit_data <- base_flow_order[gage_rows, ] %>%
    filter(is.finite(.data[[qcol]]), is.finite(DA), .data[[qcol]] > 0, DA > 0)
  
  m  <- lm(as.formula(paste0("log(", qcol, ") ~ log(DA)")), data = fit_data)
  ci <- confint(m)["log(DA)", ]
  
  tibble(
    month    = mon,
    n_gages  = nrow(fit_data),
    b        = unname(coef(m)[2]),
    se       = unname(summary(m)$coefficients["log(DA)", "Std. Error"]),
    ci_lower = unname(ci[1]),
    ci_upper = unname(ci[2]),
    r2       = summary(m)$r.squared
  )
}

monthly_exponent_fits <- map_dfr(mon_codes, fit_b_month)

## annual (statewide, all-months-pooled-via-mean) fit for reference - already
## computed in Part C as full_lm / b_point / b_ci
monthly_exponent_fits <- monthly_exponent_fits %>%
  mutate(annual_b = b_point, annual_ci_lower = b_ci[1], annual_ci_upper = b_ci[2])

write_csv(monthly_exponent_fits, file.path(out_dir, "15_monthly_exponent_fits.csv"))

cat("\nPer-month exponent fits (vs. annual b = ", sprintf("%.3f", b_point),
    ", 95% CI ", sprintf("%.3f", b_ci[1]), "-", sprintf("%.3f", b_ci[2]), "):\n", sep = "")
print(monthly_exponent_fits, n = 12)

## ---- D2. Formal test: does the log(DA)-log(Q) slope actually vary by month? ----
## Pool every month's (Q, DA) pairs into one long table, then compare a model
## that shares one exponent across all months against one that lets the
## exponent vary by month (interaction term). A significant interaction means
## the exponent is NOT well described as a single seasonally-constant value.
long_data <- map_dfr(mon_codes, function(mm) {
  qcol <- paste0("Qgage_", mm)
  base_flow_order[gage_rows, ] %>%
    filter(is.finite(.data[[qcol]]), is.finite(DA), .data[[qcol]] > 0, DA > 0) %>%
    transmute(Q = .data[[qcol]], DA = DA, month = mm)
})

model_shared_b      <- lm(log(Q) ~ log(DA) + month, data = long_data)
model_month_varies_b <- lm(log(Q) ~ log(DA) * month, data = long_data)
seasonal_b_test <- anova(model_shared_b, model_month_varies_b)

seasonal_b_test_summary <- tibble(
  df1        = seasonal_b_test$Df[2],
  df2        = seasonal_b_test$Res.Df[2],
  F_stat     = round(seasonal_b_test$F[2], 3),
  p_value    = signif(seasonal_b_test[["Pr(>F)"]][2], 4),
  n_obs      = nrow(long_data),
  conclusion = if_else(seasonal_b_test[["Pr(>F)"]][2] < 0.05,
                       "Exponent differs significantly by month (p < 0.05)",
                       "No significant evidence the exponent differs by month (p >= 0.05)")
)
write_csv(seasonal_b_test_summary, file.path(out_dir, "16_seasonal_exponent_test.csv"))

cat("\nFormal test - does log(DA) slope (the exponent) vary by month?\n")
cat(sprintf("  F(%d, %d) = %.3f, p = %s\n",
            seasonal_b_test_summary$df1, seasonal_b_test_summary$df2,
            seasonal_b_test_summary$F_stat, format.pval(seasonal_b_test_summary$p_value, digits = 4)))
cat(sprintf("  %s\n", seasonal_b_test_summary$conclusion))

## ---- Plot: monthly exponent estimates + CI, vs. the adopted annual value ---
p_seasonal <- ggplot(monthly_exponent_fits, aes(x = month, y = b, group = 1)) +
  geom_hline(yintercept = b_point, color = "grey40", linetype = "dashed") +
  geom_ribbon(aes(ymin = annual_ci_lower, ymax = annual_ci_upper), fill = "grey70", alpha = 0.25) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "#4B6CB7") +
  geom_point(size = 2.2, color = "#4B6CB7") +
  labs(
    title = "Monthly DA-ratio exponent estimates vs. the adopted annual value",
    subtitle = sprintf("Dashed line + shaded band = annual exponent (b=%.3f) and its 95%% CI. Points = per-month fits, own 95%% CI.", b_point),
    x = "Month", y = "Fitted exponent (b)"
  ) +
  theme_bw()
ggsave(file.path(out_dir, "17_seasonal_exponent_plot.png"), p_seasonal, width = 9, height = 5.5, dpi = 150)

cat(sprintf("\nPart D outputs written to ./%s/:\n", out_dir))
cat("  15_monthly_exponent_fits.csv   - b, SE, 95% CI, R^2, n per month, vs. the annual value\n")
cat("  16_seasonal_exponent_test.csv  - nested-model F-test for whether the exponent varies by month\n")
cat("  17_seasonal_exponent_plot.png  - per-month b with CI error bars, against the annual b + its CI band\n\n")
cat("Reading this for the reviewer response:\n")
cat("  - If every month's CI error bar overlaps the annual CI band, and the F-test p-value is >= 0.05,\n")
cat("    that's a defensible basis for keeping ONE fixed annual exponent - the apparent month-to-month\n")
cat("    wiggle in D1 is just sampling noise around a single true value, not a real seasonal effect.\n")
cat("  - If several months' CIs sit clearly outside the annual band (especially if they cluster in\n")
cat("    snowmelt/high-flow months like May-July) AND the F-test is significant, that's evidence the\n")
cat("    fixed annual exponent is misspecified in those months - worth flagging as a limitation, and a\n")
cat("    plausible explanation for the May-July positive bias already seen in the DAR LOOCV results\n")
cat("    (Part B) rather than something further threshold or confluence tuning would fix.\n")
cat("  - Either way, this is a check of the ASSUMPTION, not a change to the adopted method - nothing\n")
cat("    above this section is affected by what you find here.\n")
