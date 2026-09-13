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
         across(starts_with("WW_var_"),  ~coalesce(as.numeric(.x), 0)))
for (mm in mon_codes) {
  mcol <- paste0("Qgage_", mm, "_prop"); ww_col <- paste0("WW_mean_", mm)
  seed_idx <- which(flow_order$Hydroseq %in% same_seg_ww_adj_seed)
  for (i in seed_idx) flow_order[[mcol]][i] <- flow_order[[mcol]][i] + coalesce(as.numeric(flow_order[[ww_col]][i]), 0)
}
base_flow_order <- flow_order
message("Data prep complete: ", nrow(base_flow_order), " network segments.")

## Note: this reproduces your BFS-fill data prep exactly through
## `base_flow_order`. If app.R does additional confluence-summing (merging
## multiple upstream branches at a junction) AFTER the single-parent
## propagation below, that step isn't reproduced here since I don't have its
## exact code - the comparison and LOOCV below use single-parent propagation
## only. Paste that function if you want it folded in; it would apply
## identically to all three methods below, so it shouldn't change which
## method looks best, only the absolute coverage/flow numbers.

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

propagate_fast <- function(seed_rows, mon, gate_mode = c("qema", "dar"),
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

run_full_network <- function(gate_mode, b_vec = NULL, label) {
  per_month <- map_dfr(mon_codes, function(mm) {
    run <- propagate_fast(gage_rows, mm, gate_mode = gate_mode, b_vec = b_vec)
    tibble(
      month           = mm,
      n_gaged_segs    = sum(is.finite(base_flow_order[[paste0("Qgage_", mm)]][gage_rows])),
      n_propagated    = sum(is.finite(run$Q)) - length(gage_rows),
      pct_network_gaged_or_filled = round(100 * sum(is.finite(run$Q)) / n_seg, 2),
      released_flow_median = median(run$Q, na.rm = TRUE),
      released_flow_mean   = mean(run$Q, na.rm = TRUE)
    )
  })
  per_month$method <- label
  per_month
}

coverage_baseline   <- run_full_network("qema", label = "baseline")
coverage_dar_state  <- run_full_network("dar", b_vec = rep(full_fit$b_statewide, n_seg), label = "dar_statewide")
coverage_dar_region <- run_full_network("dar", b_vec = full_fit$b_vec, label = "dar_regional")

coverage_comparison <- bind_rows(coverage_baseline, coverage_dar_state, coverage_dar_region) %>%
  select(method, everything())
write_csv(coverage_comparison, file.path(out_dir, "1_network_coverage_by_method.csv"))

compute_pctWWef <- function(run, mon) {
  ww_col <- paste0("WW_mean_", mon)
  compute_pct(base_flow_order[[ww_col]], run$Q)
}

hl_baseline   <- propagate_fast(gage_rows, HIGHLIGHT_MONTH, "qema")
hl_dar_state  <- propagate_fast(gage_rows, HIGHLIGHT_MONTH, "dar", b_vec = rep(full_fit$b_statewide, n_seg))
hl_dar_region <- propagate_fast(gage_rows, HIGHLIGHT_MONTH, "dar", b_vec = full_fit$b_vec)

pctWWef_summary <- bind_rows(
  tibble(method = "baseline",      pctWWef = compute_pctWWef(hl_baseline,   HIGHLIGHT_MONTH)),
  tibble(method = "dar_statewide", pctWWef = compute_pctWWef(hl_dar_state,  HIGHLIGHT_MONTH)),
  tibble(method = "dar_regional",  pctWWef = compute_pctWWef(hl_dar_region, HIGHLIGHT_MONTH))
) %>%
  filter(is.finite(pctWWef)) %>%
  group_by(method) %>%
  summarize(
    n      = n(),
    median = round(median(pctWWef), 3),
    mean   = round(mean(pctWWef), 3),
    p90    = round(quantile(pctWWef, 0.90), 3),
    .groups = "drop"
  )
write_csv(pctWWef_summary, file.path(out_dir, sprintf("2_pctWWef_summary_month_%s.csv", HIGHLIGHT_MONTH)))

cat("\n=== PART A: full-network comparison ===\n")
print(coverage_comparison)
cat(sprintf("\n%%WWef summary (month %s):\n", HIGHLIGHT_MONTH))
print(pctWWef_summary)

# ==============================================================================
# PART B - leave-one-gage-out cross-validation (accuracy check), ALL 12 MONTHS
# ==============================================================================
# Run across every month, not just HIGHLIGHT_MONTH, since month-to-month flow
# can differ a lot and a single-month accuracy number could be an outlier.
#
# Cost-saving note: the DAR exponent is fit on ANNUAL mean flow (ann_mean_flow),
# which doesn't depend on month - so for a given held-out gage, the "exclude
# this gage" exponent fit is identical across all 12 months. fit_b_lookup() is
# therefore called once per gage (not once per gage per month), cutting the
# refitting cost 12x relative to naively nesting month inside the exclude-fit.
# The propagation itself (which does depend on month) still runs fresh for
# every gage x month x method combination, since that's the actual quantity
# being tested.

safe_pct_err <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs) & obs != 0
  abs(pred[ok] - obs[ok]) / abs(obs[ok]) * 100
}

loocv_results <- vector("list", length(gage_rows) * length(mon_codes))
r <- 0L

for (k in seq_along(gage_rows)) {
  hold_row  <- gage_rows[k]
  seed_rows <- setdiff(gage_rows, hold_row)
  
  ## exponent fit excluding this gage - same for every month, so done once here
  fit_k <- fit_b_lookup(exclude_row = hold_row)
  b_state_vec_k  <- rep(fit_k$b_statewide, n_seg)
  b_region_vec_k <- fit_k$b_vec
  
  for (mm in mon_codes) {
    obs_val <- base_flow_order[[paste0("Qgage_", mm)]][hold_row]
    
    pred_base       <- propagate_fast(seed_rows, mm, "qema")$Q[hold_row]
    pred_dar_state  <- propagate_fast(seed_rows, mm, "dar", b_vec = b_state_vec_k)$Q[hold_row]
    pred_dar_region <- propagate_fast(seed_rows, mm, "dar", b_vec = b_region_vec_k)$Q[hold_row]
    
    r <- r + 1L
    loocv_results[[r]] <- tibble(
      Hydroseq        = base_flow_order$Hydroseq[hold_row],
      site_no         = base_flow_order$site_no[hold_row],
      month           = mm,
      observed        = obs_val,
      pred_baseline   = pred_base,
      pred_dar_state  = pred_dar_state,
      pred_dar_region = pred_dar_region
    )
  }
  
  if (k %% 25 == 0) cat(sprintf("  LOOCV progress: %d / %d gages (all 12 months each)\n", k, length(gage_rows)))
}

loocv_raw <- bind_rows(loocv_results)
write_csv(loocv_raw, file.path(out_dir, "3_loocv_raw_predictions_all_months.csv"))

summarize_method <- function(df, pred_col) {
  pred <- df[[pred_col]]
  obs  <- df$observed
  reconstructed <- is.finite(pred)
  pe <- safe_pct_err(pred, obs)
  err_cfs <- (pred - obs)[reconstructed & is.finite(obs)]
  
  tibble(
    method               = pred_col,
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

## Overall summary, pooling all 12 months together
loocv_summary_overall <- bind_rows(
  summarize_method(loocv_raw, "pred_baseline"),
  summarize_method(loocv_raw, "pred_dar_state"),
  summarize_method(loocv_raw, "pred_dar_region")
)
write_csv(loocv_summary_overall, file.path(out_dir, "4_loocv_summary_overall_all_months.csv"))

## Per-month breakdown, so you can see whether accuracy varies seasonally
loocv_summary_by_month <- map_dfr(mon_codes, function(mm) {
  df_mm <- loocv_raw %>% filter(month == mm)
  bind_rows(
    summarize_method(df_mm, "pred_baseline"),
    summarize_method(df_mm, "pred_dar_state"),
    summarize_method(df_mm, "pred_dar_region")
  ) %>% mutate(month = mm, .before = 1)
})
write_csv(loocv_summary_by_month, file.path(out_dir, "5_loocv_summary_by_month.csv"))

plot_data <- loocv_raw %>%
  pivot_longer(cols = starts_with("pred_"), names_to = "method", values_to = "predicted") %>%
  filter(is.finite(predicted), is.finite(observed))

p <- ggplot(plot_data, aes(x = observed, y = predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ method) +
  labs(
    title = "Leave-one-gage-out accuracy check - all 12 months pooled",
    x = "Observed flow (cfs, log scale)",
    y = "Reconstructed flow (cfs, log scale)"
  ) +
  theme_bw()
ggsave(file.path(out_dir, "6_loocv_observed_vs_predicted_all_months.png"), p, width = 10, height = 4, dpi = 150)

p_month <- ggplot(loocv_summary_by_month, aes(x = month, y = median_abs_pct_error, color = method, group = method)) +
  geom_line() + geom_point() +
  labs(
    title = "Median absolute % error by month and method",
    x = "Month", y = "Median absolute % error"
  ) +
  theme_bw()
ggsave(file.path(out_dir, "7_loocv_error_by_month.png"), p_month, width = 8, height = 5, dpi = 150)

cat("\n=== PART B: LOOCV accuracy check (all 12 months pooled) ===\n")
print(loocv_summary_overall)
cat("\nPer-month breakdown (first 12 rows = January, all 3 methods):\n")
print(loocv_summary_by_month, n = 12)

