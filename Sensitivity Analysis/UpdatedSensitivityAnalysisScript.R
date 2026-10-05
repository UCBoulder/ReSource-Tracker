suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})
if (!requireNamespace("patchwork", quietly = TRUE))
  stop('Package "patchwork" is needed for the A/B panel figures. Install it once with: install.packages("patchwork")', call. = FALSE)

# ==============================================================================
# 0. CONFIG
# ==============================================================================
out_dir <- "dar_sensitivity_outputs"
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
n_mon     <- length(mon_codes)

# --- Must match app.R ---------------------------------------------------------
SIMILARITY_THRESHOLD <- 0.50
forced_hseq <- c("550025960", "550020503", "390005386", "390025451", "390005505")
same_seg_ww_adj_seed <- c("550022475", "550015536", "550062530", "390045364", "550025946",
                          "550019184", "760015713", "390005515", "550018746")
fill_facs <- character(0)   # DMR gap-fill; set to whatever app.R uses, e.g. c("GOLDEN CITY OF")

# --- Analysis settings ---------------------------------------------------------
HIGHLIGHT_MONTH    <- "06"
THRESHOLD_GRID     <- c(0.10, 0.25, 0.50, 0.75, 1.00, 1.50, Inf)
N_EXPONENT_STEPS   <- 7      # grid points across the exponent's 95% CI
MIN_MONTHS_ALT_FIT <- 10     # alternative exponent fit: only gages with >= this many months
COMP_TIGHT         <- 0.10   # tight gate for the gate-loosening composition figure
N_BOOT             <- 1000   # bootstrap resamples for the within-gate monthly bias CIs
DOM_PCT            <- 10     # "effluent-dominated" = %WWef >= this (percent) in the adopted run
RUN_LOOCV          <- TRUE   # the slowest section; set FALSE to skip while iterating
set.seed(42)

hl      <- match(HIGHLIGHT_MONTH, mon_codes)
hl_name <- month.name[hl]

# --- Plot style ----------------------------------------------------------------
# Matches the manuscript figures: white background, boxed panels, light grey grid,
# bold centred titles, A/B panel tags, one season colour set used everywhere.
# >>> The four season hex codes below are PLACEHOLDERS: replace them with the exact
#     colours used in your existing figures so the sets match. <<<
COL <- list(accent = "#2a78d6", context = "#c3c2b7", ink = "#0b0b0b", ink2 = "#52514e",
            muted = "#898781", grid = "grey90", band = "#d9d8d2", heat_lo = "#cde2fb", heat_hi = "#0d366b")
ORDINAL_BLUE  <- c("#86b6ef", "#3987e5", "#184f95")          # light -> dark = tighter -> looser gate
SEASON_LEVELS <- c("Spring", "Summer", "Fall", "Winter")
SEASON_COLS   <- c(Spring = "#4DAF4A", Summer = "#FF7F00", Fall = "#A65628", Winter = "#377EB8")
season_of <- function(mon) {                                  # "01".."12" or 1..12 -> season factor
  m <- as.integer(mon)
  factor(dplyr::case_when(m %in% c(12, 1, 2) ~ "Winter", m %in% 3:5 ~ "Spring",
                          m %in% 6:8 ~ "Summer", TRUE ~ "Fall"), levels = SEASON_LEVELS)
}

theme_resource <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = COL$grid, linewidth = 0.4),
      panel.border       = element_rect(colour = "grey40", fill = NA, linewidth = 0.6),
      plot.title         = element_text(face = "bold", hjust = 0.5, size = 12.5),
      axis.title         = element_text(face = "bold", size = 11),
      axis.text          = element_text(colour = "grey25", size = 10),
      legend.position    = "right",
      legend.title       = element_blank(),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold"),
      plot.margin        = margin(10, 12, 8, 10)
    )
}

# A/B panels in a row under one bold title (+ optional grey subtitle)
combine_panels <- function(plots, title, subtitle = NULL, widths = NULL, collect = FALSE) {
  patchwork::wrap_plots(plots, nrow = 1, widths = widths,
                        guides = if (collect) "collect" else "keep") +
    patchwork::plot_annotation(
      title = title, subtitle = subtitle, tag_levels = "A",
      theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
                    plot.subtitle = element_text(hjust = 0.5, colour = "grey30", size = 10.5),
                    plot.tag      = element_text(face = "bold", size = 14)))
}

# Line + dots with the value printed above every dot; the adopted setting is the filled dot.
# d needs columns: x (factor), y (numeric), adopted (logical)
dot_panel <- function(d, xlab, ylab, ymax = NA, digits = 1, suffix = "%", axis_lab = pct_lab) {
  ggplot(d, aes(x, y, group = 1)) +
    geom_line(colour = COL$accent, linewidth = 0.6) +
    geom_point(aes(fill = adopted, size = adopted), shape = 21, colour = COL$accent, stroke = 0.9) +
    geom_text(aes(label = paste0(formatC(y, format = "f", digits = digits), suffix)),
              vjust = -1.2, size = 3.5, fontface = "bold", colour = COL$ink) +
    scale_fill_manual(values = c("TRUE" = COL$accent, "FALSE" = "white"), guide = "none") +
    scale_size_manual(values = c("TRUE" = 4.4, "FALSE" = 3), guide = "none") +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.12)), labels = axis_lab) +
    labs(x = xlab, y = ylab) +
    theme_resource()
}

save_plot <- function(p, name, w = 9, h = 5) {
  ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 300, bg = "white")
}
fmt_num <- scales::label_number(big.mark = ",", drop0trailing = TRUE)
pct_lab <- scales::label_number(suffix = "%", big.mark = ",", drop0trailing = TRUE)

# ==============================================================================
# 1. DATA PREP  (mirrors app.R)
# ==============================================================================
message("Loading and preparing network data...")

# --- DMR wastewater: facility-month -> COMID-month mean + variance -------------
matched_dmr <- read_csv(
  paths$matched_dmr,
  col_types = cols(.default = col_guess(), geometry = col_skip(), `...1` = col_skip()),
  show_col_types = FALSE
)

dmr_monthly_fac <- matched_dmr %>%
  mutate(matched_COMID = as.character(matched_COMID),
         mon_code      = sprintf("%02d", as.integer(mon_code))) %>%
  filter(mon_code %in% mon_codes) %>%
  group_by(Facility, matched_COMID, mon_code) %>%
  summarise(mean_dmr_cfs_fac = mean(mean_dmr_cfs, na.rm = TRUE),
            sd_dmr_cfs_fac   = mean(sd_dmr_cfs,   na.rm = TRUE), .groups = "drop")

if (length(fill_facs) > 0) {   # fill missing months for verified continuous dischargers
  fac_fill <- dmr_monthly_fac %>%
    filter(Facility %in% fill_facs) %>%
    group_by(Facility, matched_COMID) %>%
    summarise(mean_dmr_cfs_fac = mean(mean_dmr_cfs_fac, na.rm = TRUE),
              sd_dmr_cfs_fac   = mean(sd_dmr_cfs_fac,   na.rm = TRUE), .groups = "drop") %>%
    tidyr::crossing(mon_code = mon_codes) %>%
    anti_join(dmr_monthly_fac, by = c("Facility", "matched_COMID", "mon_code"))
  dmr_monthly_fac <- bind_rows(dmr_monthly_fac, fac_fill)
}

dmr_monthly_wide <- dmr_monthly_fac %>%
  mutate(var_dmr_cfs_fac = sd_dmr_cfs_fac^2) %>%
  group_by(matched_COMID, mon_code) %>%
  summarise(WW_mean = sum(mean_dmr_cfs_fac, na.rm = TRUE),
            WW_var  = sum(var_dmr_cfs_fac,  na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = mon_code, values_from = c(WW_mean, WW_var), names_sep = "_")

# --- Flowlines (geometry is not needed here) ------------------------------------
flowlines <- st_read(paths$flowlines, quiet = TRUE) %>% st_drop_geometry()
if (!"TotDASqKM" %in% names(flowlines)) {
  stop("Flowlines are missing TotDASqKM (total drainage area), which DAR propagation requires.",
       call. = FALSE)
}
flowlines <- flowlines %>%
  mutate(COMID      = as.character(COMID),
         Hydroseq   = as.character(Hydroseq),
         DnHydroseq = as.character(DnHydroseq),
         DA         = as.numeric(TotDASqKM)) %>%
  left_join(dmr_monthly_wide, by = c("COMID" = "matched_COMID"))
for (cc in c(paste0("WW_mean_", mon_codes), paste0("WW_var_", mon_codes))) {
  if (!cc %in% names(flowlines)) flowlines[[cc]] <- NA_real_   # month with no DMR data at all
}

# --- Local WW per Hydroseq, then accumulated downstream (matrix form of app.R loop) ---
flow_map <- flowlines %>%
  select(Hydroseq, DnHydroseq, starts_with("WW_mean_"), starts_with("WW_var_")) %>%
  mutate(across(starts_with("WW_"), as.numeric)) %>%
  group_by(Hydroseq) %>%
  summarise(DnHydroseq = first(DnHydroseq),
            across(starts_with("WW_mean_"), ~sum(.x, na.rm = TRUE)),
            across(starts_with("WW_var_"),  ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(desc(as.numeric(Hydroseq)))          # upstream first (Hydroseq falls downstream)

fm_mean <- unname(as.matrix(flow_map[, paste0("WW_mean_", mon_codes)]))
fm_var  <- unname(as.matrix(flow_map[, paste0("WW_var_",  mon_codes)]))
fm_mean[is.na(fm_mean)] <- 0
fm_var[is.na(fm_var)]   <- 0
fm_dn   <- match(flow_map$DnHydroseq, flow_map$Hydroseq)

cum_mean <- fm_mean                              # WW accumulated from ALL upstream WWTPs
for (i in seq_len(nrow(flow_map))) {
  j <- fm_dn[i]
  if (!is.na(j)) cum_mean[j, ] <- cum_mean[j, ] + cum_mean[i, ]
}

# --- Gage monthly stats -> one gage per COMID ------------------------------------
monthly_flows <- read_csv(paths$monthly_flows, show_col_types = FALSE)
gage_monthlies <- monthly_flows %>%
  select(site_no, starts_with("mean_cfs_"), starts_with("sd_cfs_"), start_year, end_year, n_obs) %>%
  rename_with(~ gsub("mean_cfs_", "Qgage_",    .x), starts_with("mean_cfs_")) %>%
  rename_with(~ gsub("sd_cfs_",   "Qgage_sd_", .x), starts_with("sd_cfs_"))
for (mm in mon_codes) {
  sd_name <- paste0("Qgage_sd_", mm)
  if (sd_name %in% names(gage_monthlies)) gage_monthlies[[paste0("Qgage_var_", mm)]] <- gage_monthlies[[sd_name]]^2
}

gages_sf <- st_read(paths$gages_sf, layer = "gages_sf", quiet = TRUE)
gages_per_comid <- gages_sf %>%
  st_drop_geometry() %>%
  mutate(best_COMID = as.character(best_COMID)) %>%
  select(site_no, best_COMID) %>%
  left_join(gage_monthlies, by = "site_no") %>%
  group_by(best_COMID) %>%
  arrange(desc(n_obs), desc(end_year), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

flowlines <- flowlines %>% left_join(gages_per_comid %>% rename(COMID = best_COMID), by = "COMID")

# --- flow_order: low-flow screen exactly as in app.R -----------------------------
flow_order <- flowlines %>%
  select(Hydroseq, DnHydroseq, site_no, DA, starts_with("Qgage_")) %>%
  arrange(as.numeric(Hydroseq)) %>%
  mutate(low_flow = if_all(all_of(paste0("Qgage_", mon_codes)), ~ . < 0.5),
         site_no  = if_else(low_flow, NA_character_, site_no))
for (mm in mon_codes) {
  for (pref in c("Qgage_", "Qgage_sd_", "Qgage_var_")) {
    col <- paste0(pref, mm)
    if (col %in% names(flow_order)) flow_order[[col]] <- ifelse(flow_order$low_flow, NA_real_, flow_order[[col]])
  }
}
flow_order <- flow_order %>% select(-low_flow)

n_seg      <- nrow(flow_order)
DA_vec     <- flow_order$DA
force_vec  <- flow_order$Hydroseq %in% forced_hseq
dn_index   <- match(flow_order$DnHydroseq, flow_order$Hydroseq)

Qraw       <- unname(as.matrix(flow_order[, paste0("Qgage_",     mon_codes)]))   # observed gage monthly means
Vraw       <- unname(as.matrix(flow_order[, paste0("Qgage_var_", mon_codes)]))   # their interannual variances
idx_fm     <- match(flow_order$Hydroseq, flow_map$Hydroseq)
WWm_mat    <- fm_mean[idx_fm, , drop = FALSE]                                    # local WW mean (cfs)
WWv_mat    <- fm_var[idx_fm,  , drop = FALSE]                                    # local WW variance
WWcum_mat  <- cum_mean[idx_fm, , drop = FALSE]                                   # cumulative WW (the %WW numerator)
WWm_mat[is.na(WWm_mat)] <- 0
WWv_mat[is.na(WWv_mat)] <- 0

# Seed values: raw gage flow, with local WW mean AND variance baked in at the
# same-segment seeds (so downstream inherits the right denominator) - as app.R.
Q0 <- Qraw
V0 <- Vraw
seed_ww <- which(flow_order$Hydroseq %in% same_seg_ww_adj_seed)
Q0[seed_ww, ] <- Q0[seed_ww, ] + WWm_mat[seed_ww, ]
V0[seed_ww, ] <- V0[seed_ww, ] + WWv_mat[seed_ww, ]

# Gaged segments usable for the exponent fit, LOOCV and the pair checks
ann_mean_flow <- rowMeans(Qraw, na.rm = TRUE)
ann_mean_flow[is.nan(ann_mean_flow)] <- NA_real_
is_gage       <- !is.na(flow_order$site_no) & is.finite(ann_mean_flow)
gage_rows     <- which(is_gage)
n_months_obs  <- rowSums(is.finite(Qraw))

# Upstream-neighbour structure and topological order (app.R: pred_list / order_idx)
order_idx   <- order(as.numeric(flow_order$Hydroseq), decreasing = TRUE)   # most-upstream first
ord_pos     <- integer(n_seg); ord_pos[order_idx] <- seq_len(n_seg)
valid_edges <- which(!is.na(dn_index))
pred_split  <- split(valid_edges, dn_index[valid_edges])
pred_list   <- vector("list", n_seg)
pred_list[as.integer(names(pred_split))] <- pred_split
pred_list[sapply(pred_list, is.null)] <- list(integer(0))

message("Data prep complete: ", n_seg, " segments, ", length(gage_rows), " usable gages.")

# ==============================================================================
# 2. PROPAGATION ENGINE  (the app.R month-aware loop, parametrised)
# ==============================================================================
#   >= 2 branches carry a value in a month -> unconditional sum (mass balance)
#   exactly 1 branch carries a value        -> ordinary gated, DA-scaled hop
#   0 branches                              -> NA
# Seeds (real, non-held-out gages) are never overwritten. A forced Hydroseq
# bypasses the gate and passes flow through unscaled.
propagate_dar <- function(seed_rows, threshold = SIMILARITY_THRESHOLD, b = B_STATE,
                          order_use = order_idx) {
  Q <- matrix(NA_real_, n_seg, n_mon)
  V <- matrix(NA_real_, n_seg, n_mon)
  Q[seed_rows, ] <- Q0[seed_rows, ]
  V[seed_rows, ] <- V0[seed_rows, ]
  is_seed <- logical(n_seg); is_seed[seed_rows] <- TRUE
  has_val <- is_seed
  
  for (j in order_use) {
    if (is_seed[j]) next
    preds <- pred_list[[j]]
    if (length(preds) == 0L) next
    preds_val <- preds[has_val[preds]]
    n_val     <- length(preds_val)
    if (n_val == 0L) next
    
    if (n_val == 1L) {
      i  <- preds_val
      ai <- DA_vec[i]; aj <- DA_vec[j]
      gate_ok <- is.finite(ai) && is.finite(aj) && ai > 0 && aj > 0 &&
        abs(ai - aj) / max(aj, 1e-6) < threshold
      if (!gate_ok && !isTRUE(force_vec[j])) next
      sf <- if (gate_ok) (aj / ai)^b else 1
      Q[j, ] <- Q[i, ] * sf + WWm_mat[j, ]
      V[j, ] <- V[i, ] * sf^2 + WWv_mat[j, ]
    } else {
      fin   <- is.finite(Q[preds_val, , drop = FALSE])
      n_fin <- colSums(fin)
      Qj <- rep(NA_real_, n_mon); Vj <- rep(NA_real_, n_mon)
      m2 <- which(n_fin >= 2L)
      if (length(m2) > 0L) {
        Qj[m2] <- colSums(Q[preds_val, m2, drop = FALSE], na.rm = TRUE) + WWm_mat[j, m2]
        Vj[m2] <- colSums(V[preds_val, m2, drop = FALSE], na.rm = TRUE) + WWv_mat[j, m2]
      }
      for (m in which(n_fin == 1L)) {
        i  <- preds_val[fin[, m]]
        ai <- DA_vec[i]; aj <- DA_vec[j]
        gate_ok <- is.finite(ai) && is.finite(aj) && ai > 0 && aj > 0 &&
          abs(ai - aj) / max(aj, 1e-6) < threshold
        if (!gate_ok && !isTRUE(force_vec[j])) next
        sf <- if (gate_ok) (aj / ai)^b else 1
        Qj[m] <- Q[i, m] * sf + WWm_mat[j, m]
        Vj[m] <- V[i, m] * sf^2 + WWv_mat[j, m]
      }
      if (!any(is.finite(Qj))) next
      Q[j, ] <- Qj; V[j, ] <- Vj
    }
    has_val[j] <- TRUE
  }
  list(Q = Q, V = V)
}

# %WW exactly as app.R's compute_pct (also works on matrices)
compute_pct <- function(ww, denom) {
  out <- 100 * ww / denom
  bad <- is.na(denom) | is.na(ww) | denom <= 0
  out[bad] <- NA_real_
  out[!bad & out > 100] <- 100
  out[!bad & out < 0]   <- 0
  out
}

# Upstream watershed of a segment (including itself) - used to speed up LOOCV
get_ancestors <- function(target) {
  buf <- integer(n_seg); buf[1L] <- target; n <- 1L; h <- 1L
  while (h <= n) {
    p <- pred_list[[buf[h]]]; h <- h + 1L
    k <- length(p)
    if (k > 0L) { buf[(n + 1L):(n + k)] <- p; n <- n + k }
  }
  buf[seq_len(n)]
}

# Nearest upstream gages of gage d: walk upstream but stop at the first real gage on each branch
nearest_upstream_gages <- function(d) {
  p <- pred_list[[d]]; n <- length(p)
  if (n == 0L) return(integer(0))
  buf <- integer(n_seg); buf[seq_len(n)] <- p
  h <- 1L; out <- integer(0)
  while (h <= n) {
    node <- buf[h]; h <- h + 1L
    if (is_gage[node]) {
      out <- c(out, node)
    } else {
      p <- pred_list[[node]]; k <- length(p)
      if (k > 0L) { buf[(n + 1L):(n + k)] <- p; n <- n + k }
    }
  }
  out
}

# ==============================================================================
# 3. STATEWIDE EXPONENT
# ==============================================================================
fit_b <- function(rows) {
  q  <- ann_mean_flow[rows]; a <- DA_vec[rows]
  ok <- is.finite(q) & is.finite(a) & q > 0 & a > 0
  m  <- lm(log(q[ok]) ~ log(a[ok]))
  list(model = m, b = unname(coef(m)[2]), n = sum(ok), r2 = summary(m)$r.squared)
}
full_fit <- fit_b(gage_rows)
B_STATE  <- full_fit$b
B_CI     <- unname(confint(full_fit$model)[2, ])

# Alternative fit: drop partial-year gages
alt_rows <- gage_rows[n_months_obs[gage_rows] >= MIN_MONTHS_ALT_FIT]
alt_fit  <- fit_b(alt_rows)
alt_ci   <- unname(confint(alt_fit$model)[2, ])

# ==============================================================================
# 4. ALIGNMENT CHECKPOINT vs app.R
# ==============================================================================
message(sprintf("DA-ratio exponent: b = %.3f (n = %d gaged segments, R^2 = %.3f)",
                B_STATE, full_fit$n, full_fit$r2))

Q_adopted   <- propagate_dar(gage_rows)$Q
pct_adopted <- compute_pct(WWcum_mat, Q_adopted)
n_direct     <- sum(!is.na(flow_order$site_no))
n_propagated <- sum(is.na(flow_order$site_no) & rowSums(is.finite(Q_adopted)) > 0)
n_nhd_only   <- n_seg - n_direct - n_propagated

cat("\n=== ALIGNMENT CHECKPOINT (compare with the app after it runs) ===\n")
cat(sprintf("  flow_source:  gage_direct = %d | gage_propagated = %d | nhd_only = %d\n",
            n_direct, n_propagated, n_nhd_only))
cat(sprintf("  median %%WWef, %s, adopted settings = %.3f   (app: median(flowlines$pct_mean_%s, na.rm = TRUE))\n",
            hl_name, median(pct_adopted[, hl], na.rm = TRUE), HIGHLIGHT_MONTH))
cat(sprintf("  exponent: %.3f (95%% CI %.3f-%.3f), n = %d, R^2 = %.3f\n",
            B_STATE, B_CI[1], B_CI[2], full_fit$n, full_fit$r2))
cat(sprintf("  alt fit (gages with >= %d months): b = %.3f (95%% CI %.3f-%.3f), n = %d\n",
            MIN_MONTHS_ALT_FIT, alt_fit$b, alt_ci[1], alt_ci[2], alt_fit$n))
cat("  reference values from the earlier validated run: b = 0.719, n = 396; June median %WWef (percent, all segments with an estimate) about 0.14\n")
cat("  (June median can shift slightly if fill_facs is used or the DMR/gage inputs changed)\n")

# ==============================================================================
# 5. LEAVE-ONE-GAGE-OUT VALIDATION  (adopted method)
# ==============================================================================
# Hide each gage, refit the exponent without it, re-propagate from the remaining
# gages upstream of it, and compare the reconstructed monthly flow with what the
# gage actually recorded. Only the held-out gage's upstream watershed is
# processed (nothing else can affect its value).
summarise_loocv <- function(df) {
  has_obs <- is.finite(df$observed) & df$observed > 0
  rec     <- has_obs & is.finite(df$predicted)
  o <- df$observed[rec]; p <- df$predicted[rec]
  pos <- p > 0
  tibble(
    n_gage_months     = sum(has_obs),
    n_reconstructed   = sum(rec),
    pct_reconstructed = round(100 * sum(rec) / sum(has_obs), 1),
    median_abs_pct_error = round(median(abs(p - o) / o * 100), 1),
    median_bias_pct   = round((exp(median(log(p[pos] / o[pos]))) - 1) * 100, 1),   # + = over-prediction
    rmse_cfs          = round(sqrt(mean((p - o)^2)), 2),
    pearson_r         = round(suppressWarnings(cor(p, o)), 3),
    pearson_r_log     = round(suppressWarnings(cor(log(p[pos]), log(o[pos]))), 3)
  )
}

if (RUN_LOOCV) {
  message("Running leave-one-gage-out validation (", length(gage_rows), " gages)...")
  loocv_list <- vector("list", length(gage_rows))
  for (k in seq_along(gage_rows)) {
    hold  <- gage_rows[k]
    anc   <- get_ancestors(hold)
    seeds <- intersect(gage_rows, anc); seeds <- seeds[seeds != hold]
    ord_k <- anc[order(ord_pos[anc])]
    b_k   <- fit_b(gage_rows[gage_rows != hold])$b
    Qk    <- propagate_dar(seeds, SIMILARITY_THRESHOLD, b_k, ord_k)$Q
    loocv_list[[k]] <- tibble(Hydroseq = flow_order$Hydroseq[hold],
                              site_no  = flow_order$site_no[hold],
                              month    = mon_codes,
                              observed = Qraw[hold, ],
                              predicted = Qk[hold, ])
    if (k %% 50 == 0) message(sprintf("  LOOCV: %d / %d gages", k, length(gage_rows)))
  }
  loocv <- bind_rows(loocv_list)
  write_csv(loocv, file.path(out_dir, "tbl_loocv_predictions.csv"))
  
  loocv_overall  <- summarise_loocv(loocv)
  loocv_by_month <- loocv %>% group_by(month) %>% group_modify(~ summarise_loocv(.x)) %>% ungroup()
  write_csv(loocv_overall,  file.path(out_dir, "tbl_loocv_summary_overall.csv"))
  write_csv(loocv_by_month, file.path(out_dir, "tbl_loocv_summary_by_month.csv"))
  
  cat("\n=== LEAVE-ONE-GAGE-OUT VALIDATION (all months pooled) ===\n")
  print(as.data.frame(loocv_overall), row.names = FALSE)
  cat("  reference from the earlier validated run: median_abs_pct_error = 23.8\n")
  cat("\nBy month:\n"); print(as.data.frame(loocv_by_month), row.names = FALSE)
  
  # ---- Figure 1: LOOCV accuracy. A = every held-out gage-month vs the 1:1 line; B = error by month ----
  sc <- loocv %>%
    filter(is.finite(observed), is.finite(predicted), observed > 0, predicted > 0) %>%
    mutate(season = season_of(month))
  lims <- range(c(sc$observed, sc$predicted))
  
  f1a <- ggplot(sc, aes(observed, predicted, colour = season)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey30", linetype = "dashed", linewidth = 0.5) +
    geom_point(alpha = 0.45, size = 1.4) +
    annotate("text", x = lims[1], y = lims[2], hjust = 0, vjust = 1.3, size = 3.8, fontface = "italic",
             label = sprintf("Median error = %.1f%%\nn = %s gage-months",
                             loocv_overall$median_abs_pct_error, format(loocv_overall$n_reconstructed, big.mark = ","))) +
    scale_x_log10(limits = lims, labels = fmt_num) +
    scale_y_log10(limits = lims, labels = fmt_num) +
    scale_colour_manual(values = SEASON_COLS, drop = FALSE) +
    guides(colour = "none") +
    coord_equal() +
    labs(title = "Observed vs. reconstructed flow",
         x = "Observed gage flow (cfs)", y = "Reconstructed flow (cfs)") +
    theme_resource()
  
  mb <- loocv_by_month %>%
    mutate(month_f = factor(month, levels = mon_codes, labels = month.abb), season = season_of(month))
  f1b <- ggplot(mb, aes(month_f, median_abs_pct_error, fill = season)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%.0f%%", median_abs_pct_error)), vjust = -0.5, size = 3.4, fontface = "bold") +
    scale_fill_manual(values = SEASON_COLS, drop = FALSE) +
    scale_y_continuous(labels = pct_lab, expand = expansion(mult = c(0, 0.12))) +
    labs(title = "Median error by month", x = NULL, y = "Median absolute error") +
    theme_resource()
  
  fig1 <- combine_panels(list(f1a, f1b),
                         title = "Leave-one-gage-out check: can the method re-create a gage it has not seen?",
                         subtitle = "Each held-out gage is estimated only from other gages upstream of it.",
                         widths = c(1, 1.25), collect = TRUE)
  save_plot(fig1, "fig1_loocv_observed_vs_reconstructed.png", w = 12, h = 5.6)
}

# ==============================================================================
# 6. SENSITIVITY SWEEPS - shared helpers
# ==============================================================================
# One call = one full-network run with ALL real gages in place. Summaries are
# computed for every month. (Coverage = share of ALL segments carrying a flow estimate.)
# All sensitivity summaries are computed on WASTEWATER-IMPACTED segments only:
# segments with cumulative upstream effluent > 0 in that month. This set depends only
# on the WWTP data (not on b or the gate), so the SAME segments are compared in every
# scenario. Segments with no upstream WWTP are 0% %WWef under any setting.
#   pct_filled           = % of impacted segments that receive a flow estimate (coverage)
#   median_pct, p25...   = %WWef among impacted segments that have an estimate
#   median_abs_change_pp = per-segment |change| in %WWef vs the adopted run, in percentage
#                          points (segments with an estimate in both runs)
#   share_change_gt1pp   = % of those segments whose %WWef moves by more than 1 point
#   share_ge_dom         = % of impacted segments (with an estimate) at or above DOM_PCT %WWef
#   *_change_dom_pp      = same per-segment change, restricted to segments that are
#                          effluent-dominated (>= DOM_PCT) in the adopted run
# NOTE: the median %WWef of impacted segments sits in a gap between two groups (large rivers
# at a fraction of a percent, small effluent-dominated streams at 20%+), so it barely moves.
# The dominated-segment metrics are the informative ones.
ww_pos  <- !is.na(WWcum_mat) & WWcum_mat > 0
REF_PCT <- compute_pct(WWcum_mat, propagate_dar(gage_rows)$Q)     # adopted run = reference
cat(sprintf("Wastewater-impacted segments in %s: %s of %s\n", hl_name,
            format(sum(ww_pos[, hl]), big.mark = ","), format(n_seg, big.mark = ",")))

summarise_scenario <- function(Q) {
  pct <- compute_pct(WWcum_mat, Q)
  map_dfr(seq_len(n_mon), function(k) {
    imp <- ww_pos[, k]
    p <- pct[imp, k]; p <- p[is.finite(p)]
    dd <- pct[imp, k] - REF_PCT[imp, k]
    d  <- dd[is.finite(dd)]
    dom <- is.finite(REF_PCT[imp, k]) & REF_PCT[imp, k] >= DOM_PCT      # dominated in the ADOPTED run
    d_dom <- dd[dom & is.finite(dd)]
    tibble(month          = mon_codes[k],
           n_impacted     = sum(imp),
           pct_filled     = 100 * mean(is.finite(Q[imp, k])),
           pct_filled_all = 100 * mean(is.finite(Q[, k])),
           median_flow    = median(Q[imp, k], na.rm = TRUE),
           n_pct          = length(p),
           median_pct     = median(p),
           p25 = quantile(p, 0.25, names = FALSE),
           p75 = quantile(p, 0.75, names = FALSE),
           p90 = quantile(p, 0.90, names = FALSE),
           share_ge_dom = if (length(p)) 100 * mean(p >= DOM_PCT) else NA_real_,
           n_dominated  = sum(dom),
           median_abs_change_dom_pp = if (length(d_dom)) median(abs(d_dom)) else NA_real_,
           p90_abs_change_dom_pp    = if (length(d_dom)) quantile(abs(d_dom), 0.90, names = FALSE) else NA_real_,
           median_abs_change_pp = if (length(d)) median(abs(d)) else NA_real_,
           share_change_gt1pp   = if (length(d)) 100 * mean(abs(d) > 1) else NA_real_)
  })
}
sweep_scenarios <- function(values, run_fn) {
  map_dfr(values, function(v) summarise_scenario(run_fn(v)) %>% mutate(value = v, .before = 1))
}

# ==============================================================================
# 7. THRESHOLD SENSITIVITY  (computation; figures are drawn in section 8b)
# ==============================================================================
message("Threshold sweep...")
THRESHOLD_GRID <- sort(unique(c(THRESHOLD_GRID, SIMILARITY_THRESHOLD)))
th_labels <- ifelse(is.infinite(THRESHOLD_GRID), "No gate", sprintf("%g%%", THRESHOLD_GRID * 100))
adopted_th_label <- paste0(th_labels[THRESHOLD_GRID == SIMILARITY_THRESHOLD], "\n(adopted)")
th_labels[THRESHOLD_GRID == SIMILARITY_THRESHOLD] <- adopted_th_label

thr <- sweep_scenarios(THRESHOLD_GRID, function(th) propagate_dar(gage_rows, threshold = th)$Q) %>%
  mutate(threshold = value,
         th_f = factor(th_labels[match(value, THRESHOLD_GRID)], levels = th_labels))
write_csv(thr %>% select(-value, -th_f), file.path(out_dir, "tbl_threshold_sensitivity.csv"))

thr_hl <- thr %>% filter(month == HIGHLIGHT_MONTH)
cat(sprintf("\n=== THRESHOLD SENSITIVITY (%s) ===\n", hl_name))
print(as.data.frame(thr_hl %>% transmute(threshold, pct_filled = round(pct_filled, 1),
                                         median_pct = round(median_pct, 2), p25 = round(p25, 2),
                                         p75 = round(p75, 2), p90 = round(p90, 2),
                                         share_ge_dom = round(share_ge_dom, 1),
                                         chg_dom_med = round(median_abs_change_dom_pp, 2),
                                         chg_dom_p90 = round(p90_abs_change_dom_pp, 2))), row.names = FALSE)

# ==============================================================================
# 8. EXPONENT SENSITIVITY  (95% CI of the fitted b, plus the partial-year alternative)
# ==============================================================================
message("Exponent sweep...")
B_GRID <- sort(unique(c(seq(B_CI[1], B_CI[2], length.out = N_EXPONENT_STEPS), B_STATE)))
expo <- sweep_scenarios(B_GRID, function(b) propagate_dar(gage_rows, b = b)$Q) %>% mutate(b = value)
write_csv(expo %>% select(-value), file.path(out_dir, "tbl_exponent_sensitivity.csv"))

alt_row <- summarise_scenario(propagate_dar(gage_rows, b = alt_fit$b)$Q) %>%
  filter(month == HIGHLIGHT_MONTH) %>% mutate(b = alt_fit$b)
write_csv(alt_row, file.path(out_dir, "tbl_exponent_alt_fit_partial_year_excluded.csv"))

expo_hl <- expo %>% filter(month == HIGHLIGHT_MONTH)
cat(sprintf("\n=== EXPONENT SENSITIVITY (%s) ===\n", hl_name))
print(as.data.frame(expo_hl %>% transmute(b = round(b, 3), median_pct = round(median_pct, 2),
                                          p25 = round(p25, 2), p75 = round(p75, 2),
                                          share_ge_dom = round(share_ge_dom, 1),
                                          chg_dom_med = round(median_abs_change_dom_pp, 2),
                                          chg_dom_p90 = round(p90_abs_change_dom_pp, 2))), row.names = FALSE)
cat(sprintf("  alt fit (>= %d months), b = %.3f: median %%WWef = %.2f\n",
            MIN_MONTHS_ALT_FIT, alt_fit$b, alt_row$median_pct))

# ------------------------------------------------------------------------------
# 8b. FIGURES 2 and 3 - threshold and exponent (same y-axis so the two are comparable)
# ------------------------------------------------------------------------------
row_at_b <- function(b) expo_hl[which.min(abs(expo_hl$b - b)), ]
SENS_YMAX <- 100

# Figure 2: A = how the answer changes with the gate; B = how much of the network gets an estimate
thr_dots <- thr_hl %>% transmute(x = th_f, adopted = threshold == SIMILARITY_THRESHOLD,
                                 dom = share_ge_dom, filled = pct_filled, chg = median_abs_change_dom_pp)
th_xlab <- "Similarity threshold\n(max relative change in drainage area per hop)"
dom_lab <- sprintf("Impacted segments with %%WWef \u2265 %g%%", DOM_PCT)
pt_axis <- scales::label_number(drop0trailing = TRUE)
f2a <- dot_panel(thr_dots %>% transmute(x, y = dom, adopted), th_xlab, dom_lab, ymax = 100, digits = 1) +
  labs(title = "Result vs. threshold")
f2b <- dot_panel(thr_dots %>% transmute(x, y = filled, adopted), th_xlab,
                 "Impacted segments with a flow estimate", ymax = 100, digits = 0) +
  labs(title = "Coverage vs. threshold")
f2c <- dot_panel(thr_dots %>% transmute(x, y = chg, adopted), th_xlab,
                 sprintf("Median change in %%WWef, dominated segments (points)"), digits = 2, suffix = " pt",
                 axis_lab = pt_axis) +
  labs(title = "Change vs. adopted threshold")
fig2 <- combine_panels(list(f2a, f2b, f2c),
                       title = "Sensitivity to the drainage-area similarity threshold",
                       subtitle = sprintf("%s, wastewater-impacted segments. 'Dominated' = %%WWef \u2265 %g%% at the adopted threshold. Filled dot = adopted %g%%.",
                                          hl_name, DOM_PCT, SIMILARITY_THRESHOLD * 100))
save_plot(fig2, "fig2_threshold_sensitivity.png", w = 16, h = 5.6)

# Figure 3: median %WWef under four exponent values, ordered by b
ex_d <- tibble(
  scenario = c("Lower 95% CI", sprintf("Partial-year gages\nexcluded (<%d months)", MIN_MONTHS_ALT_FIT),
               "Fitted\n(adopted)", "Upper 95% CI"),
  b        = c(B_CI[1], alt_fit$b, B_STATE, B_CI[2]),
  y        = c(row_at_b(B_CI[1])$share_ge_dom, alt_row$share_ge_dom, row_at_b(B_STATE)$share_ge_dom, row_at_b(B_CI[2])$share_ge_dom),
  chg      = c(row_at_b(B_CI[1])$median_abs_change_dom_pp, alt_row$median_abs_change_dom_pp,
               row_at_b(B_STATE)$median_abs_change_dom_pp, row_at_b(B_CI[2])$median_abs_change_dom_pp),
  chg90    = c(row_at_b(B_CI[1])$p90_abs_change_dom_pp, alt_row$p90_abs_change_dom_pp,
               row_at_b(B_STATE)$p90_abs_change_dom_pp, row_at_b(B_CI[2])$p90_abs_change_dom_pp),
  adopted  = c(FALSE, FALSE, TRUE, FALSE)
) %>%
  arrange(b) %>%
  mutate(x_chr = sprintf("%s\nb = %.3f", scenario, b),
         x     = factor(x_chr, levels = x_chr))
write_csv(ex_d %>% select(scenario, b, share_ge_dom = y, median_abs_change_dom_pp = chg, p90_abs_change_dom_pp = chg90, adopted) %>% mutate(scenario = gsub("\n", " ", scenario)),
          file.path(out_dir, "tbl_exponent_scenarios_summary.csv"))
ex_xlab <- "Drainage-area exponent scenario"
f3a <- dot_panel(ex_d %>% transmute(x, y, adopted), ex_xlab, dom_lab, ymax = 100, digits = 1) +
  labs(title = "Result vs. exponent")
f3b <- dot_panel(ex_d %>% transmute(x, y = chg, adopted), ex_xlab,
                 "Median change in %WWef, dominated segments (points)", digits = 2, suffix = " pt",
                 axis_lab = pt_axis) +
  labs(title = "Change vs. fitted exponent")
fig3 <- combine_panels(list(f3a, f3b),
                       title = "Sensitivity to the drainage-area exponent",
                       subtitle = sprintf("%s, wastewater-impacted segments. Share of impacted segments at or above %g%% %%WWef: %.1f%% to %.1f%% across all exponents. Filled dot = adopted.",
                                          hl_name, DOM_PCT, min(ex_d$y), max(ex_d$y)))
save_plot(fig3, "fig3_exponent_sensitivity.png", w = 13, h = 5.6)

# ------------------------------------------------------------------------------
# Supplement S1: every month x every threshold (table-like heatmap)
# ------------------------------------------------------------------------------
hm <- thr %>%
  mutate(month_f = factor(month, levels = rev(mon_codes), labels = rev(month.abb)))
hm_mid <- (min(hm$share_ge_dom, na.rm = TRUE) + max(hm$share_ge_dom, na.rm = TRUE)) / 2
adopt_idx <- match(adopted_th_label, levels(hm$th_f))
figS1 <- ggplot(hm, aes(th_f, month_f, fill = share_ge_dom)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f", share_ge_dom), colour = share_ge_dom > hm_mid), size = 3.4) +
  annotate("rect", xmin = adopt_idx - 0.5, xmax = adopt_idx + 0.5, ymin = 0.5, ymax = n_mon + 0.5,
           fill = NA, colour = COL$ink, linewidth = 0.9) +
  scale_fill_gradient(low = COL$heat_lo, high = COL$heat_hi, name = sprintf("%% of impacted\nsegments with\n%%WWef \u2265 %g%%", DOM_PCT)) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = COL$ink), guide = "none") +
  labs(title = sprintf("Effluent-dominated segments (%%WWef \u2265 %g%%) by month and similarity threshold", DOM_PCT),
       subtitle = "Wastewater-impacted segments only. Numbers are % of those segments. Outlined column = adopted threshold.",
       x = th_xlab, y = NULL) +
  theme_resource() +
  theme(panel.grid = element_blank(), plot.subtitle = element_text(hjust = 0.5, colour = "grey30"))
save_plot(figS1, "figS1_threshold_by_month_heatmap.png", w = 9, h = 6.4)

# ==============================================================================
# 9. IS THE EXPONENT SEASONAL?  (monthly fits + nested-model F-test)
# ==============================================================================
message("Seasonal exponent check...")
monthly_fits <- map_dfr(seq_len(n_mon), function(k) {
  q  <- Qraw[gage_rows, k]; a <- DA_vec[gage_rows]
  ok <- is.finite(q) & is.finite(a) & q > 0 & a > 0
  m  <- lm(log(q[ok]) ~ log(a[ok])); ci <- confint(m)[2, ]
  tibble(month = mon_codes[k], n_gages = sum(ok), b = unname(coef(m)[2]),
         ci_lo = unname(ci[1]), ci_hi = unname(ci[2]), r2 = summary(m)$r.squared)
})
write_csv(monthly_fits, file.path(out_dir, "tbl_monthly_exponent_fits.csv"))

long_data <- map_dfr(seq_len(n_mon), function(k) {
  q  <- Qraw[gage_rows, k]; a <- DA_vec[gage_rows]
  ok <- is.finite(q) & is.finite(a) & q > 0 & a > 0
  tibble(Q = q[ok], A = a[ok], month = mon_codes[k])
})
ft <- anova(lm(log(Q) ~ log(A) + month, data = long_data),
            lm(log(Q) ~ log(A) * month, data = long_data))
seasonal_test <- tibble(df1 = ft$Df[2], df2 = ft$Res.Df[2], F_stat = round(ft$F[2], 3),
                        p_value = signif(ft[["Pr(>F)"]][2], 4), n_obs = nrow(long_data))
write_csv(seasonal_test, file.path(out_dir, "tbl_seasonal_exponent_test.csv"))
cat("\n=== IS THE EXPONENT SEASONAL? ===\n")
print(as.data.frame(monthly_fits %>% mutate(across(c(b, ci_lo, ci_hi, r2), ~ round(.x, 3)))), row.names = FALSE)
cat(sprintf("  Nested F-test (exponent varies by month?): F(%d, %d) = %.3f, p = %s\n",
            seasonal_test$df1, seasonal_test$df2, seasonal_test$F_stat, format.pval(seasonal_test$p_value, digits = 3)))

mf_plot <- monthly_fits %>%
  mutate(month_f = factor(month, levels = mon_codes, labels = month.abb), season = season_of(month))
fig4 <- ggplot(mf_plot, aes(month_f, b)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = B_CI[1], ymax = B_CI[2], fill = COL$band, alpha = 0.6) +
  geom_hline(yintercept = B_STATE, colour = "grey30", linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, colour = season), width = 0.2, linewidth = 0.6) +
  geom_point(aes(colour = season), size = 3.2) +
  geom_text(aes(y = ci_hi, label = sprintf("%.2f", b)), vjust = -0.8, size = 3.3, fontface = "bold") +
  annotate("text", x = 0.6, y = -Inf, hjust = 0, vjust = -0.8, size = 3.6, fontface = "italic",
           label = sprintf("Gray band = annual fit 95%% CI (line: b = %.3f)\nF-test for a month effect: F(%d, %d) = %.2f, p = %s",
                           B_STATE, seasonal_test$df1, seasonal_test$df2, seasonal_test$F_stat,
                           format.pval(seasonal_test$p_value, digits = 3))) +
  scale_colour_manual(values = SEASON_COLS, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.1))) +
  labs(title = "Is the drainage-area exponent seasonal?",
       x = NULL, y = "Fitted exponent, b") +
  theme_resource()
save_plot(fig4, "fig4_seasonal_exponent.png", w = 10, h = 5.4)

# ==============================================================================
# 10. WHAT DOES LOOSENING THE GATE ADD?  (composition of the network)
# ==============================================================================
message("Gate-loosening composition...")
Q_tight <- propagate_dar(gage_rows, threshold = COMP_TIGHT)$Q[, hl]
Q_adopt <- propagate_dar(gage_rows, threshold = SIMILARITY_THRESHOLD)$Q[, hl]
Q_loose <- propagate_dar(gage_rows, threshold = Inf)$Q[, hl]
pct_loose <- compute_pct(WWcum_mat[, hl], Q_loose)       # %WWef at the loosest gate

g1 <- sprintf("Reached at\n%g%% gate", COMP_TIGHT * 100)
g2 <- sprintf("Added from %g%%\nto %g%%", COMP_TIGHT * 100, SIMILARITY_THRESHOLD * 100)
g3 <- sprintf("Added beyond\n%g%% gate", SIMILARITY_THRESHOLD * 100)
grp <- rep(NA_character_, n_seg)
grp[is.finite(Q_loose)] <- g3
grp[is.finite(Q_adopt)] <- g2
grp[is.finite(Q_tight)] <- g1

comp <- tibble(row = which(!is.na(grp))) %>%
  mutate(group = factor(grp[row], levels = c(g1, g2, g3)), DA = DA_vec[row], pct = pct_loose[row])

summ_group <- function(d, label) {
  tibble(group = label, n_segments = nrow(d),
         median_DA_km2 = round(median(d$DA, na.rm = TRUE), 1),
         mean_DA_km2   = round(mean(d$DA, na.rm = TRUE), 1),
         median_pctWWef = round(median(d$pct, na.rm = TRUE), 3),
         mean_pctWWef   = round(mean(d$pct, na.rm = TRUE), 3),
         share_no_upstream_WW_pct = round(100 * mean(d$pct == 0, na.rm = TRUE), 1))
}
comp_tbl <- bind_rows(
  summ_group(comp %>% filter(group == g1), gsub("\n", " ", g1)),
  summ_group(comp %>% filter(group == g2), gsub("\n", " ", g2)),
  summ_group(comp %>% filter(group == g3), gsub("\n", " ", g3)),
  summ_group(comp %>% filter(group != g1), sprintf("Added beyond %g%% gate (both groups)", COMP_TIGHT * 100))
)
write_csv(comp_tbl, file.path(out_dir, "tbl_gate_loosening_composition.csv"))
cat(sprintf("\n=== WHAT LOOSENING THE GATE ADDS (%s) ===\n", hl_name))
print(as.data.frame(comp_tbl), row.names = FALSE)

grp_stats <- comp %>% group_by(group) %>%
  summarise(n = n(), median_DA = median(DA, na.rm = TRUE),
            no_ww = 100 * mean(pct == 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(x_lab = sprintf("%s\n(n = %s)", as.character(group), format(n, big.mark = ",")),
         x_lab = factor(x_lab, levels = x_lab[order(as.integer(group))]))

grp_bar <- function(d, ycol, ylab, suffix = "", digits = 0) {
  ggplot(d, aes(x_lab, .data[[ycol]], fill = group)) +
    geom_col(width = 0.62) +
    geom_text(aes(label = paste0(formatC(.data[[ycol]], format = "f", digits = digits, big.mark = ","), suffix)),
              vjust = -0.5, size = 3.6, fontface = "bold") +
    scale_fill_manual(values = setNames(ORDINAL_BLUE, c(g1, g2, g3)), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = ylab) +
    theme_resource() +
    theme(axis.text.x = element_text(lineheight = 0.95))
}
f5a <- grp_bar(grp_stats, "median_DA", "Median drainage area (km²)", digits = 0) +
  labs(title = "Segments added are smaller")
f5b <- grp_bar(grp_stats, "no_ww", "Segments with no upstream WWTP", suffix = "%", digits = 0) +
  scale_y_continuous(limits = c(0, 100), labels = pct_lab, expand = expansion(mult = c(0, 0.12))) +
  labs(title = "...and mostly have no wastewater upstream")
fig5 <- combine_panels(list(f5a, f5b),
                       title = "What does loosening the similarity gate add to the network?",
                       subtitle = sprintf("%s. Segments grouped by the first gate that reaches them (darker = looser gate).", hl_name))
save_plot(fig5, "fig5_gate_loosening_composition.png", w = 11.5, h = 5.6)

# ==============================================================================
# 11. REAL-GAGE PAIRS: DOES OBSERVED FLOW GROW LIKE (Ad/Au)^b?
# ==============================================================================
# Pairs = each gage with its NEAREST upstream gage on every branch (no other real
# gage in between): the real-data analogue of one propagation hop. Uses observed
# gage flows only; no propagation involved.
#   log_over = log(predicted ratio) - log(observed ratio)   -> positive = DAR over-predicts
#   abs_pct_err = |observed ratio / predicted ratio - 1| x 100
message("Nested real-gage pair checks...")
nested_pairs <- map_dfr(gage_rows, function(d) {
  u <- nearest_upstream_gages(d)
  if (length(u) == 0L) return(NULL)
  tibble(down_row = d, up_row = u)
})
cat(sprintf("\n%d nearest-upstream gage pairs found.\n", nrow(nested_pairs)))

pair_long <- map_dfr(seq_len(n_mon), function(k) {
  Qd <- Qraw[nested_pairs$down_row, k]; Qu <- Qraw[nested_pairs$up_row, k]
  Ad <- DA_vec[nested_pairs$down_row];  Au <- DA_vec[nested_pairs$up_row]
  ok <- is.finite(Qd) & is.finite(Qu) & Qd > 0 & Qu > 0 &
    is.finite(Ad) & is.finite(Au) & Au > 0 & Ad > Au
  tibble(month = mon_codes[k], da_ratio = (Ad / Au)[ok], obs_ratio = (Qd / Qu)[ok])
}) %>%
  mutate(pred_ratio  = da_ratio^B_STATE,
         log_over    = log(pred_ratio) - log(obs_ratio),
         abs_pct_err = abs(obs_ratio / pred_ratio - 1) * 100)

GATE_MAX_RATIO <- 1 / (1 - SIMILARITY_THRESHOLD)   # |Ai-Aj|/Aj < tau  <=>  Aj/Ai < 1/(1-tau)
cat(sprintf("Single-hop area ratio allowed by the %g%% gate: Aj/Ai < %.2f\n",
            SIMILARITY_THRESHOLD * 100, GATE_MAX_RATIO))

# --- 11a. error and bias by size of the area-ratio jump ---------------------------
bin_breaks <- sort(unique(c(1, 1.25, 1.5, GATE_MAX_RATIO, 3, 5, 10, Inf)))
bin_hi     <- bin_breaks[-1]
bin_labels <- map_chr(seq_along(bin_hi), function(i) {
  if (is.infinite(bin_hi[i])) sprintf(">%.2f", bin_breaks[i]) else sprintf("%.2f-%.2f", bin_breaks[i], bin_hi[i])
})
bin_summary <- pair_long %>%
  mutate(bin = cut(da_ratio, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE)) %>%
  filter(!is.na(bin)) %>%
  group_by(bin) %>%
  summarise(n = n(),
            median_abs_pct_err = round(median(abs_pct_err), 1),
            bias_pct = round((exp(median(log_over)) - 1) * 100, 1),
            median_log_ratio_bias_obs_minus_pred = round(-median(log_over), 3), .groups = "drop") %>%
  mutate(within_gate = bin_hi[as.integer(bin)] <= GATE_MAX_RATIO + 1e-9)
write_csv(bin_summary, file.path(out_dir, "tbl_pairs_by_area_ratio_bin.csv"))
cat("\n=== REAL-GAGE PAIRS BY AREA-RATIO JUMP (all months pooled) ===\n")
print(as.data.frame(bin_summary), row.names = FALSE)

bin_x <- bin_summary %>%
  mutate(x_lab = sprintf("%s\n(n = %s)", as.character(bin), format(n, big.mark = ",")),
         x_lab = factor(x_lab, levels = x_lab))
gate_fill <- c("TRUE" = COL$accent, "FALSE" = COL$context)
gate_lab  <- c("TRUE" = sprintf("Within the %g%% gate (ratio ≤ %.2f)", SIMILARITY_THRESHOLD * 100, GATE_MAX_RATIO),
               "FALSE" = "Beyond the gate")
pair_bar <- function(ycol, ylab, title) {
  ggplot(bin_x, aes(x_lab, .data[[ycol]], fill = within_gate)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%.0f%%", .data[[ycol]]), vjust = ifelse(.data[[ycol]] >= 0, -0.5, 1.4)),
              size = 3.5, fontface = "bold") +
    scale_fill_manual(values = gate_fill, labels = gate_lab, drop = FALSE) +
    scale_y_continuous(labels = pct_lab, expand = expansion(mult = c(0.12, 0.12))) +
    labs(title = title, x = "Downstream / upstream drainage-area ratio", y = ylab) +
    theme_resource() +
    theme(axis.text.x = element_text(size = 9, lineheight = 0.95))
}
f6a <- pair_bar("median_abs_pct_err", "Median absolute error", "Error")
f6b <- pair_bar("bias_pct", "Median bias (+ = flow over-predicted)", "Bias")
fig6 <- combine_panels(list(f6a, f6b),
                       title = "Drainage-area scaling vs. real paired gages",
                       subtitle = "Each gage compared with its nearest upstream gage (observed flows only; all months pooled).",
                       collect = TRUE)
save_plot(fig6, "fig6_pairs_error_by_area_ratio.png", w = 12.5, h = 5.6)

# --- 11b. within the gate's own range: is the small bias seasonal? -----------------
boot_median_ci <- function(x, n_boot = N_BOOT) {
  if (length(x) < 5) return(c(NA_real_, NA_real_))
  meds <- replicate(n_boot, median(sample(x, replace = TRUE)))
  unname(quantile(meds, c(0.025, 0.975)))
}
gate_month <- pair_long %>%
  filter(da_ratio <= GATE_MAX_RATIO) %>%
  group_by(month) %>%
  summarise(n = n(), med = median(log_over),
            median_abs_pct_err = round(median(abs_pct_err), 1),
            ci = list(boot_median_ci(log_over)), .groups = "drop") %>%
  mutate(bias_pct = (exp(med) - 1) * 100,
         lo_pct   = (exp(map_dbl(ci, 1)) - 1) * 100,
         hi_pct   = (exp(map_dbl(ci, 2)) - 1) * 100,
         median_log_ratio_bias_obs_minus_pred = -med) %>%
  select(-ci, -med)
write_csv(gate_month %>% mutate(across(c(bias_pct, lo_pct, hi_pct), ~ round(.x, 2))),
          file.path(out_dir, "tbl_pairs_within_gate_bias_by_month.csv"))
cat(sprintf("\n=== BIAS WITHIN THE GATE'S OWN RANGE (Aj/Ai <= %.2f), by month ===\n", GATE_MAX_RATIO))
print(as.data.frame(gate_month %>% mutate(across(c(bias_pct, lo_pct, hi_pct,
                                                   median_log_ratio_bias_obs_minus_pred), ~ round(.x, 3)))),
      row.names = FALSE)

# the finer bin x month detail, as a table only
gate_bin_month <- pair_long %>%
  filter(da_ratio <= GATE_MAX_RATIO) %>%
  mutate(bin = cut(da_ratio, breaks = sort(unique(c(1, 1.25, 1.5, GATE_MAX_RATIO))), include.lowest = TRUE)) %>%
  filter(!is.na(bin)) %>%
  group_by(month, bin) %>%
  summarise(n = n(), bias_pct = round((exp(median(log_over)) - 1) * 100, 1),
            median_abs_pct_err = round(median(abs_pct_err), 1), .groups = "drop")
write_csv(gate_bin_month, file.path(out_dir, "tbl_pairs_within_gate_bin_by_month.csv"))

# --- Figure 7: seasonal bias. A = LOOCV bias by month; B = real-pair bias inside the gate (with CI) ---
gm_plot <- gate_month %>%
  mutate(month_f = factor(month, levels = mon_codes, labels = month.abb), season = season_of(month))
bias_vals <- c(0, gm_plot$bias_pct, gm_plot$lo_pct, gm_plot$hi_pct)
if (RUN_LOOCV) bias_vals <- c(bias_vals, loocv_by_month$median_bias_pct)
bias_lim <- range(bias_vals, na.rm = TRUE); bias_lim <- bias_lim + c(-1, 1) * 0.12 * diff(bias_lim)
dir_note <- "Above 0: flow over-predicted\n(%WWef under-estimated)"

f7b <- ggplot(gm_plot, aes(month_f, bias_pct, colour = season)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lo_pct, ymax = hi_pct), width = 0.2, linewidth = 0.6) +
  geom_point(size = 3.2) +
  annotate("text", x = 0.6, y = bias_lim[2], hjust = 0, vjust = 1.2, size = 3.5, fontface = "italic",
           colour = "grey30", label = dir_note) +
  scale_colour_manual(values = SEASON_COLS, drop = FALSE) +
  scale_y_continuous(limits = bias_lim, labels = pct_lab) +
  labs(title = "Real gage pairs (inside the gate)", x = NULL, y = "Median bias (95% CI)") +
  theme_resource()

if (RUN_LOOCV) {
  lb <- loocv_by_month %>%
    mutate(month_f = factor(month, levels = mon_codes, labels = month.abb), season = season_of(month))
  f7a <- ggplot(lb, aes(month_f, median_bias_pct, colour = season)) +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.5) +
    geom_segment(aes(xend = month_f, yend = 0), linewidth = 0.6) +
    geom_point(size = 3.2) +
    annotate("text", x = 0.6, y = bias_lim[2], hjust = 0, vjust = 1.2, size = 3.5, fontface = "italic",
             colour = "grey30", label = dir_note) +
    scale_colour_manual(values = SEASON_COLS, drop = FALSE) +
    scale_y_continuous(limits = bias_lim, labels = pct_lab) +
    labs(title = "Leave-one-gage-out reconstruction", x = NULL, y = "Median bias") +
    theme_resource()
  fig7 <- combine_panels(list(f7a, f7b), title = "Is the bias seasonal?",
                         subtitle = "Same y-axis in both panels.", collect = TRUE)
  save_plot(fig7, "fig7_seasonal_bias.png", w = 12.5, h = 5.4)
} else {
  fig7 <- f7b + labs(title = "Is the bias seasonal? Real gage pairs (inside the gate)")
  save_plot(fig7, "fig7_seasonal_bias.png", w = 8, h = 5.4)
}

# ==============================================================================
# DONE
# ==============================================================================
cat(sprintf("\nAll outputs written to ./%s/\n", out_dir))
cat("  fig1_loocv_observed_vs_reconstructed.png    fig2_threshold_sensitivity.png\n")
cat("  fig3_exponent_sensitivity.png               fig4_seasonal_exponent.png\n")
cat("  fig5_gate_loosening_composition.png         fig6_pairs_error_by_area_ratio.png\n")
cat("  fig7_seasonal_bias.png                      figS1_threshold_by_month_heatmap.png\n")
cat("  tbl_*.csv: the numbers behind every figure (table view)\n")

# ---- (1) Low-flow screen: how many gages were removed? -------------------------
Qcols <- paste0("Qgage_", mon_codes)

pre <- flowlines %>%
  filter(!is.na(site_no)) %>%                                   # gaged segments, BEFORE the screen
  mutate(mean_annual_cfs = rowMeans(across(all_of(Qcols)), na.rm = TRUE),
         n_months        = rowSums(!is.na(across(all_of(Qcols)))),
         max_monthly_cfs = suppressWarnings(
           do.call(pmax, c(unname(as.list(across(all_of(Qcols)))), na.rm = TRUE))),
         has_data        = n_months > 0)
pre$mean_annual_cfs[is.nan(pre$mean_annual_cfs)] <- NA_real_
pre <- pre %>%
  mutate(low_flow_excluded = has_data & max_monthly_cfs < 0.5,   # every observed month < 0.5 cfs
         in_final          = site_no %in% flow_order$site_no[is_gage])

attrition <- tibble(
  step = c("Gages in the gage file",
           "Unique flowline segments (COMIDs) with a gage (one gage kept per COMID)",
           "Gaged segments with monthly flow statistics",
           "Removed: all observed monthly means < 0.5 cfs",
           "Gages used (exponent fit, LOOCV, propagation)"),
  n = c(nrow(gages_sf),
        dplyr::n_distinct(gages_sf$best_COMID, na.rm = TRUE),
        sum(pre$has_data),
        sum(pre$low_flow_excluded),
        length(gage_rows)))
cat("\n=== GAGE ATTRITION (low-flow screen) ===\n")
print(as.data.frame(attrition), row.names = FALSE)
cat(sprintf("  low-flow gages removed: %d of %d with data (%.1f%%)\n",
            sum(pre$low_flow_excluded), sum(pre$has_data),
            100 * sum(pre$low_flow_excluded) / sum(pre$has_data)))
if (sum(pre$has_data) - sum(pre$low_flow_excluded) != length(gage_rows))
  cat(sprintf("  NOTE: %d gaged segments with data are neither removed by the screen nor in the final set ",
              sum(pre$has_data) - sum(pre$low_flow_excluded) - length(gage_rows)),
      "(e.g., no finite annual mean) - inspect `pre` for in_final == FALSE.\n", sep = "")
write_csv(attrition, file.path(out_dir, "tbl_gage_attrition.csv"))
write_csv(pre %>% filter(low_flow_excluded) %>%
            select(site_no, DA_km2 = DA, mean_annual_cfs, max_monthly_cfs, n_months),
          file.path(out_dir, "tbl_gages_removed_low_flow.csv"))

# ---- (2) Exponent figure: how b was determined ---------------------------------
# Fit: ln(annual mean flow) = a + b * ln(drainage area), ordinary least squares,
# one point per gage. (Fit uses full_fit from the script - same gages and model.)
fit_dat <- tibble(site_no = flow_order$site_no[gage_rows],
                  A       = DA_vec[gage_rows],
                  Q       = ann_mean_flow[gage_rows]) %>%
  filter(is.finite(A), is.finite(Q), A > 0, Q > 0)

cf  <- coef(full_fit$model); V <- vcov(full_fit$model)
tq  <- qt(0.975, df = full_fit$model$df.residual)
xg  <- seq(log(min(fit_dat$A)), log(max(fit_dat$A)), length.out = 200)
band <- tibble(A  = exp(xg),
               fit = exp(cf[1] + cf[2] * xg),
               se  = sqrt(V[1, 1] + 2 * xg * V[1, 2] + xg^2 * V[2, 2]),
               lo  = exp(cf[1] + cf[2] * xg - tq * se),
               hi  = exp(cf[1] + cf[2] * xg + tq * se))

removed <- pre %>% filter(low_flow_excluded, is.finite(DA), DA > 0, is.finite(mean_annual_cfs), mean_annual_cfs > 0)

p_scatter <- ggplot() +
  geom_ribbon(data = band, aes(A, ymin = lo, ymax = hi), fill = COL$band, alpha = 0.7) +
  geom_point(data = removed, aes(DA, mean_annual_cfs), shape = 4, size = 1.8, colour = COL$muted) +
  geom_point(data = fit_dat, aes(A, Q), colour = COL$accent, alpha = 0.55, size = 1.8) +
  geom_line(data = band, aes(A, fit), colour = COL$ink, linewidth = 0.8) +
  annotate("text", x = min(fit_dat$A), y = max(fit_dat$Q), hjust = 0, vjust = 1, size = 3.9,
           fontface = "italic",
           label = sprintf("b = %.3f (95%% CI %.3f–%.3f)\nR² = %.2f, n = %d gages",
                           B_STATE, B_CI[1], B_CI[2], full_fit$r2, full_fit$n)) +
  scale_x_log10(labels = fmt_num) +
  scale_y_log10(labels = fmt_num) +
  labs(title = "Annual mean flow vs. drainage area",
       x = "Drainage area (km²)", y = "Annual mean gage flow (cfs)") +
  theme_resource()

# Panel B: is b robust? Same regression under three gage selections, with b = 1 for reference
b_row <- function(label, rows) {
  q <- ann_mean_flow[rows]; a <- DA_vec[rows]
  ok <- is.finite(q) & is.finite(a) & q > 0 & a > 0
  m  <- lm(log(q[ok]) ~ log(a[ok])); ci <- confint(m)[2, ]
  tibble(label = label, b = unname(coef(m)[2]), lo = ci[1], hi = ci[2], n = sum(ok))
}
no_screen <- {
  q <- pre$mean_annual_cfs[pre$has_data]; a <- pre$DA[pre$has_data]
  ok <- is.finite(q) & is.finite(a) & q > 0 & a > 0
  m  <- lm(log(q[ok]) ~ log(a[ok])); ci <- confint(m)[2, ]
  tibble(label = "Low-flow gages\nkept (no screen)", b = unname(coef(m)[2]), lo = ci[1], hi = ci[2], n = sum(ok))
}
b_tab <- bind_rows(
  b_row("All usable gages\n(adopted)", gage_rows),
  b_row(sprintf("Gages with\n≥%d months", MIN_MONTHS_ALT_FIT), gage_rows[n_months_obs[gage_rows] >= MIN_MONTHS_ALT_FIT]),
  no_screen) %>%
  mutate(x = factor(sprintf("%s\n(n = %d)", label, n), levels = sprintf("%s\n(n = %d)", label, n)),
         adopted = row_number() == 1)
write_csv(b_tab %>% select(label, b, lo, hi, n) %>% mutate(label = gsub("\n", " ", label)),
          file.path(out_dir, "tbl_exponent_fit_robustness.csv"))

p_robust <- ggplot(b_tab, aes(x, b)) +
  geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, colour = COL$accent, linewidth = 0.7) +
  geom_point(aes(fill = adopted, size = adopted), shape = 21, colour = COL$accent, stroke = 0.9) +
  geom_text(aes(y = hi, label = sprintf("%.3f", b)), vjust = -0.9, size = 3.6, fontface = "bold") +
  annotate("text", x = 0.55, y = 1, hjust = 0, vjust = -0.6, size = 3.4, fontface = "italic",
           colour = "grey30", label = "b = 1: flow proportional to area") +   # after the data layers
  scale_fill_manual(values = c("TRUE" = COL$accent, "FALSE" = "white"), guide = "none") +
  scale_size_manual(values = c("TRUE" = 4.4, "FALSE" = 3), guide = "none") +
  scale_y_continuous(limits = c(0.5, 1.1), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Estimated exponent under different gage selections",
       x = NULL, y = "Exponent b (95% CI)") +
  theme_resource()

fig_exp <- combine_panels(list(p_scatter, p_robust),
                          title = "How the drainage-area exponent (b) was determined",
                          subtitle = "Ordinary least squares on log-transformed values, one point per gage. Crosses = gages removed by the < 0.5 cfs screen.",
                          widths = c(1.15, 1))
save_plot(fig_exp, "fig_exponent_fit.png", w = 12.5, h = 5.8)

