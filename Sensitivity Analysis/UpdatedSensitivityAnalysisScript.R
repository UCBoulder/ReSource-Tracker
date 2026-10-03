suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

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
RUN_LOOCV          <- TRUE   # the slowest section; set FALSE to skip while iterating
set.seed(42)

hl      <- match(HIGHLIGHT_MONTH, mon_codes)
hl_name <- month.name[hl]

# --- Plot style (palette validated with the dataviz validator: accent/second
#     pass categorical checks; the ordinal blue ramp passes the ordinal checks) --
COL <- list(accent = "#2a78d6", second = "#eb6834", context = "#c3c2b7",
            band = "#cde2fb", ink = "#0b0b0b", ink2 = "#52514e", muted = "#898781",
            grid = "#e1e0d9", axis = "#c3c2b7")
ORDINAL_BLUE <- c("#86b6ef", "#3987e5", "#184f95")

theme_sens <- function() {
  theme_minimal(base_size = 12) +
    theme(
      text               = element_text(colour = COL$ink),
      plot.title         = element_text(face = "bold", size = 13, margin = margin(b = 4)),
      plot.subtitle      = element_text(size = 10, colour = COL$ink2, margin = margin(b = 10)),
      plot.caption       = element_text(size = 8.5, colour = COL$muted, hjust = 0, margin = margin(t = 10)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      axis.text          = element_text(colour = COL$muted, size = 9.5),
      axis.title         = element_text(colour = COL$ink2, size = 10),
      panel.grid.major.y = element_line(colour = COL$grid, linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.line.x        = element_line(colour = COL$axis, linewidth = 0.4),
      strip.text         = element_text(face = "bold", size = 10.5, hjust = 0),
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.text        = element_text(colour = COL$ink2, size = 9.5),
      plot.margin        = margin(12, 14, 10, 12)
    )
}
save_plot <- function(p, name, w = 9, h = 5) {
  ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 300, bg = "white")
}
fmt_num <- scales::label_number(big.mark = ",", drop0trailing = TRUE)

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
cat("  reference values from the earlier validated run: b = 0.719, n = 396, June median %WWef = 14.4\n")
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
  
  # Fig 1 - observed vs reconstructed
  sc <- loocv %>% filter(is.finite(observed), is.finite(predicted), observed > 0, predicted > 0)
  lims <- range(c(sc$observed, sc$predicted))
  fig1 <- ggplot(sc, aes(observed, predicted)) +
    geom_abline(slope = 1, intercept = 0, colour = COL$ink2, linewidth = 0.4) +
    geom_point(colour = COL$accent, alpha = 0.25, size = 1.3) +
    scale_x_log10(limits = lims, labels = fmt_num) +
    scale_y_log10(limits = lims, labels = fmt_num) +
    coord_equal() +
    labs(title = "Leave-one-gage-out check: reconstructed vs. observed monthly flow",
         subtitle = sprintf("%s gage-months, all 12 months pooled. Median absolute error %.1f%%; 1:1 line shown.",
                            format(loocv_overall$n_reconstructed, big.mark = ","), loocv_overall$median_abs_pct_error),
         x = "Observed gage flow (cfs, log scale)", y = "Reconstructed flow (cfs, log scale)",
         caption = "Each point is one held-out gage in one month; the gage's own record is never used to estimate it.") +
    theme_sens() +
    theme(panel.grid.major.x = element_line(colour = COL$grid, linewidth = 0.3))
  save_plot(fig1, "fig1_loocv_observed_vs_reconstructed.png", w = 6.8, h = 6.4)
  
  # Fig 2 - error and bias by month (two measures -> two panels, never a dual axis)
  lm_long <- bind_rows(
    loocv_by_month %>% transmute(month, panel = "Median absolute error (%)", y = median_abs_pct_error),
    loocv_by_month %>% transmute(month, panel = "Median bias (%; positive = over-prediction)", y = median_bias_pct)
  ) %>%
    mutate(month_f = factor(month, levels = mon_codes, labels = month.abb),
           panel   = factor(panel, levels = unique(panel)))
  fig2 <- ggplot(lm_long, aes(month_f, y, group = 1)) +
    geom_hline(yintercept = 0, colour = COL$axis, linewidth = 0.4) +
    geom_line(colour = COL$accent, linewidth = 0.8) +
    geom_point(shape = 21, fill = COL$accent, colour = "white", stroke = 0.8, size = 3.2) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    labs(title = "Leave-one-gage-out accuracy by month",
         subtitle = "A positive bias means reconstructed flow is too high, so %WWef would be under-estimated.",
         x = NULL, y = NULL) +
    theme_sens()
  save_plot(fig2, "fig2_loocv_error_and_bias_by_month.png", w = 8, h = 6.2)
}

# ==============================================================================
# 6. SENSITIVITY SWEEPS - shared helpers
# ==============================================================================
# One call = one full-network run with ALL real gages in place. Summaries are
# computed for every month so the figures can show whether the result is
# seasonal. (Coverage = share of ALL segments carrying a flow estimate.)
summarise_scenario <- function(Q) {
  pct <- compute_pct(WWcum_mat, Q)
  map_dfr(seq_len(n_mon), function(k) {
    p <- pct[, k]; p <- p[is.finite(p)]
    tibble(month       = mon_codes[k],
           pct_filled  = 100 * mean(is.finite(Q[, k])),
           median_flow = median(Q[, k], na.rm = TRUE),
           n_pct       = length(p),
           median_pct  = median(p),
           p25 = quantile(p, 0.25, names = FALSE),
           p75 = quantile(p, 0.75, names = FALSE),
           p90 = quantile(p, 0.90, names = FALSE))
  })
}
sweep_scenarios <- function(values, run_fn) {
  map_dfr(values, function(v) summarise_scenario(run_fn(v)) %>% mutate(value = v, .before = 1))
}
month_levels <- c(setdiff(mon_codes, HIGHLIGHT_MONTH), HIGHLIGHT_MONTH)   # highlighted month drawn last (on top)
series_cols  <- setNames(c(COL$accent, COL$context), c(hl_name, "Other months"))
series_lwd   <- setNames(c(1.1, 0.45),               c(hl_name, "Other months"))

# ==============================================================================
# 7. THRESHOLD SENSITIVITY
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
                                         p75 = round(p75, 2), p90 = round(p90, 2))), row.names = FALSE)

panel_lv <- c("Segments with a flow estimate (% of network)", "Median %WWef")
thr_long <- bind_rows(
  thr %>% transmute(th_f, month, panel = panel_lv[1], y = pct_filled),
  thr %>% transmute(th_f, month, panel = panel_lv[2], y = median_pct)
) %>%
  mutate(panel   = factor(panel, levels = panel_lv),
         month_f = factor(month, levels = month_levels),
         series  = if_else(month == HIGHLIGHT_MONTH, hl_name, "Other months"))
thr_band  <- thr_hl %>% transmute(th_f, panel = factor(panel_lv[2], levels = panel_lv), ymin = p25, ymax = p75)
thr_adopt <- thr_long %>% filter(month == HIGHLIGHT_MONTH, as.character(th_f) == adopted_th_label)

fig3 <- ggplot(thr_long, aes(th_f, y)) +
  geom_vline(xintercept = match(adopted_th_label, levels(thr_long$th_f)), colour = COL$axis, linewidth = 0.5) +
  geom_ribbon(data = thr_band, aes(x = th_f, ymin = ymin, ymax = ymax, group = 1),
              inherit.aes = FALSE, fill = COL$band) +
  geom_line(aes(group = month_f, colour = series, linewidth = series)) +
  geom_point(data = thr_adopt, shape = 21, fill = COL$accent, colour = "white", stroke = 0.9, size = 3.8) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = series_cols, name = NULL) +
  scale_linewidth_manual(values = series_lwd, name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Sensitivity to the drainage-area similarity threshold",
       subtitle = sprintf("Adopted method (statewide exponent, gated confluences). Highlighted: %s; shaded band: interquartile range of %%WWef in %s.",
                          hl_name, hl_name),
       x = "Similarity threshold, |Ai - Aj| / Aj  (relative change in drainage area per hop)", y = NULL,
       caption = "Each line is one calendar month. Looser gates reconstruct more of the network, mostly small segments with little or no upstream wastewater.") +
  theme_sens()
save_plot(fig3, "fig3_threshold_sensitivity.png", w = 11, h = 5)

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
                                          p25 = round(p25, 2), p75 = round(p75, 2))), row.names = FALSE)
cat(sprintf("  alt fit (>= %d months), b = %.3f: median %%WWef = %.2f\n",
            MIN_MONTHS_ALT_FIT, alt_fit$b, alt_row$median_pct))

expo_plot <- expo %>%
  mutate(month_f = factor(month, levels = month_levels),
         series  = if_else(month == HIGHLIGHT_MONTH, hl_name, "Other months"))
expo_band  <- expo_hl %>% transmute(b, ymin = p25, ymax = p75)
expo_adopt <- expo_hl %>% filter(abs(b - B_STATE) < 1e-12)
alt_label  <- sprintf("Fit on gages with %d+ months\n(b = %.3f)", MIN_MONTHS_ALT_FIT, alt_fit$b)

fig4 <- ggplot(expo_plot, aes(b, median_pct)) +
  geom_ribbon(data = expo_band, aes(x = b, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = COL$band) +
  geom_line(aes(group = month_f, colour = series, linewidth = series)) +
  geom_point(data = expo_adopt, shape = 21, fill = COL$accent, colour = "white", stroke = 0.9, size = 3.8) +
  geom_point(data = alt_row, aes(x = b, y = median_pct), shape = 21, fill = COL$second,
             colour = "white", stroke = 0.9, size = 3.8, inherit.aes = FALSE) +
  geom_text(data = alt_row, aes(x = b, y = median_pct, label = alt_label), inherit.aes = FALSE,
            hjust = 0, vjust = 1.7, nudge_x = 0.002, size = 3.2, colour = COL$ink2, lineheight = 0.95) +
  scale_colour_manual(values = series_cols, name = NULL) +
  scale_linewidth_manual(values = series_lwd, name = NULL) +
  scale_x_continuous(breaks = c(B_CI[1], B_STATE, B_CI[2]),
                     labels = c(sprintf("%.3f\n(95%% CI low)", B_CI[1]),
                                sprintf("%.3f\n(fitted, adopted)", B_STATE),
                                sprintf("%.3f\n(95%% CI high)", B_CI[2]))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Sensitivity of median %WWef to the fitted drainage-area exponent",
       subtitle = sprintf("Exponent varied across its 95%% CI at the adopted %g%% gate. Highlighted: %s; shaded band: interquartile range in %s.",
                          SIMILARITY_THRESHOLD * 100, hl_name, hl_name),
       x = "Drainage-area exponent, b", y = "Median %WWef",
       caption = "Orange point: exponent re-fitted after dropping partial-year gages (all other settings unchanged).") +
  theme_sens()
save_plot(fig4, "fig4_exponent_sensitivity.png", w = 8.5, h = 5.2)

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

mf_plot <- monthly_fits %>% mutate(month_f = factor(month, levels = mon_codes, labels = month.abb))
fig5 <- ggplot(mf_plot, aes(month_f, b)) +
  geom_ribbon(aes(ymin = B_CI[1], ymax = B_CI[2], group = 1), fill = COL$band) +
  geom_hline(yintercept = B_STATE, colour = COL$accent, linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.18, colour = COL$ink2, linewidth = 0.45) +
  geom_point(shape = 21, fill = COL$accent, colour = "white", stroke = 0.8, size = 3.2) +
  annotate("text", x = 0.7, y = B_CI[2], label = sprintf("Annual fit, 95%% CI (b = %.3f)", B_STATE),
           hjust = 0, vjust = -0.7, size = 3.2, colour = COL$ink2) +
  labs(title = "Is the drainage-area exponent seasonal?",
       subtitle = sprintf("Monthly fits (points, 95%% CI bars) against the single annual fit (line and shaded 95%% CI). Nested F-test for a month effect: F(%d, %d) = %.2f, p = %s.",
                          seasonal_test$df1, seasonal_test$df2, seasonal_test$F_stat,
                          format.pval(seasonal_test$p_value, digits = 3)),
       x = NULL, y = "Fitted exponent, b",
       caption = "Each month is fit separately on that month's gage flows vs. drainage area (log-log OLS).") +
  theme_sens()
save_plot(fig5, "fig5_seasonal_exponent.png", w = 9, h = 5.2)

# ==============================================================================
# 10. WHAT DOES LOOSENING THE GATE ADD?  (composition of the network)
# ==============================================================================
message("Gate-loosening composition...")
Q_tight <- propagate_dar(gage_rows, threshold = COMP_TIGHT)$Q[, hl]
Q_adopt <- propagate_dar(gage_rows, threshold = SIMILARITY_THRESHOLD)$Q[, hl]
Q_loose <- propagate_dar(gage_rows, threshold = Inf)$Q[, hl]
pct_loose <- compute_pct(WWcum_mat[, hl], Q_loose)       # %WWef at the loosest gate, as in the earlier analysis

g1 <- sprintf("Reconstructed at\n%g%% gate", COMP_TIGHT * 100)
g2 <- sprintf("Added between %g%%\nand %g%% (adopted)", COMP_TIGHT * 100, SIMILARITY_THRESHOLD * 100)
g3 <- sprintf("Added only beyond %g%%\n(up to no gate)", SIMILARITY_THRESHOLD * 100)
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

lab_tbl <- comp %>% group_by(group) %>%
  summarise(n = n(), no_ww = 100 * mean(pct == 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(x_lab = sprintf("%s\nn = %s; %.0f%% with no\nupstream WWTP", as.character(group), format(n, big.mark = ","), no_ww),
         x_lab = factor(x_lab, levels = x_lab[order(as.integer(group))]))
comp_long <- bind_rows(
  comp %>% transmute(group, panel = "Drainage area (km², log scale)", value = DA),
  comp %>% filter(is.finite(pct), pct > 0) %>%
    transmute(group, panel = "%WWef, segments with upstream WWTP (log scale)", value = pct)
) %>%
  filter(is.finite(value), value > 0) %>%
  left_join(lab_tbl %>% select(group, x_lab), by = "group") %>%
  mutate(panel = factor(panel, levels = unique(panel)))

fig6 <- ggplot(comp_long, aes(x_lab, value, fill = group)) +
  geom_boxplot(width = 0.55, linewidth = 0.4, colour = COL$ink2,
               outlier.size = 0.6, outlier.alpha = 0.2, outlier.colour = COL$muted) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_y_log10(labels = fmt_num) +
  scale_fill_manual(values = setNames(ORDINAL_BLUE, c(g1, g2, g3))) +
  guides(fill = "none") +
  labs(title = "What loosening the similarity gate adds to the reconstructed network",
       subtitle = sprintf("%s. Segments grouped by the first gate that reaches them; darker = looser gate.", hl_name),
       x = NULL, y = NULL,
       caption = "%WWef shown at the loosest (no-gate) run. Segments with no upstream WWTP (0%) cannot be shown on a log axis; their share is in the axis labels.") +
  theme_sens() +
  theme(axis.text.x = element_text(size = 8.5, lineheight = 0.95))
save_plot(fig6, "fig6_gate_loosening_composition.png", w = 11.5, h = 5.8)

# ==============================================================================
# 11. REAL-GAGE PAIRS: DOES OBSERVED FLOW GROW LIKE (Ad/Au)^b?
# ==============================================================================
# Pairs = each gage with its NEAREST upstream gage on every branch (no other real
# gage in between): the real-data analogue of one propagation hop. Uses observed
# gage flows only; no propagation involved.
#   log_over = log(predicted ratio) - log(observed ratio)   -> positive = DAR over-predicts
#   abs_pct_err = |observed ratio / predicted ratio - 1| x 100   (as in the original Part E-G tables)
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

bin_long <- bind_rows(
  bin_summary %>% transmute(bin, n, within_gate, panel = "Median absolute error (%)", y = median_abs_pct_err),
  bin_summary %>% transmute(bin, n, within_gate, panel = "Median bias (%; positive = DAR over-predicts)", y = bias_pct)
) %>%
  mutate(x_lab = factor(sprintf("%s\n(n = %s)", as.character(bin), format(n, big.mark = ",")),
                        levels = unique(sprintf("%s\n(n = %s)", as.character(bin_summary$bin),
                                                format(bin_summary$n, big.mark = ",")))),
         panel = factor(panel, levels = unique(panel)))
fig7 <- ggplot(bin_long, aes(x_lab, y, fill = within_gate)) +
  geom_hline(yintercept = 0, colour = COL$axis, linewidth = 0.4) +
  geom_col(width = 0.62) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = COL$accent, "FALSE" = COL$context),
                    labels = c("TRUE" = "Within the gate's single-hop range", "FALSE" = "Beyond the gate"),
                    name = NULL) +
  labs(title = "Drainage-area scaling vs. real paired gages, by size of the area-ratio jump",
       subtitle = sprintf("Each gage paired with its nearest upstream gage; all months pooled. The %g%% gate allows ratios up to %.2f in one hop.",
                          SIMILARITY_THRESHOLD * 100, GATE_MAX_RATIO),
       x = "Downstream / upstream drainage-area ratio", y = NULL,
       caption = "Observed flow ratio compared with (Ad/Au)^b using the statewide exponent. No propagation involved.") +
  theme_sens() +
  theme(axis.text.x = element_text(size = 8.5, lineheight = 0.95))
save_plot(fig7, "fig7_pairs_error_by_area_ratio.png", w = 11.5, h = 5.2)

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

gm_plot <- gate_month %>% mutate(month_f = factor(month, levels = mon_codes, labels = month.abb))
fig8 <- ggplot(gm_plot, aes(month_f, bias_pct)) +
  geom_hline(yintercept = 0, colour = COL$axis, linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo_pct, ymax = hi_pct), width = 0.18, colour = COL$ink2, linewidth = 0.45) +
  geom_point(shape = 21, fill = COL$accent, colour = "white", stroke = 0.8, size = 3.2) +
  labs(title = "Bias of drainage-area scaling within the gate's own single-hop range, by month",
       subtitle = sprintf("Nearest-upstream gage pairs with Aj/Ai <= %.2f only. Positive = scaled flow too high (%%WWef under-estimated). Bars: bootstrap 95%% CI of the median.",
                          GATE_MAX_RATIO),
       x = NULL, y = "Median bias (%)",
       caption = sprintf("Pairs per month: %d-%d. A bar that excludes zero indicates a bias distinguishable from sampling noise.",
                         min(gate_month$n), max(gate_month$n))) +
  theme_sens()
save_plot(fig8, "fig8_pairs_within_gate_bias_by_month.png", w = 9, h = 5.2)

# ==============================================================================
# DONE
# ==============================================================================
cat(sprintf("\nAll outputs written to ./%s/\n", out_dir))
cat("  fig1_loocv_observed_vs_reconstructed.png    fig2_loocv_error_and_bias_by_month.png\n")
cat("  fig3_threshold_sensitivity.png              fig4_exponent_sensitivity.png\n")
cat("  fig5_seasonal_exponent.png                  fig6_gate_loosening_composition.png\n")
cat("  fig7_pairs_error_by_area_ratio.png          fig8_pairs_within_gate_bias_by_month.png\n")
cat("  tbl_*.csv: the numbers behind every figure (table view)\n")

