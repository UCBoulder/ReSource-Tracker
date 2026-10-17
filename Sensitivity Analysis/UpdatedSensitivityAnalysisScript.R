#Required packages 
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})
if (!requireNamespace("patchwork", quietly = TRUE))
  stop('Package "patchwork" is needed for the A/B panel figures. Install it once with: install.packages("patchwork")', call. = FALSE)

#Configuration
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

#Alignment with main app 
SIMILARITY_THRESHOLD <- 0.50
forced_hseq <- c("550025960", "550020503", "390005386", "390025451", "390005505")
same_seg_ww_adj_seed <- c("550022475", "550015536", "550062530", "390045364", "550025946",
                          "550019184", "760015713", "390005515", "550018746")

# Analysis settings ---
HIGHLIGHT_MONTH    <- "06"
THRESHOLD_GRID     <- c(0.10, 0.25, 0.50, 0.75, Inf)   # Inf = gate removed
N_EXPONENT_STEPS   <- 7      # grid points across the exponent's 95% CI
MIN_MONTHS_ALT_FIT <- 10     # alternative exponent fit: only gages with >= this many months
COMP_TIGHT         <- 0.10   # tight gate for the gate-loosening composition figure
N_BOOT             <- 1000   # bootstrap resamples for the within-gate monthly bias CIs
RUN_LOOCV          <- TRUE   
set.seed(42)

hl      <- match(HIGHLIGHT_MONTH, mon_codes)
hl_name <- month.name[hl]

# Plot theme
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

# A/B panels in a row under one bold title 
combine_panels <- function(plots, title, subtitle = NULL, widths = NULL, collect = FALSE) {
  patchwork::wrap_plots(plots, nrow = 1, widths = widths,
                        guides = if (collect) "collect" else "keep") +
    patchwork::plot_annotation(
      title = title, subtitle = subtitle, tag_levels = "A",
      theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
                    plot.subtitle = element_text(hjust = 0.5, colour = "grey30", size = 10.5),
                    plot.tag      = element_text(face = "bold", size = 14)))
}

# Line + dots with the value printed above every dot
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

#Data prep, follows app.R logic 
message("Loading and preparing network data...")

# --- DMR wastewater: facility-month -> COMID-month mean + variance ------------
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

if (length(fill_facs) > 0) {   
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

# Flowlines
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

# --- Local WW per Hydroseq, then accumulated downstream ---
flow_map <- flowlines %>%
  select(Hydroseq, DnHydroseq, starts_with("WW_mean_"), starts_with("WW_var_")) %>%
  mutate(across(starts_with("WW_"), as.numeric)) %>%
  group_by(Hydroseq) %>%
  summarise(DnHydroseq = first(DnHydroseq),
            across(starts_with("WW_mean_"), ~sum(.x, na.rm = TRUE)),
            across(starts_with("WW_var_"),  ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(desc(as.numeric(Hydroseq)))          # upstream first 

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

# --- Gage monthly stats -> one gage per COMID ---------------------------------
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

# --- flow_order: low-flow screen --------------------------
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
# same-segment seeds .
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

# Upstream-neighbour structure and topological order 
order_idx   <- order(as.numeric(flow_order$Hydroseq), decreasing = TRUE)   # most-upstream first
ord_pos     <- integer(n_seg); ord_pos[order_idx] <- seq_len(n_seg)
valid_edges <- which(!is.na(dn_index))
pred_split  <- split(valid_edges, dn_index[valid_edges])
pred_list   <- vector("list", n_seg)
pred_list[as.integer(names(pred_split))] <- pred_split
pred_list[sapply(pred_list, is.null)] <- list(integer(0))

message("Data prep complete: ", n_seg, " segments, ", length(gage_rows), " usable gages.")

# Propagation
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

# %WW 
compute_pct <- function(ww, denom) {
  out <- 100 * ww / denom
  bad <- is.na(denom) | is.na(ww) | denom <= 0
  out[bad] <- NA_real_
  out[!bad & out > 100] <- 100
  out[!bad & out < 0]   <- 0
  out
}

# Upstream watershed of a segment 
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


# STATEWIDE EXPONENT
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

# ALIGNMENT CHECKPOINT vs full app
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


# LEAVE-ONE-GAGE-OUT VALIDATION
# Hide each gage, refit the exponent without it, re-propagate from the remaining
# gages upstream of it, and compare the reconstructed monthly flow with what the
# gage actually recorded. Only the held-out gage's upstream watershed is
# processed 
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
}

# ---- Figure 1 observed vs. reconstructed flow ----
sc <- loocv %>%
  filter(is.finite(observed), is.finite(predicted), observed > 0, predicted > 0) %>%
  mutate(season = season_of(month))
lims <- range(c(sc$observed, sc$predicted))

fig1 <- ggplot(sc, aes(observed, predicted, colour = season)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey30", linetype = "dashed", linewidth = 0.5) +
  geom_point(alpha = 0.45, size = 1.8) +
  annotate("text", x = lims[1], y = lims[2], hjust = 0, vjust = 1.3, size = 4.2, fontface = "italic",
           label = sprintf("Median error = %.1f%%\nn = %s gage-months",
                           loocv_overall$median_abs_pct_error, format(loocv_overall$n_reconstructed, big.mark = ","))) +
  scale_x_log10(limits = lims, labels = fmt_num) +
  scale_y_log10(limits = lims, labels = fmt_num) +
  scale_colour_manual(values = SEASON_COLS, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  coord_equal() +
  labs(title = "Leave-One-Gage-Out Analysis",
       x = "Observed Gage Flow (cfs)", y = "Reconstructed Flow (cfs)") +
  theme_resource() +
  theme(plot.subtitle = element_text(hjust = 0.5, colour = "grey30", size = 10))
save_plot(fig1, "fig1_loocv_observed_vs_reconstructed.png", w = 7.5, h = 6.5)

# ---- Figure error and bias by month ----
mb <- loocv_by_month %>%
  mutate(month_f = factor(month, levels = mon_codes, labels = month.abb), season = season_of(month))

s1a <- ggplot(mb, aes(month_f, median_abs_pct_error, fill = season)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.0f%%", median_abs_pct_error)), vjust = -0.5, size = 3.4, fontface = "bold") +
  scale_fill_manual(values = SEASON_COLS, drop = FALSE) +
  scale_y_continuous(labels = pct_lab, expand = expansion(mult = c(0, 0.12))) +
  guides(fill = "none") +
  labs(title = "Median Error by Month", x = NULL, y = "Median Absolute Error") +
  theme_resource()

c_lim <- range(c(0, mb$median_bias_pct), na.rm = TRUE)
c_lim <- c_lim + c(-1, 1) * 0.15 * diff(c_lim)
s1b <- ggplot(mb, aes(month_f, median_bias_pct, colour = season)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.5) +
  geom_segment(aes(xend = month_f, yend = 0), linewidth = 0.6) +
  geom_point(size = 3.2) +
  annotate("text", x = 0.6, y = c_lim[2], hjust = 0, vjust = 1.2, size = 3.5, fontface = "italic",
           colour = "grey30", label = "Above 0: flow over-predicted\n(%WWef under-estimated)") +
  scale_colour_manual(values = SEASON_COLS, drop = FALSE) +
  scale_y_continuous(limits = c_lim, labels = pct_lab) +
  labs(title = "Median Bias by Month", x = NULL, y = "Median Bias") +
  theme_resource()

figS_loocv <- combine_panels(list(s1a, s1b),
                             title = "Leave-One-Gage-Out Error and Bias by Month",
                             collect = TRUE)
save_plot(figS_loocv, "figS_loocv_error_bias_by_month.png", w = 12.5, h = 5.4)

#Leave one out analysis on thresholds 
loocv_at <- function(th) {
  map_dfr(seq_along(gage_rows), function(k) {
    hold  <- gage_rows[k]
    anc   <- get_ancestors(hold)
    seeds <- setdiff(intersect(gage_rows, anc), hold)
    ord_k <- anc[order(ord_pos[anc])]
    b_k   <- fit_b(gage_rows[gage_rows != hold])$b
    Qk    <- propagate_dar(seeds, th, b_k, ord_k)$Q
    tibble(month = mon_codes, observed = Qraw[hold, ], predicted = Qk[hold, ])
  })
}
loocv_thr <- map_dfr(c(0.10, 0.25, 0.50, 0.75, Inf), function(th)
  summarise_loocv(loocv_at(th)) %>% mutate(threshold = th, .before = 1))
write_csv(loocv_thr, file.path(out_dir, "tbl_loocv_by_threshold.csv"))
print(as.data.frame(loocv_thr %>%
                      select(threshold, pct_reconstructed, median_abs_pct_error, median_bias_pct)), row.names = FALSE)

#Leave one out analysis on b
loocv_at_b <- function(b_fixed) {
  map_dfr(seq_along(gage_rows), function(k) {
    hold  <- gage_rows[k]
    anc   <- get_ancestors(hold)
    seeds <- setdiff(intersect(gage_rows, anc), hold)
    ord_k <- anc[order(ord_pos[anc])]
    Qk    <- propagate_dar(seeds, SIMILARITY_THRESHOLD, b_fixed, ord_k)$Q
    tibble(month = mon_codes, observed = Qraw[hold, ], predicted = Qk[hold, ])
  })
}
b_vals <- c(B_CI[1], B_STATE, B_CI[2], 1)
loocv_b <- map_dfr(b_vals, function(b)
  summarise_loocv(loocv_at_b(b)) %>% mutate(b = b, .before = 1))
write_csv(loocv_b, file.path(out_dir, "tbl_loocv_by_exponent.csv"))
print(as.data.frame(loocv_b %>%
                      select(b, pct_reconstructed, median_abs_pct_error, median_bias_pct)), row.names = FALSE)


# SENSITIVITY ANALYSES (exponent b and similarity gate): helpers
# Each scenario = one full-network propagation with all real gages in place, compared
# with the adopted run (fitted b, 50% gate). Summaries use WASTEWATER-IMPACTED segments
# only (cumulative upstream effluent > 0); all others are 0% %WWef under any setting.
# "Change" = |%WWef in scenario - %WWef in adopted run|, in percentage points, for
# segments that have an estimate in BOTH runs. Segments that gain or lose an estimate
# are counted separately.
ww_pos  <- !is.na(WWcum_mat) & WWcum_mat > 0
REF_Q   <- propagate_dar(gage_rows)$Q
REF_PCT <- compute_pct(WWcum_mat, REF_Q)

summarise_scenario <- function(Q) {
  pct <- compute_pct(WWcum_mat, Q)
  map_dfr(seq_len(n_mon), function(k) {
    imp     <- ww_pos[, k]
    p       <- pct[imp, k]
    chg     <- abs(p - REF_PCT[imp, k])
    chg     <- chg[is.finite(chg)]
    has_ref <- is.finite(REF_Q[imp, k])
    has_new <- is.finite(Q[imp, k])
    tibble(month                = mon_codes[k],
           pct_with_estimate    = 100 * mean(has_new),
           pct_gained           = 100 * mean(has_new & !has_ref),
           pct_lost             = 100 * mean(!has_new & has_ref),
           median_pct           = median(p, na.rm = TRUE),            # all impacted segments with an estimate
           median_pct_gained    = if (any(has_new & !has_ref)) median(p[has_new & !has_ref], na.rm = TRUE) else NA_real_,
           median_pct_lost      = if (any(!has_new & has_ref)) median(REF_PCT[imp, k][!has_new & has_ref], na.rm = TRUE) else NA_real_,
           p90_pct              = quantile(p, 0.90, na.rm = TRUE, names = FALSE),
           median_abs_change_pp = if (length(chg)) median(chg) else NA_real_,
           pct_change_gt1pp     = if (length(chg)) 100 * mean(chg > 1) else NA_real_,
           p99_abs_change_pp    = if (length(chg)) quantile(chg, 0.99, names = FALSE) else NA_real_,
           max_abs_change_pp    = if (length(chg)) max(chg) else NA_real_)
  })
}
sweep_scenarios <- function(values, run_fn) {
  map_dfr(values, function(v) summarise_scenario(run_fn(v)) %>% mutate(value = v, .before = 1))
}
fmt_pp <- function(x) ifelse(x < 0.01, "<0.01", sprintf("%.2f", x))

# EXPONENT SENSITIVITY  
message("Exponent sweep...")
B_GRID <- seq(B_CI[1], B_CI[2], length.out = N_EXPONENT_STEPS)   # 7 values; midpoint = fitted b
expo   <- sweep_scenarios(B_GRID, function(b) propagate_dar(gage_rows, b = b)$Q) %>%
  rename(b = value)
write_csv(expo, file.path(out_dir, "tbl_exponent_sensitivity_by_month.csv"))

# Worst case across the 12 months for each exponent value
tbl_b <- expo %>%
  group_by(b) %>%
  summarise(med_chg = max(median_abs_change_pp, na.rm = TRUE),
            gt1     = max(pct_change_gt1pp,     na.rm = TRUE),
            p99     = max(p99_abs_change_pp,    na.rm = TRUE),
            mx      = max(max_abs_change_pp,    na.rm = TRUE), .groups = "drop") %>%
  transmute(`Exponent, b` = paste0(sprintf("%.3f", b),
                                   ifelse(abs(b - B_STATE) < 1e-6, " (fitted)", "")),
            `Median change (pp)`           = fmt_pp(med_chg),
            `Segments changing > 1 pp (%)` = sprintf("%.1f", gt1),
            `99th percentile change (pp)`  = fmt_pp(p99),
            `Maximum change (pp)`          = fmt_pp(mx))
write_csv(tbl_b, file.path(out_dir, "tblSX_exponent_sensitivity.csv"))
cat("\n=== TABLE SX: EXPONENT SENSITIVITY (worst month for each b) ===\n")
print(as.data.frame(tbl_b), row.names = FALSE)


# SIMILARITY-GATE SENSITIVITY  
message("Threshold sweep...")
THRESHOLD_GRID <- sort(unique(c(THRESHOLD_GRID, SIMILARITY_THRESHOLD)))
thr <- sweep_scenarios(THRESHOLD_GRID, function(th) propagate_dar(gage_rows, threshold = th)$Q) %>%
  rename(threshold = value)
write_csv(thr, file.path(out_dir, "tbl_threshold_sensitivity_by_month.csv"))

rng <- function(x, digits = 1) {                    # "lo-hi" across months; "-" if no values
  x <- x[is.finite(x)]
  if (!length(x)) return("-")
  lo <- min(x); hi <- max(x)
  if (abs(hi - lo) < 10^(-digits) / 2) sprintf(paste0("%.", digits, "f"), hi)
  else sprintf(paste0("%.", digits, "f-%.", digits, "f"), lo, hi)
}

tbl_th <- thr %>%
  group_by(threshold) %>%
  summarise(cov      = rng(pct_with_estimate, 0),
            gained   = rng(pct_gained),
            lost     = rng(pct_lost),
            med_all  = rng(median_pct, 2),
            med_gain = rng(median_pct_gained, 2),
            med_lost = rng(median_pct_lost, 2),
            gt1      = max(pct_change_gt1pp,  na.rm = TRUE),
            mx       = max(max_abs_change_pp, na.rm = TRUE), .groups = "drop") %>%
  transmute(`Gate threshold` = ifelse(is.infinite(threshold), "Removed",
                                      paste0(sprintf("%g%%", threshold * 100),
                                             ifelse(threshold == SIMILARITY_THRESHOLD, " (adopted)", ""))),
            `Max. downstream/upstream area ratio` = ifelse(is.infinite(threshold), "None",
                                                           sprintf("%.2f", 1 / (1 - threshold))),
            `Impacted segments with a gage-based flow estimate (%)` = cov,
            `Gained a gage-based estimate vs. adopted (%)`          = gained,
            `Lost a gage-based estimate vs. adopted (%)`            = lost,
            `Median %WWef, all segments with an estimate (%)`       = med_all,
            `Median %WWef, segments that gained an estimate (%)`    = med_gain,
            `Median %WWef, segments that lost an estimate (%)`      = med_lost,
            `Segments changing > 1 pp vs. adopted (%), maximum across months` = sprintf("%.1f", gt1),
            `Largest change vs. adopted (pp), maximum across months`          = fmt_pp(mx))
write_csv(tbl_th, file.path(out_dir, "tblSY_threshold_sensitivity.csv"))
cat("\n=== TABLE SY: SIMILARITY-GATE SENSITIVITY ===\n")
print(as.data.frame(tbl_th), row.names = FALSE)

# IS THE EXPONENT SEASONAL?  (monthly fits + nested-model F-test)
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
  labs(title = "Seasonality Evaluation for the Drainage-Area Exponent, b",
       x = NULL, y = "Fitted Exponent, b") +
  theme_resource()
save_plot(fig4, "fig4_seasonal_exponent.png", w = 10, h = 5.4)


# WHAT DOES LOOSENING THE GATE ADD?  (composition of the network)
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

# REAL-GAGE PAIRS: DOES OBSERVED FLOW GROW LIKE (Ad/Au)^b?
# Pairs = each gage with its NEAREST upstream gage on every branch (no other real
# gage in between) 
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

# --- error and bias by size of the area-ratio jump ---------------------------
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

short_lab <- map_chr(seq_along(bin_hi), function(i) {
  lo <- format(bin_breaks[i], drop0trailing = TRUE)
  if (is.infinite(bin_hi[i])) paste0(">", lo)
  else paste0(lo, "\u2013", format(bin_hi[i], drop0trailing = TRUE))
})
bin_x <- bin_summary %>%
  mutate(x_lab = factor(short_lab[as.integer(bin)], levels = short_lab))

gate_fill <- c("TRUE" = COL$accent, "FALSE" = COL$context)
gate_lab  <- c("TRUE"  = sprintf("Within %g%% gate (ratio \u2264 %g)", SIMILARITY_THRESHOLD * 100, GATE_MAX_RATIO),
               "FALSE" = "Beyond gate")

pair_bar <- function(ycol, ylab, title) {
  ggplot(bin_x, aes(x_lab, .data[[ycol]], fill = within_gate)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%.0f%%", .data[[ycol]]),
                  vjust = ifelse(.data[[ycol]] >= 0, -0.5, 1.4)),
              size = 3.5, fontface = "bold") +
    scale_fill_manual(values = gate_fill, labels = gate_lab, drop = FALSE) +
    scale_y_continuous(labels = pct_lab, expand = expansion(mult = c(0.12, 0.12))) +
    labs(title = title, x = "Downstream / Upstream Drainage-Area Ratio", y = ylab) +
    theme_resource()
}
f6a <- pair_bar("median_abs_pct_err", "Median Absolute Error", "Error")
f6b <- pair_bar("bias_pct", "Median Bias (+ = over-predicted)", "Bias")
fig6 <- combine_panels(list(f6a, f6b),
                       title = "Scaling Error Grows with Drainage-Area Ratio",
                       collect = TRUE) & theme(legend.position = "bottom")
save_plot(fig6, "fig6_pairs_error_by_area_ratio.png", w = 11, h = 5.4)

# --- Is the bias seasonal? (real gage pairs inside the gate) ---
gm_plot <- gate_month %>%
  mutate(month_f = factor(month, levels = mon_codes, labels = month.abb), season = season_of(month))
bias_lim <- range(c(0, gm_plot$lo_pct, gm_plot$hi_pct), na.rm = TRUE)
bias_lim <- bias_lim + c(-1, 1) * 0.12 * diff(bias_lim)

fig7 <- ggplot(gm_plot, aes(month_f, bias_pct, colour = season)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lo_pct, ymax = hi_pct), width = 0.2, linewidth = 0.6) +
  geom_point(size = 3.2) +
  annotate("text", x = 0.6, y = bias_lim[2], hjust = 0, vjust = 1.2, size = 3.5,
           fontface = "italic", colour = "grey30", label = "Above 0: flow over-predicted") +
  scale_colour_manual(values = SEASON_COLS, drop = FALSE) +
  scale_y_continuous(limits = bias_lim, labels = pct_lab) +
  labs(title = "Bias by Month, Paired Gages Within the Gate",
       x = NULL, y = "Median Bias (95% CI)") +
  theme_resource()
save_plot(fig7, "fig7_seasonal_bias.png", w = 8, h = 5.4)

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

#Low-flow screen: how many gages were removed? 
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

#Exponent Determination 
# Fit: ln(annual mean flow) = a + b * ln(drainage area), ordinary least squares,
fit_dat <- tibble(site_no = flow_order$site_no[gage_rows],
                  A       = DA_vec[gage_rows],
                  Q       = ann_mean_flow[gage_rows]) %>%
  filter(is.finite(A), is.finite(Q), A > 0, Q > 0)

cf  <- coef(full_fit$model); V <- vcov(full_fit$model)
tq  <- qt(0.975, df = full_fit$model$df.residual)
xg  <- seq(log(min(fit_dat$A)), log(max(fit_dat$A)), length.out = 200)
band <- tibble(A   = exp(xg),
               fit = exp(cf[1] + cf[2] * xg),
               se  = sqrt(V[1, 1] + 2 * xg * V[1, 2] + xg^2 * V[2, 2]),
               lo  = exp(cf[1] + cf[2] * xg - tq * se),
               hi  = exp(cf[1] + cf[2] * xg + tq * se))

removed <- pre %>% filter(low_flow_excluded, is.finite(DA), DA > 0,
                          is.finite(mean_annual_cfs), mean_annual_cfs > 0)

fig_exp <- ggplot() +
  geom_ribbon(data = band, aes(A, ymin = lo, ymax = hi), fill = COL$band, alpha = 0.7) +
  geom_point(data = removed, aes(DA, mean_annual_cfs), shape = 4, size = 1.8, colour = COL$muted) +
  geom_point(data = fit_dat, aes(A, Q), colour = COL$accent, alpha = 0.55, size = 1.8) +
  geom_line(data = band, aes(A, fit), colour = COL$ink, linewidth = 0.8) +
  annotate("text", x = min(fit_dat$A), y = max(fit_dat$Q), hjust = 0, vjust = 1, size = 3.9,
           fontface = "italic",
           label = sprintf("b = %.3f (95%% CI %.3f\u2013%.3f)\nR\u00b2 = %.2f, n = %d gages",
                           B_STATE, B_CI[1], B_CI[2], full_fit$r2, full_fit$n)) +
  scale_x_log10(labels = fmt_num) +
  scale_y_log10(labels = fmt_num) +
  labs(title    = "Determination of the Drainage-Area Exponent (b)",
       x = "Drainage Area (km\u00b2)", y = "Annual Mean Fage Glow (cfs)") +
  theme_resource() +
  theme(plot.subtitle = element_text(hjust = 0.5, colour = "grey30", size = 10))

save_plot(fig_exp, "fig_exponent_fit.png", w = 8, h = 5.8)

# ---- How often does %WWef reach the 100% cap? ----
raw_pct <- 100 * WWcum_mat / Q_adopted                  # uncapped %WWef
imp     <- ww_pos & is.finite(Q_adopted) & Q_adopted > 0 # impacted segments with a gage-based estimate
at_cap  <- imp & raw_pct >= 100

cap_by_month <- map_dfr(seq_len(n_mon), function(k) {
  i_all <- imp[, k]; i_gag <- i_all & is_gage
  tibble(month           = mon_codes[k],
         n_impacted      = sum(i_all),
         n_at_cap        = sum(at_cap[, k]),
         pct_at_cap      = 100 * sum(at_cap[, k]) / sum(i_all),
         pct_at_cap_gage = if (any(i_gag)) 100 * sum(at_cap[i_gag, k]) / sum(i_gag) else NA_real_,
         median_uncapped = if (any(at_cap[, k])) median(raw_pct[at_cap[, k], k]) else NA_real_)
})
write_csv(cap_by_month, file.path(out_dir, "tbl_pct_wwef_at_cap.csv"))
print(as.data.frame(cap_by_month %>% mutate(across(where(is.double), ~ round(.x, 1)))), row.names = FALSE)

ever_cap <- rowSums(at_cap) > 0
cat(sprintf("\nPooled: %.1f%% of impacted segment-months are at the cap (%s of %s)\n",
            100 * sum(at_cap) / sum(imp), format(sum(at_cap), big.mark = ","), format(sum(imp), big.mark = ",")))
cat(sprintf("%.1f%% of impacted segments are at the cap in at least one month\n",
            100 * mean(ever_cap[rowSums(imp) > 0])))

cap_any <- which(rowSums(at_cap[, c(5, 7)]) > 0)          # capped in May or July
cat("June %WWef at those segments:\n"); print(summary(raw_pct[cap_any, 6]))
cat("\nEffluent (cfs) at those segments, by month:\n"); print(round(colSums(WWcum_mat[cap_any, ], na.rm = TRUE), 1))
cat("\nMean estimated flow (cfs) at those segments, by month:\n"); print(round(colMeans(Q_adopted[cap_any, ], na.rm = TRUE), 1))
