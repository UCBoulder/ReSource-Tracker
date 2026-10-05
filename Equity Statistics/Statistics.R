# ==============================================================================
# equity_stats_offline.R  -  run ONCE on your machine; writes 4 spreadsheets (CSV)
# that the app can read later. Does not touch app.R. Covers Gage-based AND NHDPlus-based.
# Not run on real data - parse-checked only.
#
# REVIEWER COMMENT                                   -> OUTPUT FILE
#  1. Stream hydrology confounds the EJ gradients     -> equity_associations.csv  (unadjusted vs adjusted columns)
#  2. Neighbouring tracts are spatially correlated    -> moran_check.csv  + CIs/p-values that resample whole watersheds
#  3. Seasonal differences need effect sizes / CIs    -> seasonality_effects.csv
#  4. Rural vs urban needs OR / CI / rank-biserial    -> rural_urban_effects.csv
#
# Run in the same R session where the app's objects exist (fl_cejst, cejst, cejst_lgl_cols, mon_codes).
# Packages: dplyr tidyr tibble readr sf sandwich spdep rlang
# ==============================================================================
list2env(readRDS(file.path(path.expand("~"), "equity_inputs.rds")), envir = .GlobalEnv)
sapply(c("fl_cejst", "cejst", "cejst_lgl_cols", "mon_codes", "tracts_sf"), exists)
head(grep("pct_mean", names(fl_cejst), value = TRUE))
grep("DA|Order|HUC", names(fl_cejst), value = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(sf); library(sandwich); library(spdep)
})
stopifnot(exists("fl_cejst"), exists("cejst"), exists("cejst_lgl_cols"), exists("mon_codes"))

# ---- CONFIG (the only part you should need to edit) -------------------------------
SOURCES  <- c("Gage-based", "NHDPlus-based")
MONTHS   <- mon_codes                    # e.g. c("03","06","09","12") to run faster
DA_COL   <- "TotDASqKM"                  # drainage-area column in fl_cejst
ORD_COL  <- "StreamOrde"                 # stream-order column (set NA to skip)
HUC_COL  <- NA_character_                # e.g. "HUC8"; NA -> use county (GEOID digits 1-5)
B_BOOT   <- 2000                         # resamples for rural/urban CIs
B_SEASON <- 500                          # resamples for seasonality CIs
out_dir  <- "equity_stats_outputs"; dir.create(out_dir, showWarnings = FALSE)
set.seed(42)

# ---- 0. Prepare data ----------------------------------------------------------------
if (!DA_COL %in% names(fl_cejst)) stop("DA_COL not found in fl_cejst: ", DA_COL)
has_ord <- !is.na(ORD_COL) && ORD_COL %in% names(fl_cejst)
dat0 <- fl_cejst[!is.na(fl_cejst$GEOID), ]
tract_xy <- NULL                          # tract location (for Moran's I) = centroid of the tract polygon
if (exists("tracts_sf")) {
  cen <- sf::st_coordinates(suppressWarnings(sf::st_centroid(sf::st_transform(sf::st_geometry(tracts_sf), 5070))))
  tract_xy <- tibble(GEOID = tracts_sf$GEOID, X = cen[, 1], Y = cen[, 2]) %>% distinct(GEOID, .keep_all = TRUE)
} else message("tracts_sf not found: Moran's I will be skipped.")
dat <- if (inherits(dat0, "sf")) sf::st_drop_geometry(dat0) else dat0
dat$.cl  <- if (!is.na(HUC_COL)) substr(as.character(dat[[HUC_COL]]), 1, 8) else substr(as.character(dat$GEOID), 1, 5)
dat$.lDA <- log10(as.numeric(dat[[DA_COL]]))
dat$.lDA[!is.finite(dat$.lDA)] <- NA_real_      # DA of 0 or negative gives -Inf/NaN, which breaks the models
dat$.ord <- if (has_ord) as.numeric(dat[[ORD_COL]]) else NA_real_
cat(sprintf("Clusters: %s | stream order adjustment: %s\n",
            if (!is.na(HUC_COL)) HUC_COL else "county", if (has_ord) "yes" else "no"))

ind_cols <- setdiff(intersect(c(names(cejst)[sapply(cejst, is.numeric)], cejst_lgl_cols), names(dat)),
                    c("GEOID", DA_COL, ORD_COL))
tract_ind <- dat %>% group_by(GEOID) %>%
  summarise(rurality = first(rurality), across(all_of(ind_cols), ~ first(.x)), .groups = "drop") %>%
  mutate(across(all_of(intersect(cejst_lgl_cols, ind_cols)), as.numeric),
         rurality_rural = as.numeric(rurality == "Rural"))
ind_all <- c("rurality_rural", ind_cols)

# one row per tract for a given source + month (same definitions as the app)
build_tract <- function(mm, src) {
  pc <- if (src == "Gage-based") paste0("pct_mean_", mm) else paste0("pct_mean_", mm, "_nhd")
  wc <- paste0("WW_cum_mean_", mm)
  if (!all(c(pc, wc) %in% names(dat))) return(NULL)
  tibble(GEOID = dat$GEOID, cl = dat$.cl, lDA = dat$.lDA, ord = dat$.ord,
         ww = as.numeric(dat[[wc]]), pct = as.numeric(dat[[pc]])) %>%
    mutate(imp_ww = is.finite(ww) & ww > 0, imp_pct = is.finite(pct) & pct > 0) %>%
    group_by(GEOID) %>%
    summarise(cl = first(cl), n_seg = n(), any_imp = as.integer(any(imp_ww)),
              lDA_all = median(lDA, na.rm = TRUE), ord_all = median(ord, na.rm = TRUE),
              sev     = if (any(imp_pct)) median(pct[imp_pct]) else NA_real_,
              lDA_imp = if (any(imp_pct)) median(lDA[imp_pct], na.rm = TRUE) else NA_real_,
              ord_imp = if (any(imp_pct)) median(ord[imp_pct], na.rm = TRUE) else NA_real_,
              .groups = "drop") %>%
    mutate(lnseg = log(n_seg),
           across(c(lDA_all, lDA_imp, ord_all, ord_imp), ~ ifelse(is.finite(.x), .x, NA_real_))) %>%
    left_join(tract_ind, by = "GEOID")
}

# ---- helpers --------------------------------------------------------------------------
zs <- function(x) as.numeric(scale(x));  zr <- function(x) as.numeric(scale(rank(x)))

robust_coef <- function(m, term, cl, odds = FALSE) {       # CI that accounts for clustering of tracts in watersheds
  G <- length(unique(cl)); V <- sandwich::vcovCL(m, cluster = cl, type = "HC1")
  b <- unname(coef(m)[term]); se <- sqrt(V[term, term]); tc <- qt(0.975, G - 1)
  r <- c(est = b, lo = b - tc * se, hi = b + tc * se); if (odds) r <- exp(r)
  c(r, p = 2 * pt(-abs(b / se), G - 1))
}
cluster_boot <- function(cl, B, f) {                         # resample whole watersheds; f(index) -> named numbers
  idx <- split(seq_along(cl), cl); G <- length(idx); pt <- f(seq_along(cl))
  bs <- matrix(vapply(seq_len(B), function(i) f(unlist(idx[sample.int(G, G, TRUE)], use.names = FALSE)),
                      numeric(length(pt))), nrow = length(pt))
  tibble(stat = names(pt), est = unname(pt), lo = apply(bs, 1, quantile, 0.025, na.rm = TRUE),
         hi = apply(bs, 1, quantile, 0.975, na.rm = TRUE))
}
cliff <- function(x, g) {                                    # rank-biserial: P(Rural > Urban) - P(Rural < Urban)
  a <- x[g == 1]; b <- x[g == 0]; if (!length(a) || !length(b)) return(NA_real_)
  2 * (sum(rank(c(a, b))[seq_along(a)]) - length(a) * (length(a) + 1) / 2) / (length(a) * length(b)) - 1
}
.lw <- new.env()
moran_of <- function(v, geoid) {                             # Moran's I with 8 nearest tracts
  if (is.null(tract_xy)) return(c(I = NA_real_, p = NA_real_, n = NA_real_))
  keep <- is.finite(v) & geoid %in% tract_xy$GEOID
  key <- rlang::hash(geoid[keep])
  if (!exists(key, envir = .lw, inherits = FALSE)) {
    co <- as.matrix(tract_xy[match(geoid[keep], tract_xy$GEOID), c("X", "Y")]); co <- co + matrix(rnorm(length(co)), ncol = 2)
    assign(key, nb2listw(knn2nb(knearneigh(co, k = 8)), style = "W", zero.policy = TRUE), envir = .lw)
  }
  mt <- moran.test(v[keep], get(key, envir = .lw), zero.policy = TRUE)
  c(I = unname(mt$estimate[1]), p = mt$p.value, n = sum(keep))
}

# ---- 1+2. Equity associations (unadjusted vs adjusted for stream size) + Moran's I ------
# Severity:   estimate = rank correlation (rho) between the indicator and %WW among impacted tracts
# Prevalence: estimate = odds ratio of tract having any impacted segment
#             (per +1 SD for continuous indicators; Yes vs No for yes/no indicators)
# "Adjusted" = also controls for drainage area (+ stream order; + number of segments for prevalence).
fit_one <- function(d, ind, covs, mode) {
  covs <- covs[vapply(covs, function(v) dplyr::n_distinct(d[[v]][is.finite(d[[v]])]) > 1, logical(1))]
  dd <- bind_cols(tibble(cl = d$cl, y = d$y, x = d[[ind]]), d[covs])
  numc <- vapply(dd, is.numeric, logical(1))
  dd <- dd[complete.cases(dd) & Reduce(`&`, lapply(dd[numc], is.finite)), ]
  G <- n_distinct(dd$cl)
  if (nrow(dd) < 50 || G < 15 || n_distinct(dd$x) < 2) return(NULL)
  is_bin <- all(dd$x %in% c(0, 1))
  if (mode == "Severity") {
    dd$ys <- zr(dd$y); dd$xs <- if (is_bin) zs(dd$x) else zr(dd$x); for (v in covs) dd[[v]] <- zr(dd[[v]])
    m0 <- lm(ys ~ xs, dd); m1 <- lm(reformulate(c("xs", covs), "ys"), dd); odds <- FALSE
  } else {
    dd$xs <- if (is_bin) dd$x else zs(dd$x); for (v in covs) dd[[v]] <- zs(dd[[v]])
    m0 <- suppressWarnings(glm(y ~ xs, binomial, dd)); m1 <- suppressWarnings(glm(reformulate(c("xs", covs), "y"), binomial, dd)); odds <- TRUE
  }
  a <- robust_coef(m0, "xs", dd$cl, odds); b <- robust_coef(m1, "xs", dd$cl, odds)
  tibble(indicator = ind, type = if (is_bin) "yes/no" else "continuous", n_tracts = nrow(dd), n_watersheds = G,
         effect_unadj = a["est"], lo_unadj = a["lo"], hi_unadj = a["hi"], p_unadj = a["p"],
         effect_adj = b["est"], lo_adj = b["lo"], hi_adj = b["hi"], p_adj = b["p"])
}

assoc <- list(); moran_rows <- list(); tr_store <- list()
for (src in SOURCES) for (mm in MONTHS) {
  tr <- build_tract(mm, src); if (is.null(tr)) next
  tr_store[[paste(src, mm)]] <- tr
  for (mode in c("Prevalence", "Severity")) {
    d <- if (mode == "Severity") filter(tr, is.finite(sev)) else tr
    d$y  <- if (mode == "Severity") d$sev else d$any_imp
    covs <- if (mode == "Severity") c("lDA_imp", if (has_ord) "ord_imp") else c("lDA_all", "lnseg", if (has_ord) "ord_all")
    # Moran's I of the outcome, and of what is left after removing stream-size effects
    ok <- complete.cases(d[, c("y", covs)]) & Reduce(`&`, lapply(d[, c("y", covs)], is.finite)); dm <- d[ok, ]
    res <- if (mode == "Severity") residuals(lm(reformulate(covs, "ysr"), transform(dm, ysr = rank(y), lDA_imp = dm$lDA_imp)))
    else residuals(suppressWarnings(glm(reformulate(covs, "y"), binomial, dm)), type = "pearson")
    mo <- moran_of(if (mode == "Severity") rank(dm$y) else dm$y, dm$GEOID); mr <- moran_of(res, dm$GEOID)
    moran_rows[[length(moran_rows) + 1]] <- tibble(source = src, month = mm, mode = mode, n_tracts = nrow(dm),
                                                   n_watersheds = n_distinct(dm$cl), moran_I_outcome = mo["I"], p_outcome = mo["p"],
                                                   moran_I_after_stream_size = mr["I"], p_after_stream_size = mr["p"])
    rows <- bind_rows(lapply(ind_all, function(ind) fit_one(d, ind, covs, mode)))
    if (nrow(rows)) {
      rows$q_unadj <- p.adjust(rows$p_unadj, "BH"); rows$q_adj <- p.adjust(rows$p_adj, "BH")   # BH across indicators
      assoc[[length(assoc) + 1]] <- mutate(rows, source = src, month = mm, mode = mode, .before = 1)
    }
    message(sprintf("done: %s | %s | %s", src, mm, mode))
  }
}

assoc <- bind_rows(assoc); moran_tab <- bind_rows(moran_rows)
write_csv(assoc, file.path(out_dir, "equity_associations.csv"))
write_csv(moran_tab, file.path(out_dir, "moran_check.csv"))
cat("\n--- Share of indicators significant after BH (q < 0.05): before vs after adjusting for stream size ---\n")
print(assoc %>% group_by(source, mode) %>%
        summarise(unadjusted = mean(q_unadj < 0.05, na.rm = TRUE), adjusted = mean(q_adj < 0.05, na.rm = TRUE), .groups = "drop"))
cat("\n--- Moran's I (positive = neighbouring tracts are similar) ---\n"); print(moran_tab %>% head(12))

# ---- 4. Rural vs urban: rank-biserial, median/risk difference, odds ratio, all with CIs ----
ru <- list()
for (key in names(tr_store)) {
  tr <- tr_store[[key]]; src <- sub(" .*", "", key); mm <- sub(".* ", "", key)
  d <- filter(tr, !is.na(rurality_rural), !is.na(cl)); s <- filter(d, is.finite(sev))
  bs <- cluster_boot(s$cl, B_BOOT, function(i) c(rank_biserial = cliff(s$sev[i], s$rurality_rural[i]),
                                                 median_diff_pp = median(s$sev[i][s$rurality_rural[i] == 1]) - median(s$sev[i][s$rurality_rural[i] == 0])))
  bp <- cluster_boot(d$cl, B_BOOT, function(i) c(risk_diff_pp = 100 * (mean(d$any_imp[i][d$rurality_rural[i] == 1]) -
                                                                         mean(d$any_imp[i][d$rurality_rural[i] == 0]))))
  a <- filter(assoc, source == src, month == mm, indicator == "rurality_rural")
  ru[[key]] <- bind_rows(
    mutate(bs, mode = "Severity"), mutate(bp, mode = "Prevalence"),
    a %>% filter(mode == "Prevalence") %>% transmute(stat = "odds_ratio_unadj", est = effect_unadj, lo = lo_unadj, hi = hi_unadj, mode),
    a %>% filter(mode == "Prevalence") %>% transmute(stat = "odds_ratio_adj",   est = effect_adj,   lo = lo_adj,   hi = hi_adj,   mode),
    a %>% filter(mode == "Severity")   %>% transmute(stat = "rank_corr_unadj",  est = effect_unadj, lo = lo_unadj, hi = hi_unadj, mode),
    a %>% filter(mode == "Severity")   %>% transmute(stat = "rank_corr_adj",    est = effect_adj,   lo = lo_adj,   hi = hi_adj,   mode)
  ) %>% mutate(source = src, month = mm, .before = 1)
}
ru <- bind_rows(ru); write_csv(ru, file.path(out_dir, "rural_urban_effects.csv"))
cat("\n--- Rural vs urban (direction: Rural relative to Urban; first month shown) ---\n")
print(ru %>% filter(month == MONTHS[1]), n = 24)

# ---- 3. Seasonality: Kendall's W + paired month differences with CIs ---------------------
season_effects <- function(src, units) {
  use_tracts <- identical(units, "tracts")
  long <- bind_rows(lapply(MONTHS, function(mm) {
    col <- if (src == "Gage-based") paste0("pct_mean_", mm) else paste0("pct_mean_", mm, "_nhd")
    if (!col %in% names(dat)) return(NULL)
    tibble(GEOID = dat$GEOID, Hydroseq = dat$Hydroseq, cl = dat$.cl, mon = mm, pct = suppressWarnings(as.numeric(dat[[col]])))
  })) %>% mutate(unit = if (use_tracts) GEOID else as.character(Hydroseq)) %>%
    group_by(unit, mon) %>% summarise(cl = first(cl), pct = suppressWarnings(median(pct, na.rm = TRUE)), .groups = "drop") %>%
    mutate(pct = if_else(is.finite(pct), pct, NA_real_))
  ever <- long %>% group_by(unit) %>% summarise(ever = any(!is.na(pct) & pct > 0), .groups = "drop")
  wide <- long %>% inner_join(ever, by = "unit") %>% filter(ever) %>% select(-ever) %>%
    pivot_wider(id_cols = c(unit, cl), names_from = mon, values_from = pct) %>% drop_na()
  if (nrow(wide) < 30) return(NULL)
  X <- as.matrix(wide[, MONTHS]); n <- nrow(X); k <- ncol(X); fr <- friedman.test(X)
  overall <- tibble(source = src, units = units, n_complete_units = n, n_watersheds = n_distinct(wide$cl),
                    friedman_p = fr$p.value, kendalls_W = unname(fr$statistic) / (n * (k - 1)))
  pairs <- combn(MONTHS, 2, simplify = FALSE); np <- length(pairs)
  nm <- vapply(pairs, paste, "", collapse = "|")
  pair_stats <- function(i) {
    Xs <- X[i, , drop = FALSE]; out <- numeric(2 * np)
    for (j in seq_len(np)) {
      d <- Xs[, pairs[[j]][1]] - Xs[, pairs[[j]][2]]; out[j] <- median(d); nz <- d[d != 0]
      out[np + j] <- if (length(nz)) { rk <- rank(abs(nz)); (sum(rk[nz > 0]) - sum(rk[nz < 0])) / sum(rk) } else 0
    }
    names(out) <- c(paste0("med|", nm), paste0("rb|", nm)); out
  }
  bt <- cluster_boot(wide$cl, B_SEASON, pair_stats) %>%
    separate(stat, c("measure", "m1", "m2"), sep = "\\|") %>% pivot_wider(names_from = measure, values_from = c(est, lo, hi))
  pw <- bind_rows(lapply(pairs, function(p) tibble(m1 = p[1], m2 = p[2],
                                                   p = wilcox.test(X[, p[1]], X[, p[2]], paired = TRUE, exact = FALSE)$p.value))) %>% mutate(q_BH = p.adjust(p, "BH"))
  list(overall = overall, pairs = left_join(bt, pw, by = c("m1", "m2")) %>%
         transmute(source = src, units = units, month1 = month.abb[as.integer(m1)], month2 = month.abb[as.integer(m2)],
                   median_diff_pp = est_med, md_lo = lo_med, md_hi = hi_med,
                   rank_biserial = est_rb, rb_lo = lo_rb, rb_hi = hi_rb, q_BH))
}
se <- list(); for (src in SOURCES) for (u in c("segments", "tracts")) { r <- season_effects(src, u); if (!is.null(r)) se[[paste(src, u)]] <- r }
write_csv(bind_rows(lapply(se, `[[`, "pairs")), file.path(out_dir, "seasonality_effects.csv"))
write_csv(bind_rows(lapply(se, `[[`, "overall")), file.path(out_dir, "seasonality_overall.csv"))
cat("\n--- Seasonality: Kendall's W (0 = no seasonal pattern, 1 = identical month ranking everywhere) ---\n")
print(bind_rows(lapply(se, `[[`, "overall")))
cat(sprintf("\nDone. Spreadsheets are in ./%s\n", out_dir))

# =====================================================================
# equity_figures.R
# Three simple figures from equity_stats_outputs/equity_associations.csv
# Runs for BOTH sources (Gage-based and NHDPlus-based). Nothing here
# touches app.R. Run it in a fresh R session.
#
#   Fig 1  Forest plot       unadjusted vs stream-size-adjusted effect (one month)
#   Fig 2  Attenuation plot  how much each effect moves after adjusting (one month)
#   Fig 3  Consistency map   do the top indicators hold up in all 12 months?
#
# Each figure also writes a CSV of exactly what was plotted (the "table twin").
# =====================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr); library(patchwork)
})
has_repel <- requireNamespace("ggrepel", quietly = TRUE)

# ---- Config (the only thing you should need to edit) -----------------
in_csv        <- "equity_stats_outputs/equity_associations.csv"
out_dir       <- "equity_figures"
FEATURE_MONTH <- NULL      # e.g. "Jul" or "07" (must match the 'month' column); NULL = first month in the file
TOP_N         <- 15        # indicators shown per panel (ranked by size of adjusted effect)
ALPHA         <- 0.05      # BH q threshold for "significant"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

COL <- list(accent = "#2a78d6", unadj = "#8a8f98", up = "#2a78d6", down = "#e07b39",
            ns = "#e5e7eb", grid = "#e6e8eb", ink = "#1f2328")

# ---- Read + tidy ------------------------------------------------------
a <- read_csv(in_csv, show_col_types = FALSE,
              col_types = cols(.default = col_guess(), month = col_character()))
a <- a %>% filter(!grepl("percentile", indicator, ignore.case = TRUE))   # drop duplicate "(percentile)" versions
# wastewater variables are exposure-related, so they are circular next to the equity indicators; drop them here
a <- a %>% filter(!grepl("wastewater", indicator, ignore.case = TRUE))
months <- unique(a$month)
a$month <- factor(a$month, levels = months)
if (is.null(FEATURE_MONTH)) FEATURE_MONTH <- months[1]
stopifnot(FEATURE_MONTH %in% months)
message("Featured month: ", FEATURE_MONTH)

# Prevalence effects are odds ratios (plot on a log scale); severity effects are rank correlations (linear).
sc   <- function(x, mode) if (mode == "Prevalence") log(pmax(x, 1e-6)) else x
wrap <- function(s, w = 30) vapply(s, function(z) paste(strwrap(z, w), collapse = "\n"), "")

theme_resource <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          panel.border     = element_rect(colour = "#b8bdc4", fill = NA, linewidth = 0.5),
          panel.grid.major = element_line(colour = COL$grid, linewidth = 0.3),
          panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold", hjust = 0.5, size = 12),
          axis.title       = element_text(colour = COL$ink),
          axis.text        = element_text(colour = COL$ink),
          axis.text.y      = element_text(size = 8.5, lineheight = 0.85),
          legend.position  = "bottom")
}

# pick the TOP_N indicators for one source/mode, ranked by |adjusted effect| in the featured month
pick_top <- function(src, md) {
  a %>% filter(source == src, mode == md, month == FEATURE_MONTH, is.finite(effect_adj)) %>%
    mutate(rk = abs(sc(effect_adj, md))) %>% slice_max(rk, n = TOP_N, with_ties = FALSE) %>%
    arrange(sc(effect_adj, md)) %>% pull(indicator)
}

# ---- Fig 1: forest plot ----------------------------------------------
forest_panel <- function(src, md) {
  keep <- pick_top(src, md)
  d <- a %>% filter(source == src, mode == md, month == FEATURE_MONTH, indicator %in% keep)
  long <- bind_rows(
    d %>% transmute(indicator, model = "Unadjusted",              est = effect_unadj, lo = lo_unadj, hi = hi_unadj, q = q_unadj),
    d %>% transmute(indicator, model = "Adjusted for stream size", est = effect_adj,   lo = lo_adj,   hi = hi_adj,   q = q_adj)) %>%
    mutate(indicator = factor(wrap(indicator), levels = wrap(keep)),
           model = factor(model, levels = c("Unadjusted", "Adjusted for stream size")),
           sig = factor(ifelse(!is.na(q) & q < ALPHA, "sig", "ns"), levels = c("sig", "ns")))
  p <- ggplot(long, aes(est, indicator, colour = model, group = model)) +
    geom_vline(xintercept = if (md == "Prevalence") 1 else 0, colour = "#6b7280", linewidth = 0.4) +
    geom_linerange(aes(xmin = lo, xmax = hi), orientation = "y", linewidth = 0.6,
                   position = position_dodge(width = 0.65)) +
    geom_point(aes(shape = sig), size = 2.8, stroke = 1,
               position = position_dodge(width = 0.65)) +
    scale_shape_manual(values = c(sig = 16, ns = 1), name = NULL, drop = FALSE,
                       labels = c(sig = paste0("q < ", ALPHA), ns = paste0("q ≥ ", ALPHA))) +
    scale_colour_manual(values = c("Unadjusted" = "#9aa0a8", "Adjusted for stream size" = COL$accent), name = NULL) +
    guides(colour = guide_legend(order = 1, override.aes = list(shape = 16)),
           shape  = guide_legend(order = 2, override.aes = list(colour = "#1f2328"))) +
    labs(y = NULL, title = md,
         x = if (md == "Prevalence") "Odds ratio of any wastewater impact (log scale)" else "Rank correlation with median %WWef") +
    theme_resource()
  if (md == "Prevalence") p <- p + scale_x_log10()
  list(plot = p, data = long %>% mutate(source = src, mode = md, month = FEATURE_MONTH))
}

fig_forest <- function(src) {
  pp <- forest_panel(src, "Prevalence"); ps <- forest_panel(src, "Severity")
  fig <- (pp$plot | ps$plot) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A",
                    title = paste0(src, ", ", FEATURE_MONTH, ": association with wastewater exposure"),
                    caption = paste0("Filled = BH q < ", ALPHA, "; open = not significant. Bars: 95% cluster-robust CI. Top ", TOP_N, " indicators per panel (wastewater variables excluded)."),
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5),
                                  plot.caption = element_text(hjust = 0))) &
    theme(legend.position = "bottom")
  tag <- gsub("[^A-Za-z]", "", src)
  ggsave(file.path(out_dir, paste0("fig1_forest_", tag, ".png")), fig, width = 12, height = 9, dpi = 300, bg = "white")
  write_csv(bind_rows(pp$data, ps$data), file.path(out_dir, paste0("fig1_forest_", tag, "_data.csv")))
  invisible(fig)
}

# ---- Fig 2: attenuation (unadjusted -> adjusted) -----------------------
slope_panel <- function(src, md) {
  d <- a %>% filter(source == src, mode == md, month == FEATURE_MONTH, is.finite(effect_adj), is.finite(effect_unadj)) %>%
    transmute(indicator, unadj = sc(effect_unadj, md), adj = sc(effect_adj, md), q_adj) %>%
    mutate(status = ifelse(!is.na(q_adj) & q_adj < ALPHA, "Still significant after adjustment", "Not significant after adjustment"))
  long <- d %>% pivot_longer(c(unadj, adj), names_to = "step", values_to = "y") %>%
    mutate(step = factor(step, levels = c("unadj", "adj"), labels = c("Unadjusted", "Adjusted for\nstream size")))
  lab <- long %>% filter(step == "Adjusted for\nstream size", status == "Still significant after adjustment") %>%
    slice_max(abs(y), n = 8, with_ties = FALSE)
  p <- ggplot(long, aes(step, y, group = indicator, colour = status)) +
    geom_hline(yintercept = 0, colour = "#6b7280", linewidth = 0.4) +
    geom_line(linewidth = 0.5, alpha = 0.8) + geom_point(size = 1.8) +
    scale_colour_manual(values = c("Still significant after adjustment" = COL$accent,
                                   "Not significant after adjustment" = COL$unadj), name = NULL) +
    scale_x_discrete(expand = expansion(add = c(0.3, 1.6))) +
    labs(x = NULL, title = md,
         y = if (md == "Prevalence") "log odds ratio" else "Rank correlation (rho)") +
    theme_resource()
  p <- if (has_repel)
    p + ggrepel::geom_text_repel(data = lab, aes(label = indicator), size = 2.6, hjust = 0, direction = "y",
                                 nudge_x = 0.1, segment.size = 0.2, show.legend = FALSE, colour = COL$ink)
  else p + geom_text(data = lab, aes(label = indicator), size = 2.6, hjust = -0.05, show.legend = FALSE, colour = COL$ink)
  list(plot = p, data = d %>% mutate(source = src, mode = md, month = FEATURE_MONTH))
}

fig_slope <- function(src) {
  pp <- slope_panel(src, "Prevalence"); ps <- slope_panel(src, "Severity")
  fig <- (pp$plot | ps$plot) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A",
                    title = paste0(src, ", ", FEATURE_MONTH, ": how much does stream size explain?"),
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))) &
    theme(legend.position = "bottom")
  tag <- gsub("[^A-Za-z]", "", src)
  ggsave(file.path(out_dir, paste0("fig2_attenuation_", tag, ".png")), fig, width = 12, height = 6.5, dpi = 300, bg = "white")
  write_csv(bind_rows(pp$data, ps$data), file.path(out_dir, paste0("fig2_attenuation_", tag, "_data.csv")))
  invisible(fig)
}

# ---- Fig 3: consistency across months ----------------------------------
heat_panel <- function(src, md) {
  keep <- pick_top(src, md)
  d <- a %>% filter(source == src, mode == md, indicator %in% keep) %>%
    mutate(dir = sc(effect_adj, md),
           status = case_when(is.na(q_adj) | q_adj >= ALPHA ~ "Not significant",
                              dir > 0 ~ "Higher exposure", TRUE ~ "Lower exposure"),
           status = factor(status, levels = c("Higher exposure", "Lower exposure", "Not significant")),
           indicator = factor(wrap(indicator), levels = wrap(keep)))
  lim <- max(abs(d$dir), na.rm = TRUE)
  p <- ggplot(d, aes(month, indicator, fill = dir)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_point(data = filter(d, status != "Not significant"), shape = 21, fill = "white",
               colour = COL$ink, size = 1.6, stroke = 0.4) +
    scale_fill_gradient2(low = COL$down, mid = "#f1f2f4", high = COL$up, midpoint = 0, limits = c(-lim, lim),
                         name = if (md == "Prevalence") "log odds ratio" else "rank correlation (rho)") +
    labs(x = NULL, y = NULL, title = md) + theme_resource() +
    theme(panel.grid = element_blank(), legend.key.width = unit(1.4, "cm"))
  list(plot = p, data = d %>% select(source, mode, month, indicator, effect_adj, lo_adj, hi_adj, q_adj, status))
}

fig_heat <- function(src) {
  pp <- heat_panel(src, "Prevalence"); ps <- heat_panel(src, "Severity")
  fig <- (pp$plot | ps$plot) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A",
                    title = paste0(src, ": adjusted associations by month"),
                    caption = paste0("Adjusted for stream size. Blue = higher exposure, orange = lower. White dot = BH q < ", ALPHA, ". Indicators are the top ", TOP_N, " in ", FEATURE_MONTH, "."),
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5),
                                  plot.caption = element_text(hjust = 0))) &
    theme(legend.position = "bottom")
  tag <- gsub("[^A-Za-z]", "", src)
  ggsave(file.path(out_dir, paste0("fig3_months_", tag, ".png")), fig, width = 13, height = 6.5, dpi = 300, bg = "white")
  write_csv(bind_rows(pp$data, ps$data) %>% mutate(indicator = gsub("\n", " ", indicator)),
            file.path(out_dir, paste0("fig3_months_", tag, "_data.csv")))
  invisible(fig)
}

# ---- Run for both sources -------------------------------------------------
for (src in unique(a$source)) {
  message("Source: ", src)
  fig_forest(src); fig_slope(src); fig_heat(src)
}
message("Done. Figures and table-twin CSVs are in: ", normalizePath(out_dir))

