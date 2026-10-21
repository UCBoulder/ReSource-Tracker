# ==============================================================================
# equity_stats_offline.R  -  run ONCE; writes 4 spreadsheets (CSV)
# Run in the same R session where the app's objects exist (fl_cejst, cejst, cejst_lgl_cols, mon_codes)
# ==============================================================================
list2env(readRDS(file.path(path.expand("~"), "equity_inputs.rds")), envir = .GlobalEnv)
sapply(c("fl_cejst", "cejst", "cejst_lgl_cols", "mon_codes", "tracts_sf"), exists)
head(grep("pct_mean", names(fl_cejst), value = TRUE))
grep("DA|Order|HUC", names(fl_cejst), value = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(sf); library(sandwich); library(spdep)
})
stopifnot(exists("fl_cejst"), exists("cejst"), exists("cejst_lgl_cols"), exists("mon_codes"))

# ---- CONFIG  -------------------------------
SOURCES  <- c("Gage-based", "NHDPlus-based")
MONTHS   <- mon_codes                    #
DA_COL   <- "TotDASqKM"                  # drainage-area column in fl_cejst 
#severity: drainage area; prevalence: drainage area + number of segments)
HUC_COL  <- NA_character_                # e.g. "HUC8"; NA -> use county (GEOID digits 1-5)
B_BOOT   <- 2000                         # resamples for rural/urban CIs
B_SEASON <- 500                          # resamples for seasonality CIs
out_dir  <- "equity_stats_outputs"; dir.create(out_dir, showWarnings = FALSE)
set.seed(42)

# ---- Prepare data ----------------------------------------------------------------
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
dat$.lDA[!is.finite(dat$.lDA)] <- NA_real_      # DA of 0 or negative gives -Inf/NaN 
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

# One row per tract for a given source + month
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

# ----Equity associations (unadjusted vs adjusted for stream size) + Moran's I ------
# Severity:   estimate = rank correlation (rho) between the indicator and %WW among impacted tracts
# Prevalence: estimate = odds ratio of tract having any impacted segment
#             (per +1 SD for continuous indicators; Yes vs No for yes/no indicators)
# "Adjusted" = also controls for drainage area (+ number of segments for prevalence).
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

# ---- Rural vs urban: rank-biserial, median/risk difference, odds ratio, all with CIs ----
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

# ---- Seasonality: Kendall's W + paired month differences with CIs ---------------------
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