# ============================================================
# ReSource Tracker plots for publication

# Manuscript plots: seasonality + field vs. tool comparisons
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(purrr)
library(cowplot)

#--------------------------------------------------------------
# 0. Paths, shared settings, helpers
#--------------------------------------------------------------

# Combined data workbook (seasonality sheets, nitrate/conductivity sheet)
xlsx_path <- "/Users/kylieboenischoakes/Library/CloudStorage/OneDrive-UCB-O365/General/Boenisch-Oakes, Kylie/Manuscripts/De Facto Reuse Model + DBPs/Data_ManuscriptModel/CombinedData_ManuscriptPlots.xlsx"

# Field vs. tool workbook
model_path <- "/Users/kylieboenischoakes/Library/CloudStorage/OneDrive-UCB-O365/General/Boenisch-Oakes, Kylie/Research/Chapter 4/ModelVToolData.xlsx"

# Month ordering
month_order_full <- month.name      # "January"..."December"
month_order_abbr <- month.abb       # "Jan"..."Dec"

# Output folder
out_dir <- "/Users/kylieboenischoakes/Library/CloudStorage/OneDrive-UCB-O365/General/Boenisch-Oakes, Kylie/Manuscripts/De Facto Reuse Model + DBPs/Figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#Save figures
save_fig <- function(plot, name, w, h) {
  png_dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(file.path(out_dir, paste0(name, ".png")), plot,
         width = w, height = h, units = "in", dpi = 600, bg = "white",
         device = png_dev)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = w, height = h, units = "in", device = grDevices::pdf)
  message("Saved: ", name)
}

# Text sizes
TXT_TITLE <- 20; TXT_SUB <- 15; TXT_AXIS_TITLE <- 17; TXT_AXIS <- 15; TXT_LEG <- 14
GEOM_TXT  <- 6        
PT_SIZE   <- 4.5      

# Manuscript plot theme
theme_ms <- function() {
  theme_bw(base_size = TXT_AXIS) +
    theme(
      plot.title    = element_text(face = "bold", size = TXT_TITLE, hjust = 0.5),
      plot.subtitle = element_text(size = TXT_SUB, hjust = 0.5),
      axis.title.x  = element_text(face = "bold", size = TXT_AXIS_TITLE, margin = margin(t = 10)),
      axis.title.y  = element_text(face = "bold", size = TXT_AXIS_TITLE, margin = margin(r = 10)),
      axis.text.x   = element_text(size = TXT_AXIS, colour = "black"),
      axis.text.y   = element_text(size = TXT_AXIS, colour = "black"),
      legend.title  = element_text(face = "bold", size = TXT_LEG + 1),
      legend.text   = element_text(size = TXT_LEG),
      legend.key.size = unit(1.1, "lines"),
      panel.grid.major.x = element_blank(),
      plot.margin   = margin(10, 14, 10, 10)
    )
}

#--------------------------------------------------------------
# 1. Seasonality: line + ribbon plots
#--------------------------------------------------------------

seasonality_gage <- read_excel(xlsx_path, sheet = "Seasonality_Gage")
seasonality_nhd  <- read_excel(xlsx_path, sheet = "Seasonality_NHD")

prep_seasonality <- function(df) {
  df %>%
    mutate(
      month      = factor(month, levels = month_order_full),
      month_abbr = factor(month.abb[mon_code], levels = month_order_abbr)
    )
}

seasonality_gage <- prep_seasonality(seasonality_gage)
seasonality_nhd  <- prep_seasonality(seasonality_nhd)

plot_seasonality <- function(df, denom_label) {
  ggplot(df, aes(x = month_abbr, y = med_pct, group = 1)) +
    geom_ribbon(
      aes(ymin = lo_iqr, ymax = hi_iqr),
      fill   = "grey85",
      alpha  = 0.6,
      colour = NA
    ) +
    geom_line(linewidth = 1.6) +
    geom_point(size = 3.8) +
    labs(
      x = "Month",
      y = "Median % Wastewater by Census Tract",
      title = denom_label,
      subtitle = NULL
    ) +
    theme_ms() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
}

p_season_gage <- plot_seasonality(seasonality_gage, "Gage-based")
p_season_nhd  <- plot_seasonality(seasonality_nhd,  "NHDPlus-based")

# Show plots side by side
p_season_final <- (p_season_gage | p_season_nhd) +
  plot_annotation(
    tag_levels = "A",
    caption = "Line = median; shaded band = interquartile range",
    theme = theme(plot.caption = element_text(size = TXT_SUB, hjust = 0.5),
                  plot.tag = element_text(size = TXT_TITLE, face = "bold"))
  )
print(p_season_final)
save_fig(p_season_final, "fig_seasonality", w = 14, h = 6.5)

# ============================================================
# 2. Field vs. tool comparison (gage-based & NHDPlus)
# ============================================================

# Read and organize input data 
rdata <- read_excel(model_path, sheet = "RData") %>%
  rename(
    location    = `Location`,
    field_mean  = `%WW Field`,
    field_sd    = `%WW Field Sd`,
    gage_pctww  = `%WW Gage`,
    gage_max    = `Max Gage%`,
    gage_min    = `Min Gage %`,
    nhd_pctww   = `%WW NHDPlus`
  ) %>%
  mutate(
    
    month_chr = str_trim(as.character(Month)),
    
    month_str = case_when(
    
      str_detect(month_chr, "^[0-9]+(\\.[0-9]+)?$") ~
        format(as.Date(as.numeric(month_chr), origin = "1899-12-30"), "%b-%y"),
      
      
      str_detect(month_chr, "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}$|^[0-9]{4}-[0-9]{2}-[0-9]{2}$") ~
        format(lubridate::parse_date_time(month_chr, orders = c("mdy", "ymd")), "%b-%y"),
      
      
      TRUE ~ month_chr
    ),
    
    # convert fractional columns to percent
    across(
      c(field_mean, field_sd, gage_pctww, gage_max, gage_min, nhd_pctww),
      ~ as.numeric(.x) * 100
    ),
    
    month_idx = match(substr(month_str, 1, 3), month.abb),
    site_date = paste0("Site ", location, " - ", month_str)
  ) %>%
  arrange(location, month_idx) %>%
  mutate(
    site_date   = factor(site_date, levels = unique(site_date)),
    location    = factor(location),
    field_lower = pmax(0,   field_mean - field_sd),
    field_upper = pmin(100, field_mean + field_sd)
  )


# Transform data into long format 
tool_long <- rdata %>%
  select(
    site_date,
    location,
    month_str,
    gage_pctww,
    gage_min,
    gage_max,
    nhd_pctww
  ) %>%
  pivot_longer(
    cols      = c(gage_pctww, nhd_pctww),
    names_to  = "model",
    values_to = "tool_pctww"
  ) %>%
  mutate(
    model = case_match(
      model,
      "gage_pctww" ~ "Tool - gage-based",
      "nhd_pctww"  ~ "Tool - NHDPlus"
    ),
    raw_min = ifelse(model == "Tool - gage-based", gage_min, NA_real_),
    raw_max = ifelse(model == "Tool - gage-based", gage_max, NA_real_),
    gage_lower = dplyr::case_when(
      !is.na(raw_min) & !is.na(raw_max) ~ raw_min,
      !is.na(raw_min) &  is.na(raw_max) ~ raw_min,
      is.na(raw_min)  & !is.na(raw_max) ~ tool_pctww,
      TRUE                               ~ NA_real_
    ),
    gage_upper = dplyr::case_when(
      !is.na(raw_min) & !is.na(raw_max) ~ raw_max,
      !is.na(raw_min) &  is.na(raw_max) ~ tool_pctww,
      is.na(raw_min)  & !is.na(raw_max) ~ raw_max,
      TRUE                               ~ NA_real_
    ),
    gage_lower = ifelse(model == "Tool - NHDPlus", tool_pctww, gage_lower),
    gage_upper = ifelse(model == "Tool - NHDPlus", tool_pctww, gage_upper)
  )

# Field bars + tool points
pd <- position_dodge(width = 0.5)

# Gage min-max ranges, clipped to 0-100
gage_range <- tool_long %>%
  filter(model == "Tool - gage-based", !is.na(gage_lower), !is.na(gage_upper)) %>%
  mutate(lo     = pmax(0,   gage_lower),
         hi     = pmin(100, gage_upper),
         capped = hi >= 99.99)

p_site_date <- ggplot() +
  geom_col(
    data  = rdata,
    aes(x = site_date, y = field_mean),
    width  = 0.7,
    fill   = "grey75",
    colour = "grey40"
  ) +
  geom_errorbar(
    data = rdata,
    aes(x = site_date, ymin = field_lower, ymax = field_upper),
    width  = 0.2,
    colour = "black"
  ) +
  # Gage min-max range
  # solid with caps = full range below the 100% cap
  geom_errorbar(
    data = gage_range %>% filter(!capped),
    aes(x = site_date, ymin = lo, ymax = hi, colour = model),
    width     = 0.15,
    linewidth = 0.7,
    alpha     = 0.55,
    position  = position_nudge(x = -0.125)
  ) +
  # bottom cap tick for the capped ranges (matches the solid ranges)
  geom_errorbar(
    data = gage_range %>% filter(capped),
    aes(x = site_date, ymin = lo, ymax = lo, colour = model),
    width     = 0.15,
    linewidth = 0.7,
    alpha     = 0.55,
    position  = position_nudge(x = -0.125)
  ) +
  # dashed + arrow = range reaches the 100% cap 
  geom_segment(
    data = gage_range %>% filter(capped),
    aes(x = site_date, xend = site_date, y = lo, yend = hi, colour = model),
    linewidth = 0.7,
    alpha     = 0.55,
    linetype  = "dashed",
    arrow     = arrow(length = unit(0.12, "in"), type = "open", ends = "last"),
    position  = position_nudge(x = -0.125)
  ) +
  geom_point(
    data = tool_long,
    aes(
      x      = site_date,
      y      = tool_pctww,
      shape  = model,
      colour = model
    ),
    size     = 4.2,
    position = pd
  ) +
  scale_colour_manual(
    values = c(
      "Tool - gage-based" = "#1b9e77",
      "Tool - NHDPlus"    = "#7570b3"
    ),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(
      "Tool - gage-based" = 16,
      "Tool - NHDPlus"    = 17
    ),
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x        = "Site ID and Sample Date",
    y        = "% Wastewater in Stream (%WWef)",
    title    = "Field-Measured vs Tool-Estimated Wastewater Contribution",
    subtitle = paste0("Grey bars = field mean (±1 SD); green = gage-based with min-max range; purple = NHDPlus\n",
                      "Dashed arrow = gage range reaches the 100% cap")
  ) +
  theme_ms() +
  theme(
    axis.text.x = element_text(
      angle = 60, hjust = 1, vjust = 1, size = 12
    ),
    legend.position = "top",
    legend.direction = "horizontal"
  )

print(p_site_date)
save_fig(p_site_date, "fig_field_vs_tool_bars", w = 16, h = 9)

# 1:1 plots with season
axis_max <- 100

# Fix Excel serial dates and add season
all_locations <- sort(unique(rdata$location))

rdata <- rdata %>%
  mutate(
    Month = ifelse(grepl("^[0-9]+$", Month),
                   format(as.Date(as.numeric(Month), origin = "1899-12-30"), "%b"),
                   Month),
    season = case_when(
      grepl("^Dec|^Jan|^Feb", Month, ignore.case = TRUE) ~ "Winter",
      grepl("^Mar|^Apr|^May", Month, ignore.case = TRUE) ~ "Spring",
      grepl("^Jun|^Jul|^Aug", Month, ignore.case = TRUE) ~ "Summer",
      grepl("^Sep|^Oct|^Nov", Month, ignore.case = TRUE) ~ "Fall"
    ),
    season   = factor(season, levels = c("Spring", "Summer", "Fall", "Winter")),
    location = factor(location, levels = all_locations)
  )

# Color and shape scales
season_colors <- c(
  "Spring" = "#4daf4a",
  "Summer" = "#ff7f00",
  "Fall"   = "#a65628",
  "Winter" = "#377eb8"
)

shape_vals <- c(16, 17, 15, 18, 8, 10, 11, 12)[seq_along(all_locations)]
names(shape_vals) <- all_locations

# Shared scales as functions
shared_colour <- function() scale_colour_manual(values = season_colors, name = "Season",   drop = FALSE)
shared_shape  <- function() scale_shape_manual(values  = shape_vals,    name = "Location", drop = FALSE)

# build one 1:1 panel
make_1to1_plot <- function(df, x_col, y_col,
                           x_label, y_label, title_str,
                           through_zero = TRUE) {
  
  df$location <- factor(df$location, levels = all_locations)
  df$season   <- factor(df$season,   levels = c("Spring", "Summer", "Fall", "Winter"))
  
  if (through_zero) {
    fit      <- lm(df[[y_col]] ~ df[[x_col]] + 0)
    lm_label <- sprintf("y = %.2fx  (through 0)\nR² = %.2f",
                        coef(fit)[1], summary(fit)$r.squared)
    smooth_formula <- y ~ x + 0
  } else {
    fit      <- lm(df[[y_col]] ~ df[[x_col]])
    lm_label <- sprintf("y = %.2fx %s %.1f\nR² = %.2f",
                        coef(fit)[2], ifelse(coef(fit)[1] >= 0, "+", "-"),
                        abs(coef(fit)[1]), summary(fit)$r.squared)
    smooth_formula <- y ~ x
  }
  
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 colour = season, shape = location)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey40", linewidth = 1.1) +
    geom_smooth(method      = "lm",
                formula     = smooth_formula,
                se          = TRUE,
                colour      = "grey30",
                fill        = "grey70",
                alpha       = 0.15,
                linewidth   = 1.4,
                inherit.aes = FALSE,
                aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_errorbarh(
      aes(xmin = field_lower, xmax = field_upper),
      height    = 0,
      linewidth = 0.8,
      alpha     = 0.6
    ) +
    geom_point(size = PT_SIZE, alpha = 0.9) +
    annotate("text",
             x = 0.02 * axis_max, y = 0.98 * axis_max,
             label    = lm_label,
             hjust    = 0, vjust = 1,
             size     = GEOM_TXT,
             fontface = "italic", colour = "black") +
    shared_colour() +
    shared_shape() +
    scale_x_continuous(limits = c(0, axis_max),
                       breaks = seq(0, axis_max, 20),
                       labels = function(x) paste0(x, "%")) +
    scale_y_continuous(limits = c(0, axis_max),
                       breaks = seq(0, axis_max, 20),
                       labels = function(x) paste0(x, "%")) +
    coord_fixed() +
    labs(x = x_label, y = y_label, title = title_str) +
    theme_ms()
}

build_1to1_fig <- function(through_zero, fit_type, file_name) {
  # Panel A: Gage-based vs Field
  p_1to1_gage <- make_1to1_plot(
    df        = rdata %>% filter(!is.na(field_mean) & !is.na(gage_pctww)),
    x_col     = "field_mean",
    y_col     = "gage_pctww",
    x_label   = "Field-Measured %WW (mean ± 1 SD)",
    y_label   = "Tool-Estimated %WW (Gage-based)",
    title_str = "A: Gage-based vs. Field",
    through_zero = through_zero
  )
  
  # Panel B: NHDPlus vs Field
  p_1to1_nhd <- make_1to1_plot(
    df        = rdata %>% filter(!is.na(field_mean) & !is.na(nhd_pctww)),
    x_col     = "field_mean",
    y_col     = "nhd_pctww",
    x_label   = "Field-Measured %WW (mean ± 1 SD)",
    y_label   = "Tool-Estimated %WW (NHDPlus)",
    title_str = "B: NHDPlus vs. Field",
    through_zero = through_zero
  )
  
  # Extract legend from one panel only
  legend <- get_legend(
    p_1to1_nhd + theme(legend.position = "right")
  )
  
  # Strip legends from both panels
  p_gage_noleg <- p_1to1_gage + theme(legend.position = "none")
  p_nhd_noleg  <- p_1to1_nhd  + theme(legend.position = "none")
  
  # Combine with cowplot
  p_1to1 <- plot_grid(
    plot_grid(p_gage_noleg, p_nhd_noleg, nrow = 1, labels = NULL),
    legend,
    nrow = 1,
    rel_widths = c(1, 0.15)
  )
  
  # Add title
  title <- ggdraw() +
    draw_label(
      paste0("Field-Measured vs. Tool-Estimated Wastewater Contribution (", fit_type, ")"),
      fontface = "bold", size = TXT_TITLE, hjust = 0.5
    )
  subtitle <- ggdraw() +
    draw_label(
      "Horizontal error bars = ±1 SD of field mean; dashed line = 1:1 agreement",
      size = TXT_SUB, hjust = 0.5
    )
  
  p_1to1_final <- plot_grid(title, subtitle, p_1to1,
                            ncol = 1, rel_heights = c(0.07, 0.05, 1))
  
  print(p_1to1_final)
  save_fig(p_1to1_final, file_name, w = 17, h = 8.5)
  invisible(p_1to1_final)
}

# Constrained (through 0) and unconstrained (free intercept) versions
p_1to1_zero <- build_1to1_fig(TRUE,  "regression through 0",  "fig_field_vs_tool_1to1_through0")
p_1to1_free <- build_1to1_fig(FALSE, "regression with intercept", "fig_field_vs_tool_1to1_free")

# Deviations from 1:1
make_resid_plot <- function(df, x_col, y_col,
                            x_label, y_label, title_str, y_lim) {
  
  df$location <- factor(df$location, levels = all_locations)
  df$season   <- factor(df$season,   levels = c("Spring", "Summer", "Fall", "Winter"))
  
  df$resid <- df[[y_col]] - df[[x_col]]  # deviation from 1:1 line
  
  ggplot(df, aes(x = .data[[x_col]], y = resid,
                 colour = season, shape = location)) +
    geom_hline(yintercept = 0,
               linetype = "dashed", colour = "grey40", linewidth = 1.1) +
    geom_point(size = PT_SIZE, alpha = 0.9) +
    shared_colour() +
    shared_shape() +
    scale_x_continuous(limits = c(0, axis_max),
                       breaks = seq(0, axis_max, 20),
                       labels = function(x) paste0(x, "%")) +
    scale_y_continuous(
      limits = c(-y_lim, y_lim),
      breaks = seq(-y_lim, y_lim, by = 10),
      labels = function(x) paste0(x, "%")
    ) +
    labs(x = x_label, y = y_label, title = title_str) +
    theme_ms()
}

# Shared y-axis limit from 1:1 deviations
gage_dev  <- rdata %>% filter(!is.na(field_mean) & !is.na(gage_pctww)) %>% mutate(dev = gage_pctww - field_mean) %>% pull(dev)
nhd_dev   <- rdata %>% filter(!is.na(field_mean) & !is.na(nhd_pctww))  %>% mutate(dev = nhd_pctww  - field_mean) %>% pull(dev)
resid_lim <- ceiling(max(abs(c(gage_dev, nhd_dev))) / 10) * 10

# Panel A: Gage-based deviations
p_resid_gage <- make_resid_plot(
  df        = rdata %>% filter(!is.na(field_mean) & !is.na(gage_pctww)),
  x_col     = "field_mean",
  y_col     = "gage_pctww",
  x_label   = "Field-Measured %WW (mean ± 1 SD)",
  y_label   = "Deviation from 1:1 (%WW)",
  title_str = "A: Gage-based Deviations from 1:1",
  y_lim     = resid_lim
)

# Panel B: NHDPlus deviations
p_resid_nhd <- make_resid_plot(
  df        = rdata %>% filter(!is.na(field_mean) & !is.na(nhd_pctww)),
  x_col     = "field_mean",
  y_col     = "nhd_pctww",
  x_label   = "Field-Measured %WW (mean ± 1 SD)",
  y_label   = "Deviation from 1:1 (%WW)",
  title_str = "B: NHDPlus Deviations from 1:1",
  y_lim     = resid_lim
)

# Legend from panel B 
legend_resid <- get_legend(
  p_resid_nhd + theme(legend.position = "right")
)

p_resid_gage_noleg <- p_resid_gage + theme(legend.position = "none")
p_resid_nhd_noleg  <- p_resid_nhd  + theme(legend.position = "none")

p_resid <- plot_grid(
  plot_grid(p_resid_gage_noleg, p_resid_nhd_noleg, nrow = 1, labels = NULL),
  legend_resid,
  nrow = 1,
  rel_widths = c(1, 0.15)
)

title_resid <- ggdraw() +
  draw_label(
    "Deviation from 1:1 Agreement: Tool vs. Field-Measured %WW",
    fontface = "bold", size = TXT_TITLE, hjust = 0.5
  )
subtitle_resid <- ggdraw() +
  draw_label(
    "Deviation = tool estimate minus field measurement; dashed line = perfect agreement",
    size = TXT_SUB, hjust = 0.5
  )

p_resid_final <- plot_grid(title_resid, subtitle_resid, p_resid,
                           ncol = 1, rel_heights = c(0.07, 0.05, 1))

print(p_resid_final)
save_fig(p_resid_final, "fig_field_vs_tool_deviations", w = 17, h = 7.5)


# Agreement metrics (tool vs. field)
agreement_metrics <- function(df, pred_col, label) {
  d   <- df %>% filter(!is.na(field_mean), !is.na(.data[[pred_col]]))
  obs <- d$field_mean
  sim <- d[[pred_col]]
  fit_free <- lm(sim ~ obs)
  fit_zero <- lm(sim ~ obs + 0)
  tibble(
    tool           = label,
    n              = nrow(d),
    MAE_pp         = mean(abs(sim - obs)),
    RMSE_pp        = sqrt(mean((sim - obs)^2)),
    mean_error_pp  = mean(sim - obs),
    PBIAS_pct      = 100 * sum(sim - obs) / sum(obs),
    slope_free     = unname(coef(fit_free)[2]),
    intercept_free = unname(coef(fit_free)[1]),
    R2_free        = summary(fit_free)$r.squared,
    slope_zero     = unname(coef(fit_zero)[1]),
    R2_zero        = summary(fit_zero)$r.squared,
    n_with_sd      = sum(!is.na(d$field_lower) & !is.na(d$field_upper)),
    within_1SD     = sum(sim >= d$field_lower & sim <= d$field_upper, na.rm = TRUE)
  )
}

metrics <- bind_rows(
  agreement_metrics(rdata, "gage_pctww", "Gage-based"),
  agreement_metrics(rdata, "nhd_pctww",  "NHDPlus-based")
) %>% mutate(across(where(is.numeric), ~ round(.x, 2)))
print(as.data.frame(metrics))

# Bias by season (is there a seasonal pattern in the errors?)
season_bias <- bind_rows(
  rdata %>% filter(!is.na(field_mean), !is.na(gage_pctww)) %>%
    group_by(season) %>%
    summarise(tool = "Gage-based", n = n(),
              mean_error_pp = mean(gage_pctww - field_mean),
              MAE_pp = mean(abs(gage_pctww - field_mean)), .groups = "drop"),
  rdata %>% filter(!is.na(field_mean), !is.na(nhd_pctww)) %>%
    group_by(season) %>%
    summarise(tool = "NHDPlus-based", n = n(),
              mean_error_pp = mean(nhd_pctww - field_mean),
              MAE_pp = mean(abs(nhd_pctww - field_mean)), .groups = "drop")
) %>% mutate(across(where(is.numeric), ~ round(.x, 1)))
print(as.data.frame(season_bias))

# Gage range checks: does the gage min-max range contain the field mean? how many hit the 100% cap?
gage_range_check <- rdata %>%
  filter(!is.na(field_mean), !is.na(gage_pctww)) %>%
  summarise(n = n(),
            field_in_gage_range = sum(field_mean >= gage_min & field_mean <= gage_max, na.rm = TRUE),
            gage_range_capped   = sum(gage_max >= 99.99, na.rm = TRUE),
            gage_range_missing  = sum(is.na(gage_min) | is.na(gage_max)))
print(as.data.frame(gage_range_check))

write.csv(metrics,      file.path(out_dir, "tool_vs_field_metrics.csv"),      row.names = FALSE)
write.csv(season_bias,  file.path(out_dir, "tool_vs_field_season_bias.csv"),  row.names = FALSE)

# Nitrate & conductivity vs. average %WW (1:1 style plots)

nitrate_cond <- read_excel(xlsx_path, sheet = "Nitrate_CondVEstimate") %>%
  rename(
    site        = Site,
    month_year  = `Mon. + Year`,
    cond_pct    = `Cond.`,
    nitrate_pct = Nitrate,
    avg_pct     = `Average (%)`
  ) %>%
  mutate(
    # Excel % cells come in as decimals (0.23); multiply by 100
    across(c(cond_pct, nitrate_pct, avg_pct),
           ~ {
             x <- as.numeric(.x)
             # If all non-NA values are <= 1, assume 0-1 scale and convert
             if (all(x[!is.na(x)] <= 1)) x * 100 else x
           }),
    site = factor(site)
  )

# Shared axis max
nc_axis_max <- 100

# Fit stats 
fit_label <- function(df, x_col, y_col) {
  fit   <- lm(df[[y_col]] ~ df[[x_col]])
  slope <- coef(fit)[2]
  inter <- coef(fit)[1]
  r2    <- summary(fit)$r.squared
  sign_str <- ifelse(inter >= 0, "+", "-")
  sprintf("y = %.2fx %s %.1f\nR² = %.2f", slope, sign_str, abs(inter), r2)
}

# Site color/shape scales
site_levels <- levels(nitrate_cond$site)
nc_shapes   <- c(16, 17, 15, 18, 8, 10)[seq_along(site_levels)]
names(nc_shapes) <- site_levels

nc_colors <- RColorBrewer::brewer.pal(max(3, length(site_levels)), "Dark2")[seq_along(site_levels)]
names(nc_colors) <- site_levels

# Panel A: Conductivity vs Average
lab_cond <- fit_label(nitrate_cond, "cond_pct", "avg_pct")

p_cond_vs_avg <- ggplot(nitrate_cond,
                        aes(x = cond_pct, y = avg_pct,
                            colour = site, shape = site)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40", linewidth = 1.1) +
  geom_smooth(method = "lm", formula = y ~ x,
              se = TRUE, colour = "grey30", fill = "grey70",
              alpha = 0.15, linewidth = 1.4,
              inherit.aes = FALSE,
              aes(x = cond_pct, y = avg_pct),
              data = nitrate_cond) +
  geom_point(size = PT_SIZE, alpha = 0.9) +
  annotate("text",
           x = 0.02 * nc_axis_max, y = 0.98 * nc_axis_max,
           label = lab_cond, hjust = 0, vjust = 1,
           size = GEOM_TXT, fontface = "italic", colour = "black") +
  scale_colour_manual(values = nc_colors, name = "Site") +
  scale_shape_manual(values  = nc_shapes, name = "Site") +
  scale_x_continuous(limits = c(0, nc_axis_max),
                     breaks = seq(0, nc_axis_max, 20),
                     labels = function(x) paste0(x, "%")) +
  scale_y_continuous(limits = c(0, nc_axis_max),
                     breaks = seq(0, nc_axis_max, 20),
                     labels = function(x) paste0(x, "%")) +
  coord_fixed() +
  labs(
    x = "Conductivity-Estimated %WW",
    y = "Multi-Constituent Average %WW",
    title = "A: Conductivity vs. Average %WW"
  ) +
  theme_ms()

# Panel B: Nitrate vs Average
lab_nit <- fit_label(nitrate_cond, "nitrate_pct", "avg_pct")

p_nit_vs_avg <- ggplot(nitrate_cond,
                       aes(x = nitrate_pct, y = avg_pct,
                           colour = site, shape = site)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40", linewidth = 1.1) +
  geom_smooth(method = "lm", formula = y ~ x,
              se = TRUE, colour = "grey30", fill = "grey70",
              alpha = 0.15, linewidth = 1.4,
              inherit.aes = FALSE,
              aes(x = nitrate_pct, y = avg_pct),
              data = nitrate_cond) +
  geom_point(size = PT_SIZE, alpha = 0.9) +
  annotate("text",
           x = 0.02 * nc_axis_max, y = 0.98 * nc_axis_max,
           label = lab_nit, hjust = 0, vjust = 1,
           size = GEOM_TXT, fontface = "italic", colour = "black") +
  scale_colour_manual(values = nc_colors, name = "Site") +
  scale_shape_manual(values  = nc_shapes, name = "Site") +
  scale_x_continuous(limits = c(0, nc_axis_max),
                     breaks = seq(0, nc_axis_max, 20),
                     labels = function(x) paste0(x, "%")) +
  scale_y_continuous(limits = c(0, nc_axis_max),
                     breaks = seq(0, nc_axis_max, 20),
                     labels = function(x) paste0(x, "%")) +
  coord_fixed() +
  labs(
    x = "Nitrate-Estimated %WW",
    y = "Multi-Constituent Average %WW",
    title = "B: Nitrate vs. Average %WW"
  ) +
  theme_ms()

# Combine with shared legend
legend_nc <- get_legend(p_nit_vs_avg + theme(legend.position = "right"))

p_nc_final <- plot_grid(
  plot_grid(
    p_cond_vs_avg + theme(legend.position = "none"),
    p_nit_vs_avg  + theme(legend.position = "none"),
    nrow = 1
  ),
  legend_nc,
  nrow = 1,
  rel_widths = c(1, 0.12)
)

title_nc <- ggdraw() +
  draw_label(
    "Single-Constituent Estimates vs. Multi-Constituent Average %WW",
    fontface = "bold", size = TXT_TITLE, hjust = 0.5
  )
subtitle_nc <- ggdraw() +
  draw_label(
    "Dashed line = 1:1 agreement; shaded band = 95% CI of regression",
    size = TXT_SUB, hjust = 0.5
  )

p_nc_final <- plot_grid(title_nc, subtitle_nc, p_nc_final,
                        ncol = 1, rel_heights = c(0.07, 0.05, 1))

print(p_nc_final)
save_fig(p_nc_final, "fig_nitrate_conductivity_vs_avg", w = 17, h = 8.5)


# Rural vs. urban: prevalence (shown once) + severity across 12 months
effects_sheet <- "RuralUrban_Effects"   

if (!effects_sheet %in% excel_sheets(xlsx_path)) {
  stop("Sheet '", effects_sheet, "' not found. Sheets in the workbook: ",
       paste(excel_sheets(xlsx_path), collapse = ", "),
       "\nSet effects_sheet (Section 4) to the right name.")
}
eff <- read_excel(xlsx_path, sheet = effects_sheet) %>%
  mutate(month_f = factor(as.integer(month), levels = 1:12, labels = month.abb))

src_cols <- c("Gage-based" = "#1b9e77", "NHDPlus-based" = "#7570b3")
dodge    <- position_dodge(width = 0.55)

# Local theme 
theme_eff <- function() {
  theme_ms() + theme(plot.title    = element_text(face = "bold", size = TXT_TITLE, hjust = 0),
                     plot.title.position = "plot",     
                     plot.subtitle = element_text(size = TXT_SUB, hjust = 0),
                     legend.title  = element_blank())
}

# Prevalence odds ratio, unadjusted vs. adjusted (shown once) 
prev <- eff %>%
  filter(source == "Gage-based", stat %in% c("odds_ratio_unadj", "odds_ratio_adj")) %>%
  group_by(stat) %>%
  summarise(est = median(est), lo = median(lo), hi = median(hi), .groups = "drop") %>%
  mutate(model = factor(ifelse(stat == "odds_ratio_unadj",
                               "Unadjusted", "Adjusted for\nDrainage Area"),
                        levels = c("Adjusted for\nDrainage Area", "Unadjusted")),
         label = sprintf("OR = %.2f\n(%.2f-%.2f)", est, lo, hi))

p_prev <- ggplot(prev, aes(est, model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 1) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 1.2) +
  geom_point(size = 5) +
  geom_text(aes(label = label), vjust = -0.9, size = 5.2, lineheight = 0.95) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8),
                labels = c("0.25", "0.5", "1", "2", "4", "8")) +
  coord_cartesian(ylim = c(0.6, 2.6)) +
  labs(title = "A  Prevalence",
       subtitle = "Rural vs. Urban Odds of an\nImpacted Segment",
       x = "Odds Ratio (Log Scale)\n>1 = More in Rural", y = NULL) +
  theme_eff() + theme(panel.grid.major.y = element_blank())

# Severity rank correlation, unadjusted vs. adjusted 
# Median across the 12 monthly estimates (and CI limits), gage-based rows.
sev_sum <- eff %>%
  filter(source == "Gage-based", stat %in% c("rank_corr_unadj", "rank_corr_adj")) %>%
  group_by(stat) %>%
  summarise(est = median(est), lo = median(lo), hi = median(hi), .groups = "drop") %>%
  mutate(model = factor(ifelse(stat == "rank_corr_unadj", "Unadjusted", "Adjusted for\nDrainage Area"),
                        levels = c("Adjusted for\nDrainage Area", "Unadjusted")),
         label = sprintf("rho = %.2f\n(%.2f to %.2f)", est, lo, hi))

p_sev_sum <- ggplot(sev_sum, aes(est, model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 1) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 1.2) +
  geom_point(size = 5) +
  geom_text(aes(label = label), vjust = -0.9, size = 5.2, lineheight = 0.95) +
  coord_cartesian(xlim = c(-0.6, 0.1), ylim = c(0.6, 2.6)) +
  labs(title = "B  Severity",
       subtitle = "Rank Correlation of Rurality With %WWef\nAmong Impacted Segments",
       x = "Rank Correlation (Median Across Months)\n<0 = Lower in Rural", y = NULL) +
  theme_eff() + theme(panel.grid.major.y = element_blank())

p_rural_urban <- (p_prev | p_sev_sum) +
  plot_annotation(
    title = "Rural vs. Urban Wastewater Effluent Prevalence and Severity",
    theme = theme(plot.title = element_text(face = "bold", size = TXT_TITLE + 2))
  )
print(p_rural_urban)
save_fig(p_rural_urban, "fig_rural_urban_prevalence_severity", w = 16, h = 7.5)

# Severity by month 
# rank correlation, unadjusted vs adjusted. B, C: probability urban is higher and
# median difference 
model_cols <- c("Unadjusted" = "grey35", "Adjusted for Drainage Area" = "#2a78d6")

make_severity_si <- function(src, file_name) {
  src_col <- unname(src_cols[src])
  theme_si <- function() theme_eff() +
    theme(plot.title    = element_text(face = "bold", size = TXT_TITLE - 3, hjust = 0),
          plot.subtitle = element_text(size = TXT_SUB - 2, hjust = 0),
          plot.margin   = margin(10, 30, 10, 14))
  d <- eff %>% filter(source == src)
  
  p_rc <- d %>%
    filter(stat %in% c("rank_corr_unadj", "rank_corr_adj")) %>%
    mutate(model = factor(ifelse(stat == "rank_corr_unadj", "Unadjusted", "Adjusted for Drainage Area"),
                          levels = names(model_cols))) %>%
    ggplot(aes(month_f, est, colour = model)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.9) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, linewidth = 0.9,
                  alpha = 0.75, position = dodge) +
    geom_point(size = 4, position = dodge) +
    scale_colour_manual(values = model_cols, name = NULL) +
    labs(title = "A  Rank Correlation, Unadjusted vs. Adjusted",
         subtitle = "Negative = Lower %WWef in Rural Tracts",
         x = NULL, y = "Rank Correlation") +
    theme_si() + theme(legend.position = "top")
  
  # rank-biserial r = P(rural > urban) - P(rural < urban)
  # P(urban higher) = (1 - r)/2. 
  p_prob <- d %>%
    filter(stat == "rank_biserial") %>%
    mutate(prob = 100 * (1 - est) / 2, p_lo = 100 * (1 - hi) / 2, p_hi = 100 * (1 - lo) / 2) %>%
    ggplot(aes(month_f, prob)) +
    geom_hline(yintercept = 50, colour = "grey40", linewidth = 0.9) +
    geom_errorbar(aes(ymin = p_lo, ymax = p_hi), width = 0.25, linewidth = 0.9,
                  alpha = 0.75, colour = src_col) +
    geom_point(size = 4, colour = src_col) +
    scale_y_continuous(limits = c(40, 90), breaks = seq(40, 90, 10),
                       labels = function(x) paste0(x, "%")) +
    labs(title = "B  How Often Urban Is Higher (Unadjusted)",
         subtitle = "Chance a Random Urban Tract Has Higher %WWef\nThan a Rural One (50% = No Difference)",
         x = NULL, y = "Probability Urban Is Higher") +
    theme_si()
  
  p_md <- d %>%
    filter(stat == "median_diff_pp") %>%
    ggplot(aes(month_f, est)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.9) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, linewidth = 0.9,
                  alpha = 0.75, colour = src_col) +
    geom_point(size = 4, colour = src_col) +
    labs(title = "C  Size of Difference (Unadjusted)",
         subtitle = "Rural Minus Urban Median %WWef\nAmong Impacted Segments",
         x = NULL, y = "Median Difference (% Points)") +
    theme_si()
  
  fig <- (p_rc / p_prob / p_md) +
    plot_annotation(
      title = paste0("Rural vs. Urban Wastewater Effluent\nSeverity by Month (", src, ")"),
      theme = theme(plot.title = element_text(face = "bold", size = TXT_TITLE, hjust = 0),
                    plot.margin = margin(10, 10, 10, 10)))
  print(fig)
  save_fig(fig, file_name, w = 12, h = 17)
  invisible(fig)
}
make_severity_si("Gage-based",    "figS_rural_urban_severity_Gage")
make_severity_si("NHDPlus-based", "figS_rural_urban_severity_NHDPlus")

# Prevalence odds ratio by month 
prev_month <- eff %>%
  filter(source == "Gage-based", stat %in% c("odds_ratio_unadj", "odds_ratio_adj")) %>%
  mutate(model = factor(ifelse(stat == "odds_ratio_unadj", "Unadjusted", "Adjusted for Drainage Area"),
                        levels = c("Unadjusted", "Adjusted for Drainage Area")))

p_prev_month <- ggplot(prev_month, aes(month_f, est, colour = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 1) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, linewidth = 0.9,
                alpha = 0.7, position = dodge) +
  geom_point(size = 4, position = dodge) +
  scale_colour_manual(values = c("Unadjusted" = "grey35", "Adjusted for Drainage Area" = "#2a78d6"),
                      name = NULL) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8),
                labels = c("0.25", "0.5", "1", "2", "4", "8")) +
  labs(title = "Prevalence of Wastewater Effluent in\nRural vs. Urban Tracts, by Month",
       subtitle = "Odds a rural tract contains an impacted segment,\nrelative to an urban tract (>1 = more in rural)",
       x = NULL, y = "Odds Ratio (Log Scale)") +
  theme_ms() +
  theme(plot.title    = element_text(face = "bold", size = TXT_TITLE - 2, hjust = 0),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = TXT_SUB - 2, hjust = 0),
        legend.position = "top")

print(p_prev_month)
save_fig(p_prev_month, "figS_prevalence_by_month", w = 12, h = 7.5)


# Equity: severity associations 

eq_sheet <- "Equity_Associations"   # equity_associations.csv
sb_sheet <- "SpatialBlock_Check"    # spatial_block_check.csv

MIN_MONTHS      <- 9          # "most months" = at least this many of 12
TOP_N           <- 15         # indicators shown (ranked by |median adjusted association|, gage-based)
ALPHA           <- 0.05       # BH q threshold
EXCLUDE_PATTERN <- "wastewater|percentile|rurality"  

for (s in c(eq_sheet, sb_sheet)) {
  if (!s %in% excel_sheets(xlsx_path)) {
    stop("Sheet '", s, "' not found. Sheets in the workbook: ",
         paste(excel_sheets(xlsx_path), collapse = ", "))
  }
}

eqa <- read_excel(xlsx_path, sheet = eq_sheet) %>%
  mutate(month = as.integer(month)) %>%
  filter(mode == "Severity",
         !grepl(EXCLUDE_PATTERN, indicator, ignore.case = TRUE),
         is.finite(effect_adj))

sbc <- read_excel(xlsx_path, sheet = sb_sheet) %>%
  mutate(month = as.integer(month), still_sig = as.logical(still_sig)) %>%
  filter(mode == "Severity",
         !grepl(EXCLUDE_PATTERN, indicator, ignore.case = TRUE))

# One row per source x indicator, summarised across months
eq_sum <- eqa %>%
  group_by(source, indicator) %>%
  summarise(
    med_adj    = median(effect_adj),
    min_adj    = min(effect_adj),
    max_adj    = max(effect_adj),
    med_unadj  = median(effect_unadj, na.rm = TRUE),
    sig_months = sum(q_adj < ALPHA & sign(effect_adj) == sign(median(effect_adj)), na.rm = TRUE),
    .groups = "drop"
  )

# Robust to the spatial checks: still significant in most months under BOTH
# the block-fixed-effects and the Denver-excluded runs 
sb_sum <- sbc %>%
  group_by(source, indicator, spec) %>%
  summarise(ok = sum(still_sig, na.rm = TRUE), .groups = "drop") %>%
  group_by(source, indicator) %>%
  summarise(spatial_checks_hold = n_distinct(spec) == 2 && all(ok >= MIN_MONTHS), .groups = "drop")

# Theme 
theme_levels <- c("Air pollution & traffic",
                  "Facilities & land use",
                  "Natural hazard risk",
                  "Demographics",
                  "Other burdens")

assign_theme <- function(x) {
  case_when(
    grepl("diesel|pm2|pm 2|particulate|ozone|traffic|air toxics", x, ignore.case = TRUE) ~ theme_levels[1],
    grepl("flood|wildfire|expected .*loss|natural hazard", x, ignore.case = TRUE)         ~ theme_levels[3],
    grepl("hazardous|risk management|rmp|superfund|impervious|cropland|underground|lead paint|facilit|proximity", x, ignore.case = TRUE) ~ theme_levels[2],
    grepl("^percent|^share of (people|persons|residents)|age |race|hispanic|asian|white|black|language|education", x, ignore.case = TRUE) ~ theme_levels[4],
    TRUE ~ theme_levels[5]
  )
}

eq_sum <- eq_sum %>%
  left_join(sb_sum, by = c("source", "indicator")) %>%
  mutate(spatial_checks_hold = coalesce(spatial_checks_hold, FALSE),
         consistent = sig_months >= MIN_MONTHS,
         theme = factor(assign_theme(indicator), levels = theme_levels))

# All indicators, both denominators 
si_tbl <- eq_sum %>%
  arrange(source, theme, desc(abs(med_adj))) %>%
  transmute(source, theme, indicator,
            median_rho_adj = round(med_adj, 3),
            min_rho_adj    = round(min_adj, 3),
            max_rho_adj    = round(max_adj, 3),
            median_rho_unadj = round(med_unadj, 3),
            months_significant = sig_months,
            holds_in_spatial_checks = spatial_checks_hold)
write.csv(si_tbl, file.path(out_dir, "tableS_equity_associations_all_indicators.csv"), row.names = FALSE)

cat("\nThemes assigned (gage-based, shown indicators):\n")

# Figure: gage-based only, top indicators, grouped by theme 
wrap_lab <- function(x, w = 38) vapply(x, function(z) paste(strwrap(z, w), collapse = "\n"), "")

eq_gage <- eq_sum %>% filter(source == "Gage-based")
keep <- eq_gage %>% slice_max(abs(med_adj), n = TOP_N, with_ties = FALSE) %>% pull(indicator)

eq_plot <- eq_gage %>%
  filter(indicator %in% keep) %>%
  droplevels() %>%
  group_by(theme) %>%
  arrange(med_adj, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(indicator = factor(wrap_lab(indicator), levels = unique(wrap_lab(indicator))),
         consistent_lab = factor(ifelse(consistent,
                                        "Significant in most months",
                                        "Not consistently significant"),
                                 levels = c("Significant in most months", "Not consistently significant")))

print(eq_plot %>% select(theme, indicator, med_adj, sig_months, spatial_checks_hold) %>% arrange(theme, med_adj), n = Inf)

gage_col <- unname(src_cols["Gage-based"])

p_equity <- ggplot(eq_plot, aes(med_adj, indicator)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.9) +
  geom_linerange(aes(xmin = min_adj, xmax = max_adj), orientation = "y",
                 colour = gage_col, linewidth = 1.2, alpha = 0.6) +
  geom_point(aes(shape = consistent_lab), colour = gage_col, fill = "white",
             size = 5, stroke = 1.5) +
  scale_shape_manual(values = c("Significant in most months" = 16,
                                "Not consistently significant" = 1),
                     limits = c("Significant in most months", "Not consistently significant"),
                     name = NULL, drop = FALSE) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    title    = "Screening-Level Associations Between Tract Indicators\nand Wastewater Effluent Severity",
    x = "Association With %WWef (Rank Correlation, Adjusted for Drainage Area)",
    y = NULL
  ) +
  theme_ms() +
  theme(plot.title    = element_text(face = "bold", size = TXT_TITLE - 2, hjust = 0),
        plot.title.position = "plot",     
        axis.text.y   = element_text(size = TXT_AXIS - 2, lineheight = 0.9),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = TXT_AXIS - 2),
        strip.background  = element_blank(),
        panel.spacing.y   = unit(0.6, "lines"),
        # legend only needed if some plotted indicators are not consistently significant
        legend.position = if (all(eq_plot$consistent)) "none" else "bottom",
        legend.title = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92"))

print(p_equity)
save_fig(p_equity, "fig_equity_severity_associations", w = 15, h = 11)