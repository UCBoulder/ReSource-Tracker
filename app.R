# --- Libraries ---------------------------------------------------------------
library(tidyverse)
library(readxl)
library(leaflet)
library(sf)
library(shiny)
library(tigris)
library(viridisLite)
library(mapgl)

setwd('/Users/kylieboenischoakes/Documents/ReSource-Tracker')
# --- Relative Paths ----------------------------------------------------------
data_dir <- "data-raw"

paths <- list(
  dwtp          = file.path(data_dir, "filtereddwtpCOFinal.csv"),
  cejst         = file.path(data_dir, "ColoradoCEJSTData.xlsx"),
  matched_dmr   = file.path(data_dir, "matchedDMRwwtpupdated6.csv"),
  wwtp          = file.path(data_dir, "PermittedPOTW.xlsx"),
  flowlines     = file.path(data_dir, "Flowlines_Colorado_Filtered.gpkg"),
  monthly_flows = file.path(data_dir, "CO_Seasonal_Flow_Summary_Feb2026.csv"),
  gages_sf      = file.path(data_dir, "processed_gages_sf_combinedV4.gpkg"), 
  ur10          = file.path(data_dir, "CO_2010_URBANRURAL_TRACTS.csv")
)

# ---- Required file checks ---------------------------------------------------
check_required_file <- function(path, label, url = NULL) {
  if (!file.exists(path)) {
    msg <- paste0(
      "Missing required data file: ", label, "\n",
      "Expected at:\n  ", normalizePath(path, winslash = "/", mustWork = FALSE), "\n\n"
    )
    if (!is.null(url)) {
      msg <- paste0(
        msg,
        "Download it from:\n  ", url, "\n\n",
        "Save it at the path above and then rerun the app.\n"
      )
    }
    stop(msg, call. = FALSE)
  }
}

# Flowlines geopackage: too large to store on repo, provided via GitHub release
check_required_file(
  paths$flowlines,
  "Flowlines_Colorado_Filtered.gpkg",
  "https://github.com/sheldonmasters/Masters-Research-Group/releases/tag/v0.1.0"
)

# --- Constants ---------------------------------------------------------------
mon_codes    <- sprintf("%02d", 1:12)          # "01"–"12"
month_labels <- setNames(month.name, mon_codes)  # names "01".."12" → "January"...
month_colors <- setNames(
  c(
    "#3B4F8F", # 01 Jan  
    "#4B6CB7", # 02 Feb  
    "#3C9D82", # 03 Mar  
    "#54B887", # 04 Apr  
    "#7BC96F", # 05 May  
    "#F4D35E", # 06 Jun  
    "#F6B35E", # 07 Jul  
    "#F28E52", # 08 Aug 
    "#E76F51", # 09 Sep  
    "#D95F5F", # 10 Oct  
    "#A8558C", # 11 Nov  
    "#6D4C9B"  # 12 Dec  
  ),
  mon_codes
)
# Consistent colors for rurality categories (Urban / Rural)
rurality_colors <- c(
  "Urban"    = "#1b9e77",
  "Rural" = "#7570b3"
)

# Columns that need to be adjusted to whole number from decimal percent
proportion_cols <- c(
  "Percent Black or African American alone",
  "Percent American Indian / Alaska Native",
  "Percent Asian",
  "Percent Native Hawaiian or Pacific",
  "Percent two or more races",
  "Percent White",
  "Percent Hispanic or Latino",
  "Percent other races",
  "Percent age under 10",
  "Percent age 10 to 64",
  "Percent age over 64",
  "Adjusted percent of individuals below 200% Federal Poverty Line",
  "Share of homes with no kitchen or indoor plumbing (percent)")

# Similarity gate for gage propagation, propagation is only allowed when
# |Ai -Aj|/Aj < 0.5, accounting for differences in drainage area 
SIMILARITY_THRESHOLD = 0.5 

# --- Load Datasets -----------------------------------------------------------

# Drinking water treatment plants
geocoded_dwtp <- read.csv(paths$dwtp) %>%
  rename(
    `PWS Name`          = PWS.Name,
    `Population Served` = Population.Served,
    `Primary Source`    = Primary.Source
  )

# CEJST data
cejst <- read_excel(paths$cejst) %>%
  rename(
    GEOID          = `Census tract 2010 ID`,
    wastewater_pct = `Wastewater discharge (percentile)`,
    low_income     = `Is low income?`
  ) %>%
  mutate(
    GEOID  = sprintf("%011.0f", GEOID),
    status = case_when(
      wastewater_pct > 90 & low_income ~ "Both",
      wastewater_pct > 90              ~ "High Wastewater",
      low_income                       ~ "Low Income",
      TRUE                             ~ NA_character_
    )
  )

# ------------------ 2010 TRACTS + UR10-BASED RURALITY -----------------------
options(tigris_use_cache = TRUE)

# 2010 tracts for Colorado (geometry + GEOID10 key)
tracts10 <- tracts(state = "08", year = 2010, class = "sf") %>%
  st_make_valid() %>%
  st_transform(5070) %>%
  mutate(
    GEOID10 = sprintf("%011s", as.character(GEOID10))
  )

# UR10 rurality table (precomputed)
ur10 <- readr::read_csv(paths$ur10, show_col_types = FALSE) %>%
  mutate(
    GEOID10  = sprintf("%011s", as.character(GEOID)),
    rurality = as.character(rurality)
  ) %>%
  select(GEOID10, rurality)

# Attach Census-based rurality back to CEJST table
cejst <- cejst %>%
  mutate(GEOID10 = sprintf("%011s", as.character(GEOID))) %>%
  left_join(ur10, by = "GEOID10") %>%
  select(-GEOID10)

# 2010 tracts + CEJST (including rurality), in WGS84 for the app
tracts_sf <- tracts10 %>%
  st_transform(4326) %>%
  left_join(
    cejst %>%
      mutate(GEOID10 = sprintf("%011s", as.character(GEOID))) %>%
      select(GEOID10, GEOID, everything()),
    by = "GEOID10"
  )

# ------------------------- DMR (WW) mean ± SD (monthly) ----------------------
# DMR monthly flows matched to NHDPlus COMIDs
matched_dmr <- read_csv(
  paths$matched_dmr,
  col_types = cols(.default = col_guess(), geometry = col_skip(), `...1` = col_skip()),
  show_col_types = FALSE
)

dmr_monthly <- matched_dmr %>%
  mutate(
    matched_COMID = as.character(matched_COMID),
    mon_code      = sprintf("%02d", as.integer(mon_code))
  ) %>%
  filter(mon_code %in% mon_codes)

# 1) Average within Facility–COMID–month to collapse duplicate outfalls
# (mean across outfalls; combine SDs conservatively by averaging SD field)
dmr_monthly_fac <- dmr_monthly %>%
  group_by(Facility, matched_COMID, mon_code) %>%
  summarise(
    mean_dmr_cfs_fac = mean(mean_dmr_cfs, na.rm = TRUE),
    sd_dmr_cfs_fac   = mean(sd_dmr_cfs,   na.rm = TRUE),
    .groups = "drop"
  )

# 2) Sum across facilities per COMID-month
# carry variance (sd^2) to accumulate uncertainty downstream
dmr_monthly_comid <- dmr_monthly_fac %>%
  mutate(var_dmr_cfs_fac = sd_dmr_cfs_fac^2) %>%
  group_by(matched_COMID, mon_code) %>%
  summarise(
    WW_mean = sum(mean_dmr_cfs_fac, na.rm = TRUE),
    WW_var  = sum(var_dmr_cfs_fac,  na.rm = TRUE),
    .groups = "drop"
  )

# Wide format: keep mean + variance 
dmr_monthly_wide <- dmr_monthly_comid %>%
  pivot_wider(
    names_from  = mon_code,
    values_from = c(WW_mean, WW_var),
    names_sep   = "_"
  )

# --- WWTP design data (for popups, markers) ----------------------------------
wwtp_df <- read_excel(paths$wwtp) %>%
  filter(!is.na(LATITUDE83), !is.na(LONGITUDE83)) %>%
  rename(
    Facility   = FACILITY_NAME,
    Latitude   = LATITUDE83,
    Longitude  = LONGITUDE83,
    DesignFlow = TOTAL_DESIGN_FLOW_NMBR
  )

# Flowlines
flowlines <- st_read(paths$flowlines)

# --- WWTP to Flowline Mapping ------------------------------------------------
wwtp_pts <- st_as_sf(wwtp_df, coords = c("Longitude", "Latitude"), crs = 4326)
wwtp_pts <- st_transform(wwtp_pts, st_crs(flowlines))
wwtp_pts$matched_COMID <- flowlines$COMID[st_nearest_feature(wwtp_pts, flowlines)]

# Join:
#   (1) DMR monthly WW mean + variance (COMID-level)
#   (2) a single representative Facility name per COMID (for contributor labels)
flowlines <- flowlines %>%
  mutate(COMID = as.character(COMID)) %>%
  left_join(dmr_monthly_wide, by = c("COMID" = "matched_COMID")) %>%
  left_join(
    wwtp_pts %>%
      st_drop_geometry() %>%
      distinct(matched_COMID, .keep_all = TRUE) %>%
      select(Facility, matched_COMID) %>%
      mutate(matched_COMID = as.character(matched_COMID)),
    by = c("COMID" = "matched_COMID")
  )

# --- Downstream WW Flow & Contributor Accumulation (monthly) -----------------
flowlines$Hydroseq   <- as.character(flowlines$Hydroseq)
flowlines$DnHydroseq <- as.character(flowlines$DnHydroseq)

# Reduced routing table: Hydroseq, DnHydroseq, and WW mean/var by month
flow_map <- flowlines %>%
  st_drop_geometry() %>%
  select(
    Hydroseq, DnHydroseq,
    starts_with("WW_mean_"),
    starts_with("WW_var_")
  ) %>%
  mutate(across(starts_with("WW_"), as.numeric)) %>%
  arrange(desc(as.numeric(Hydroseq)))

# Aggregate any duplicated Hydroseq rows
# (sum mean and sum variance across duplicates)
flow_map <- flow_map %>%
  group_by(Hydroseq) %>%
  summarise(
    DnHydroseq = first(DnHydroseq),
    across(starts_with("WW_mean_"), ~sum(.x, na.rm = TRUE)),
    across(starts_with("WW_var_"),  ~sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(desc(as.numeric(Hydroseq)))

# Check: Hydroseq should decrease downstream
check_hseq <- flow_map %>%
  mutate(
    Hydroseq_num   = as.numeric(Hydroseq),
    DnHydroseq_num = as.numeric(DnHydroseq)
  ) %>%
  filter(!is.na(DnHydroseq_num))

if (any(check_hseq$DnHydroseq_num >= check_hseq$Hydroseq_num)) {
  warning("Some DnHydroseq have Hydroseq >= upstream segment; verify network ordering.")
}

# ------------------- CONTRIBUTORS (DMR-seeded) + MONTH-SPECIFIC HTML -------------------

# COMID -> Hydroseq lookup
comid_hseq <- flowlines %>%
  st_drop_geometry() %>%
  transmute(COMID = as.character(COMID), Hydroseq = as.character(Hydroseq)) %>%
  distinct()

# Facility "source" discharge table from DMR:
# each facility enters the network at its matched_COMID's Hydroseq.
facility_src <- dmr_monthly_fac %>%
  mutate(matched_COMID = as.character(matched_COMID)) %>%
  left_join(comid_hseq, by = c("matched_COMID" = "COMID")) %>%
  filter(!is.na(Hydroseq), mon_code %in% mon_codes) %>%
  group_by(Hydroseq, Facility, mon_code) %>%
  summarise(WW_fac = sum(mean_dmr_cfs_fac, na.rm = TRUE), .groups = "drop")

# Hydroseq -> facilities list (seed contributors only where DMR actually has flow)
fac_list_by_hseq <- facility_src %>%
  filter(is.finite(WW_fac), WW_fac > 0) %>%
  group_by(Hydroseq) %>%
  summarise(facs = list(unique(Facility)), .groups = "drop")

fac_lookup <- setNames(fac_list_by_hseq$facs, fac_list_by_hseq$Hydroseq)

# Initialize contributors list-column on flow_map
flow_map$contributors <- lapply(flow_map$Hydroseq, function(hseq) {
  if (!is.na(hseq) && hseq %in% names(fac_lookup)) fac_lookup[[hseq]] else character(0)
})

# Hydroseq → row index, and downstream neighbor index vector
idx    <- setNames(seq_len(nrow(flow_map)), flow_map$Hydroseq)
dn_idx <- idx[flow_map$DnHydroseq]

# Initialize cumulative mean + variance
for (mm in mon_codes) {
  flow_map[[paste0("WW_cum_mean_", mm)]] <- dplyr::coalesce(flow_map[[paste0("WW_mean_", mm)]], 0)
  flow_map[[paste0("WW_cum_var_",  mm)]] <- dplyr::coalesce(flow_map[[paste0("WW_var_",  mm)]],  0)
}

# Accumulate wastewater downstream:
# for each segment, add its cumulative mean+var to its downstream neighbor
for (i in seq_len(nrow(flow_map))) {
  j <- dn_idx[i]
  if (!is.na(j)) {
    for (mm in mon_codes) {
      flow_map[[paste0("WW_cum_mean_", mm)]][j] <- flow_map[[paste0("WW_cum_mean_", mm)]][j] +
        flow_map[[paste0("WW_cum_mean_", mm)]][i]
      
      flow_map[[paste0("WW_cum_var_", mm)]][j] <- flow_map[[paste0("WW_cum_var_", mm)]][j] +
        flow_map[[paste0("WW_cum_var_", mm)]][i]
    }
    # Union contributor facility names along the downstream path
    flow_map$contributors[[j]] <- union(flow_map$contributors[[j]], flow_map$contributors[[i]])
  }
}

# Now compute cumulative SD, low, high downstream (derived from mean+var)
for (mm in mon_codes) {
  sd_col <- paste0("WW_cum_sd_", mm)
  flow_map[[sd_col]] <- sqrt(pmax(flow_map[[paste0("WW_cum_var_", mm)]], 0))
  
  flow_map[[paste0("WW_cum_low_",  mm)]] <- pmax(flow_map[[paste0("WW_cum_mean_", mm)]] - flow_map[[sd_col]], 0)
  flow_map[[paste0("WW_cum_high_", mm)]] <-       flow_map[[paste0("WW_cum_mean_", mm)]] + flow_map[[sd_col]]
}

# Build a fast lookup for facility monthly discharge for all months
facility_month_lookup <- facility_src %>%
  filter(is.finite(WW_fac), WW_fac > 0) %>%
  group_by(mon_code, Facility) %>%
  summarise(WW_fac = sum(WW_fac, na.rm = TRUE), .groups = "drop")

# Create month-specific contributor HTML for each Hydroseq
for (mm in mon_codes) {
  ww_tot_col <- paste0("WW_cum_mean_", mm)
  out_col    <- paste0("contributors_html_", mm)
  
  m_lut <- facility_month_lookup %>%
    filter(mon_code == mm) %>%
    { setNames(.$WW_fac, .$Facility) }
  
  flow_map[[out_col]] <- vapply(seq_len(nrow(flow_map)), function(i) {
    ww_total <- flow_map[[ww_tot_col]][i]
    if (!is.finite(ww_total) || ww_total <= 0) return("")
    
    contribs <- flow_map$contributors[[i]]
    if (length(contribs) == 0) return("")
    
    ww_fac <- unname(m_lut[contribs])
    ww_fac[is.na(ww_fac)] <- 0
    
    keep <- is.finite(ww_fac) & ww_fac > 0
    if (!any(keep)) return("")
    
    facs <- contribs[keep]
    ww   <- ww_fac[keep]
    
    ord <- order(ww, decreasing = TRUE)
    facs <- facs[ord]
    ww   <- ww[ord]
    
    pct <- round(100 * ww / ww_total, 1)
    
    top_n <- min(3, length(facs))
    lines <- character(length(facs))
    
    lines[seq_len(top_n)] <- paste0("<b>", facs[seq_len(top_n)], ": ", pct[seq_len(top_n)], "%</b>")
    if (length(facs) > top_n) {
      lines[(top_n + 1):length(facs)] <- paste0(facs[(top_n + 1):length(facs)], ": ", pct[(top_n + 1):length(facs)], "%")
    }
    
    paste(lines, collapse = "<br>")
  }, FUN.VALUE = character(1))
}

# Create a simple comma-separated list of all upstream contributors
flow_map$contributors_all_str <- vapply(flow_map$contributors, function(x) {
  paste(unique(na.omit(x)), collapse = ", ")
}, FUN.VALUE = character(1))

# Join cumulative WW bands + contributors back to flowlines
flowlines <- flowlines %>%
  left_join(
    flow_map %>%
      select(
        Hydroseq,
        contributors_all_str,
        starts_with("contributors_html_"),
        starts_with("WW_cum_mean_"),
        starts_with("WW_cum_sd_"),
        starts_with("WW_cum_low_"),
        starts_with("WW_cum_high_")
      ) %>%
      distinct(),
    by = "Hydroseq"
  )

# Also keep a "contributors_str" column so existing helper doesn't break
if (!"contributors_str" %in% names(flowlines)) flowlines$contributors_str <- ""
flowlines$contributors_str <- flowlines$contributors_all_str

# --- Load and Process Gage Flow Data (monthly) -------------------------------
# Monthly gage stats for each USGS site
monthly_flows <- read_csv(
  paths$monthly_flows,
  show_col_types = FALSE
)

# Rename columns to Qgage_* pattern for easier handling
gage_monthlies <- monthly_flows %>%
  select(
    site_no,
    starts_with("mean_cfs_"),
    starts_with("sd_cfs_"),
    start_year, end_year, n_obs
  ) %>%
  # rename means and sds
  rename_with(~ gsub("mean_cfs_", "Qgage_",    .x), starts_with("mean_cfs_")) %>%
  rename_with(~ gsub("sd_cfs_",   "Qgage_sd_", .x), starts_with("sd_cfs_"))

# Build Qgage_var_* 
for (mm in mon_codes) {
  sd_name  <- paste0("Qgage_sd_",  mm)
  var_name <- paste0("Qgage_var_", mm)
  if (sd_name %in% names(gage_monthlies)) {
    gage_monthlies[[var_name]] <- gage_monthlies[[sd_name]]^2
  }
}


# Load gage-to-COMID spatial mapping
gages_sf <- st_read(
  paths$gages_sf,
  layer = "gages_sf"
)

gages_sf$best_COMID <- as.character(gages_sf$best_COMID)

# Ensure flowlines IDs are character
flowlines$Hydroseq   <- as.character(flowlines$Hydroseq)
flowlines$DnHydroseq <- as.character(flowlines$DnHydroseq)
flowlines$COMID      <- as.character(flowlines$COMID)

# Join gage monthly stats to each COMID (via best_COMID)
gages_joined <- gages_sf %>%
  st_drop_geometry() %>%
  select(site_no, best_COMID) %>%
  left_join(gage_monthlies, by = "site_no")

# De-duplicate: pick 1 best gage per COMID (prioritize longer and more recent records)
gages_per_comid <- gages_joined %>%
  group_by(best_COMID) %>%
  arrange(desc(n_obs), desc(end_year), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Attach gage stats to flowlines
flowlines <- flowlines %>%
  left_join(
    gages_per_comid %>% rename(COMID = best_COMID),
    by = "COMID"
  )

# --- Build flow order data for propagation (monthly) -------------------------
# Ensure draingage area is present in dataset 
if (!"TotDASqKM" %in% names(flowlines)) {
  stop("Flowlines are missing the TotDASqKM (total drainage area) column, ",
       "which the drainage-area-ratio propagation requires.", call. = FALSE)
}
# Set up a flow_order table used for gage propagation down the network
flow_order <- flowlines %>%
  st_drop_geometry() %>%
  mutate(DA = as.numeric(TotDASqKM)) %>%
  select(
    Hydroseq, DnHydroseq, site_no, QE_MA, DA,
    starts_with("Qgage_")
  ) %>%
  arrange(as.numeric(Hydroseq))

# Low-flow screen: ignore gages where all monthly flows are less than 0.5 cfs
flow_order <- flow_order %>%
  mutate(
    low_flow = if_all(all_of(paste0("Qgage_", mon_codes)), ~ . < 0.5), #changed to if_all on 1/28/26
    site_no  = if_else(low_flow, NA_character_, site_no)
  )

# Blank out Qgage_* where low_flow is TRUE
for (mm in mon_codes) {
  mcol  <- paste0("Qgage_",      mm)
  sdcol <- paste0("Qgage_sd_",   mm)
  vcol  <- paste0("Qgage_var_",  mm)
  
  if (mcol  %in% names(flow_order)) flow_order[[mcol]]  <- ifelse(flow_order$low_flow, NA_real_, flow_order[[mcol]])
  if (sdcol %in% names(flow_order)) flow_order[[sdcol]] <- ifelse(flow_order$low_flow, NA_real_, flow_order[[sdcol]])
  if (vcol  %in% names(flow_order)) flow_order[[vcol]]  <- ifelse(flow_order$low_flow, NA_real_, flow_order[[vcol]])
}
flow_order <- flow_order %>% select(-low_flow)

# --- Manual overrides for segments where NHDPlus is clearly wrong ----------
# These Hydroseqs will always accept propagated gage flows (force_gage = TRUE)
forced_hseq <- c(
  "550025960",
  "550020503",
  "390005386",
  "390025451",
  "390005505"
) # add more Hydroseq IDs here if needed

# ── Define seed segments and initialize needs_ww_adj flag ──────────
same_seg_ww_adj_seed <- c(
  "550022475", #Loveland 
  "550015536", #Robert H White, gage upstream
  "550062530", #Louisville 
  "390045364", #Tri-lakes
  "550025946", #Boulder
  "550019184", #Mulberry 
  "760015713", #Aspen
  "390005515", #Las Animas
  "550018746"  #BoxElder
  )

# One row per Hydroseq
flow_order <- flow_order %>%
  mutate(
    force_gage   = Hydroseq %in% forced_hseq
  )

# Donor gage Hydroseq: segments with a direct gage know their own Hydroseq
flow_order <- flow_order %>%
  mutate(
    gage_hseq = if_else(!is.na(site_no), Hydroseq, NA_character_)
  )

# Precompute downstream index in flow_order for fast propagation
dn_index <- match(flow_order$DnHydroseq, flow_order$Hydroseq)

# Initialize propagation columns
flow_order$gage_used <- flow_order$site_no
for (mm in mon_codes) {
  mcol <- paste0("Qgage_",     mm)
  vcol <- paste0("Qgage_var_", mm)
  
  flow_order[[paste0(mcol, "_prop")]] <- as.numeric(flow_order[[mcol]])
  flow_order[[paste0(vcol, "_prop")]] <- as.numeric(flow_order[[vcol]])
}


# Mark confluences: any Hydroseq that appears as DnHydroseq more than once
flow_order$is_confluence <- flow_order$Hydroseq %in% flow_order$DnHydroseq[duplicated(flow_order$DnHydroseq)]

# Join local WW onto flow_order so propagation can add it to Qgage_prop
ww_local_wide <- flow_map %>%
  select(Hydroseq, starts_with("WW_mean_"), starts_with("WW_var_"))
flow_order <- flow_order %>%
  left_join(ww_local_wide, by = "Hydroseq") %>%
  mutate(
    across(starts_with("WW_mean_"), ~coalesce(as.numeric(.x), 0)),
    across(starts_with("WW_var_"),  ~coalesce(as.numeric(.x), 0))
  )

# Bake local WW into Qgage_prop at same-seg seeds so downstream inherits correct denominator
for (mm in mon_codes) {
  mcol     <- paste0("Qgage_", mm, "_prop")
  ww_col   <- paste0("WW_mean_", mm)
  seed_idx <- which(flow_order$Hydroseq %in% same_seg_ww_adj_seed)
  for (i in seed_idx) {
    ww_val <- coalesce(as.numeric(flow_order[[ww_col]][i]), 0)
    flow_order[[mcol]][i] <- flow_order[[mcol]][i] + ww_val
  }
}

# Statewide drainage-area-ratio exponent 
raw_q_cols <- paste0("Qgage_", mon_codes)
flow_order$ann_mean_flow <- rowMeans(as.matrix(flow_order[, raw_q_cols]), na.rm = TRUE)
flow_order$ann_mean_flow[is.nan(flow_order$ann_mean_flow)] <- NA_real_

da_fit_data <- flow_order %>%
  filter(!is.na(site_no), is.finite(ann_mean_flow), is.finite(DA),
         ann_mean_flow > 0, DA > 0)
da_fit <- lm(log(ann_mean_flow) ~ log(DA), data = da_fit_data)
b_da   <- unname(coef(da_fit)[2])

if (!is.finite(b_da)) stop("DA-ratio exponent fit failed; check TotDASqKM and the gage flow columns.", call. = FALSE)
if (b_da < 0.5 || b_da > 1.0) warning(sprintf("Fitted DA-ratio exponent b = %.3f is outside the expected 0.5-1.0 range; verify inputs.", b_da))
message(sprintf("DA-ratio exponent: b = %.3f (n = %d gaged segments, R^2 = %.3f)",
                b_da, nrow(da_fit_data), summary(da_fit)$r.squared))

# Network structure for a single upstream-to-downstream pass 
n_seg     <- nrow(flow_order)
DA_vec    <- flow_order$DA
force_vec <- flow_order$force_gage
order_idx <- order(as.numeric(flow_order$Hydroseq), decreasing = TRUE)

valid_edges <- which(!is.na(dn_index))
pred_split  <- split(valid_edges, dn_index[valid_edges])
pred_list   <- vector("list", n_seg)                 # pred_list[[j]] = direct upstream rows of j
pred_list[as.integer(names(pred_split))] <- pred_split
pred_list[sapply(pred_list, is.null)] <- list(integer(0))

# State carried through the pass (12 months as matrix columns) 
known     <- which(rowSums(is.finite(as.matrix(flow_order[, paste0("Qgage_", mon_codes)]))) > 0)
has_val   <- seq_len(n_seg) %in% known                # carries a flow estimate (direct or propagated)
protected <- has_val | !is.na(flow_order$site_no)     # real gages are never overwritten

Qprop   <- as.matrix(flow_order[, paste0("Qgage_",     mon_codes, "_prop")])
Vprop   <- as.matrix(flow_order[, paste0("Qgage_var_", mon_codes, "_prop")])
WWm_mat <- as.matrix(flow_order[, paste0("WW_mean_",   mon_codes)])
WWv_mat <- as.matrix(flow_order[, paste0("WW_var_",    mon_codes)])

gage_used_vec <- flow_order$gage_used
gage_hseq_vec <- flow_order$gage_hseq

# Propagation math
for (j in order_idx) {
  if (protected[j]) next
  
  preds <- pred_list[[j]]
  if (length(preds) == 0L) next
  preds_val <- preds[has_val[preds]]
  n_val     <- length(preds_val)
  if (n_val == 0L) next
  
  if (n_val == 1L) {
    # ----  single-parent hop: gated, drainage-area-ratio scaled ----
    i  <- preds_val
    ai <- DA_vec[i]; aj <- DA_vec[j]
    gate_ok <- is.finite(ai) && is.finite(aj) && ai > 0 && aj > 0 &&
      abs(ai - aj) / max(aj, 1e-6) < SIMILARITY_THRESHOLD
    if (!gate_ok && !isTRUE(force_vec[j])) next        # blocked; forced Hydroseqs still pass
    
    sf <- if (gate_ok) (aj / ai)^b_da else 1           # forced hop passes through unscaled
    
    # local WW is added at every newly filled segment 
    Qprop[j, ] <- Qprop[i, ] * sf + WWm_mat[j, ]
    Vprop[j, ] <- Vprop[i, ] * sf^2 + WWv_mat[j, ]
    
    gage_used_vec[j] <- gage_used_vec[i]
    gage_hseq_vec[j] <- gage_hseq_vec[i]
    
  } else {
    # Confluence logic 
    fin   <- is.finite(Qprop[preds_val, , drop = FALSE])   # branches x months
    n_fin <- colSums(fin)
    
    Qj <- rep(NA_real_, length(mon_codes))
    Vj <- rep(NA_real_, length(mon_codes))
    
    m2 <- which(n_fin >= 2L)
    if (length(m2) > 0L) {
      Qj[m2] <- colSums(Qprop[preds_val, m2, drop = FALSE], na.rm = TRUE) + WWm_mat[j, m2]
      Vj[m2] <- colSums(Vprop[preds_val, m2, drop = FALSE], na.rm = TRUE) + WWv_mat[j, m2]
    }
    
    for (m in which(n_fin == 1L)) {
      i  <- preds_val[fin[, m]]
      ai <- DA_vec[i]; aj <- DA_vec[j]
      gate_ok <- is.finite(ai) && is.finite(aj) && ai > 0 && aj > 0 &&
        abs(ai - aj) / max(aj, 1e-6) < SIMILARITY_THRESHOLD
      if (!gate_ok && !isTRUE(force_vec[j])) next          # skips to next month
      sf <- if (gate_ok) (aj / ai)^b_da else 1
      Qj[m] <- Qprop[i, m] * sf + WWm_mat[j, m]
      Vj[m] <- Vprop[i, m] * sf^2 + WWv_mat[j, m]
    }
    
    if (!any(is.finite(Qj))) next                            # nothing filled in any month
    Qprop[j, ] <- Qj
    Vprop[j, ] <- Vj
    
    gage_used_vec[j] <- paste(gage_used_vec[preds_val], collapse = "+")
    gage_hseq_vec[j] <- NA_character_
  }
  
  has_val[j] <- TRUE
}

# Write results back to flow_order
for (k in seq_along(mon_codes)) {
  mm <- mon_codes[k]
  flow_order[[paste0("Qgage_",     mm, "_prop")]] <- Qprop[, k]
  flow_order[[paste0("Qgage_var_", mm, "_prop")]] <- Vprop[, k]
}
flow_order$gage_used <- gage_used_vec
flow_order$gage_hseq <- gage_hseq_vec


# Derive low/high bands from propagated mean + variance (mean ± 1 SD)
for (mm in mon_codes) {
  mcol  <- paste0("Qgage_",     mm, "_prop")
  vcol  <- paste0("Qgage_var_", mm, "_prop")
  
  sd_vec <- sqrt(pmax(flow_order[[vcol]], 0))
  
  flow_order[[paste0("Qgage_low_",  mm, "_prop")]]  <- pmax(flow_order[[mcol]] - sd_vec, 0)
  flow_order[[paste0("Qgage_high_", mm, "_prop")]]  <- flow_order[[mcol]] + sd_vec
}

# Finalize columns:
# QE_mm_final / QE_mm_low_final / QE_mm_high_final are the gage-based denominators
for (mm in mon_codes) {
  flow_order[[paste0("QE_", mm, "_final")]]      <- flow_order[[paste0("Qgage_", mm, "_prop")]]
  flow_order[[paste0("QE_", mm, "_low_final")]]  <- flow_order[[paste0("Qgage_low_", mm, "_prop")]]
  flow_order[[paste0("QE_", mm, "_high_final")]] <- flow_order[[paste0("Qgage_high_", mm, "_prop")]]
}

# Join assigned gage flows & IDs to flowlines
flowlines <- flowlines %>%
  select(-site_no, -starts_with("Qgage_")) %>%
  left_join(
    flow_order %>%
      select(
        Hydroseq, site_no, gage_used, gage_hseq,
        matches("^QE_\\d{2}_final$"),
        matches("^QE_\\d{2}_low_final$"),
        matches("^QE_\\d{2}_high_final$")
      ),
    by = "Hydroseq"
  ) %>%
  mutate(
    flow_source = case_when(
      !is.na(site_no)                     ~ "gage_direct",
      is.na(site_no) & !is.na(gage_used)  ~ "gage_propagated",
      TRUE                                ~ "nhd_only"
    )
  )

# --- % Wastewater using monthly bands (gage-based + NHDPlus-based) ----------
# Helper: compute %WW, but return NA if denom or ww is NA; cap at 100%.
compute_pct <- function(ww, denom) {
  ifelse(
    is.na(denom) | is.na(ww) | denom <= 0,
    NA_real_,
    pmin(100, pmax(0, 100 * ww / denom))
  )
}

# NHDPlus monthly denominators 
# QE_mm_nhd is the original NHDPlus monthly flow
for (mm in mon_codes) {
  flowlines[[paste0("QE_", mm, "_nhd")]] <- as.numeric(flowlines[[paste0("QE_", mm)]])
}

# Gage-based %WW for each month:
#   - - pct_mean_mm: WW_cum / Qgage_prop (WW already baked into denominator via propagation)
#   - pct_min_mm, pct_max_mm: bands using WW_low/WW_high and Q_high/Q_low
for (mm in mon_codes) {
  ww_mean <- paste0("WW_cum_mean_", mm)
  ww_low  <- paste0("WW_cum_low_",  mm)
  ww_high <- paste0("WW_cum_high_", mm)
  
  q_mean  <- paste0("QE_", mm, "_final")
  q_low   <- paste0("QE_", mm, "_low_final")
  q_high  <- paste0("QE_", mm, "_high_final")
  
  pct_mean_col      <- paste0("pct_mean_",      mm)
  pct_min_col       <- paste0("pct_min_",       mm)
  pct_max_col       <- paste0("pct_max_",       mm)
  pct_mean_base_col <- paste0("pct_mean_", mm, "_base")
  
  WW_seg_mean <- flowlines[[ww_mean]]
  WW_seg_low  <- flowlines[[ww_low]]
  WW_seg_high <- flowlines[[ww_high]]
  
  Q_mean <- flowlines[[q_mean]]
  Q_low  <- flowlines[[q_low]]
  Q_high <- flowlines[[q_high]]
  
  flowlines[[pct_mean_base_col]] <- compute_pct(WW_seg_mean, Q_mean)
  flowlines[[pct_mean_col]]      <- compute_pct(WW_seg_mean, Q_mean)
  flowlines[[pct_min_col]]       <- compute_pct(WW_seg_low,  Q_high)
  flowlines[[pct_max_col]]       <- compute_pct(WW_seg_high, Q_low)
}

# NHDPlus-based %WW (parallel to the gage-based %WW, but using NHDPlus denominators)
for (mm in mon_codes) {
  ww_mean <- paste0("WW_cum_mean_", mm)
  ww_low  <- paste0("WW_cum_low_",  mm)
  ww_high <- paste0("WW_cum_high_", mm)
  
  q_nhd   <- paste0("QE_", mm, "_nhd")
  
  flowlines[[paste0("pct_mean_", mm, "_nhd")]] <- compute_pct(flowlines[[ww_mean]], flowlines[[q_nhd]])
  flowlines[[paste0("pct_min_",  mm, "_nhd")]] <- compute_pct(flowlines[[ww_low]],  flowlines[[q_nhd]])
  flowlines[[paste0("pct_max_",  mm, "_nhd")]] <- compute_pct(flowlines[[ww_high]], flowlines[[q_nhd]])
}

# Drop Z/M dimensions from geometry (if present)
flowlines <- st_zm(flowlines, drop = TRUE, what = "ZM")
  
# --- UI ----------------------------------------------------------------------
ui <- navbarPage(
  "ReSource Tracker",
  
  # ========================= TAB 1: MAP ======================================
  tabPanel(
    "Map",
    fluidPage(
      fluidRow(
        column(
          width = 3,
          
          h4("Map Controls"),
          helpText(
            "Explore how much of streamflow is contributed by upstream ",
            "wastewater discharges. Choose the denominator source and month."
          ),
          radioButtons(
            "src", "Denominator source",
            choices  = c("Gage-based", "NHDPlus-based"),
            selected = "Gage-based", inline = FALSE
          ),
          selectInput(
            "month", "Month (flow conditions)",
            choices  = setNames(mon_codes, month_labels),
            selected = "06"
          ),
          
          actionButton(
            "reset_map",
            label = "Reset Map View",
            icon  = icon("refresh"),
            width = "100%"
          ),
          
          tags$hr(),
          
          h4("Quick guide"),
          tags$small(
            "Hover over a stream and click to see:",
            tags$br(),
            "• Selected month wastewater flow (cfs)",
            tags$br(),
            "• % wastewater for the chosen denominator",
            tags$br(),
            "• Upstream WWTP contributors and gage info",
            tags$br(), tags$br(),
            "For equity analysis using Climate and Economic Justice Screening Tool (CEJST) indicators, use the 'Equity' tab."
          ), 
         
          tags$hr(),
          h4("Selected census tract"),
          verbatimTextOutput("tract_info")
        ),
        
        column(
          width = 9,
          maplibreOutput("map", height = "900px")
        )
      )
    )
  ),
  
  # ====================== TAB 2: HOW TO USE =================================
  tabPanel(
    "How to Use",
    fluidPage(
      fluidRow(
        column(
          width = 8,
          h3("How to use ReSource Tracker"),
          p(
            "This tool estimates the fraction of streamflow composed of ",
            "treated wastewater effluent for each NHDPlus flowline in Colorado."
          ),
          h4("1. Choose a denominator source (Map tab)"),
          tags$ul(
            tags$li(
              strong("Gage-based: "),
              "USGS stream gage data are processed into monthly means ",
              "and propagated downstream along segments of similar drainage area, ",
              "with flow scaled by the drainage-area ratio between segments. ",
              "This is the preferred denominator where available."
            ),
            tags$li(
              strong("NHDPlus-based: "),
              "Uses NHDPlus mean monthly flows (original NHDPlus volumes). ",
              "These provide statewide coverage but can overestimate flows in ",
              "regulated or drought-impacted streams."
            )
          ),
          h4("2. Select a month (Map tab)"),
          p(
            "Use the Month selector to switch between January–December. ",
            "All wastewater flows and percent wastewater values are computed ",
            "for the selected month."
          ),
          h4("3. Interpret the stream colors (Map tab)"),
          p(
            "Stream colors show the estimated percent wastewater for the ",
            "selected month and denominator source."
          ),
          tags$ul(
            tags$li("Dark colors ≈ low % wastewater (near 0%)."),
            tags$li("Bright colors ≈ high % wastewater (approaching 100%).")
          ),
          h4("4. Explore equity impacts (Equity tab)"),
          p(
            "The Equity tab links the stream network to Climate and Economic Justice Screening Tool (CEJST) census tracts. ",
            "For each tract, it summarizes tract-level wastewater prevalence (whether any intersecting segment is impacted) ",
            "and severity (how large % wastewater is among impacted segments), ",
            "and compares those exposures across bins of tract-level CEJST indicators."
          ),
          h4("5. Reading the popups (Map tab)"),
          p("Click on a stream segment to see a popup with:"),
          tags$ul(
            tags$li("The NHDPlus Hydroseq ID (unique segment identifier)."),
            tags$li("Monthly wastewater flow (mean, low, high) in cfs."),
            tags$li("% wastewater (mean/min/max) for the selected month."),
            tags$li("Which gage(s) are used for the denominator."),
            tags$li("Whether flow comes from a direct gage, propagated gage, or NHDPlus only."),
            tags$li("Top contributing wastewater treatment plants on that segment.")
          )
        )
      )
    )
  ),
  
  # ==================== TAB 3: EQUITY =======================================
  tabPanel(
    "Equity",
    fluidPage(
      fluidRow(
        column(
          width = 3,
          h4("CEJST & Exposure Settings"),
          helpText(
            "Two modes:",
            tags$br(),
            tags$ul(
              tags$li(strong("Prevalence:"), " tract-level (does the tract have ANY impacted segment?)."),
              tags$li(strong("Severity:"), " among impacted segments (%WW > 0), how large the contribution is.")
            )
          ),
          
          radioButtons(
            "cejst_mode", "Analysis mode",
            choices  = c("Prevalence", "Severity"),
            selected = "Prevalence",
            inline   = FALSE
          ),
          
          conditionalPanel(
            condition = "input.cejst_mode == 'Severity'",
            radioButtons(
              "cejst_src", "Denominator source (severity only)",
              choices  = c("Gage-based", "NHDPlus-based"),
              selected = "Gage-based",
              inline   = FALSE
            )
          ),
          
          selectInput(
            "cejst_month", "Month (flow conditions)",
            choices  = setNames(mon_codes, month_labels),
            selected = "06"
          ),
          
          selectInput(
            "rurality",
            "Rurality category",
            choices  = c("All", "Urban", "Rural"),
            selected = "All"
          ),
          
          ## <-- Download data block 
          tags$hr(),
          h4("Download"),
          helpText("Export tract-level data used in the current equity plot."),
          downloadButton("download_cejst_data", "Download equity CSV"),
          
          tags$hr(),
          h4("Grouping variable (x-axis bins)"),
          helpText(
            "Choose a CEJST numeric indicator (binned) or Rurality (categorical).",
            tags$br(), tags$br(),
            strong("To compare Urban/Suburban/Rural:"),
            tags$br(),
            HTML("&bull; Set <b>Rurality category</b> to <b>All</b>"),
            tags$br(),
            HTML("&bull; Set <b>Grouping variable</b> to <b>Rurality</b>")
          ),
          
          selectInput(
            "cejst_var_box",
            "CEJST variable",
            choices = NULL
          ),
          
          conditionalPanel(
            condition = "input.cejst_var_box != 'rurality'",
            sliderInput(
              "cejst_nbins",
              "Number of bins (numeric variables only)",
              min   = 2,
              max   = 6,
              value = 4,
              step  = 1,
              round = TRUE,
              ticks = FALSE
            )
          )
        ),
        
        column(
          width = 9,
          h3("Exposure Plot: By CEJST Indicator or Rurality"),
          plotOutput("cejst_plot", height = "520px"),
          tags$hr(),
          h4("Statistical comparison between groups"),
          verbatimTextOutput("cejst_stats")
        )
      )
    )
  ),
  
  #==================  TAB 4: SEASONALITY ===================================
  tabPanel(
    "Seasonality",
    fluidPage(
      fluidRow(
        column(
          width = 3,
          h4("Seasonality settings"),
          helpText("Monthly summaries across all tracts (or filtered)."),
          
          radioButtons(
            "season_src", "Denominator source (%WW)",
            choices = c("Gage-based", "NHDPlus-based"),
            selected = "Gage-based",
            inline = FALSE
          ),
          
          selectInput(
            "season_rurality",
            "Rurality category",
            choices  = c("All", "Urban", "Rural"),
            selected = "All"
          ),
          
          checkboxInput(
            "season_use_tracts",
            "Summarize at tract-level (instead of segment-level)",
            TRUE
          ),
          helpText("Tract-level collapses many segments into one value per tract per month."),
          
          tags$hr(),
          
          h4("Download"),
          helpText("Export the data used to draw the seasonality plot."),
          downloadButton("download_season_data", "Download seasonality CSV")
        ),
        
        column(
          width = 9,
          h3("Seasonal patterns in wastewater influence"),
          plotOutput("season_plot_median", height = "420px"),
          tags$hr(),
          verbatimTextOutput("season_text"),
          
          tags$hr(),
          h4("Month-to-month statistics"),
          verbatimTextOutput("season_stats_text")
        )
      )
    )
  ),
  
  # ================== TAB 5: DATA SOURCES ===================================
  tabPanel(
    "Data Sources",
    fluidPage(
      fluidRow(
        column(
          width = 8,
          h3("Core input datasets"),
          p(
            "This tab documents the main input datasets used to build the ",
            "ReSource Tracker. Where possible, links point to ",
            "public, externally maintained sources. Processed versions (e.g., ",
            "Colorado-only subsets, network-joined tables) are distributed via ",
            "the accompanying manuscript / repository."
          ),
          
          h4("1. Stream network and base hydrology"),
          tags$ul(
            tags$li(
              strong("NHDPlus flowlines: "),
              "National Hydrography Dataset Plus (NHDPlus) flowlines and ",
              "associated attributes (e.g., COMID, Hydroseq, mean annual flow). ",
              "Used as the base river network and for routing wastewater and ",
              "streamflow downstream."
            ),
            tags$li(
              strong("NHDPlus mean monthly flows: "),
              "NHDPlus-supplied mean monthly flows are used as an alternative ",
              "denominator when USGS gage data are unavailable."
            )
          ),
          
          h4("2. USGS stream gage records"),
          tags$ul(
            tags$li(
              strong("USGS NWIS daily streamflow: "),
              "Daily discharge (cfs) time series for Colorado gages, retrieved ",
              "from the USGS National Water Information System (NWIS). These ",
              "records are aggregated to monthly means and propagated along the ",
              "NHDPlus network to estimate contemporary streamflow conditions."
            ),
            tags$li(
              strong("Gage–COMID linkage: "),
              "Each gage is spatially matched to its nearest NHDPlus COMID; ",
              "only one gage per COMID is retained, favoring the longest and ",
              "most recent record."
            )
          ),
          
          h4("3. Wastewater treatment plants and discharges"),
          tags$ul(
            tags$li(
              strong("WWTP locations: "),
              "Facility coordinates and design flows for municipal wastewater ",
              "treatment plants discharging to surface waters in Colorado. ",
              "Used to place outfalls on the map and assign them to NHDPlus ",
              "flowlines."
            ),
            tags$li(
              strong("Discharge Monitoring Report (DMR) data: "),
              "Monthly-average effluent flows are compiled from publicly ",
              "reported DMR data. Flows are converted to cfs, aggregated by ",
              "facility–COMID–month, and then accumulated downstream along the ",
              "river network."
            )
          ),
          
          h4("4. Drinking water systems and intakes"),
          tags$ul(
            tags$li(
              strong("Drinking water treatment plants (DWTPs): "),
              "Geocoded locations of community drinking water systems, including ",
              "population served and primary source type (e.g., surface water, ",
              "groundwater, mixed). These points are displayed as map markers."
            )
          ),
          
          h4("5. Climate and Economic Justice Screening Tool (CEJST)"),
          tags$ul(
            tags$li(
              strong("CEJST tract-level indicators: "),
              "Census-tract–level CEJST data, including demographic, economic, ",
              "health, climate, and pollution burden indicators. These fields ",
              "are joined to each tract and then intersected with the stream ",
              "network to summarize wastewater exposure by community."
            )
          ),
          
          h4("6. Processed intermediate tables"),
          p(
            "The following intermediate tables are generated within the analysis ",
            "workflow and are not standalone external datasets, but are important ",
            "for reproducibility:"
          ),
          tags$ul(
            tags$li(
              strong("Seasonal gage summary: "),
              "Per-gage monthly mean, standard deviation, and low / high flow ",
              "estimates used in the map popups."
            ),
            tags$li(
              strong("WWTP monthly flows: "),
              "Facility–COMID–month tables of effluent discharge (cfs) before ",
              "downstream accumulation."
            ),
            tags$li(
              strong("Tract-level exposure summary: "),
              "For each tract, the median, 90th percentile, and maximum % ",
              "wastewater across all intersecting stream segments for each ",
              "month and denominator source."
            )
          ),
          
          h4("Accessing the inputs"),
          p(
            "Public links and permanent DOIs for the processed inputs (e.g., ",
            "Colorado-only flowlines, gage summaries, CEJST-joined tables) will ",
            "be provided in the associated manuscript, data repository, or ",
            "project website. This tab serves as a living data dictionary for ",
            "interpreting those files."
          )
        )
      )
    )
  ),
  
  # ================== TAB 6: DOWNLOADS & METHODS =============================
  tabPanel(
    "Downloads & Methods",
    fluidPage(
      fluidRow(
        column(
          width = 8,
          h3("Data & Methods"),
          p(
            "This tab will be updated to allow downloading of data for users to conduct their own analysis."
          ),
          h4("Flow estimation and propagation (summary)"),
          tags$ul(
            tags$li(
              "USGS daily streamflow records (post–2001) are aggregated to ",
              "monthly means, with variability summarized using the interannual ",
              "standard deviation and corresponding mean ± 1 SD bands (low/high) ",
              "for each calendar month."
            ),
            tags$li(
              "Each gage is mapped to its closest NHDPlus flowline (COMID); ",
              "only one gage per COMID is retained, favoring the longest and ",
              "most recent record."
            ),
            tags$li(
              "Monthly flows are propagated downstream one segment at a time, ",
              "scaled by the drainage-area ratio between segments ",
              "(Q_downstream = Q_upstream × (A_downstream / A_upstream)^b, where b is ",
              "a statewide exponent fitted to the gage network). Propagation continues ",
              "only while drainage area changes by less than 50% (relative to the ",
              "downstream segment) from one segment to the next, and stops at larger ",
              "changes or at new gages."
            ),
            tags$li(
              "At confluences, flows from multiple upstream branches that carry ",
              "gage-based estimates are summed (mass balance). Segments downstream of ",
              "a confluence are then subject to the same drainage-area gate as any ",
              "other segment."
            )
          ),
          h4("Wastewater accumulation (summary)"),
          tags$ul(
            tags$li(
              "DMR-reported effluent flows are converted to cfs and averaged ",
              "by facility, COMID, and month, collapsing multiple outfalls per ",
              "facility–COMID combination."
            ),
            tags$li(
              "Wastewater flows are summed for all WWTPs mapped to the same ",
              "COMID and accumulated downstream using the NHDPlus flow network."
            ),
            tags$li(
              "For each segment and month, % wastewater is computed using ",
              "either the gage-based or NHDPlus-based denominator."
            )
          )
        )
      )
    )
  )
)

# --- SERVER ------------------------------------------------------------------
server <- function(input, output, session) {
  
  #helper
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  selected_tract_id <- reactiveVal(NULL)
  # ---- GAGES (markers; monthly stats) ---------------------------------------
  usgs_gages <- readr::read_csv(
    paths$monthly_flows,
    show_col_types = FALSE
  )
  
  safe_round <- function(x) ifelse(is.finite(x), round(x, 1), "NA")
  pal        <- leaflet::colorNumeric(palette = "plasma", domain = c(0, 100), na.color = "transparent") ## Last use of leaflet function
  
  # Build popup string per gage with monthly mean±sd and low/high
  usgs_gages$popup_str <- sapply(seq_len(nrow(usgs_gages)), function(i) {
    lines <- sapply(mon_codes, function(mm) {
      paste0(
        month_labels[mm], " mean±sd: ",
        safe_round(usgs_gages[[paste0("mean_cfs_", mm)]][i]), " cfs ",
        "(sd ",  safe_round(usgs_gages[[paste0("sd_cfs_",   mm)]][i]), "; ",
        "low ", safe_round(usgs_gages[[paste0("low_cfs_",  mm)]][i]), ", ",
        "high ",safe_round(usgs_gages[[paste0("high_cfs_", mm)]][i]), ")"
      )
    })
    paste0(
      "<b>", usgs_gages$station_nm[i], "</b> (", usgs_gages$site_no[i], ")<br>",
      paste(lines, collapse = "<br>"), "<br><br>",
      "<i>Years: ", usgs_gages$start_year[i], "–", usgs_gages$end_year[i],
      ", N = ", usgs_gages$n_obs[i], "</i>"
    )
  })
  
  usgs_gages_sf <- sf::st_as_sf(usgs_gages, coords = c("dec_long_va", "dec_lat_va"), crs = 4326)
  
  #For looking at tract polygons. Tracts is called later in mapgl function
  observeEvent(input$map_feature_click, {
    click <- input$map_feature_click
    if (is.null(click)) return()
    if (!identical(click$layer, "tracts")) return()
    
    ## mapgl generates a numeric feature id, reading GEOID from properties
    gid <- click$properties$GEOID
    if (is.null(gid)) return()
    
    # make sure it exists
    if (!gid %in% tracts_sf$GEOID) return()
    
    selected_tract_id(gid)
  })
  
  # Ensure contributors_str exists
  if (!"contributors_str" %in% names(flowlines)) flowlines$contributors_str <- ""
  
  contrib_block <- function(x) {
    ifelse(!is.na(x) & nzchar(x),
           paste0(
             "<details><summary style='cursor:pointer; color:#1E90FF; font-weight:bold; text-decoration:underline;'>",
             "Show WWTP Contributors</summary>",
             "<div style='margin-left:10px;'>", x, "</div></details><br>"
           ),
           ""
    )
  }
  
  # Column selector for map: which %WW column to use based on denominator source
  cols_for <- reactive({
    m <- input$month
    if (input$src == "Gage-based") {
      list(pct = paste0("pct_mean_", m))
    } else {
      list(pct = paste0("pct_mean_", m, "_nhd"))
    }
  })
  
  filtered_flowlines <- reactive({ flowlines })
  
  # ---- Shared ggplot theme for all app plots --------------------------------
  theme_app <- function() {
    theme_minimal(base_size = 15) +
      theme(
        plot.title    = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 13),
        axis.title.x  = element_text(size = 15, face = "bold"),
        axis.title.y  = element_text(size = 15, face = "bold"),
        axis.text.x   = element_text(size = 13, face = "bold"),
        axis.text.y   = element_text(size = 13, face = "bold"),
        legend.position = "none"
      )
  }
  # ---- Precompute flowlines–CEJST join for equity tab ----------------------
  flowlines_4326 <- sf::st_transform(flowlines, 4326)
  fl_cejst_sf <- sf::st_join(flowlines_4326, tracts_sf, left = TRUE)
  fl_cejst <- fl_cejst_sf %>% sf::st_drop_geometry()
  
  # Populate CEJST variable choices (numeric, non-percentile) + rurality
  observe({
    cejst_num_all <- names(cejst)[sapply(cejst, is.numeric)]
    
    cejst_num_raw <- cejst_num_all[!grepl("percentile", cejst_num_all, ignore.case = TRUE) &
                                     !grepl("pct_", cejst_num_all, ignore.case = TRUE) &
                                     !grepl("_pct$", cejst_num_all, ignore.case = TRUE)]
    cejst_num_raw <- setdiff(cejst_num_raw, "wastewater_pct")
    if (length(cejst_num_raw) == 0) cejst_num_raw <- cejst_num_all
    
    if ("rurality" %in% names(cejst)) {
      var_values <- c("rurality", cejst_num_raw)
      var_labels <- c("Rurality", cejst_num_raw)
    } else {
      var_values <- cejst_num_raw
      var_labels <- cejst_num_raw
    }
    
    updateSelectInput(
      session, "cejst_var_box",
      choices  = stats::setNames(var_values, var_labels),
      selected = var_values[1]
    )
  })
  
  # ---- BASE MAP (Map tab) ----------------------------------------------------
  ## This is the bulk of the changes to mapgl. 
  ## mapgl functions reference the entire 'sf' object when sent to functions. 
  ## This helper function trims the tracts to what the map actually calls to reduce some comptuaiton. 

  ## Trimming flowlines data in tracts_sf
  tracts_map <- tracts_sf %>%
    dplyr::mutate(
      popup_html = paste0(
        "<b>Census tract: </b>", GEOID, "<br>",
        "<b>Rurality: </b>", rurality, "<br>", 
        "<b> Low income tract? </b>",
        ifelse(low_income, "Yes", "No"), "<br>", 
        ifelse(!is.na(status), 
                paste0("<b>Status: </b>", status, "<br>"), 
                ""), 
        "<i>Click to see tract-level wastewater metrics in left panel.</i>"
      )
    )%>%
    dplyr::select(GEOID, popup_html)

## Generates static object that is developed once. 
  dwtp_map <- geocoded_dwtp
  if (!inherits(dwtp_map, "sf")) {
    dwtp_map <- dwtp_map %>%
      dplyr::filter(!is.na(Longitude), !is.na(Latitude)) %>%
      sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)
  }
## Transmute drops trims columns not called explictly. Repeated for dwtp and wwtp
  dwtp_map <- dwtp_map %>%
    dplyr::transmute(
      pws_name   = `PWS Name`,
      popup_html = paste0("<b>", `PWS Name`, "</b><br>Pop: ", `Population Served`,
                          "<br>Source: ", `Primary Source`)
    )
  wwtp_map <- wwtp_pts %>%
    dplyr::transmute(
      facility = Facility,
      popup_html = paste0("<b>", Facility, "</b><br>Design Flow: ", DesignFlow, "MGD")
    )
  ## Uses select to filter gage sf
  gages_map <- usgs_gages_sf %>%
    dplyr::select(station_nm, popup_str)


  # ---- STREAMS: colours + popups for the current month / denominator ---------
  ## Updating color logic for heatmap with mapgl functions. 
  ## Note steams_layer_data is reactive; depends on month and source to determine color, 
  ## and thus are recomputed when dropdowns are changed. 
  streams_layer_data <- reactive({
    spec <- cols_for()
    dat  <- filtered_flowlines()
    m    <- input$month
    mlab <- month_labels[m]

    vals     <- dat[[spec$pct]]
    ww_color <- ifelse(is.na(vals), "#d9d9d9", pal(vals))

    if (input$src == "Gage-based") {
      ww_mean_col <- paste0("WW_cum_mean_", m)
      ww_low_col  <- paste0("WW_cum_low_",  m)
      ww_high_col <- paste0("WW_cum_high_", m)
      sd_col      <- paste0("WW_cum_sd_",   m)

      pct_mean_col <- paste0("pct_mean_", m)
      pct_min_col  <- paste0("pct_min_",  m)
      pct_max_col  <- paste0("pct_max_",  m)

      popup_html <- paste0(
        "<b>Hydroseq:</b> ", dat$Hydroseq, "<br><br>",
        "<b>", mlab, " WW Flow (mean±sd):</b> ",
        safe_round(dat[[ww_mean_col]]), " cfs",
        " (sd ", safe_round(dat[[sd_col]]),
        "; low ",  safe_round(dat[[ww_low_col]]),
        ", high ", safe_round(dat[[ww_high_col]]), ")<br><br>",
        "<b>", mlab, " % Wastewater (Gage-based):</b><br>",
        "mean/min/max: ",
        safe_round(dat[[pct_mean_col]]), "% / ",
        safe_round(dat[[pct_min_col]]),  "% / ",
        safe_round(dat[[pct_max_col]]),  "%<br><br>",
        ifelse(!is.na(dat$start_year),
               paste0("<i>Gage years: ", dat$start_year, "–", dat$end_year,
                      ", N = ", dat$n_obs, "</i><br>"), ""),
        ifelse(!is.na(dat$gage_used),
               paste0("<b>Gage used:</b> ", dat$gage_used, "<br>"), ""),
        ifelse(!is.na(dat$flow_source),
               paste0("<b>Flow source:</b> ", dat$flow_source, "<br>"), ""),
        contrib_block(dat[[paste0("contributors_html_", m)]])
      )
    } else {
      q_nhd_col    <- paste0("QE_", m, "_nhd")
      pct_mean_nhd <- paste0("pct_mean_", m, "_nhd")

      popup_html <- paste0(
        "<b>Hydroseq:</b> ", dat$Hydroseq, "<br><br>",
        "<b>", mlab, " NHDPlus flow:</b> ",
        safe_round(dat[[q_nhd_col]]), " cfs<br><br>",
        "<b>", mlab, " % Wastewater (NHDPlus-based):</b><br>",
        "mean: ", safe_round(dat[[pct_mean_nhd]]), "%<br><br>",
        contrib_block(dat[[paste0("contributors_html_", m)]])
      )
    }

    sf::st_sf(
      data.frame(
        hydroseq   = dat$Hydroseq,
        pct        = vals,
        ww_color   = ww_color,
        line_w     = scales::rescale(log1p(dat$QE_MA), to = c(0.2, 4)),
        popup_html = popup_html,
        stringsAsFactors = FALSE
      ),
      geometry = sf::st_geometry(dat)
    )
  })

  pct_legend_title <- function(src, m) {
    paste0("% Wastewater (", src, ", ", month_labels[m], ")")
  }


    # Custom HTML legend with circle markers
    symbol_legend_html <-
      "<div style='background: white; padding: 6px 8px; font-size: 12px;
               box-shadow: 0 0 10px rgba(0,0,0,0.2);'>
     <b>Map symbols</b><br>
     <div>
       <span style='display:inline-block; width:10px; height:10px;
                    border-radius:50%; background:blue; margin-right:6px;'></span>
       Drinking water treatment plant (DWTP)
     </div>
     <div>
       <span style='display:inline-block; width:10px; height:10px;
                    border-radius:50%; background:brown; margin-right:6px;'></span>
       Wastewater treatment plant (WWTP)
     </div>
     <div>
       <span style='display:inline-block; width:10px; height:10px;
                    border-radius:50%; background:black; margin-right:6px;'></span>
       USGS gage
     </div>
     <div style='margin-top:4px;'>
       <span style='display:inline-block; width:20px; height:3px;
                    background:#d9d9d9; margin-right:6px;'></span>
       Stream with no %WW data
     </div>
   </div>"

  # ---- BASE MAP ----------------------------------------------
  # MapLibre draws layers inthe order they are called.
  output$map <- renderMaplibre({
    mapgl::maplibre(
      style = mapgl::openfreemap_style("positron") # This matches Carto theme best I think 
    ) %>%
      mapgl::fit_bounds(flowlines, animate = FALSE) %>%
      mapgl::add_navigation_control(position = "top-right") %>%
      mapgl::add_scale_control(position = "bottom-left") %>%

      ## --- Census Tract Polygons ------------------------------
      ## Those drawn first are on the 'bottom' of the map
      ## The color options here can be changed based on readability etc
      mapgl::add_fill_layer(
        id = "tracts", 
        source = tracts_map, 
        fill_color = "#CCCCCC", 
        fill_opacity = 0.15, 
        fill_outline_color = "#666666", 
        hover_options = list(fill_color = "#000000", fill_opacity = 0.35), 
        popup = "popup_html"
      ) %>%
      
      ## --- Streams ---------------------------------------------------
      ## source calls isolate() to not recompute on changes to month/source dropdown
      mapgl::add_line_layer(
        id = "streams", 
        source = isolate(streams_layer_data()), 
        line_color = mapgl::get_column("ww_color"), 
        line_width = mapgl::get_column("line_w"), 
        line_opacity = 0.85, 
        popup = "popup_html"
      ) %>%
    
    ## --- DWTP Markers -------------------------------------------------
    mapgl::add_circle_layer(
      id = "dwtp", 
      source = dwtp_map, 
      circle_color = "blue", 
      circle_radius = 6, 
      circle_opacity = 0.8, 
      popup = "popup_html", 
      tooltip = "pws_name"
    ) %>%

    ## --- WWTP Markers -------------------------------------------------
    mapgl::add_circle_layer(
      id = "Wwtp", 
      source = wwtp_map, 
      circle_color = "brown", 
      circle_radius = 6, 
      circle_opacity = 0.8, 
      popup = "popup_html", 
      tooltip = "facility"
    ) %>%

    ## --- USGS Gages -------------------------------------------------
    mapgl::add_circle_layer(
      id = "gages", 
      source = gages_map, 
      circle_color = "black", 
      circle_radius = 5, 
      circle_opacity = 1, 
      circle_stroke_color = "#ffffff", 
      circle_stroke_width = 1,
      popup = "popup_html", 
      tooltip = "station_nm"
    ) %>%

    ## --- Controls -------------------------------------------------
      mapgl::add_layers_control(
        position = "top-left", 
        layers = list(
          "Streams" = "streams", 
          "DWTP" = "dwtp", 
          "WWTP" = "wwtp", 
          "USGS Gages" = "gages", 
          "Census tracts" = "tracts"
        ), 
        collapsible = TRUE
      )%>%
    mapgl::add_legend(
      legend_title = pct_legend_title(isolate(input$src), isolate(input$month)), 
      values = c("0%", "100%"), 
      colors = pal(seq(0, 100, length.out = 11)), ##  uses same pallete from leaflet
      type = "continuous",
      position = "bottom-right", 
      unique_id = "pctLegend"
    ) %>%
    mapgl::add_control(
      html = symbol_legend_html, 
      position = "bottom-left", 
      id = "symbolLegend"
    )
  })
  
  observeEvent(input$reset_map, {
    mapgl::maplibre_proxy("map") %>%
      mapgl::fit_bounds(flowlines, animate = TRUE)
  })

  # --- Recolor and repopulate streams when dropdown changes --------------
  observeEvent(list(input$src, input$month), {
    mapgl::maplibre_proxy("map") %>%
      mapgl::clear_layer("streams") %>%
      mapgl::add_line_layer(
        id           = "streams",
        source       = streams_layer_data(),
        line_color   = mapgl::get_column("ww_color"),
        line_width   = mapgl::get_column("line_w"),
        line_opacity = 0.85,
        popup        = "popup_html",
        before_id    = "dwtp"   # keep streams under the point layers
      ) %>%
      mapgl::clear_legend("pctLegend") %>%
      mapgl::add_legend(
        legend_title = pct_legend_title(input$src, input$month),
        values       = c("0%", "100%"),
        colors       = pal(seq(0, 100, length.out = 11)),
        type         = "continuous",
        position     = "bottom-right",
        unique_id    = "pctLegend"
      )
  }, ignoreInit = TRUE)   # initial draw happens during renderMaplibre()

  # ================= SEASONALITY TAB LOGIC ====================
  season_cols_for <- reactive({
    src <- input$season_src %||% "Gage-based"
    list(
      pct = function(mm) if (identical(src, "Gage-based")) paste0("pct_mean_", mm) else paste0("pct_mean_", mm, "_nhd"),
      ww  = function(mm) paste0("WW_cum_mean_", mm)
    )
  })
  
  season_base <- reactive({
    dat <- fl_cejst
    if (!is.null(input$season_rurality) && input$season_rurality != "All") {
      if ("rurality" %in% names(dat)) dat <- dplyr::filter(dat, rurality == input$season_rurality) else return(NULL)
    }
    dplyr::filter(dat, !is.na(GEOID), !is.na(Hydroseq))
  })
  
  season_monthly <- reactive({
    dat <- season_base()
    req(dat)
    
    cols <- season_cols_for()
    
    use_tracts <- isTRUE(input$season_use_tracts)
    unit_lab   <- ifelse(use_tracts, "tracts", "segments")
    
    # long: keep BOTH GEOID and Hydroseq to switch units safely
    long <- dplyr::bind_rows(lapply(mon_codes, function(mm) {
      pct_col <- cols$pct(mm)
      ww_col  <- cols$ww(mm)
      if (!pct_col %in% names(dat) || !ww_col %in% names(dat)) return(NULL)
      
      tibble::tibble(
        GEOID    = dat$GEOID,
        Hydroseq = dat$Hydroseq,
        rurality = dat$rurality,
        mon_code = mm,
        month    = month_labels[mm],
        pct      = suppressWarnings(as.numeric(dat[[pct_col]])),
        ww       = suppressWarnings(as.numeric(dat[[ww_col]]))
      )
    }))
    req(nrow(long) > 0)
    
    # choose unit id
    long <- long %>%
      dplyr::mutate(unit_id = if (use_tracts) GEOID else as.character(Hydroseq))
    
    # collapse to ONE value per unit_id x month
    long_unit <- long %>%
      dplyr::group_by(unit_id, mon_code, month) %>%
      dplyr::summarise(
        pct = stats::median(pct, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(pct = dplyr::if_else(is.finite(pct), pct, NA_real_))
    
    # ever-impacted membership across months (units that have >0 at least once)
    ever_imp <- long_unit %>%
      dplyr::group_by(unit_id) %>%
      dplyr::summarise(
        ever_impacted = any(!is.na(pct) & pct > 0, na.rm = TRUE),
        .groups = "drop"
      )
    
    long_imp <- long_unit %>%
      dplyr::inner_join(ever_imp, by = "unit_id") %>%
      dplyr::filter(ever_impacted)
    
    # monthly summary (median + IQR) among ever-impacted units, zeros retained
    out <- long_imp %>%
      dplyr::group_by(mon_code, month) %>%
      dplyr::summarise(
        n_units_total = dplyr::n_distinct(unit_id),
        n_units       = sum(!is.na(pct)),
        med_pct = stats::median(pct, na.rm = TRUE),
        lo_iqr  = stats::quantile(pct, 0.25, na.rm = TRUE, names = FALSE),
        hi_iqr  = stats::quantile(pct, 0.75, na.rm = TRUE, names = FALSE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        mon_code = factor(mon_code, levels = mon_codes),
        month    = factor(month, levels = unname(month_labels[mon_codes]))
      )
    
    out
  })
  
  season_stats <- reactive({
    dat <- season_base()
    req(dat)
    
    cols <- season_cols_for()
    use_tracts <- isTRUE(input$season_use_tracts)
    
    long <- dplyr::bind_rows(lapply(mon_codes, function(mm) {
      pct_col <- cols$pct(mm)
      if (!pct_col %in% names(dat)) return(NULL)
      tibble::tibble(
        GEOID    = dat$GEOID,
        Hydroseq = dat$Hydroseq,
        mon_code = mm,
        pct      = suppressWarnings(as.numeric(dat[[pct_col]]))
      )
    }))
    req(nrow(long) > 0)
    
    long <- long %>%
      dplyr::mutate(unit_id = if (use_tracts) GEOID else as.character(Hydroseq)) %>%
      dplyr::group_by(unit_id, mon_code) %>%
      dplyr::summarise(
        pct = stats::median(pct, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(pct = dplyr::if_else(is.finite(pct), pct, NA_real_))
    
    # ever-impacted membership
    ever_imp <- long %>%
      dplyr::group_by(unit_id) %>%
      dplyr::summarise(ever_impacted = any(!is.na(pct) & pct > 0, na.rm = TRUE), .groups = "drop")
    
    long2 <- long %>%
      dplyr::inner_join(ever_imp, by = "unit_id") %>%
      dplyr::filter(ever_impacted)
    
    wide <- long2 %>% tidyr::pivot_wider(names_from = mon_code, values_from = pct)
    wide_cc <- wide %>% tidyr::drop_na()
    
    if (nrow(wide_cc) < 10) {
      return(list(n = nrow(wide_cc), friedman = NULL, pairwise = NULL))
    }
    
    X <- as.matrix(wide_cc[, mon_codes, drop = FALSE])
    fried <- stats::friedman.test(X)
    
    pairs <- combn(mon_codes, 2, simplify = FALSE)
    pw <- lapply(pairs, function(p) {
      a <- wide_cc[[p[1]]]; b <- wide_cc[[p[2]]]
      wt <- stats::wilcox.test(a, b, paired = TRUE, exact = FALSE)
      data.frame(m1 = p[1], m2 = p[2],
                 med_diff = stats::median(a - b, na.rm = TRUE),
                 p_value = wt$p.value,
                 stringsAsFactors = FALSE)
    }) %>% dplyr::bind_rows()
    
    pw$p_adj_BH <- p.adjust(pw$p_value, method = "BH")
    
    month_meds <- tibble::tibble(
      mon_code = mon_codes,
      month    = unname(month_labels[mon_codes]),
      med      = apply(X, 2, stats::median, na.rm = TRUE)
    ) %>% dplyr::arrange(dplyr::desc(med))
    
    list(n = nrow(wide_cc), friedman = fried, pairwise = pw, month_meds = month_meds)
  })
  
  output$season_stats_text <- renderPrint({
    s <- season_stats()
    req(s)
    
    unit_lab <- ifelse(isTRUE(input$season_use_tracts), "tracts", "segments")
    
    if (is.null(s$friedman) || is.null(s$pairwise) || is.null(s$month_meds)) {
      cat("Not enough complete units across all 12 months to run paired month-to-month tests.\n")
      cat("Tip: switch to tract-level, or loosen missingness (but paired tests need complete cases).\n")
      return(invisible(NULL))
    }
    
    # 1) Overall test
    cat(sprintf("Complete paired units (n_complete): %d\n", s$n))
    cat("Overall seasonality (Friedman):\n")
    cat(sprintf("  chi-squared = %.3f, df = %d, p = %s\n\n",
                unname(s$friedman$statistic),
                unname(s$friedman$parameter),
                format.pval(s$friedman$p.value, digits = 2, eps = 1e-100)))
    
    # 2) Ranked months by median
    top_n <- 4
    bot_n <- 4
    cat("Highest months by median %WW :\n")
    print(utils::head(s$month_meds, top_n), row.names = FALSE)
    cat("\nLowest months by median %WW:\n")
    print(utils::tail(s$month_meds, bot_n), row.names = FALSE)
    
    # 3) Only show “meaningful” significant pairwise differences
    pw <- s$pairwise %>%
      dplyr::mutate(
        sig = !is.na(p_adj_BH) & p_adj_BH < 0.05,
        abs_med_diff = abs(med_diff)
      ) %>%
      dplyr::filter(sig) %>%
      dplyr::arrange(p_adj_BH, dplyr::desc(abs_med_diff))
    
    cat("\n")
    
    if (nrow(pw) == 0) {
      cat("No month-to-month pairs are significant after BH correction (alpha = 0.05).\n")
      return(invisible(NULL))
    }
    
    show_n <- min(8, nrow(pw))
    
    cat(sprintf(
      "Top significant month-to-month differences (BH-adjusted p < 0.05; showing %d):\n",
      show_n
    ))
    
    # Prettify month labels (01 -> Jan)
    mm_to_abb <- function(mm) month.abb[as.integer(mm)]
    
    out <- pw %>%
      dplyr::mutate(
        month1   = mm_to_abb(m1),
        month2   = mm_to_abb(m2),
        med_diff = round(med_diff, 2),
        p_adj_BH_fmt = ifelse(
          p_adj_BH < 1e-100,
          "< 1e-100",
          formatC(p_adj_BH, format = "e", digits = 2)
        )
      ) %>%
      dplyr::select(month1, month2, med_diff, p_adj_BH = p_adj_BH_fmt) %>%
      utils::head(show_n)
    
    print(out, row.names = FALSE)
    
    cat("\nNotes:\n")
    cat("  med_diff = median(month1 - month2) across the same units (paired).\n")
    cat("  BH-adjusted p-values control false discovery rate across many pairwise tests.\n")
  })
  
  output$season_plot_median <- renderPlot({
    df <- season_monthly()
    req(df)
    
    src_lab  <- ifelse((input$season_src %||% "Gage-based") == "Gage-based", "gage-based", "NHDPlus-based")
    unit_lab <- ifelse(isTRUE(input$season_use_tracts), "tracts", "segments")
    
    # y-limit with a bit of headroom
    ymax <- max(df$hi_iqr, na.rm = TRUE)
    
    ggplot(df, aes(x = mon_code, y = med_pct, group = 1)) +
      # IQR ribbon: one neutral fill for all months 
      geom_ribbon(
        aes(ymin = lo_iqr, ymax = hi_iqr),
        fill  = "grey85",
        alpha = 0.4,
        color = NA
      ) +
      # Seasonal colors for line + points
      geom_line(aes(color = mon_code), linewidth = 1.1) +
      geom_point(aes(color = mon_code), size = 2.8) +
      
      scale_x_discrete(labels = setNames(month.abb, mon_codes)) +
      scale_color_manual(values = month_colors, guide = "none") +
      
      coord_cartesian(ylim = c(0, ymax * 1.12)) +
      
      labs(
        x = "Month",
        y = paste0("Median %WW ", unit_lab, "; zeros retained)"),
        title = paste0("Seasonality of wastewater influence (", src_lab, " denominator)"),
        subtitle = paste0("Shaded band: IQR (25–75%) ", unit_lab)
      ) +
      theme_app() +
      theme(
        axis.text.x = element_text(
          angle = 0,
          hjust = 0.5
        )
      )
  })
  
  
  output$season_text <- renderPrint({
    df <- season_monthly()
    req(df)
    
    unit_lab <- ifelse(isTRUE(input$season_use_tracts), "tracts", "segments")
    
    peak_med <- df %>%
      dplyr::filter(med_pct == max(med_pct, na.rm = TRUE)) %>%
      dplyr::slice(1)
    
    cat("Monthly summary (ever-impacted units only; zeros retained):\n")
    print(df %>% dplyr::select(mon_code, month, n_units, med_pct, lo_iqr, hi_iqr))
    
    cat("\nInterpretation notes:\n")
    cat(sprintf("  n_units = number of ever-impacted %s with a value that month.\n", unit_lab))
    cat("  Error bars show IQR (25–75%) of %WW across those units.\n\n")
    
    cat("Peak median month:\n")
    cat(sprintf("  %s (%s): median %.1f%% (IQR %.1f–%.1f%%)\n",
                as.character(peak_med$month), as.character(peak_med$mon_code),
                peak_med$med_pct, peak_med$lo_iqr, peak_med$hi_iqr))
   
  })
  
  # ---- Download: Seasonality data -----------------
  output$download_season_data <- downloadHandler(
    filename = function() {
      paste0("ReSourceTracker_seasonality_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- season_monthly()
      req(df)
      
      # Add context columns (denominator, rurality filter, unit type)
      df_out <- df %>%
        mutate(
          denominator_source = input$season_src %||% "Gage-based",
          rurality_filter    = input$season_rurality %||% "All",
          unit_type          = if (isTRUE(input$season_use_tracts)) "tract" else "segment"
        ) %>%
        # move context columns to front
        select(denominator_source, rurality_filter, unit_type, everything())
      
      readr::write_csv(df_out, file)
    }
  )
  
  # ================= CEJST & EQUITY TAB LOGIC =================
  
  cejst_pct_col <- reactive({
    req(input$cejst_month)
    m <- input$cejst_month
    
    # Only used for Severity
    src <- input$cejst_src %||% "Gage-based"
    
    if (identical(src, "Gage-based")) {
      paste0("pct_mean_", m)
    } else {
      paste0("pct_mean_", m, "_nhd")
    }
  })
  
  cejst_summary <- reactive({
    req(input$cejst_mode, input$rurality, input$cejst_month)
    
    m       <- input$cejst_month
    ww_col  <- paste0("WW_cum_mean_", m)   # prevalence uses numerator-only
    pct_col <- cejst_pct_col()            # severity uses denominator-specific %WW
    
    dat <- fl_cejst
    
    # Filter by rurality if requested
    if (!is.null(input$rurality) && input$rurality != "All") {
      if ("rurality" %in% names(dat)) {
        dat <- dplyr::filter(dat, rurality == input$rurality)
      } else {
        return(NULL)
      }
    }
    
    dat <- dplyr::filter(dat, !is.na(GEOID))
    
    # numeric CEJST columns from original cejst object
    cejst_num_cols <- names(cejst)[sapply(cejst, is.numeric)]
    
    # ---- PREVALENCE: tract has ANY impacted segment (0/1) -------------------
    if (identical(input$cejst_mode, "Prevalence")) {
      
      if (!ww_col %in% names(dat)) return(NULL)
      
      dat_prev <- dat %>%
        dplyr::mutate(
          ww = .data[[ww_col]],
          impacted = is.finite(ww) & !is.na(ww) & ww > 0
        )
      
      out <- dat_prev %>%
        dplyr::group_by(GEOID) %>%
        dplyr::summarise(
          any_impacted = as.integer(any(impacted, na.rm = TRUE)),
          metric_value = any_impacted, # 0/1 per tract
          dplyr::across(
            dplyr::all_of(cejst_num_cols),
            ~ dplyr::first(.x),
            .names = "{.col}"
          ),
          rurality = dplyr::first(rurality),
          .groups = "drop"
        )
      
      return(out)
    }
    
    # ---- SEVERITY: conditional on %WW > 0 -----------------------------------
    if (!pct_col %in% names(dat)) {
      warning("Missing severity %WW column: ", pct_col)
      return(NULL)
    }
    
    dat_sev <- dat %>%
      dplyr::mutate(
        pct = .data[[pct_col]],
        impacted = is.finite(pct) & !is.na(pct) & pct > 0
      )
    
    tract_base <- dat_sev %>%
      dplyr::group_by(GEOID) %>%
      dplyr::summarise(
        n_segments_total = dplyr::n(),
        n_impacted       = sum(impacted, na.rm = TRUE),
        dplyr::across(
          dplyr::all_of(cejst_num_cols),
          ~ dplyr::first(.x),
          .names = "{.col}"
        ),
        rurality = dplyr::first(rurality),
        .groups = "drop"
      )
    
    dat_imp <- dat_sev %>%
      dplyr::filter(impacted)
    
    sev_vals <- dat_imp %>%
      dplyr::group_by(GEOID) %>%
      dplyr::summarise(
        metric_value = stats::median(pct, na.rm = TRUE),
        .groups = "drop"
      )
    
    out <- tract_base %>%
      dplyr::inner_join(sev_vals, by = "GEOID")
    
    out
  })
  
  # Create bins (shared for prevalence/severity)
  cejst_binned <- reactive({
    tab <- cejst_summary()
    req(tab)
    
    var <- input$cejst_var_box
    req(var)
    
    # --- Rurality path (categorical bins) ---
    if (identical(var, "rurality")) {
      req("rurality" %in% names(tab))
      
      df <- tab %>%
        dplyr::select(GEOID, metric_value, rurality) %>%
        dplyr::filter(is.finite(metric_value), !is.na(rurality))
      
      df$var_bin <- factor(df$rurality, levels = c("Urban", "Rural"))
      return(df)
    }
    
    # --- Numeric CEJST variable path ---
    req(var %in% names(tab))
    
    df <- tab %>%
      dplyr::transmute(
        GEOID,
        metric_value,
        var_value = suppressWarnings(as.numeric(.data[[var]]))
      ) %>%
      dplyr::filter(is.finite(metric_value), is.finite(var_value))
    
    req(nrow(df) > 0)
    
    # Clamp bins 
    cur_nb <- input$cejst_nbins %||% 4L
    nb     <- max(2L, min(6L, as.integer(cur_nb)))
    
    n_uniq <- dplyr::n_distinct(df$var_value)
    if (n_uniq <= 1) {
      df$var_bin <- factor("All tracts")
      return(df)
    }
    
    nb_eff <- min(nb, n_uniq)
    
    # bin by rank -> stable bins even with ties
    df <- df %>%
      dplyr::mutate(
        .r      = dplyr::dense_rank(var_value),
        .bin_id = dplyr::ntile(.r, nb_eff)
      )
    
    # label bins by actual value ranges
    lab_tbl <- df %>%
      dplyr::group_by(.bin_id) %>%
      dplyr::summarise(
        lo = min(var_value, na.rm = TRUE),
        hi = max(var_value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(.bin_id) %>%
      dplyr::mutate(
        label = if (var %in% proportion_cols) {
          paste0("[", round(lo * 100, 0), ", ", round(hi * 100, 0), "%]")
        } else {
          paste0("[", round(lo, 1), ", ", round(hi, 1), "]")
        },
        label = dplyr::if_else(.bin_id == 1, paste0(label, " (lowest)"), label),
        label = dplyr::if_else(.bin_id == dplyr::n(), paste0(label, " (highest)"), label)
      )
    
    # Make sure levels are unique 
    lab_levels <- unique(lab_tbl$label)
    
    df <- df %>%
      dplyr::left_join(
        lab_tbl %>% dplyr::select(.bin_id, label),
        by = ".bin_id"
      ) %>%
      dplyr::mutate(
        var_bin = factor(label, levels = lab_levels)
      ) %>%
      dplyr::select(GEOID, metric_value, var_bin)
    
    df
  })
  
  # ---- Dynamic max bins for slider  -----------------------------------------
  observeEvent({
    list(
      input$cejst_var_box,
      input$cejst_mode,
      input$cejst_month,
      input$rurality,
      input$cejst_src
    )
  }, {
    var <- input$cejst_var_box
    req(var)
    
    # If rurality, slider is hidden
    if (identical(var, "rurality")) {
      updateSliderInput(
        session, "cejst_nbins",
        min   = 2,
        max   = 2,
        value = 2
      )
      return()
    }
    
    tab <- cejst_summary()
    req(tab)
    req("metric_value" %in% names(tab))
    req(var %in% names(tab))
    
    # Same numeric-filtering logic as cejst_binned()
    df <- tab %>%
      dplyr::transmute(
        GEOID,
        metric_value,
        var_value = suppressWarnings(as.numeric(.data[[var]]))
      ) %>%
      dplyr::filter(
        is.finite(metric_value),
        is.finite(var_value)
      )
    
    if (!nrow(df)) return()
    
    # Helper: compute the largest k (<= 6) that yields k DISTINCT labels
    max_feasible_bins <- function(v) {
      n_uniq <- dplyr::n_distinct(v)
      if (n_uniq <= 1) return(1L)
      
      max_try <- min(6L, n_uniq)
      
      for (k in seq(from = max_try, to = 2L, by = -1L)) {
        df_k <- tibble::tibble(var_value = v) %>%
          dplyr::mutate(
            .r      = dplyr::dense_rank(var_value),
            .bin_id = dplyr::ntile(.r, k)
          )
        
        lab_tbl <- df_k %>%
          dplyr::group_by(.bin_id) %>%
          dplyr::summarise(
            lo = min(var_value, na.rm = TRUE),
            hi = max(var_value, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          dplyr::arrange(.bin_id)
        
        # Key based on numeric range only
        key <- paste0(round(lab_tbl$lo, 6), "_", round(lab_tbl$hi, 6))
        
        if (dplyr::n_distinct(key) == k) {
          return(k)
        }
      }
      
      1L
    }
    
    mx_eff <- max_feasible_bins(df$var_value)
    mx     <- max(2L, min(6L, as.integer(mx_eff)))
    
    updateSliderInput(
      session,
      "cejst_nbins",
      min   = 2,
      max   = mx,
      value = mx  # default to the largest feasible number of bins
    )
  }, ignoreInit = FALSE)
  
  # ---- Plotting: Prevalence / Severity ---------------------------------
  output$cejst_plot <- renderPlot({
    df <- cejst_binned()
    req(df)
    df$var_bin <- droplevels(df$var_bin)
    
    var  <- input$cejst_var_box
    mlab <- month_labels[input$cejst_month]
    
    mode <- input$cejst_mode %||% "Prevalence"
    
    if (identical(mode, "Prevalence")) {
      
      prev_sum <- df %>%
        dplyr::group_by(var_bin) %>%
        dplyr::summarise(
          n    = dplyr::n(),
          k    = sum(metric_value == 1, na.rm = TRUE),
          prop = ifelse(n > 0, k / n, NA_real_),
          .groups = "drop"
        )
      
      z <- stats::qnorm(0.975)
      
      prev_sum <- prev_sum %>%
        dplyr::mutate(
          denom  = 1 + z^2 / n,
          center = (prop + z^2 / (2 * n)) / denom,
          half   = (z * sqrt((prop * (1 - prop) / n) + (z^2 / (4 * n^2)))) / denom,
          lo     = pmax(0, center - half),
          hi     = pmin(1, center + half)
        )
      
      title_txt <- paste0(
        "Prevalence by ",
        ifelse(identical(var, "rurality"), "Rurality", var),
        " (", mlab, "; prevalence uses WW>0 and is denominator-independent)"
      )
      
      ggplot(prev_sum, aes(x = var_bin, y = prop, fill = var_bin)) +
        geom_col(width = 0.75, color = "white", linewidth = 0.2) +
        geom_text(
          aes(label = paste0(round(100 * prop, 1), "% (", k, "/", n, ")")),
          vjust = -0.4,
          size = 4.2,
          fontface = "bold"
        ) +
        {
          if (identical(var, "rurality")) {
            scale_fill_manual(values = rurality_colors, drop = FALSE)
          } else {
            scale_fill_viridis_d(option = "C", end = 0.9)
          }
        } +
        labs(
          x = if (identical(var, "rurality")) "Rurality category" else paste0("Bins of ", var),
          y = "Share of tracts with any wastewater-impacted segment",
          title = title_txt
        ) +
        theme_app() +
        theme(
          axis.text.x = element_text(
            angle = if (identical(var, "rurality")) 0 else 25,
            hjust = if (identical(var, "rurality")) 0.5 else 1
          )
        )
      
    } else {
      
      src       <- input$cejst_src %||% "Gage-based"
      denom_lab <- ifelse(src == "Gage-based", "gage-based", "NHDPlus-based")
      
      ymax <- max(df$metric_value, na.rm = TRUE)
      
      sev_sum <- df %>%
        dplyr::group_by(var_bin) %>%
        dplyr::summarise(
          n   = dplyr::n(),
          med = stats::median(metric_value, na.rm = TRUE),
          mx  = max(metric_value, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          y_lab = pmin(mx + 0.06 * ymax, ymax * 1.08),
          label = paste0("Median: ", round(med, 1), "%\n(n=", n, ")")
        )
      
      ggplot(df, aes(x = var_bin, y = metric_value, fill = var_bin)) +
        geom_boxplot(outlier.alpha = 0.5) +
        geom_text(
          data = sev_sum,
          aes(x = var_bin, y = y_lab, label = label),
          inherit.aes = FALSE,
          size = 4.2,
          fontface = "bold"
        ) +
        {
          if (identical(var, "rurality")) {
            scale_fill_manual(values = rurality_colors, drop = FALSE)
          } else {
            scale_fill_viridis_d(option = "C", end = 0.9)
          }
        } +
        coord_cartesian(ylim = c(0, ymax * 1.12)) +
        labs(
          x = if (identical(var, "rurality")) "Rurality category" else paste0("Bins of ", var),
          y = "Severity: %WW among impacted segments",
          title = paste0(
            "Severity by ",
            ifelse(identical(var, "rurality"), "Rurality", var),
            " (", mlab, ", ", denom_lab, " denominator)"
          )
        ) +
        theme_app() +
        theme(
          axis.text.x = element_text(
            angle = if (identical(var, "rurality")) 0 else 25,
            hjust = if (identical(var, "rurality")) 0.5 else 1
          )
        )
    }
  })
  
  output$cejst_stats <- renderPrint({
    df <- cejst_binned()
    req(df)
    df$var_bin <- droplevels(df$var_bin)
    
    if (length(levels(df$var_bin)) < 2) {
      cat("Only one bin contains data; cannot compare bins statistically.\n")
      return(invisible(NULL))
    }
    
    cat("Tracts per bin:\n")
    print(table(df$var_bin))
    cat("\n")
    
    if (identical(input$cejst_mode, "Prevalence")) {
      
      # ---------- PREVALENCE STATS ----------
      tab <- table(df$var_bin, df$metric_value)
      cat("Contingency table (bin x any_impacted 0/1):\n")
      print(tab)
      cat("\n")
      
      chi <- try(suppressWarnings(chisq.test(tab)), silent = TRUE)
      if (!inherits(chi, "try-error")) {
        cat("Chi-squared test (overall differences in prevalence across bins):\n")
        cat(sprintf(
          "  X-squared = %.3f, df = %d, p-value = %.3e\n\n",
          unname(chi$statistic),
          unname(chi$parameter),
          chi$p.value
        ))
      } else {
        cat("Could not compute chi-squared test.\n\n")
      }
      
      prev_sum <- df %>%
        dplyr::group_by(var_bin) %>%
        dplyr::summarise(
          n    = dplyr::n(),
          k    = sum(metric_value == 1, na.rm = TRUE),
          prop = k / n,
          .groups = "drop"
        )
      
      cat("Prevalence by bin (k/n):\n")
      print(prev_sum)
      cat("\n")
      
      bins <- as.character(prev_sum$var_bin)
      if (length(bins) >= 2) {
        pairs <- combn(bins, 2, simplify = FALSE)
        pw <- lapply(pairs, function(p) {
          a <- prev_sum[prev_sum$var_bin == p[1], ]
          b <- prev_sum[prev_sum$var_bin == p[2], ]
          tt <- prop.test(
            x = c(a$k, b$k),
            n = c(a$n, b$n),
            correct = FALSE
          )
          data.frame(
            bin1      = p[1],
            bin2      = p[2],
            p_value   = tt$p.value,
            diff_prop = (a$k / a$n) - (b$k / b$n),
            stringsAsFactors = FALSE
          )
        })
        pw <- dplyr::bind_rows(pw)
        pw$p_adj_BH <- p.adjust(pw$p_value, method = "BH")
        
        cat("Pairwise proportion tests (BH-adjusted):\n")
        pw_print <- pw %>%
          dplyr::mutate(
            p_value  = formatC(p_value,  format = "e", digits = 2),
            p_adj_BH = formatC(p_adj_BH, format = "e", digits = 2),
            diff_prop = round(diff_prop, 3)
          )
        print(pw_print, row.names = FALSE)
      }
      
      cat("\nNotes:\n")
      cat("  Prevalence is tract-level: 1 if ANY intersecting segment has WW flow > 0, else 0.\n")
      cat("  Plot shows % of tracts impacted per bin.\n")
      cat("  Denominator selection does not affect prevalence (it only affects severity).\n")
      
    } else {
      
      # ---------- SEVERITY STATS ----------
      kw <- try(stats::kruskal.test(metric_value ~ var_bin, data = df), silent = TRUE)
      if (inherits(kw, "try-error")) {
        cat("Could not compute Kruskal–Wallis test.\n")
        return(invisible(NULL))
      }
      
      cat("Kruskal–Wallis test for differences between bins:\n")
      cat(sprintf(
        "  chi-squared = %.3f, df = %d, p-value = %.3e\n\n",
        kw$statistic,
        kw$parameter,
        kw$p.value
      ))
      
      if (kw$p.value < 0.05) {
        cat("Pairwise comparisons (Wilcoxon, BH-adjusted p-values):\n\n")
        pw <- try(
          stats::pairwise.wilcox.test(df$metric_value, df$var_bin, p.adjust.method = "BH"),
          silent = TRUE
        )
        if (!inherits(pw, "try-error")) {
          pmat     <- pw$p.value
          fmt_pmat <- pmat
          if (!is.null(pmat)) {
            fmt_pmat[] <- ifelse(
              is.na(pmat),
              NA,
              formatC(pmat, format = "e", digits = 2)
            )
          }
          print(fmt_pmat, quote = FALSE)
        } else {
          cat("  Could not compute pairwise Wilcoxon tests.\n")
        }
      } else {
        cat("Result: p ≥ 0.05 → no strong evidence bins differ.\n")
      }
      
      cat("\nNotes:\n")
      cat("  Severity uses %WW for the chosen denominator and is conditional on %WW > 0 (impacted segments only).\n")
      
      # ---------- flow-by-rurality summary  ----------------------------------
      if (identical(input$cejst_var_box, "rurality")) {
        m   <- input$cejst_month
        src <- input$cejst_src %||% "Gage-based"
        
        ww_col <- paste0("WW_cum_mean_", m)
        q_col  <- if (identical(src, "Gage-based")) {
          paste0("QE_", m, "_final")
        } else {
          paste0("QE_", m, "_nhd")
        }
        
        # Use the global fl_cejst table (flowlines + tracts + rurality)
        dat_flow <- fl_cejst %>%
          dplyr::filter(!is.na(GEOID), !is.na(rurality)) %>%
          dplyr::mutate(
            ww       = .data[[ww_col]],
            flow_cfs = .data[[q_col]],
            impacted = is.finite(ww) & ww > 0
          ) %>%
          dplyr::filter(impacted, is.finite(flow_cfs))
        
        if (nrow(dat_flow) >= 1) {
          flow_summary <- dat_flow %>%
            dplyr::group_by(rurality) %>%
            dplyr::summarise(
              n_segments   = dplyr::n(),
              med_flow_cfs = stats::median(flow_cfs, na.rm = TRUE),
              iqr_low_cfs  = stats::quantile(flow_cfs, 0.25, na.rm = TRUE),
              iqr_high_cfs = stats::quantile(flow_cfs, 0.75, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            dplyr::arrange(dplyr::desc(med_flow_cfs))
          
          cat("\nDenominator streamflow by rurality (impacted segments only):\n")
          print(flow_summary, row.names = FALSE)
          cat(
            "\nFlows are in cfs for", month_labels[m],
            "and use the", src, "denominator.\n"
          )
        } else {
          cat("\nNo impacted segments with finite flow for this month / filter.\n")
        }
      }
    }
  })
  
  # ---- Download: Equity data (tract-level values + bins) --------------------
  output$download_cejst_data <- downloadHandler(
    filename = function() {
      paste0("ReSourceTracker_equity_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- cejst_binned()
      req(df)
      
      mode_now    <- input$cejst_mode %||% "Prevalence"
      denom_now   <- if (identical(mode_now, "Severity")) (input$cejst_src %||% "Gage-based") else NA_character_
      month_code  <- input$cejst_month
      month_label <- month_labels[month_code]
      rural_filt  <- input$rurality %||% "All"
      var_name    <- input$cejst_var_box %||% NA_character_
      
      df_out <- df %>%
        dplyr::mutate(
          analysis_mode      = mode_now,
          denominator_source = denom_now,
          month_code         = month_code,
          month_label        = unname(month_label),
          cejst_variable     = var_name,
          rurality_filter    = rural_filt
        ) %>%
        dplyr::select(
          analysis_mode, denominator_source,
          month_code, month_label,
          cejst_variable, rurality_filter,
          GEOID, var_bin, metric_value,
          dplyr::everything()
        )
      
      readr::write_csv(df_out, file)
    }
  )
  
  # ---- Census tract info panel (Map tab) ------------------------------------
  output$tract_info <- renderPrint({
    gid <- selected_tract_id()
    if (is.null(gid)) {
      cat("Click a census tract polygon to see \ntract-level wastewater metrics here.\n")
      return(invisible(NULL))
    }
    
    # Get tract-level CEJST attributes (rurality, low_income, status, etc.)
    tr_row <- tracts_sf[tracts_sf$GEOID == gid, ]
    if (nrow(tr_row) == 0) {
      cat("Selected tract not found in tract table.\n")
      return(invisible(NULL))
    }
    
    # Which month and denominator are currently selected in the MAP tab?
    m   <- input$month %||% "06"
    src <- input$src   %||% "Gage-based"
    
    pct_col <- if (identical(src, "Gage-based")) {
      paste0("pct_mean_", m)
    } else {
      paste0("pct_mean_", m, "_nhd")
    }
    ww_col <- paste0("WW_cum_mean_", m)
    
    # All segments in this tract
    seg <- fl_cejst %>%
      dplyr::filter(GEOID == gid)
    
    if (!nrow(seg)) {
      cat("No stream segments intersect this tract.\n")
      return(invisible(NULL))
    }
    
    if (!all(c(pct_col, ww_col) %in% names(seg))) {
      cat("Missing wastewater columns for this month / denominator.\n")
      return(invisible(NULL))
    }
    
    seg <- seg %>%
      dplyr::mutate(
        pct = suppressWarnings(as.numeric(.data[[pct_col]])),
        ww  = suppressWarnings(as.numeric(.data[[ww_col]])),
        impacted = is.finite(ww) & ww > 0   # prevalence uses WW>0
      )
    
    n_seg_total <- nrow(seg)
    n_imp       <- sum(seg$impacted, na.rm = TRUE)
    prevalence  <- n_imp > 0
    
    # Severity: use %WW among impacted segments 
    seg_imp <- seg %>% dplyr::filter(impacted, is.finite(pct))
    if (nrow(seg_imp) > 0) {
      sev_med <- stats::median(seg_imp$pct, na.rm = TRUE)
      sev_p90 <- stats::quantile(seg_imp$pct, 0.9, na.rm = TRUE)
      sev_max <- max(seg_imp$pct, na.rm = TRUE)
    }
    
    # Formatted output
    cat("Census tract:", gid, "\n")
    cat("Rurality:    ", tr_row$rurality[1] %||% "NA", "\n")
    cat("Low income?: ", ifelse(isTRUE(tr_row$low_income[1]), "Yes",
                                ifelse(isFALSE(tr_row$low_income[1]), "No", "NA")), "\n")
    
    cat("\nMonth:        ", month_labels[m], "\n")
    cat("Denominator:  ", ifelse(identical(src, "Gage-based"),
                                 "Gage-based streamflow",
                                 "NHDPlus-based streamflow"), "\n\n")
    
    cat("Segments intersecting tract: ", n_seg_total, "\n")
    cat("Segments with WW>0 (impacted): ", n_imp, "\n")
    cat("Prevalence this month: ",
        ifelse(prevalence, "Impacted", "Not impacted"),
        "\n\n")
    
    if (nrow(seg_imp) == 0) {
      cat("No impacted segments to compute severity (%WW) \nfor this tract.\n")
    } else {
      cat("Severity among impacted segments (", nrow(seg_imp), " segments):\n", sep = "")
      cat("  Median %WW: ", round(sev_med, 1), "%\n", sep = "")
      cat("  90th pct %WW: ", round(sev_p90, 1), "%\n", sep = "")
      cat("  Max %WW: ", round(sev_max, 1), "%\n", sep = "")
    }
  })
}

# --- Launch App -------------------------------------------------------------
shinyApp(ui = ui, server = server)

 

