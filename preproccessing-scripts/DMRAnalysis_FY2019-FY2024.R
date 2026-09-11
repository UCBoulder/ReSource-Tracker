# --- Libraries ---------------------------------------------------------------
library(lubridate)
library(dplyr)
library(readr)

# --- Column spec for raw DMR files ------------------------------------------
dmr_headers <- c(
  "ACTIVITY_ID","EXTERNAL_PERMIT_NMBR","VERSION_NMBR","PERM_FEATURE_ID","PERM_FEATURE_NMBR",
  "PERM_FEATURE_TYPE_CODE","LIMIT_SET_ID","LIMIT_SET_DESIGNATOR","LIMIT_SET_SCHEDULE_ID",
  "LIMIT_ID","LIMIT_BEGIN_DATE","LIMIT_END_DATE","NMBR_OF_SUBMISSION","NMBR_OF_REPORT",
  "PARAMETER_CODE","PARAMETER_DESC","MONITORING_LOCATION_CODE","STAY_TYPE_CODE",
  "LIMIT_VALUE_ID","LIMIT_VALUE_TYPE_CODE","LIMIT_VALUE_NMBR","LIMIT_UNIT_CODE",
  "LIMIT_UNIT_DESC","STANDARD_UNIT_CODE","STANDARD_UNIT_DESC","LIMIT_VALUE_STANDARD_UNITS",
  "STATISTICAL_BASE_CODE","STATISTICAL_BASE_TYPE_CODE","LIMIT_VALUE_QUALIFIER_CODE",
  "OPTIONAL_MONITORING_FLAG","LIMIT_SAMPLE_TYPE_CODE","LIMIT_FREQ_OF_ANALYSIS_CODE",
  "STAY_VALUE_NMBR","LIMIT_TYPE_CODE","DMR_EVENT_ID","MONITORING_PERIOD_END_DATE",
  "DMR_SAMPLE_TYPE_CODE","DMR_FREQ_OF_ANALYSIS_CODE","REPORTED_EXCURSION_NMBR",
  "DMR_FORM_VALUE_ID","VALUE_TYPE_CODE","DMR_VALUE_ID","DMR_VALUE_NMBR","DMR_UNIT_CODE",
  "DMR_UNIT_DESC","DMR_VALUE_STANDARD_UNITS","DMR_VALUE_QUALIFIER_CODE","VALUE_RECEIVED_DATE",
  "DAYS_LATE","NODI_CODE","EXCEEDENCE_PCT","NPDES_VIOLATION_ID","VIOLATION_CODE",
  "RNC_DETECTION_CODE","RNC_DETECTION_DATE","RNC_RESOLUTION_CODE","RNC_RESOLUTION_DATE"
)

read_dmr_file <- function(path) {
  read.csv(path, header = FALSE, col.names = dmr_headers, colClasses = "character")
}

# --- Read & combine ----------------------------------------------------------
dmr_2019 <- read_dmr_file('/Users/kylieboenischoakes/Library/CloudStorage/OneDrive-UCB-O365/General/Boenisch-Oakes, Kylie/Research/Chapter 4/DMR_CO_2019.csv')
dmr_2020 <- read_dmr_file('/Users/kylieboenischoakes/Library/CloudStorage/OneDrive-UCB-O365/General/Boenisch-Oakes, Kylie/Research/Chapter 4/DMR_CO_2020.csv')
dmr_2021 <- read_dmr_file('/Users/kylieboenischoakes/Library/CloudStorage/OneDrive-UCB-O365/General/Boenisch-Oakes, Kylie/Research/Chapter 4/DMR_CO_2021.csv')
dmr_2022 <- read_dmr_file('/Users/kylieboenischoakes/Library/CloudStorage/OneDrive-UCB-O365/General/Boenisch-Oakes, Kylie/Research/Chapter 4/DMR_CO_2022.csv')
dmr_2023 <- read_dmr_file('/Users/kylieboenischoakes/Library/CloudStorage/OneDrive-UCB-O365/General/Boenisch-Oakes, Kylie/Research/Chapter 4/DMR_CO_2023.csv')
dmr_2024 <- read.csv('/Users/kylieboenischoakes/Library/CloudStorage/OneDrive-UCB-O365/General/Boenisch-Oakes, Kylie/Research/Chapter 4/DMR_CO_2024.csv',
                     colClasses = "character") # already has headers

dmr_all <- bind_rows(dmr_2019, dmr_2020, dmr_2021, dmr_2022, dmr_2023, dmr_2024)

# --- Normalize IDs & dates ---------------------------------------------------
dmr_all$EXTERNAL_PERMIT_NMBR <- trimws(as.character(dmr_all$EXTERNAL_PERMIT_NMBR))
dmr_all <- dmr_all %>% filter(!is.na(MONITORING_PERIOD_END_DATE))

# NOTE: wwtp_pts must already exist in workspace with DesignFlow (MGD) & matched_COMID (from main code file)
wwtp_pts$EXTERNAL_PERMIT_NMBR <- trimws(as.character(wwtp_pts$EXTERNAL_PERMIT_NMBR))

dmr_all$MONITORING_PERIOD_END_DATE <- as.Date(dmr_all$MONITORING_PERIOD_END_DATE, format = "%m/%d/%Y")

# --- Filter to flow, effluent gross, valid data ------------------------------
dmr_effluent <- dmr_all %>%
  filter(PARAMETER_CODE == "50050",
         PERM_FEATURE_TYPE_CODE == "EXO",          # external outfalls only
         MONITORING_LOCATION_CODE == "1",          # effluent gross only
         !nzchar(NODI_CODE)) %>%                   # drop no-discharge codes
  mutate(DMR_VALUE_STANDARD_UNITS = suppressWarnings(as.numeric(DMR_VALUE_STANDARD_UNITS))) %>%
  filter(is.finite(DMR_VALUE_STANDARD_UNITS), DMR_VALUE_STANDARD_UNITS > 0) %>%
  filter(grepl("AVG", toupper(STATISTICAL_BASE_TYPE_CODE)),   # keep averages
         !grepl("MAX", toupper(STATISTICAL_BASE_TYPE_CODE)))  # drop maxima

# --- Monthly facility flows (MGD) -------------------------------------------
dmr_monthly_facility <- dmr_effluent %>%
  mutate(period_month = lubridate::floor_date(MONITORING_PERIOD_END_DATE, "month")) %>%
  group_by(EXTERNAL_PERMIT_NMBR, period_month) %>%
  summarize(flow_mgd = mean(DMR_VALUE_STANDARD_UNITS, na.rm = TRUE), .groups = "drop")

# Assign month code "01"–"12"
dmr_monthly_facility <- dmr_monthly_facility %>%
  mutate(
    mon_code = sprintf("%02d", lubridate::month(period_month))
  )

# --- Per-permit monthly stats (mean ± sd) in cfs ----------------------------
# 1 MGD = 1.547 cfs
dmr_monthly_stats <- dmr_monthly_facility %>%
  group_by(EXTERNAL_PERMIT_NMBR, mon_code) %>%
  summarize(
    mean_dmr_cfs = mean(flow_mgd, na.rm = TRUE) * 1.547,
    sd_dmr_cfs   = sd  (flow_mgd, na.rm = TRUE) * 1.547,
    n_obs        = n(),
    .groups = "drop"
  ) %>%
  mutate(
    low_dmr_cfs  = pmax(mean_dmr_cfs - sd_dmr_cfs, 0),
    high_dmr_cfs = mean_dmr_cfs + sd_dmr_cfs
  )

# --- Join to WWTP metadata (DesignFlow in MGD -> cfs) -----------------------
matched_outfalls <- wwtp_pts %>%
  filter(EXTERNAL_PERMIT_NMBR %in% dmr_monthly_stats$EXTERNAL_PERMIT_NMBR) %>%
  mutate(
    Flow_cfs = as.numeric(DesignFlow) * 1.547
  )

wwtp_dmr_monthly <- matched_outfalls %>%
  select(EXTERNAL_PERMIT_NMBR, Facility, DesignFlow, Flow_cfs, matched_COMID) %>%
  inner_join(dmr_monthly_stats, by = "EXTERNAL_PERMIT_NMBR", relationship = "many-to-many") %>%
  mutate(ratio_dmr_to_design = mean_dmr_cfs / Flow_cfs)

options(scipen = 999)

# --- Export -----------------------------------------------------------------
write_csv(
  wwtp_dmr_monthly,
  "/Users/kylieboenischoakes/Library/CloudStorage/OneDrive-UCB-O365/General/Boenisch-Oakes, Kylie/Research/Chapter 4/matchedDMRwwtpupdated6.csv"
)
