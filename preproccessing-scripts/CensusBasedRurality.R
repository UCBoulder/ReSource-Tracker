# ==============================================================================
# ReSource Tracker census-based rurality

# This script defines the urban and rural census tracts for a given region.
# It currently defaults to Colorado, although other states can be identified. 

# Preprocessing data is only required if updating/revising a data source. This script is only
# included for completeness and is not required for current functionality
# ==============================================================================

# --- Libraries ---------------------------------------------------------------
library(tigris)
library(tidycensus)


urban10 <- tigris::urban_areas(state= 08, year = 2010, class = "sf") %>%
  st_make_valid() %>%
  st_transform(st_crs(tracts10))

sf1_2010 <- load_variables(2010, "sf1", cache = TRUE)
ur10_raw <- get_decennial(
  geography = "tract",
  state     = "08",   # This corresponds to Colorado. User can change to include other areas. 
  year      = 2010,
  variables = c(
    total = "P002001",
    urban = "P002002",
    rural = "P002005"
  ),
  sumfile  = "sf1"
)

ur10 <- ur10_raw %>%
  select(GEOID, variable, value) %>%
  tidyr::pivot_wider(names_from = variable, values_from = value) %>%
  mutate(
    total       = as.numeric(total),
    urban       = as.numeric(urban),
    rural       = as.numeric(rural),
    urban_share = if_else(total > 0, urban / total, NA_real_),
    rurality    = if_else(urban_share >= 0.5, "Urban", "Rural")
  ) %>%
  select(GEOID, rurality)

write.csv(ur10, here("data-raw", "CO_2010_URBANRURAL_TRACTS.csv"))
