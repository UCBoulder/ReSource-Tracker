# ReSource-Tracker
*ReSource Tracker: An Open-Source Approach to Tracking Treated Wastewater Effluent in Streams with Applications to Environmental Justice*

This is the repository for ReSource Tracker. This tool includes an associated publication, which can be accessed at `[add paper DOI when published]`. 

## Abstract   
*The abstract for ReSource Tracker will be updated at a later date*   

## Key Features
- Calculates monthly percent wastewater contribution in surface waters
- Supports watershed flow inputs from historical flow data and current gage data
- Incorporates actual wastewater discharge data
- Supports census-tract level analysis
- Displays results using an interactive, clickable basemap and summary plots with statistical analysis
- Provides a unique, environmental justice screening methodology utilizing CEJST variables
- Includes equity-focused views for disadvantaged communites 
- Is modular, allowing users to change the region of interest *(The current version is built around Colorado)*

## How it works
ReSource Tracker uses watershed flow data, wastewater treatment discharge data, and environmental justice data to estimate wastewater influence across watersheds. The app calculates the percentage of wastewater in each watershed on a monthly basis and displays the results through interactive maps and summary plots. 

### Usage
Upon downloading the repository, run the Shiny app in app.R. Running the app will download a required geopackage from Zenodo *(this is expanded on in **Data Dependencies**)*. The app calculates the underlying math prior to launching the user interface and requires a few moments for computation. The app allows the user to explore the base map and all equity screening tools. 

## Data Dependencies
The data used for the app is stored in 'data-raw'. This folder contains all required data EXCEPT 'Flowlines_Colorado_Filtered.gpkg' This data is downloaded from Zenodo upon launching the app. The DOI is: [https://doi.org/10.5281/zenodo.23444014](https://doi.org/10.5281/zenodo.23444014). 

### Data Inputs
The data files included in `data-raw` are as follows:   

| Data Type | File Name | Description | Expected Format | Source | Notes |
|---|---|---|---|---|---|
| Rural/Urban Census Tracts | `CO_2010_URBANRURAL_TRACTS.csv` | Data specificying rural/urban designation for Colorado tracts. | CSV | `[Add source]` | 2010 data |
| Seasonal Gage Flow | `CO_Seasonal_Flow_Summary.csv` | Monthly summary for flows from gages in Colorado. Includes avg, stdev, max, min, n observations, and location per gage. | USGS NWIS | USGS NWIS | Last updated February 2026| 
| CEJST Screening | `ColoradoCEJSTData.xlsx` | CEJST data for Colorado, binned by census ID. Categories include percentages and booleans. | XLSX | CEQ | 2010 data |
| Water Treatment Plant layer | `filteredwtpCO.csv` | Water Treatment Plant ID, service size, and location estimate. | CSV | USEPA SDWIS | Last update?|
| Colorado Flowlines | `Flowlines_Colorado_Filtered.gpkg` | Geopackage with flowlines, basin size for Colorado. | GPKG | USGS NHDPlus HR | Last update? | 
| Wastewater Treatment Locations and Discharge | `matchedDMRwwtp_19-24.csv` | DMR data for WWTPs from 2019 to 2024. Includes location and discharge volume by month. | CSV | EPA DMR | Built in `[preprocessing-scripts]` | 
| POTW for Colorado | `PermittedPOTW.xlsx` | Facility name, NPDES permit ID, longitude/latitude, and design flow, publicly-owned WWTPs. Matched to DMR data | XLSX | USEPA ICIS-NPDES | Last updated February 2026 |
| Geopackage for streamflow gages | `processed_gages_sf_combinedV4.gpkg` | Geopackage simple feature containing gages snapped to nearest flowline | GPKG | USGS | Last updated February 2026 | 

### Updating DMR Data
This repository ships with DMR data from 2019 to 2024 precalculated and matched to WWTPs (`matchedDMRwwtp_19-24.csv`). This data can be rebiult (for inclusion of newer data or different regions) utilizing the `preprocessing-scripts` pipeline. `DMRAnalysis_FY2019-FY2024.R` can be used to reconstruct this output, based on raw DMR inputs and POTW data. Notably, this current version ships with a [Zenodo](https://doi.org/10.5281/zenodo.23444014) integration, which can pull the DMR data for Colorado for 2019-2024 (reliant on function outlined in `DMR_download.R`). This integration can be modified by the user to apply their own DMR data, if needed. 

### Updating Rurality Data
This repository ships with ruraliy data from the 2010 Colorado census (`CO_2010_URBANRURAL_TRACTS.csv`). This data can be rebuilt  utilizing the `preproccessing-scripts/CensusBasedRurality.R` file. 

## Rurality Screening
The app can conduct analysis regarding the impact of rurality (e.g., urban or rural communities) on the prevalence and severity of treated wastewater impact to surface waters. This is conducted utilizing a basin-size adjustment factor to account for the impact of stream size. **Critically, this integration does not imply causality. The results from this analysis must be interpreted as high-level screening, not explicit causation. Please use with caution.** 

## Environmental Justice Screening
The app includes native support for environmental justice screening, based on criteria from the Climate and Economic Justice Screening Tool (CEJST). The app overlays census-level data from the CEJST and can conduct comparisons based on wastewater influence. **Critically, this integration does not imply causality. The results from this analysis must be interpreted as high-level screening, not explicit causation. Please use with caution.** The current vesion supports 49 unique CEJST indicators. 

## Installation
The app is built in R, and depends on several packages, including:

 - `shiny`
 - `tidyverse`
 - `readxl`
 - `leaflet`
 - `sf`
 - `mapgl`
 - `tigris`
 - `vidisisLite`
 - `here`
 - `jsonlite`

 Additionally, `preprocessing-scripts` requires the following:

 - `lubridate`
 - `dplyr`
 - `readr`
 - `tidycensus`

## Limitations
*The limitations for ReSource Tracker will be updated at a later date*   

## Acknowledgements
*The acknowledgements for ReSource Tracker will be updated at a later date*   

## License
This project is licensed under the MIT license, Copyright (c) 2026 Kylie Boenisch-Oakes. See `LICENSE` for additional information. Please cite the associated publication as `[add paper citation when published]`. 

## Citations and Data Availability. 
The data used in this tool is available at a national level, and can be accessed from the resources below:   

- Climate and Economic Justice Screening Tool. Climate and Economic Justice Screening Tool. https://screeningtool.geoplatform.gov   
- NHDPlus High Resolution | U.S. Geological Survey. https://www.usgs.gov/national-hydrography/nhdplus-high-resolution
- TNM Download v2. https://apps.nationalmap.gov/downloader/
- USEPA. Data Downloads | ECHO | US EPA. https://echo.epa.gov/tools/data-downloads
- Custom Search | ECHO | US EPA. https://echo.epa.gov/trends/loading-tool/get-data/custom-search
- US EPA, O. Safe Drinking Water Information System (SDWIS) Federal Reporting Services. https://www.epa.gov/ground-water-and-drinking-water/safe-drinking-water-information-system-sdwis-federal-reporting 

Additionally, the streamflow data from USGS gages can be automatted using the dataRetrieval package in R. Additional details regarding this package can be found at:   

- De Cicco, L. A.; Hirsch, R. M.; Lorenz, D.; Watkins, W. D.; Johnson, M. dataRetrieval: R Packages for Discovering and Retrieving Water Data Available from Federal Hydrologic Web Services, 2025. doi:10.5066/P9X4L3GE



