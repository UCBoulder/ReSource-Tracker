# ==============================================================================
# ReSource Tracker DMR Zenodo download helper

# Download a zip file from a Zenodo record, unzip it, and read every CSV inside.
# Supply either: record_id - the Zenodo record number
# OR url - a direct download link to the zip
#   reader    - function that takes a file path and returns a data frame
# This script does not download the zip or CSVs permanently. 

# Preprocessing data is only required if updating/revising a data source. This script is only
# included for completeness and is not required for current functionality
# ==============================================================================


read_zenodo_zip <- function(record_id = NULL,
                            url = NULL,
                            zip_name = NULL,
                            cache_dir = file.path(tempdir(), "zenodo"),
                            overwrite = FALSE,
                            reader = function(path) readr::read_csv(path, show_col_types = FALSE)) {

  if (is.null(record_id) && is.null(url)) {
    stop("Provide either `record_id` or `url`.")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  checksum <- NULL

  # 1. Look up the zip's download link via the Zenodo API ---------------------
  if (is.null(url)) {
    meta  <- jsonlite::fromJSON(paste0("https://zenodo.org/api/records/", record_id))
    files <- meta$files
    zips  <- files[grepl("\\.zip$", files$key, ignore.case = TRUE), ]
    if (!is.null(zip_name)) zips <- zips[zips$key == zip_name, ]

    if (nrow(zips) == 0) stop("No matching .zip file found in record ", record_id)
    if (nrow(zips) > 1) {
      stop("Multiple zips in record; set `zip_name` to one of: ",
           paste(zips$key, collapse = ", "))
    }

    url      <- zips$links$self
    zip_name <- zips$key
    checksum <- sub("^md5:", "", zips$checksum)
  }
  if (is.null(zip_name)) zip_name <- basename(sub("/content$", "", url))

  zip_path   <- file.path(cache_dir, zip_name)
  unzip_dir  <- file.path(cache_dir, tools::file_path_sans_ext(zip_name))

  # 2. Download (skip if already cached) --------------------------------------
  if (overwrite || !file.exists(zip_path)) {
    old <- options(timeout = max(600, getOption("timeout")))
    on.exit(options(old), add = TRUE)
    message("Downloading ", zip_name, " ...")
    utils::download.file(url, zip_path, mode = "wb", quiet = TRUE)
  }

  # Verify integrity when Zenodo gave us a checksum
  if (!is.null(checksum) && unname(tools::md5sum(zip_path)) != checksum) {
    file.remove(zip_path)
    stop("Checksum mismatch: download may be corrupted. Try again.")
  }

  # 3. Unzip -------------------------------------------------------------------
  if (overwrite || !dir.exists(unzip_dir)) {
    utils::unzip(zip_path, exdir = unzip_dir)
  }

  # 4. Read every CSV (searching subfolders, ignoring macOS junk) -------------
  csvs <- list.files(unzip_dir, pattern = "\\.csv$", recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  csvs <- csvs[!grepl("__MACOSX|/\\._", csvs)]
  if (length(csvs) == 0) stop("No CSV files found in ", zip_name)

  out <- lapply(csvs, reader)
  names(out) <- tools::file_path_sans_ext(basename(csvs))
  message("Read ", length(out), " CSV(s): ", paste(names(out), collapse = ", "))
  out
}