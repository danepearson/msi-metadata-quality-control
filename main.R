# Notes:
# - Parses 1.9 GB of imzML metadata for two datasets in ~30 seconds end to end (parse, checks, and figures), with no binary access.
# - Switched to using dplyr's bind_rows() to stack the spectrum_df rows over rbind to better handle more datasets as it scales.
# - Added error handling and validation before and after parseImzML()

library(CardinalIO)
library(dplyr)
library(ggplot2)
library(jsonlite)

source("R/parse.R")
source("R/qc-checks.R")
source("R/plots.R")

# Setup:

# Dataset registry: the single place to add or remove datasets.
DATASETS <- list(
  b07f13aa = list(
    imzml = "data/b07f13aa/VAN0103-LK-4-27_region2_roi-all_mz-fix-raw-TIC.imzML",
    metadata = "data/b07f13aa/metadata.json"
  ),
  b84bb453 = list(
    imzml = "data/b84bb453/VAN0031-LK-2-80-IMS_lipids_pos.imzML",
    metadata = "data/b84bb453/metadata.json"
  )
)

# Get working directory (should be where this script lives)
script_dir <- getwd()

# Make figures/ folder
figures_path <- "figures"
if (!dir.exists(figures_path)) {
  dir.create(figures_path)
}

# ------------------------------------------

# Main function which runs imzML data summarization, running related sub-functions:
summarize_imzML_data <- function(dataset_ID, paths) {
  message("Parsing ", dataset_ID)

  parsed_data <- extract_imzML_data(file = paths$imzml)
  spectrum_df <- create_spectrum_summary(parsed_data, dataset_ID = dataset_ID)

  # Attributes CardinalIO does not expose, read from the file header
  header_attrs <- extract_header_attrs(paths$imzml)

  # Scan pattern is needed by the detector, so pull it before the run summary
  scan_pattern <- extract_scan_pattern(parsed_data)

  # Detect regions before the run summary, which reports region count
  # The detector returns the table sorted by Index with Region attached;
  # use that table downstream, not the one passed in.
  detected <- detect_regions(spectrum_df, dataset_ID, scan_pattern)
  spectrum_df <- detected$spectrum_table
  region_df <- detected$region_summary

  # The declared grid should match the extent of a single region
  scan_settings <- parsed_data$scanSettingsList[[1]]
  declared_x <- as.numeric(extract_cv_value(scan_settings, "IMS:1000042"))
  declared_y <- as.numeric(extract_cv_value(scan_settings, "IMS:1000043"))

  grid_matches_region <- any(
    region_df$X_max == declared_x & region_df$Y_max == declared_y,
    na.rm = TRUE
  )

  mz_window_df <- get_json_mz_window(paths$metadata, dataset_ID)

  run_info_df <- create_run_summary(
    parsed_data,
    dataset_ID,
    run_id = header_attrs$run_id,
    data_processing_ref = header_attrs$data_processing_ref,
    n_regions = nrow(region_df),
    grid_matches_region = grid_matches_region
  )

  return(list(
    spectrum_df = spectrum_df,
    run_info_df = run_info_df,
    region_df = region_df,
    mz_window_df = mz_window_df
  ))
}

# Call summarize_imzML_data() on the datasets defined in DATASETS at top
# Iterates over names and paths together, so each call gets both its ID and its file path
summaries <- Map(
  f = summarize_imzML_data,
  dataset_ID = names(DATASETS),
  paths = DATASETS
)

# Stack the ...$spectrum_df dataframe rows into one long table
all_spectra <- bind_rows(lapply(summaries, function(dataset_result) {
  dataset_result$spectrum_df
}))
rownames(all_spectra) <- NULL

# Stack the run_info_df dataframe rows into one long table
all_run_info <- bind_rows(lapply(summaries, function(dataset_result) {
  dataset_result$run_info_df
}))
rownames(all_run_info) <- NULL

# Stack the region_df dataframe rows into one long table
all_regions <- bind_rows(lapply(summaries, function(dataset_result) {
  dataset_result$region_df
}))
rownames(all_regions) <- NULL

# Stack the mz_window_df dataframe rows into one long table
all_mz_windows <- bind_rows(lapply(summaries, function(dataset_result) {
  dataset_result$mz_window_df
}))
rownames(all_mz_windows) <- NULL

# Add the declared-grid coverage number to the region table
all_regions <- compute_coverage(all_regions, all_run_info)

# Compare each region's observed m/z range against its declared window
mz_range_checks <- check_mz_range(all_spectra, all_mz_windows)

# Run the order checks per region. The region boundary in b07f13aa made three
# of four checks fail when run per dataset; within a single region they should pass
order_checks <- Map(
  f = function(id, region) {
    verify_acquisition_quality(
      spectrum_table = filter(all_spectra, Dataset_ID == id, Region == region),
      scan_pattern_string = all_run_info$Scan_pattern[
        all_run_info$Dataset_ID == id
      ]
    )
  },
  id = all_regions$Dataset_ID,
  region = all_regions$Region
)
names(order_checks) <- paste0(all_regions$Dataset_ID, "_r", all_regions$Region)

# ------
# Results

cat("\n=== Datasets parsed ===\n")
cat(
  nrow(all_spectra),
  "spectra across",
  nrow(all_run_info),
  "datasets and",
  nrow(all_regions),
  "regions\n"
)

cat("\n=== Run summary ===\n")
# Run_ID and Full_UUID are long and make the table unreadable at the console;
# Run_ID is printed separately below.
print(select(all_run_info, -Run_ID, -Full_UUID))

cat("\n=== Run IDs ===\n")
# The run id is a filesystem path written by the acquisition software. It is
# the only place laterality appears in b07f13aa, where it disagrees with the
# organ field registered in metadata.json.
for (i in seq_len(nrow(all_run_info))) {
  cat(
    all_run_info$Dataset_ID[i],
    ":\n  ",
    all_run_info$Run_ID[i],
    "\n",
    sep = ""
  )
}

cat("\n=== Regions ===\n")
print(as.data.frame(all_regions))

cat("\n=== M/Z Windows ===\n")
all_mz_windows

cat("\n=== m/z range checks (per region) ===\n")
print(as.data.frame(mz_range_checks))

cat("\n=== Acquisition order checks (per region) ===\n")
# verify_acquisition_quality() returns a list per region; collapse the four
# verdicts into one row each so the result is readable at a glance.
order_check_summary <- bind_rows(lapply(
  names(order_checks),
  function(region_name) {
    checks <- order_checks[[region_name]]
    data.frame(
      Region = region_name,
      Rows_contiguous = checks$rows_contiguous,
      Y_non_decreasing = checks$y_non_decreasing,
      X_strictly_increasing = checks$x_strictly_increasing,
      No_duplicate_pixels = checks$no_duplicate_pixels,
      stringsAsFactors = FALSE
    )
  }
))
print(order_check_summary)

# All four checks passing in every region is what validates the region
# detector: the same checks fail on b07f13aa when run across the whole file,
# because the region boundary sits inside the data being checked.
if (all(order_check_summary[, -1])) {
  cat("\nAll acquisition order checks passed in all regions.\n")
} else {
  cat("\nSome acquisition order checks failed; see table above.\n")
}

# Registered organ per dataset, read from metadata.json
all_organs <- bind_rows(Map(
  f = function(id, paths) get_registered_organ(paths$metadata, id),
  id = names(DATASETS),
  paths = DATASETS
))

# Compare the run id's implied laterality against the registered organ
identifier_checks <- check_identifier_consistency(all_run_info, all_organs)

cat("\n=== Organ laterality checks summary ===\n")
print(as.data.frame(identifier_checks))

# Which x-coordinates (scan-columns) the unmeasured positions sit in, per region
missing_summary <- summarise_missing_positions(all_spectra, all_regions)

cat("\n=== Unmeasured positions (per region) ===\n")
print(as.data.frame(missing_summary))

# ------------
# Plots

p_occupancy <- plot_occupancy_map(all_spectra, all_regions, all_run_info)
ggsave(
  file.path(figures_path, "occupancy_map.png"),
  plot = p_occupancy,
  width = 10,
  height = 5,
  dpi = 150
)

p_missing <- plot_missing_positions(all_spectra, all_regions)
ggsave(
  file.path(figures_path, "missing_positions.png"),
  plot = p_missing,
  width = 10,
  height = 6,
  dpi = 150
)

# TIC and ArrayLength per scan line, the inputs for the drift plots and the trend metrics
summarized_scanned_lines <- summarize_scan_lines(all_spectra)

# One row per region: the ratio of the last window of scan lines to the first
tic_trend <- compute_tic_trend(summarized_scanned_lines)

cat("\n=== TIC trend (per region) ===\n")
print(as.data.frame(tic_trend))

p_tic_drift <- plot_TIC_drift(summarized_scanned_lines)
ggsave(
  file.path(figures_path, "tic_drift.png"),
  plot = p_tic_drift,
  width = 12,
  height = 6,
  dpi = 200
)

# One row per region: the ratio of the last window of scan lines to the first
arrayLength_trend <- compute_arrayLength_trend(summarized_scanned_lines)

cat("\n=== Array length trend (per region) ===\n")
print(as.data.frame(arrayLength_trend))

p_arrayLength_drift <- plot_arrayLength_drift(summarized_scanned_lines)
ggsave(
  file.path(figures_path, "peaks_per_pixel_drift.png"),
  plot = p_arrayLength_drift,
  width = 12,
  height = 6,
  dpi = 200
)

p_tic_spatial <- plot_TIC_spatial(all_spectra)
ggsave(
  file.path(figures_path, "tic_spatial.png"),
  plot = p_tic_spatial,
  width = 12,
  height = 4.5,
  dpi = 200
)

p_arrayLength_spatial <- plot_arrayLength_spatial(all_spectra)
ggsave(
  file.path(figures_path, "peaks_per_pixel_spatial.png"),
  plot = p_arrayLength_spatial,
  width = 12,
  height = 4.5,
  dpi = 200
)

# ------------
# Join all tables together

# First, add a region column separate from dataset ID in the order check summary table
order_check_summary <- bind_rows(lapply(
  seq_along(order_checks),
  function(i) {
    checks <- order_checks[[i]]
    data.frame(
      Dataset_ID = all_regions$Dataset_ID[i],
      Region = all_regions$Region[i],
      Rows_contiguous = checks$rows_contiguous,
      Y_non_decreasing = checks$y_non_decreasing,
      X_strictly_increasing = checks$x_strictly_increasing,
      No_duplicate_pixels = checks$no_duplicate_pixels,
      stringsAsFactors = FALSE
    )
  }
))

combined_data <- all_regions |>
  left_join(
    select(mz_range_checks, -N_spectra),
    by = c("Dataset_ID", "Region")
  ) |>
  left_join(
    select(tic_trend, Dataset_ID, Region, TIC_last_over_first),
    by = c("Dataset_ID", "Region")
  ) |>
  left_join(
    select(arrayLength_trend, Dataset_ID, Region, ArrayLength_last_over_first),
    by = c("Dataset_ID", "Region")
  ) |>
  left_join(order_check_summary, by = c("Dataset_ID", "Region")) |>
  left_join(
    select(
      all_run_info,
      Dataset_ID,
      N_regions,
      Grid_matches_region,
      Scan_pattern
    ),
    by = "Dataset_ID"
  )

print(as.data.frame(combined_data))
