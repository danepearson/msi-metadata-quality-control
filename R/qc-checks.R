library(dplyr)

# Compute number of regions in dataset, and their bounds
detect_regions <- function(
  spectrum_table,
  dataset_ID,
  scan_pattern_string = NA_character_
) {
  # Check that dataset_ID is a single, non-empty string
  if (
    !is.character(dataset_ID) || length(dataset_ID) != 1L || !nzchar(dataset_ID)
  ) {
    stop("'dataset_ID' must be a single non-empty string.", call. = FALSE)
  }
  # Required columns, same pattern as verify_acquisition_quality()
  required_cols <- c("Index", "X", "Y")
  missing_cols <- setdiff(required_cols, names(spectrum_table))
  if (length(missing_cols) > 0L) {
    stop(
      "spectrum_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  # Check that spectrum table is not empty
  if (nrow(spectrum_table) == 0) {
    stop("'spectrum_table' table is empty.", call. = FALSE)
  }

  # Warn if the declared pattern isn’t top-down, since the detection rule wouldn’t apply
  if (
    !is.na(scan_pattern_string) &&
      !grepl("top down", scan_pattern_string, fixed = TRUE)
  ) {
    warning(
      "Scan pattern for ",
      dataset_ID,
      " is not top-down ('",
      scan_pattern_string,
      "'); region detection assumes a top-down raster and may be unreliable.",
      call. = FALSE
    )
  }

  ### Detection for top-down raster. Y decreasing means new region
  spectrum_sorted_df <- spectrum_table |> arrange(Index)
  y_diffs <- diff(spectrum_sorted_df$Y)

  # Detect the new region bound
  y_drop_index <- which(y_diffs < 0)
  row_of_drop <- y_drop_index + 1L

  if (length(y_drop_index) > 0L) {
    message(
      length(y_drop_index) + 1L,
      " regions detected in ",
      dataset_ID,
      "; boundaries at index: ",
      paste(spectrum_sorted_df$Index[row_of_drop], collapse = ", ")
    )
  } else {
    message("Only one region detected in ", dataset_ID, ".\n")
  }

  # Create a region vector (like 1 1 1 2 2 2 2 ...)
  is_region_start <- logical(nrow(spectrum_sorted_df))
  is_region_start[row_of_drop] <- TRUE
  region_vector <- cumsum(is_region_start) + 1L

  # Attach to the sorted table, so the column and the rows don't fall out of alignment
  # Other functions should use this table rather than the one passed in
  spectrum_w_region <- spectrum_sorted_df |>
    mutate(Region = region_vector)

  # Coordinate validity, checked per region.
  coord_checks <- spectrum_w_region |>
    group_by(Region) |>
    summarise(
      min_x = min(X, na.rm = TRUE),
      min_y = min(Y, na.rm = TRUE),
      n_negative = sum(X < 0 | Y < 0, na.rm = TRUE),
      .groups = "drop"
    )

  # Negative coordinates are invalid under either convention
  n_negative_total <- sum(coord_checks$n_negative)
  if (n_negative_total > 0) {
    warning(
      n_negative_total,
      " spectra have negative pixel coordinates in ",
      dataset_ID,
      ".",
      call. = FALSE
    )
  }

  # A region whose X and Y disagree on base (one starts at 0, the other at 1)
  # suggests the coordinates were not written consistently
  mixed_base <- coord_checks |>
    filter(min_x <= 1, min_y <= 1, min_x != min_y)
  if (nrow(mixed_base) > 0) {
    warning(
      nrow(mixed_base),
      " region(s) in ",
      dataset_ID,
      " have inconsistent coordinate bases between X and Y.",
      call. = FALSE
    )
  }

  # One row per region, so region summaries stack across datasets the same way
  # spectrum and run summaries do.
  region_summary <- spectrum_w_region |>
    group_by(Region) |>
    summarize(
      Dataset_ID = dataset_ID,
      Spectrum_count = n(),
      Index_min = min(Index),
      Index_max = max(Index),
      X_min = min(X),
      X_max = max(X),
      Y_min = min(Y),
      Y_max = max(Y),
      # The region's own bounding box (also the denominator used to calculate coverage)
      Bounding_box_size = (max(X) - min(X) + 1L) * (max(Y) - min(Y) + 1L),
      .groups = "drop"
    ) |>
    mutate(
      Bounding_box_occupancy = Spectrum_count / Bounding_box_size,
      # Coordinate base is only meaningful for a region that starts at the origin
      Coordinate_base = if_else(
        X_min <= 1 & Y_min <= 1,
        pmin(X_min, Y_min),
        NA_real_
      )
    ) |>
    relocate(Dataset_ID, .before = Region)

  # Return a list: the per-region summary, and the region vector
  return(list(
    spectrum_table = spectrum_w_region,
    region_summary = region_summary
  ))
}

# Add the declared grid-coverage number to the region table
compute_coverage <- function(region_table, run_table) {
  region_table |>
    left_join(
      select(run_table, Dataset_ID, Declared_grid_size),
      by = "Dataset_ID"
    ) |>
    mutate(Declared_grid_coverage = Spectrum_count / Declared_grid_size)
}

# ------
# Within-run checks:

# Verify acquisition quality
verify_acquisition_quality <- function(spectrum_table, scan_pattern_string) {
  # Safety checks
  required_cols <- c("Index", "X", "Y")
  missing_cols <- setdiff(required_cols, names(spectrum_table))
  if (length(missing_cols) > 0L) {
    stop(
      "spectrum_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(spectrum_table) == 0L) {
    stop("spectrum_table has no rows.", call. = FALSE)
  }

  # Check 1: Are rows contiguous in index? (testing single-pass raster assumptions)
  row_spans <- spectrum_table |>
    group_by(Y) |>
    summarize(
      n_spectra = n(),
      index_span = max(Index) - min(Index) + 1L,
      .groups = "drop"
    ) |>
    mutate(is_contiguous = index_span == n_spectra)

  # Count number of non contiguous rows
  n_broken_rows <- sum(!row_spans$is_contiguous)
  rows_contiguous <- n_broken_rows == 0L

  if (!rows_contiguous) {
    warning(
      n_broken_rows,
      " of ",
      nrow(row_spans),
      " rows are not contiguous in index; the file may not be written in acquisition order.",
      call. = FALSE
    )
  }

  # Check 2: is Y non-decreasing across the run?
  spectrum_sorted_df <- spectrum_table |> arrange(Index)
  y_diffs <- diff(spectrum_sorted_df$Y)

  y_decreases <- sum(y_diffs < 0, na.rm = TRUE)
  y_non_decreasing <- y_decreases == 0L

  if (!y_non_decreasing) {
    warning(
      y_decreases,
      " positions where Y decreases; the raster may not be top-down.",
      call. = FALSE
    )
  }

  # Check 3: is X strictly increasing in each row?
  x_ordered <- spectrum_table |>
    arrange(Index) |>
    group_by(Y) |>
    summarise(
      x_increasing = all(diff(X) > 0),
      .groups = "drop"
    )

  # Count number of rows where x not strictly increasing
  n_unordered_rows <- sum(!x_ordered$x_increasing)
  x_strictly_increasing <- n_unordered_rows == 0L

  if (!x_strictly_increasing) {
    warning(
      n_unordered_rows,
      " of ",
      nrow(x_ordered),
      " rows do not have strictly increasing X; the raster may contain duplicate pixels.",
      call. = FALSE
    )
  }

  # --- Check 4: are there duplicate (X, Y) coordinate pairs?
  n_unique_pairs <- spectrum_table |>
    distinct(X, Y) |>
    nrow()

  n_duplicate_rows <- nrow(spectrum_table) - n_unique_pairs
  no_duplicate_pixels <- n_duplicate_rows == 0L

  if (!no_duplicate_pixels) {
    warning(
      n_duplicate_rows,
      " rows of ",
      nrow(spectrum_table),
      " are duplicate X, Y pairs.",
      call. = FALSE
    )
  }

  return(list(
    rows_contiguous = rows_contiguous,
    n_broken_rows = n_broken_rows,
    row_spans = row_spans,
    y_non_decreasing = y_non_decreasing,
    y_decreases = y_decreases,
    x_strictly_increasing = x_strictly_increasing,
    n_unordered_rows = n_unordered_rows,
    x_ordered = x_ordered,
    no_duplicate_pixels = no_duplicate_pixels,
    n_duplicate_rows = n_duplicate_rows
  ))
}

# Pull one field out of a parsed metadata.json and coerce it to numeric.
# JSON equivalent of extract_cv_value(): a missing key returns NULL, and
# as.numeric(NULL) gives numeric(0), which flows on silently and breaks
# something unrelated several steps later.
extract_json_numeric <- function(parsed_json, field, required = TRUE) {
  value <- parsed_json$metadata[[field]]

  if (is.null(value)) {
    if (required) {
      stop(
        "Field 'metadata$",
        field,
        "' is missing from this metadata.json.",
        call. = FALSE
      )
    }
    return(NA_real_)
  }

  # JSON stores these as strings
  numeric_value <- suppressWarnings(as.numeric(value))

  if (is.na(numeric_value)) {
    warning(
      "Field 'metadata$",
      field,
      "' could not be converted to numeric (value: '",
      value,
      "').",
      call. = FALSE
    )
  }

  return(numeric_value)
}

# Function that reads a dataset’s metadata.json and returns its declared low and high m/z, producing a two-row table you join on Dataset_ID
get_json_mz_window <- function(json_path, dataset_ID) {
  # Validate dataset_ID
  if (
    !is.character(dataset_ID) || length(dataset_ID) != 1L || !nzchar(dataset_ID)
  ) {
    stop("'dataset_ID' must be a single non-empty string.", call. = FALSE)
  }

  parsed_json <- parse_json(json_path)

  # Validate metadata.json
  if (is.null(parsed_json$metadata)) {
    stop(
      "No 'metadata' element found in ",
      basename(json_path),
      call. = FALSE
    )
  }

  # Extract and validate declared_mz_low, declared_mz_high foind in the metadata json
  declared_mz_low <- extract_json_numeric(
    parsed_json,
    "mass_to_charge_range_low_value"
  )
  declared_mz_high <- extract_json_numeric(
    parsed_json,
    "mass_to_charge_range_high_value"
  )

  # Check that the low end of the declared m/z ratio is the below the high end
  if (!is.na(declared_mz_low) && !is.na(declared_mz_high)) {
    if (declared_mz_low >= declared_mz_high) {
      warning(
        "Declared m/z window for ",
        dataset_ID,
        " has low (",
        declared_mz_low,
        ") at or above high (",
        declared_mz_high,
        ").",
        call. = FALSE
      )
    }
  }

  # Save a dataframe where each row is a dataset's declared mz ratio low end and its high end
  mz_window_df <- data.frame(
    Dataset_ID = dataset_ID,
    Declared_mz_low = declared_mz_low,
    Declared_mz_high = declared_mz_high,
    stringsAsFactors = FALSE
  )

  return(mz_window_df)
}

# Compare each region's observed m/z range against the window declared in its metadata.json
# The imzML itself declares no scan window (MS:1000500/1000501 are absent),
# so the declaration comes from the consortium metadata rather
# than the instrument's own output, making it a two-source check
check_mz_range <- function(
  spectrum_table,
  mz_window_table,
  warn_above_fraction = 0.01, # warn when this percent of spectra are violations
  mz_violation_tolerance = 0.01 # warn if if a spectrum (pixel) goes outside declared window by more than this amount
) {
  # Validation checks
  required_spectrum <- c("Dataset_ID", "Region", "Lowest_mz", "Highest_mz")
  missing_cols <- setdiff(required_spectrum, names(spectrum_table))
  if (length(missing_cols) > 0L) {
    stop(
      "spectrum_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  required_window <- c("Dataset_ID", "Declared_mz_low", "Declared_mz_high")
  missing_cols <- setdiff(required_window, names(mz_window_table))
  if (length(missing_cols) > 0L) {
    stop(
      "mz_window_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Join both dataframes by Dataset_ID
  spectra_with_mz_df <- left_join(
    spectrum_table,
    mz_window_table,
    by = "Dataset_ID"
  )

  # Check for any missing m/z windows by looking for NAs in the Declared_mz_low col
  n_missing_window <- sum(is.na(spectra_with_mz_df$Declared_mz_low))
  if (n_missing_window > 0L) {
    unchecked <- spectra_with_mz_df |>
      filter(is.na(Declared_mz_low)) |>
      distinct(Dataset_ID) |>
      pull(Dataset_ID)
    warning(
      n_missing_window,
      " spectra have no declared m/z window and were not checked (dataset(s): ",
      paste(unchecked, collapse = ", "),
      ").",
      call. = FALSE
    )
  }

  # count where Lowest_mz < Declared_mz_low
  # count where Highest_mz > Declared_mz_high
  mz_range_df <- spectra_with_mz_df |>
    group_by(Dataset_ID, Region) |>
    summarise(
      Declared_mz_low = first(Declared_mz_low),
      Declared_mz_high = first(Declared_mz_high),
      Observed_mz_low = min(Lowest_mz, na.rm = TRUE),
      Observed_mz_high = max(Highest_mz, na.rm = TRUE),
      N_spectra = n(),
      N_below_declared_low = sum(
        Lowest_mz < (Declared_mz_low - mz_violation_tolerance),
        na.rm = TRUE
      ),
      N_above_declared_high = sum(
        Highest_mz > (Declared_mz_high + mz_violation_tolerance),
        na.rm = TRUE
      ),
      # Counted as its own condition rather than summed from the two above: a
      # spectrum can break both ends, and summing would double-count it.
      N_outside_declared = sum(
        Lowest_mz < (Declared_mz_low - mz_violation_tolerance) |
          Highest_mz > (Declared_mz_high + mz_violation_tolerance),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    mutate(
      # Overlap of the observed and declared intervals, over the declared span.
      # Using the raw observed span would exceed 1 when a run observes outside
      # its declared window, which b84bb453 does on the low end.
      Declared_window_used = (pmin(Observed_mz_high, Declared_mz_high) -
        pmax(Observed_mz_low, Declared_mz_low)) /
        (Declared_mz_high - Declared_mz_low),
      Fraction_outside = N_outside_declared / N_spectra
    )

  # Track how many spectra are violating the declared window
  # A handful of violating spectra is an edge case
  # Most of a region violating the decalred window means the declaration does not describe the acquisition
  substantial <- mz_range_df |> filter(Fraction_outside > warn_above_fraction)
  if (nrow(substantial) > 0L) {
    for (i in seq_len(nrow(substantial))) {
      row <- substantial[i, ]
      warning(
        round(100 * row$Fraction_outside, 1),
        "% of spectra in ",
        row$Dataset_ID,
        " region ",
        row$Region,
        " fall outside the declared m/z window (",
        row$Declared_mz_low,
        "-",
        row$Declared_mz_high,
        "; observed ",
        round(row$Observed_mz_low, 1),
        "-",
        round(row$Observed_mz_high, 1),
        ").",
        call. = FALSE
      )
    }
  }

  return(mz_range_df)
}

summarize_scan_lines <- function(spectrum_table) {
  # Validation checks
  required_cols <- c("Dataset_ID", "Region", "Y", "TIC", "ArrayLength")
  missing_cols <- setdiff(required_cols, names(spectrum_table))
  if (length(missing_cols) > 0L) {
    stop(
      "Spectrum_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(spectrum_table) == 0L) {
    stop("spectrum_table has no rows.", call. = FALSE)
  }

  # Summarize the df's TIC and peaks per pixel (ArrayLength) properties
  summarized_df <- spectrum_table |>
    group_by(Dataset_ID, Region, Y) |>
    summarize(
      Median_TIC = median(TIC, na.rm = TRUE),
      TIC_25_perc = quantile(TIC, probs = 0.25, na.rm = TRUE),
      TIC_75_perc = quantile(TIC, probs = 0.75, na.rm = TRUE),
      Median_ArrayLength = median(ArrayLength, na.rm = TRUE),
      ArrayLength_25_perc = quantile(ArrayLength, probs = 0.25, na.rm = TRUE),
      ArrayLength_75_perc = quantile(ArrayLength, probs = 0.75, na.rm = TRUE),
      N_pixels = n(),
      .groups = "drop"
    )

  return(summarized_df)
}


compute_tic_trend <- function(
  summarized_scanned_lines_table,
  trend_window_fraction = 0.05
) {
  # Validation
  required_cols <- c("Dataset_ID", "Region", "Y", "Median_TIC")
  missing_cols <- setdiff(required_cols, names(summarized_scanned_lines_table))
  if (length(missing_cols) > 0L) {
    stop(
      "summarized_scanned_lines_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(summarized_scanned_lines_table) == 0L) {
    stop("summarized_scanned_lines_table has no rows.", call. = FALSE)
  }

  # Ratio of the last fraction of scan lines to the first, rather than a fitted
  # slope: the trajectories rise steeply then plateau, so a linear fit would
  # average two areas and describe neither
  tic_table_grouped <- summarized_scanned_lines_table |>
    arrange(Dataset_ID, Region, Y) |>
    group_by(Dataset_ID, Region) |>
    summarize(
      N_lines = n(),
      Window_size = round(n() * trend_window_fraction),
      Median_first_window = median(head(Median_TIC, Window_size)),
      Median_last_window = median(tail(Median_TIC, Window_size)),
      .groups = "drop"
    ) |>
    mutate(
      TIC_last_over_first = Median_last_window / Median_first_window
    )

  return(tic_table_grouped)
}

### Interpretation of tic_table_grouped:
# The ratio is sensitive to the window fraction, since the opening transient and the plateau
# are two distinct regimes and any fixed window averages across both. At a 10% window b84bb453
# reports 1.05; at 5% it reports 1.09; the figure shows the full transient to be around 22%.
# The metric is used here to rank regions and flag that a trajectory is not flat, not to quantify
# the transient, which would need a shape-based measure such as the number of lines to reach the plateau.

compute_arrayLength_trend <- function(
  summarized_scanned_lines_table,
  trend_window_fraction = 0.05
) {
  # Validation checks
  required_cols <- c("Dataset_ID", "Region", "Y", "Median_ArrayLength")
  missing_cols <- setdiff(required_cols, names(summarized_scanned_lines_table))
  if (length(missing_cols) > 0L) {
    stop(
      "summarized_scanned_lines_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(summarized_scanned_lines_table) == 0L) {
    stop("summarized_scanned_lines_table has no rows.", call. = FALSE)
  }

  # Ratio of the last fraction of scan lines to the first, rather than a fitted
  # slope: the trajectories rise steeply then plateau, so a linear fit would
  # average two areas and describe neither
  pixel_peaks_table_grouped <- summarized_scanned_lines_table |>
    arrange(Dataset_ID, Region, Y) |>
    group_by(Dataset_ID, Region) |>
    summarize(
      N_lines = n(),
      Window_size = round(n() * trend_window_fraction),
      Median_first_window = median(head(Median_ArrayLength, Window_size)),
      Median_last_window = median(tail(Median_ArrayLength, Window_size)),
      .groups = "drop"
    ) |>
    mutate(
      ArrayLength_last_over_first = Median_last_window / Median_first_window
    )

  return(pixel_peaks_table_grouped)
}


# Read the registered organ code from a dataset's metadata.json. HuBMAP
# registers this on sample records in the provenance chain, so a dataset may
# have several. The field's position varies between files, so walk the whole
# structure recursively and record the path to each value
get_registered_organ <- function(json_path, dataset_ID) {
  if (
    !is.character(dataset_ID) || length(dataset_ID) != 1L || !nzchar(dataset_ID)
  ) {
    stop("'dataset_ID' must be a single non-empty string.", call. = FALSE)
  }

  parsed_json <- parse_json(json_path)

  # Fields that state a laterality for this dataset. "organ" is the registered
  # value on sample records; "intended_organ" is stated on the upload record
  # and can differ from what was ultimately registered
  organ_field_names <- c("organ", "intended_organ")

  # Walk the nested list, collecting every element named "organ" along with the path that reaches it
  find_organ_fields <- function(node, path = character(0)) {
    if (!is.list(node)) {
      return(NULL)
    }

    found <- list()
    node_names <- names(node)

    for (i in seq_along(node)) {
      child <- node[[i]]
      child_name <- if (is.null(node_names) || !nzchar(node_names[i])) {
        paste0("[[", i, "]]")
      } else {
        node_names[i]
      }
      child_path <- c(path, child_name)

      if (
        child_name %in%
          organ_field_names &&
          is.character(child) &&
          length(child) == 1L &&
          nzchar(child)
      ) {
        found[[length(found) + 1L]] <- data.frame(
          Path = paste(child_path, collapse = "$"),
          Organ = child,
          stringsAsFactors = FALSE
        )
      }

      found <- c(found, find_organ_fields(child, child_path))
    }

    return(found)
  }

  organ_hits <- bind_rows(find_organ_fields(parsed_json))

  if (nrow(organ_hits) == 0L) {
    warning(
      "No 'organ' field found in ",
      basename(json_path),
      ".",
      call. = FALSE
    )
    return(data.frame(
      Dataset_ID = dataset_ID,
      Registered_organ = NA_character_,
      N_organ_fields = 0L,
      Organ_values = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  distinct_organs <- unique(organ_hits$Organ)

  # Several records naming the same organ is normal in a provenance chain.
  # Several naming *different* organs means the file does not state a single
  # answer, so return NA rather than picking one
  if (length(distinct_organs) > 1L) {
    warning(
      "Conflicting 'organ' values in ",
      basename(json_path),
      ": ",
      paste(
        paste0(organ_hits$Organ, " at ", organ_hits$Path),
        collapse = "; "
      ),
      ". Returning NA rather than choosing one.",
      call. = FALSE
    )
    registered_organ <- NA_character_
  } else {
    registered_organ <- distinct_organs[1]
  }

  return(data.frame(
    Dataset_ID = dataset_ID,
    Registered_organ = registered_organ,
    N_organ_fields = nrow(organ_hits),
    Organ_values = paste(
      paste0(organ_hits$Organ, " (", organ_hits$Path, ")"),
      collapse = "; "
    ),
    stringsAsFactors = FALSE
  ))
}

# Compare the file's internal run identifier against the organ registered for its sample
check_identifier_consistency <- function(run_table, organ_table) {
  # Validation checks
  required_run <- c("Dataset_ID", "Run_ID")
  missing_cols <- setdiff(required_run, names(run_table))
  if (length(missing_cols) > 0L) {
    stop(
      "run_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  required_organ <- c("Dataset_ID", "Registered_organ")
  missing_cols <- setdiff(required_organ, names(organ_table))
  if (length(missing_cols) > 0L) {
    stop(
      "organ_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  joined <- left_join(
    select(run_table, Dataset_ID, Run_ID),
    organ_table,
    by = "Dataset_ID"
  )

  # Match LK/RK only as delimited tokens. HuBMAP IDs contain random
  # alphanumeric components, so HBM449.RKBB.627 contains "RK" without meaning
  # right kidney; a bare substring scan produces false positives.
  extract_laterality <- function(run_id) {
    if (is.na(run_id)) {
      return(NA_character_)
    }
    # Look for LK or RK matches where no letters or nums precede or come after it
    hits <- regmatches(
      run_id,
      gregexpr("(?<![A-Za-z0-9])[LR]K(?![A-Za-z0-9])", run_id, perl = TRUE)
    )[[1]]
    hits <- unique(hits)
    if (length(hits) == 0L) {
      return(NA_character_)
    }
    if (length(hits) > 1L) {
      return("AMBIGUOUS")
    }
    return(hits)
  }

  # Create a dataframe to summarize the laterality and verdict (organ not registered, more than one hit, no match, match, or mismatch)
  identifier_checks <- joined |>
    mutate(
      Run_id_laterality = vapply(Run_ID, extract_laterality, character(1)),
      Verdict = case_when(
        is.na(Registered_organ) &
          N_organ_fields > 0 ~ "conflicting registered organ",
        is.na(Registered_organ) ~ "no registered organ",
        Run_id_laterality == "AMBIGUOUS" ~ "ambiguous run id",
        is.na(Run_id_laterality) ~ "not determinable",
        Run_id_laterality == Registered_organ ~ "match",
        TRUE ~ "mismatch"
      )
    ) |>
    select(
      Dataset_ID,
      Registered_organ,
      Run_id_laterality,
      Organ_values,
      Verdict
    )

  # Identify and alert of mismatches
  mismatches <- identifier_checks |> filter(Verdict == "mismatch")
  if (nrow(mismatches) > 0L) {
    for (i in seq_len(nrow(mismatches))) {
      row <- mismatches[i, ]
      warning(
        "Run id in ",
        row$Dataset_ID,
        " implies '",
        row$Run_id_laterality,
        "' but the registered organ is '",
        row$Registered_organ,
        "'.",
        call. = FALSE
      )
    }
  }

  return(identifier_checks)
}


# Positions inside each region's bounding box that have no spectrum. Built by
# constructing every position the box contains and subtracting the measured
# ones, since at ~99.8% occupancy the gaps are invisible on an occupancy map.
# Known limitation: materialises ~805k rows before the anti-join discards
# nearly all of them, which would not scale to many large regions
find_missing_positions <- function(spectrum_table, region_table) {
  required_spectrum <- c("Dataset_ID", "Region", "X", "Y")
  missing_cols <- setdiff(required_spectrum, names(spectrum_table))
  if (length(missing_cols) > 0L) {
    stop(
      "spectrum_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  required_region <- c(
    "Dataset_ID",
    "Region",
    "X_min",
    "X_max",
    "Y_min",
    "Y_max"
  )
  missing_cols <- setdiff(required_region, names(region_table))
  if (length(missing_cols) > 0L) {
    stop(
      "region_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # One region at a time, since each has its own extent
  full_grid <- bind_rows(lapply(seq_len(nrow(region_table)), function(i) {
    region_row <- region_table[i, ]
    expand.grid(
      X = region_row$X_min:region_row$X_max,
      Y = region_row$Y_min:region_row$Y_max
    ) |>
      mutate(
        Dataset_ID = region_row$Dataset_ID,
        Region = region_row$Region
      )
  }))

  missing_positions <- anti_join(
    full_grid,
    spectrum_table,
    by = c("Dataset_ID", "Region", "X", "Y")
  )

  return(missing_positions)
}

# One row per region: how many positions are unmeasured and which columns they
# sit in. The column list is the point: if every affected column is at or near
# a region's edge, the gaps are a raster edge effect rather than dropout.
summarise_missing_positions <- function(spectrum_table, region_table) {
  missing_positions <- find_missing_positions(spectrum_table, region_table)

  if (nrow(missing_positions) == 0L) {
    message("No unmeasured positions in any region.")
    return(NULL)
  }

  missing_summary <- missing_positions |>
    group_by(Dataset_ID, Region) |>
    summarise(
      N_missing = n(),
      # n_distinct() counts unique values; sort() so the list reads in order
      N_columns_affected = n_distinct(X),
      Columns_affected = paste(sort(unique(X)), collapse = ", "),
      .groups = "drop"
    ) |>
    # Scan line counts and extents come from the region table, not the missing
    # positions: a region could in principle have a line with no gaps at all.
    left_join(
      region_table |>
        select(Dataset_ID, Region, X_min, X_max, Y_min, Y_max) |>
        mutate(N_lines = Y_max - Y_min + 1L),
      by = c("Dataset_ID", "Region")
    ) |>
    select(
      Dataset_ID,
      Region,
      N_lines,
      N_missing,
      X_min,
      X_max,
      N_columns_affected,
      Columns_affected
    )

  return(missing_summary)
}
