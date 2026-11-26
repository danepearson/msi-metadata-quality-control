library(jsonlite)

EXTRA <- c(
  TIC = "MS:1000285",
  base_peak_mz = "MS:1000504",
  base_peak_int = "MS:1000505",
  lowest_mz = "MS:1000528",
  highest_mz = "MS:1000527"
)

# Function to extract metadata from the imzML file and return it as parsed data
extract_imzML_data <- function(file) {
  # Fail before parsing rather than after, since parsing a large header is slow
  if (!is.character(file) || length(file) != 1L) {
    stop("'file' must be a single file path.", call. = FALSE)
  }
  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }
  # grepl() with ignore.case catches both .imzML and .imzml
  if (!grepl("\\.imzML$", file, ignore.case = TRUE)) {
    stop("Expected a '.imzML' file, got: ", basename(file), call. = FALSE)
  }

  parsed_data <- parseImzML(
    file = file,
    ibd = FALSE,
    extra = EXTRA
  )
  return(parsed_data)
}


# Pull one field out of the parsed object, check it, and convert it to numeric.
# "container" = The list holding the field (like 'extra' or 'positions' element)
# "field" = The name within "container" (like 'external array length')
# "accession" is only used to make the error message informative.
extract_field <- function(container, field, expected_length, accession = NA) {
  metric_values <- container[[field]]

  # A missing cvParam returns NULL, which as.numeric() would turn into numeric(0)
  if (is.null(metric_values)) {
    stop(
      "Field '",
      field,
      "'",
      if (!is.na(accession)) paste0(" (", accession, ")") else "",
      " is missing from this file.",
      call. = FALSE
    )
  }

  # Verify number of values for each metric vector matches expected number of them
  if (length(metric_values) != expected_length) {
    stop(
      "Field '",
      field,
      "' has ",
      length(metric_values),
      " values but ",
      expected_length,
      " spectra were found.",
      call. = FALSE
    )
  }

  # Convert the metrics to numeric
  numeric_values <- suppressWarnings(as.numeric(metric_values))

  # Track how many values could not be converted to numeric
  introduced_nas <- sum(is.na(numeric_values))
  if (introduced_nas > 0) {
    warning(
      introduced_nas,
      " value(s) in '",
      field,
      "' could not be converted to numeric.",
      call. = FALSE
    )
  }

  return(numeric_values)
}


# Function to create a one-row-per-spectrum summary table from parsed imzML file
create_spectrum_summary <- function(parsed_data, dataset_ID) {
  # Safety checks:
  if (
    !is.character(dataset_ID) || length(dataset_ID) != 1L || !nzchar(dataset_ID)
  ) {
    stop("'dataset_ID' must be a single non-empty string.", call. = FALSE)
  }

  spectrum_list <- parsed_data$run$spectrumList
  if (is.null(spectrum_list)) {
    stop(
      "No spectrumList found. Make sure you parsed with parseImzML()",
      call. = FALSE
    )
  }

  positions <- spectrum_list$positions
  if (is.null(positions) || nrow(positions) == 0L) {
    stop("No spectrum positions found.", call. = FALSE)
  }

  # Define counts to check other fields against
  n_spectra <- nrow(positions)
  index <- seq_len(n_spectra)

  # Define and validate built-in metrics using the function made earlier
  x <- extract_field(positions, "position x", n_spectra)
  y <- extract_field(positions, "position y", n_spectra)
  array_length <- extract_field(
    spectrum_list$mzArrays,
    "external array length",
    n_spectra
  )

  # Define custom metrics, named in the EXTRA constant at the top
  extras <- spectrum_list$extra
  if (is.null(extras)) {
    stop(
      "No 'extra' fields found. make sure the file was parsed with the EXTRA vector.",
      call. = FALSE
    )
  }
  tic <- extract_field(extras, "TIC", n_spectra, EXTRA[["TIC"]])
  base_peak_mz <- extract_field(
    extras,
    "base_peak_mz",
    n_spectra,
    EXTRA[["base_peak_mz"]]
  )
  base_peak_intensity <- extract_field(
    extras,
    "base_peak_int",
    n_spectra,
    EXTRA[["base_peak_int"]]
  )
  lowest_observed_mz <- extract_field(
    extras,
    "lowest_mz",
    n_spectra,
    EXTRA[["lowest_mz"]]
  )
  highest_observed_mz <- extract_field(
    extras,
    "highest_mz",
    n_spectra,
    EXTRA[["highest_mz"]]
  )

  # Create dataframe with all the validated metrics
  spectrum_summary_df <- data.frame(
    Index = index,
    X = x,
    Y = y,
    TIC = tic,
    BasePeak_mz = base_peak_mz,
    BasePeak_intensity = base_peak_intensity,
    Lowest_mz = lowest_observed_mz,
    Highest_mz = highest_observed_mz,
    ArrayLength = array_length,
    Dataset_ID = dataset_ID
  )

  # Check for missing coordinates before it gets to a plot or summary table
  missing_coords <- sum(
    is.na(spectrum_summary_df$X) | is.na(spectrum_summary_df$Y)
  )
  if (missing_coords > 0) {
    warning(
      missing_coords,
      " spectra have missing pixel coordinates.",
      call. = FALSE
    )
  }

  # Non-integer coordinates would mean the position params hold something other than a pixel index
  non_integer_coords <- sum(
    spectrum_summary_df$X %% 1 != 0 | spectrum_summary_df$Y %% 1 != 0,
    na.rm = TRUE
  )
  if (non_integer_coords > 0) {
    warning(
      non_integer_coords,
      " spectra have non-integer pixel coordinates.",
      call. = FALSE
    )
  }

  # Highest should exceed lowest in every spectrum; if not, the fields are probaby swapped or one is unpopulated
  low_exceeding_high <- sum(
    spectrum_summary_df$Highest_mz < spectrum_summary_df$Lowest_mz,
    na.rm = TRUE
  )
  if (low_exceeding_high > 0) {
    warning(
      low_exceeding_high,
      " spectra have highest m/z below lowest m/z.",
      call. = FALSE
    )
  }

  return(spectrum_summary_df)
}


# CardinalIO does not expose XML attributes on the <run> and <spectrumList>
# elements (it surfaces only the parsed child content), so read them directly
# from the file. Both tags sit after the preamble but before the spectrum body,
# so a bounded read from the top of the file is enough.
extract_header_attrs <- function(file, max_bytes = 5e6) {
  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }

  # Read as raw bytes rather than lines: imzML has no guaranteed line breaks,
  # and readLines() on a file with very long lines is slow and memory-hungry.
  con <- file(file, open = "rb")
  on.exit(close(con), add = TRUE)
  bytes <- readBin(con, what = "raw", n = max_bytes)

  header <- rawToChar(bytes)
  # Mark as bytes to avoid encoding conversion on what may be an invalid UTF-8
  # sequence at the truncation point.
  Encoding(header) <- "bytes"

  # Pull one attribute from one element. [^>]* skips other attributes that may
  # appear before the one we want inside the same tag.
  get_attr <- function(element, attribute) {
    pattern <- paste0("<", element, "[^>]*\\s", attribute, '="[^"]*"')
    match <- regmatches(header, regexpr(pattern, header, useBytes = TRUE))
    if (length(match) == 0L) {
      return(NA_character_)
    }
    sub(
      paste0('.*\\s', attribute, '="([^"]*)".*'),
      "\\1",
      match[1],
      useBytes = TRUE
    )
  }

  run_id <- get_attr("run", "id")
  data_processing_ref <- get_attr("spectrumList", "defaultDataProcessingRef")

  if (is.na(run_id)) {
    warning(
      "No <run> id found in the first ",
      format(max_bytes, big.mark = ","),
      " bytes of ",
      basename(file),
      ". Try increasing 'max_bytes'.",
      call. = FALSE
    )
  }

  return(list(run_id = run_id, data_processing_ref = data_processing_ref))
}

# Pull a single cvParam value from parsed container using its accession
# Run-level equivalent of extract_field() to catch errors too
extract_cv_value <- function(container, accession, required = TRUE) {
  param <- container[[accession]]

  if (is.null(param)) {
    if (required) {
      stop(
        "cvParam '",
        accession,
        "' is missing from this file.",
        call. = FALSE
      )
    }
    return(NA_character_)
  }

  value <- param[["value"]]
  if (is.null(value)) {
    return(NA_character_)
  }

  return(unname(value))
}

# Pull the declared raster pattern from a parsed imzML file. Shared by
# detect_regions(), which needs it to check its top-down assumption, and
# create_run_summary(), which records it
extract_scan_pattern <- function(parsed_data) {
  scan_settings_list <- parsed_data$scanSettingsList
  if (is.null(scan_settings_list) || length(scan_settings_list) == 0L) {
    return(NA_character_)
  }
  scan_settings <- scan_settings_list[[1]]

  # Tags for top down, one way, flyback, horizontal line scan, line scan left right
  scan_pattern_tags <- c(
    "IMS:1000401",
    "IMS:1000411",
    "IMS:1000413",
    "IMS:1000480",
    "IMS:1000491"
  )

  present_tags <- scan_settings[names(scan_settings) %in% scan_pattern_tags]
  if (length(present_tags) == 0L) {
    return(NA_character_)
  }

  # These params carry their meaning in the name and have no value, so take the name field
  paste(
    vapply(
      present_tags,
      function(param) unname(param[["name"]]),
      character(1)
    ),
    collapse = "; "
  )
}


# Function to create a one-row-per-run summary from a parsed imzML file,
# combining declared acquisition settings with observed file facts
create_run_summary <- function(
  parsed_data,
  dataset_ID,
  run_id = NA_character_,
  data_processing_ref = NA_character_,
  n_regions = NA_integer_,
  grid_matches_region = NA
) {
  # Safety checks
  if (
    !is.character(dataset_ID) || length(dataset_ID) != 1L || !nzchar(dataset_ID)
  ) {
    stop("'dataset_ID' must be a single non-empty string.", call. = FALSE)
  }

  # Scan setting configurations (accession and name pairs, like "IMS:1000401 - top down") (note to self)
  scan_settings_list <- parsed_data$scanSettingsList
  if (is.null(scan_settings_list) || length(scan_settings_list) == 0L) {
    stop("No scanSettingsList found in this file.", call. = FALSE)
  }
  if (length(scan_settings_list) > 1L) {
    warning(
      length(scan_settings_list),
      " scan settings found; using the first only.",
      call. = FALSE
    )
  }
  scan_settings <- scan_settings_list[[1]]

  # Define and ensure existence of x, y, z positions from run
  positions <- parsed_data$run$spectrumList$positions
  if (is.null(positions) || nrow(positions) == 0L) {
    stop("No spectrum positions found.", call. = FALSE)
  }

  # Declared grid dimensions
  max_pixels_x <- as.numeric(extract_cv_value(scan_settings, "IMS:1000042"))
  max_pixels_y <- as.numeric(extract_cv_value(scan_settings, "IMS:1000043"))

  # Detect which scan pattern was used
  scan_pattern <- extract_scan_pattern(parsed_data)

  # Number of pixels (spectrums) recorded during the run
  spectrum_count <- nrow(positions)

  # Binary storage format, read from the first spectrum's m/z array. These are
  # run-level in practice but recorded per array, so take the first and flag
  # if they are not uniform.
  mz_arrays <- parsed_data$run$spectrumList$mzArrays
  binary_type <- NA_character_
  compression <- NA_character_
  if (!is.null(mz_arrays) && nrow(mz_arrays) > 0L) {
    binary_type <- as.character(mz_arrays[["binary data type"]][1])
    compression <- as.character(mz_arrays[["binary data compression type"]][1])

    if (length(unique(mz_arrays[["binary data type"]])) > 1L) {
      warning("Mixed binary data types across spectra.", call. = FALSE)
    }
    # Positions and arrays should describe the same spectra
    if (nrow(mz_arrays) != spectrum_count) {
      warning(
        "mzArrays has ",
        nrow(mz_arrays),
        " rows but ",
        spectrum_count,
        " spectra were found.",
        call. = FALSE
      )
    }
  }

  full_uuid <- extract_cv_value(
    parsed_data$fileDescription$fileContent,
    "IMS:1000080",
    required = FALSE
  )

  declared_grid_size <- max_pixels_x * max_pixels_y

  # One row, so run summaries can be stacked the same way spectrum summaries are
  run_info_df <- data.frame(
    Dataset_ID = dataset_ID,
    Run_ID = run_id,
    Full_UUID = full_uuid,
    Spectrum_count = spectrum_count,
    N_regions = n_regions,
    Max_pixels_x = max_pixels_x,
    Max_pixels_y = max_pixels_y,
    Declared_grid_size = declared_grid_size,
    Grid_matches_region = grid_matches_region,
    Scan_pattern = scan_pattern,
    Binary_data_type = binary_type,
    Compression = compression,
    Data_processing_ref = data_processing_ref,
    stringsAsFactors = FALSE
  )

  # Sanity check: more spectra than declared grid positions means the grid dimensions and the acquisition disagree.
  if (!is.na(max_pixels_x) && !is.na(max_pixels_y)) {
    if (spectrum_count > declared_grid_size) {
      warning(
        "Spectrum count (",
        spectrum_count,
        ") exceeds declared grid size (",
        max_pixels_x * max_pixels_y,
        ").",
        call. = FALSE
      )
    }
  }
  if (!is.na(n_regions) && n_regions > 1L) {
    warning(
      n_regions,
      " acquisition regions detected in ",
      dataset_ID,
      "; spatial metrics must be computed per region.",
      call. = FALSE
    )
  }
  # the declared grid should match the extent of at least one single region, and comparing against the union is useless when a file has more than one
  if (!is.na(grid_matches_region) && !grid_matches_region) {
    warning(
      "Declared grid (",
      max_pixels_x,
      " x ",
      max_pixels_y,
      ") matches no single region's extent in ",
      dataset_ID,
      "; it may describe the union of multiple regions.",
      call. = FALSE
    )
  }

  return(run_info_df)
}


parse_json <- function(json_path) {
  # Validation checks
  if (!is.character(json_path) || length(json_path) != 1L) {
    stop("'json_path' must be a single file path.", call. = FALSE)
  }
  if (!file.exists(json_path)) {
    stop("File not found: ", json_path, call. = FALSE)
  }
  if (tolower(tools::file_ext(json_path)) != "json") {
    stop(
      "Expected a '.json' file, got: ",
      basename(json_path),
      call. = FALSE
    )
  }
  # simplifyVector = FALSE keeps the result as a plain nested list rather than
  # collapsing arrays into vectors and data frames, which is easier to walk
  # when flattening the metadata for comparison
  parsed_json <- fromJSON(json_path, simplifyVector = FALSE)
  return(parsed_json)
}
