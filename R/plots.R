library(ggplot2)
library(grid)

# Figures for the MSI QC pipeline. These functions take tables and return
# ggplot objects; saving is the caller's job, so the report can control size
# and placement without re-running the computation.

# Every position measured in each region, drawn at true scale and position.
# Takes run_table as well so each panel can show its own declared grid.
plot_occupancy_map <- function(spectrum_table, region_table, run_table) {
  required_spectrum <- c("Dataset_ID", "Region", "X", "Y")
  missing_cols <- setdiff(required_spectrum, names(spectrum_table))
  if (length(missing_cols) > 0L) {
    stop(
      "spectrum_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  required_region <- c("Dataset_ID", "Region")
  missing_cols <- setdiff(required_region, names(region_table))
  if (length(missing_cols) > 0L) {
    stop(
      "region_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  required_run <- c("Dataset_ID", "Max_pixels_x", "Max_pixels_y")
  missing_cols <- setdiff(required_run, names(run_table))
  if (length(missing_cols) > 0L) {
    stop(
      "run_table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # One readable facet label. Built at display time rather than stored, so the
  # underlying tables keep Dataset_ID and Region as separate keys.
  spectrum_table <- spectrum_table |>
    mutate(Panel = paste0(Dataset_ID, " region ", Region))

  # The declared grid per region, drawn as a reference rectangle. Facets share
  # one axis range, so empty space does not always mean uncovered grid:
  # b84bb453's panel is wider than its own grid because b07f13aa's region 1
  # reaches x=1322. The rectangle makes each region's own reference explicit.
  # A multi-region file has one declared grid across all its regions, so the
  # right_join duplicates it per region, which is the point: both regions are
  # measured against a grid that describes neither.
  declared_boxes <- run_table |>
    select(Dataset_ID, Max_pixels_x, Max_pixels_y) |>
    right_join(
      select(region_table, Dataset_ID, Region),
      by = "Dataset_ID"
    ) |>
    mutate(
      Panel = paste0(Dataset_ID, " region ", Region),
      xmin = 0,
      xmax = Max_pixels_x,
      ymin = 0,
      ymax = Max_pixels_y
    )

  plot <- ggplot(spectrum_table, aes(x = X, y = Y)) +
    # inherit.aes = FALSE because declared_boxes has no X/Y columns
    geom_rect(
      data = declared_boxes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = NA,
      colour = "grey60",
      linetype = "dashed",
      linewidth = 0.4
    ) +
    # geom_tile() rather than geom_raster(): coverage is not perfectly regular,
    # and geom_raster() shifts cells to force a uniform grid.
    geom_tile() +
    # Fixed scales here, unlike the missing-positions plot: the point is to
    # show the three regions at true relative scale and position.
    facet_wrap(~Panel) +
    coord_fixed() +
    # imzML raster is top-down, so Y=1 is the top row
    scale_y_reverse() +
    labs(
      title = "Occupancy map",
      subtitle = "Dashed outline is each region's declared grid",
      x = "X coordinate",
      y = "Y coordinate"
    ) +
    theme_minimal() +
    theme(panel.spacing = unit(1, "lines"))

  return(plot)
}

# Positions inside a region's bounding box that were never measured. At ~99.8%
# occupancy the gaps are invisible on the occupancy map, so they are found by
# constructing the full grid and subtracting what exists.
plot_missing_positions <- function(spectrum_table, region_table) {
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

  # Build every position each region's bounding box contains, one region at a
  # time since each has its own extent. Known limitation: this materialises
  # ~805k rows before the anti-join discards nearly all of them, which would
  # not scale to many large regions.
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

  # Keep grid positions that have no matching spectrum
  missing_positions <- anti_join(
    full_grid,
    spectrum_table,
    by = c("Dataset_ID", "Region", "X", "Y")
  )

  message(
    nrow(missing_positions),
    " unmeasured positions across ",
    nrow(region_table),
    " regions."
  )

  if (nrow(missing_positions) == 0L) {
    warning("No unmeasured positions found; nothing to plot.", call. = FALSE)
    return(NULL)
  }

  missing_positions <- missing_positions |>
    mutate(Panel = paste0(Dataset_ID, " region ", Region))

  # Region 2's gaps all sit at a single x value, so free scales would collapse
  # its panel to a fraction of a coordinate unit with meaningless decimal axis
  # labels. Seeding each panel with its region's bounding box forces the axis
  # to span the real extent without drawing anything.
  panel_extents <- region_table |>
    mutate(Panel = paste0(Dataset_ID, " region ", Region)) |>
    select(Panel, X_min, X_max, Y_min, Y_max) |>
    tidyr::pivot_longer(cols = c(X_min, X_max), values_to = "X") |>
    mutate(Y = Y_min)

  # Free x here, unlike the occupancy map: the gaps are a few hundred points
  # on an edge, and a shared 0-1322 axis leaves two panels nearly empty.
  # coord_fixed() is omitted because it cannot coexist with free scales.
  plot <- ggplot(missing_positions, aes(x = X, y = Y)) +
    geom_blank(data = panel_extents) +
    geom_point(size = 0.3, alpha = 0.6) +
    facet_wrap(~Panel, scales = "free_x") +
    scale_y_reverse() +
    labs(
      title = "Unmeasured positions within each region's bounding box",
      subtitle = "X axis spans each region's own extent",
      x = "X coordinate",
      y = "Y coordinate"
    ) +
    theme_minimal() +
    theme(
      panel.spacing = unit(1, "lines"),
      panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.4)
    )

  return(plot)
}


plot_TIC_drift <- function(summarized_spectrum_table) {
  # Validation
  required_cols <- c(
    "Median_TIC",
    "Region",
    "Y",
    "TIC_25_perc",
    "TIC_75_perc"
  )
  missing_cols <- setdiff(required_cols, names(summarized_spectrum_table))
  if (length(missing_cols) > 0L) {
    stop(
      "Summarized spectrum table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(summarized_spectrum_table) == 0L) {
    stop("Summarized spectrum table has no rows.", call. = FALSE)
  }

  # One panel per region
  # Initially I was faceting on Region alone but realized it merges datasets, since region numbers repeat:
  # b07f13aa region 1 and b84bb453 region 1 were overlayed in the first panel
  summarized_spectrum_table <- summarized_spectrum_table |>
    mutate(Panel = paste0(Dataset_ID, " region ", Region))

  plot <- ggplot(summarized_spectrum_table, aes(x = Y, y = Median_TIC)) +
    geom_ribbon(
      aes(ymin = TIC_25_perc, ymax = TIC_75_perc),
      alpha = 0.2
    ) +
    geom_line(linewidth = 0.4, color = "#2b2b2b") +
    facet_wrap(~Panel, nrow = 1, scales = "free") +
    labs(
      title = "Median TIC per scan line (faceted by region)",
      subtitle = paste(
        "Ribbon is the interquartile range of pixel TIC within each line.",
        "X axis is scan line order; spatial position and acquisition time",
        "are confounded by the raster."
      ),
      x = "Y coordinate (scan line)",
      y = "Median TIC",
    ) +
    theme_minimal() +
    theme(
      panel.spacing = unit(1, "lines"),
      panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.4)
    )

  return(plot)
}

### Interpretation of Median TIC per scan line drift plot:
# Design: Median TIC within each scan line against scan line
# order, with the interquartile range of pixel TIC as a ribbon. Faceted per
# region with free scales, since the regions differ in both line count
# (193, 215, 1041) and absolute TIC level (~88k vs ~275k).
#
# All three regions show the same shape: median TIC rises steeply over the opening
# lines and then plateaus. The rise is roughly 3-4% in b07f13aa's two regions
# and ~20% in b84bb453 (225k to 275k), and in all three cases it is concentrated in
# the first 10-20% of lines. After the plateau the signal is stable to the end
# of each acquisition, with no little-to-no decline in TIC.
#
# This could be a result of a warm-up or stabilization period at the start of
# each acquisition rather than drift. The same shape appearing in three independent
# acquisitions across two instruments, two donors, and separate sessions is
# the main evidence for an instrument or source effect rather than issues with
# any one sample.
#
# The IQR ribbon widens alongside the median rather than independently, so
# this is a shift in level across the whole pixel distribution rather than a
# change in variability.

plot_arrayLength_drift <- function(summarized_spectrum_table) {
  # Validation checks
  required_cols <- c(
    "Median_ArrayLength",
    "Region",
    "Y",
    "ArrayLength_25_perc",
    "ArrayLength_75_perc"
  )
  missing_cols <- setdiff(required_cols, names(summarized_spectrum_table))
  if (length(missing_cols) > 0L) {
    stop(
      "Summarized spectrum table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(summarized_spectrum_table) == 0L) {
    stop("Summarized spectrum table has no rows.", call. = FALSE)
  }

  # One panel per region
  # Initially I was faceting on Region alone but realized it merges datasets, since region numbers repeat:
  # b07f13aa region 1 and b84bb453 region 1 were overlayed in the first panel
  summarized_spectrum_table <- summarized_spectrum_table |>
    mutate(Panel = paste0(Dataset_ID, " region ", Region))

  plot <- ggplot(
    summarized_spectrum_table,
    aes(x = Y, y = Median_ArrayLength)
  ) +
    geom_ribbon(
      aes(ymin = ArrayLength_25_perc, ymax = ArrayLength_75_perc),
      alpha = 0.2
    ) +
    geom_line(linewidth = 0.4, color = "#2b2b2b") +
    facet_wrap(~Panel, nrow = 1, scales = "free") +
    labs(
      title = "Median peaks (array length) per scan line (faceted by region)",
      subtitle = paste(
        "Ribbon is the interquartile range of pixel array length within each line.",
        "X axis is scan line order; spatial position and acquisition time",
        "are confounded by the raster."
      ),
      x = "Y coordinate (scan line)",
      y = "Median array length",
    ) +
    theme_minimal() +
    theme(
      panel.spacing = unit(1, "lines"),
      panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.4)
    )

  return(plot)
}


plot_TIC_spatial <- function(spectrum_table) {
  # Validation
  required_cols <- c(
    "X",
    "Y",
    "TIC",
    "Dataset_ID",
    "Region"
  )
  missing_cols <- setdiff(required_cols, names(spectrum_table))
  if (length(missing_cols) > 0L) {
    stop(
      "Summarized spectrum table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(spectrum_table) == 0L) {
    stop("Summarized spectrum table has no rows.", call. = FALSE)
  }

  spectrum_table <- spectrum_table |>
    mutate(Panel = paste0(Dataset_ID, " region ", Region)) |>
    group_by(Dataset_ID, Region) |>
    # The regions differ a lot in TIC, so raw values limit b07f13aa at the
    # bottom of the shared colour range.
    # Dividing by each region's own median puts all panels on the same scale
    mutate(TIC_relative = TIC / median(TIC, na.rm = TRUE)) |>
    ungroup()

  plot <- ggplot(spectrum_table, aes(x = X, y = Y, fill = TIC_relative)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0.5, 1.5), oob = scales::squish) +
    facet_wrap(~Panel, nrow = 1) +
    # Keep the x and y axes proportional
    coord_fixed() +
    # imzML raster is top-down, so Y=1 is the top row
    scale_y_reverse() +
    labs(
      title = "Relative TIC per pixel (faceted by region)",
      x = "Pixel X coordinate",
      y = "Pixel Y coordinate",
      fill = "TIC / region median"
    )

  return(plot)
}

### Interpretation of Spatial TIC plot:
# Plot design: TIC per pixel, normalized to each region's own median.
# Normalization is needed because ggplot cannot give facets separate fill
# scales and the regions differ ~300% in absolute TIC; raw values put b07f13aa's
# two regions in the bottom sixth of a shared colour range and they render
# as flat purple. Limits are clipped to 0.5-1.5 with squish, since a few
# extreme pixels otherwise stretch the scale and hide the texture.
#
# This plot is the control for the per-scan-line drift plot. Row-by-row
# rastering confounds spatial position with acquisition order, so a trend in
# TIC against scan line could be either instrument drift or tissue structure.
# Comparing the two figures is the only way to tell them apart from this data.
#
# What it shows: b84bb453 showsisolated single pixels and small clusters with elevated
# relative TIC, plus one larger coherent dark patch near (350, 650) with a brighter rim.
#
# Interpretation is limited by the fact that TIC is just the sum of all intensities
# in a spectrum, so higher TIC does not tell us anything about the specific ions.
# Potential causes of TIC variation includes uneven matrix deposition, debris,
#  defects such as folds or tears, and compositional variation.

plot_arrayLength_spatial <- function(spectrum_table) {
  # Validation
  required_cols <- c(
    "X",
    "Y",
    "ArrayLength",
    "Dataset_ID",
    "Region"
  )
  # Validation checks
  missing_cols <- setdiff(required_cols, names(spectrum_table))
  if (length(missing_cols) > 0L) {
    stop(
      "Summarized spectrum table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(spectrum_table) == 0L) {
    stop("Summarized spectrum table has no rows.", call. = FALSE)
  }

  spectrum_table <- spectrum_table |>
    mutate(Panel = paste0(Dataset_ID, " region ", Region)) |>
    group_by(Dataset_ID, Region) |>
    # The regions might differ a lot in ArrayLength, so raw values could limit a dataset at the
    # bottom of the shared colour range.
    # Dividing by each region's own median puts all panels on the same scale
    mutate(
      ArrayLegth_relative = ArrayLength / median(ArrayLength, na.rm = TRUE)
    ) |>
    ungroup()

  plot <- ggplot(
    spectrum_table,
    aes(x = X, y = Y, fill = ArrayLegth_relative)
  ) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0.5, 1.5), oob = scales::squish) +
    facet_wrap(~Panel, nrow = 1) +
    # Keep the x and y axes proportional
    coord_fixed() +
    # imzML raster is top-down, so Y=1 is the top row
    scale_y_reverse() +
    labs(
      title = "Relative Peaks per pixel (array length, faceted by region)",
      x = "Pixel X coordinate",
      y = "Pixel Y coordinate",
      fill = "Array Length /\nregion median"
    )

  return(plot)
}
