# Shared helpers for mapping/*.Rmd
# Suggested order: 01_SpatialDataProcessing.Rmd -> 02_SpeciesRichnessMaps.Rmd ->
#   03_SkateRayRichnessPublication.Rmd (or 03b_SkateRayRichnessPublicationReversed.Rmd);
#   optional 04_RichnessLatitude.Rmd, 04b_RichnessLatitudeReversed.Rmd, 06_SharkRichness.Rmd.

calculate_hex_metrics2 <- function(ranges, hex) {
  # Ensure CRS match
  ranges <- st_transform(ranges, st_crs(hex))

  # Calculate intersection between ranges and hex grid
  range_intersect_hex <- st_intersects(ranges, hex)

  # Calculate range size per species (number of hex cells occupied)
  ranges$hex_range_size <- purrr::map_int(range_intersect_hex, length)

  # Setup parallel processing
  future::plan(multisession, workers = future::availableCores() - 1)

  #### Species Richness ####
  hex$species_richness <- furrr::future_map_int(1:nrow(hex), function(i) {
    sum(purrr::map_lgl(range_intersect_hex, ~ i %in% .x))
  })

  #### Threatened Species Richness ####
  thr_indices <- which(ranges$category %in% c("CR", "EN", "VU"))
  ranges_thr <- ranges[thr_indices, ]
  range_intersect_hex_thr <- range_intersect_hex[thr_indices]

  hex$thr_species_richness <- furrr::future_map_int(1:nrow(hex), function(i) {
    sum(purrr::map_lgl(range_intersect_hex_thr, ~ i %in% .x))
  })

  #### Range-Weighted Richness of Threatened Species ####
  weights <- purrr::map_dbl(ranges_thr$hex_range_size, ~ 1 / .x)
  hex$rwr_thr <- furrr::future_map_dbl(1:nrow(hex), function(i) {
    presences <- purrr::map_lgl(range_intersect_hex_thr, ~ i %in% .x)
    sum(as.numeric(presences) * weights, na.rm = TRUE)
  })

  # Reset plan
  future::plan(multisession)

  return(hex)
}

plot_species_richness_map <- function(hex, map_title = "Species Richness", ylim_lat_deg = NULL) {
  # Filter out hex cells with no species presence
  hex_sel <- hex %>%
    filter(species_richness > 0)

  # Load and transform world map to match hex CRS
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  world <- st_transform(world, st_crs(hex))

  # Base map template (y scale from apply_lat_axis_labels(); lon breaks in ° with default_crs below)
  template_map <- ggplot() +
    geom_sf(data = world, fill = "grey80", color = "grey90", size = 0.2) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0), breaks = c(-120, -60, 0, 60, 120)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.05, 0.5),
      legend.title = element_blank(),
      legend.margin = margin(c(6, 6, 6, 6)),
      axis.text = element_text(size = 8),
      plot.margin = margin(l = 0),
      plot.tag = element_text(),
      plot.tag.position = c(0.02, 0.95),
      plot.title.position = "plot", # position title relative to plot
      plot.title = element_text(hjust = 0, size = 12, face = "bold") # left-align
    ) +
    labs(title = map_title)

  # Final species richness map: axis annotation in ° (hex is projected metres)
  crs_hex <- st_crs(hex)
  coord_sf_args <- list(
    crs = crs_hex,
    default_crs = sf::st_crs(4326),
    expand = FALSE
  )
  if (!is.null(ylim_lat_deg)) {
    coord_sf_args$ylim <- ylim_lat_deg
  }
  richness_map <- template_map +
    geom_sf(data = hex_sel, aes(fill = species_richness), color = NA) +
    scale_fill_viridis_c(direction = -1) +
    do.call(ggplot2::coord_sf, coord_sf_args)

  return(richness_map)
}

# coord_sf fixup_graticule_labels() passes graticule N-line *degree* values (length can differ from
# scale breaks). Fixed label vectors always break; waiver() can become a length-5 character vector
# during scale training when breaks are set, which still mismatches. A label function must return
# exactly length(x) strings for whatever x ggplot/sf pass (maps + SpeciesRichnessMaps scatters).
lat_break_labels <- function(x) {
  n <- length(x)
  if (!n) {
    return(character())
  }
  out <- character(n)
  for (i in seq_len(n)) {
    xi <- suppressWarnings(as.numeric(x[[i]]))
    if (length(xi) != 1L || !is.finite(xi)) {
      out[i] <- ""
    } else if (abs(xi + 90) < 0.5) {
      out[i] <- "90°S"
    } else if (abs(xi + 60) < 0.5) {
      out[i] <- "60°S"
    } else if (abs(xi + 30) < 0.5) {
      out[i] <- "30°S"
    } else if (abs(xi) < 0.5) {
      out[i] <- "0°"
    } else if (abs(xi - 30) < 0.5) {
      out[i] <- "30°N"
    } else if (abs(xi - 60) < 0.5) {
      out[i] <- "60°N"
    } else if (abs(xi - 90) < 0.5) {
      out[i] <- "90°N"
    } else {
      out[i] <- as.character(xi)
    }
  }
  out
}

apply_lat_axis_labels <- function(plot_obj) {
  plot_obj +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = seq(-90, 90, by = 30),
      labels = lat_break_labels
    ) +
    theme(
      axis.text.y = element_text(),
      axis.ticks.y = element_line()
    )
}
