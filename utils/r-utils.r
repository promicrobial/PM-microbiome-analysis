#######################################
#                                     #
#         R helper functions          #
#                                     #
#######################################

##############################################################################
# Author: Nathaniel Cole (nc564@cornell.edu)                                 #
# Date:                                                                      #
# License: MIT                                                               #
# Version: 1.0                                                               #
#                                                                            #
# Last updated:                                                              #
##############################################################################

##############################################################################
#  Configuration                                                             #
##############################################################################

suppressWarnings(suppressMessages({
  library("phyloseq")
  library("tidyverse")
  library("ggpubr")
  library("RColorBrewer")
  library("kableExtra")
  library("vegan")
  library("rstatix")
  library("gtsummary")
  library("gt")
  conflicted::conflicts_prefer(
    dplyr::filter,
    dplyr::select,
    dplyr::mutate,
    base::load,
    base::attr,
    stats::update,
    purrr::modify
  )
}))

options(max.print = 50, scipen = 999, max.width = 100, digits = 3, prType='html')
set.seed(1234)

Sys.setlocale("LC_CTYPE", "en_GB.UTF-8")
##############################################################################
#  #start GGplot Settings for theming                                             #
##############################################################################

fontsize <- 14

theme_set(
  theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 70,
        hjust = 1,
        vjust = 1,
        size = fontsize,
        face = 'bold',
        color = 'black'
      ),
      axis.text.y = element_text(
        size = fontsize,
        face = 'bold',
        color = 'black'
      ),

      axis.title.x = element_text(
        size = fontsize,
        face = 'bold',
        color = 'black'
      ),
      axis.title.y = element_text(
        size = fontsize,
        face = 'bold',
        color = 'black'
      ),

      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),

      legend.title = element_blank(),
      legend.position = 'bottom',
      legend.location = "plot",
      legend.text = element_text(
        size = rel(1),
        face = 'bold',
        color = 'black'
      ),
      legend.justification = "center",
      
      strip.text = element_text(
        size = fontsize,
        face = 'bold',
        color = 'black'
      ),
      strip.background.x = element_rect(fill = "white"),
      strip.background.y = element_rect(fill = "lightgrey"),
      
      plot.margin = margin(1, 1, 1, 1, "cm"),
      plot.caption = element_text(size = (fontsize - 3)),

      panel.spacing = unit(3, "mm")

    )
)

#end
##############################################################################
#  #start gtsummary Theme Settings                                                  #
##############################################################################

suppressMessages(set_gtsummary_theme(
  theme_gtsummary_journal("lancet")
))

pvalue_format <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "",
    x < 0.0001 ~ "<0.0001",
    x < 0.001 ~ "<0.001",
    x < 0.01 ~ "<0.01",
    x < 0.05 ~ "<0.05",
    TRUE ~ as.character(round(x, 3))
  )
}

#' Format P-values According to Scientific Convention
#'
#' @description
#' Formats p-values following standard scientific reporting conventions, converting
#' numeric p-values to character strings with appropriate thresholds and decimal places.
#'
#' @param x Numeric vector of p-values to format
#' @param digits Integer. Number of decimal places for p-values above 0.05 (default: J)
#' @param notation Character. Output notation: "inequality" for < symbols or
#'        "scientific" for scientific notation (default: "inequality")
#' @param thresholds Numeric vector. Custom thresholds for p-value formatting
#'        (default: c(0.0001, 0.001, 0.01, 0.05))
#' @param na_string Character. String to use for NA values (default: "")
#'
#' @return Character vector of formatted p-values
#'
#' @examples
#' # Basic usage
#' pvalue_format(c(0.0032, 0.047, 0.00001, 0.5))
#'
#' @export
#BUG fails to format
pvalue_format <- function(
  x,
  digits = 3,
  notation = c("inequality", "scientific"),
  na_string = ""
) {
  # Input validation
  if (!is.numeric(x) && !all(is.na(x))) {
    stop("Input must be numeric or NA")
  }

  if (!is.numeric(digits) || digits < 0) {
    stop("'digits' must be a non-negative integer")
  }

  notation <- match.arg(notation)

  # Create the formatted output using case_when
  formatted <- dplyr::case_when(
    is.na(x) ~ na_string,
    x == 0 ~ "<0.0001",
    x < 0.0001 ~ "<0.0001",
    x < 0.001 ~ "<0.001",
    x < 0.01 ~ "<0.01",
    x < 0.05 ~ "<0.05",
    x == 1 ~ ">0.999",
    TRUE ~ as.character(round(x, digits))
  )

  if (notation == "scientific") {
    # Convert inequality notation to scientific notation
    formatted <- gsub("<0\\.", "1e-", formatted)
  }

  # Add class for potential method dispatch
  class(formatted) <- c("formatted_pvalue", "character")

  return(formatted)
}

#' Print method for formatted p-values
#' @export
print.formatted_pvalue <- function(x, ...) {
  print(unclass(x), ...)
}

#' Knit print method for formatted p-values
#' @export
knit_print.formatted_pvalue <- function(x, ...) {
  knitr::knit_print(unclass(x), ...)
}

#' Create Formatted Summary Statistics Table with GT
#'
#' @description
#' Creates a formatted summary statistics table using gtsummary and gt packages,
#' with customized formatting for continuous and categorical variables, p-values,
#' and q-values. Designed for easy integration with Quarto/RMarkdown documents.
#'
#' @param data A data frame containing the variables to summarize
#' @param by Character or formula. The grouping variable for comparison
#' @param continuous_stat Character. Format for continuous variables statistics
#'        (default: "{median} ({IQR}; {min}-{max})")
#' @param categorical_stat Character. Format for categorical variables statistics
#'        (default: "{n} / {N} ({p}%)")
#' @param digits Integer. Number of decimal places (default: 2)
#' @param missing Character. How to handle missing values: "no", "ifany", or "always"
#'        (default: "no")
#' @param include_q Logical. Whether to include q-values (default: TRUE)
#'
#' @return A gt table object containing formatted summary statistics
#'
#' @details
#' The function applies the following formatting:
#' \itemize{
#'   \item Continuous variables: Median, IQR, min-max
#'   \item Categorical variables: Counts and percentages
#'   \item P-values: Formatted with thresholds (<0.05, <0.01, etc.)
#'   \item Q-values: Included if specified
#'   \item Bold significant results
#' }
#'
#' @examples
#' # Basic usage
#' gt_summary(mtcars, by = "am")
#'
#' # With custom statistics format
#' gt_summary(
#'   mtcars,
#'   by = "am",
#'   continuous_stat = "{mean} ± {sd}",
#'   categorical_stat = "{n} ({p}%)"
#' )
#'
#' # Without q-values
#' gt_summary(mtcars, by = "am", include_q = FALSE)
#'
#' @section Quarto/RMarkdown Usage:
#' ```{r}
#' #| label: summary-table
#' #| tbl-cap: "Summary Statistics by Group"
#'
#' gt_summary(your_data, by = group)
#' ```
#'
#' @importFrom gt gt fmt_markdown fmt_number
#' @importFrom gtsummary tbl_summary add_p add_q bold_p
#' @export
gt_summary <- function(
  data,
  by,
  continuous_stat = "{median}  \n({IQR}; {min}-{max})",
  categorical_stat = "{n} / {N} ({p}%)",
  digits = 2,
  missing = "no",
  include_q = TRUE
) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (missing(by)) {
    stop("'by' argument must be specified")
  }

  if (!missing %in% c("no", "ifany", "always")) {
    stop("'missing' must be one of: 'no', 'ifany', 'always'")
  }

  # Check for required packages
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' needed. Please install it.")
  } else {
    require(gt)
  }

  if (!requireNamespace("gtsummary", quietly = TRUE)) {
    stop("Package 'gtsummary' needed. Please install it.")
  } else {
    require(gtsummary)
  }

  # Create summary table
  summary_tbl <- tbl_summary(
    data = data,
    by = by,
    missing = missing,
    statistic = list(
      all_continuous() ~ continuous_stat,
      all_categorical() ~ categorical_stat
    )
  ) %>%
    add_p(pvalue_fun = pvalue_format)

  # Add q-values if requested
  if (include_q) {
    summary_tbl <- summary_tbl %>%
      add_q(pvalue_fun = pvalue_format) %>%
      bold_p(q = TRUE)
  } else {
    summary_tbl <- summary_tbl %>%
      bold_p()
  }

  # Convert to gt and apply formatting
  formatted_tbl <- summary_tbl %>%
    as_gt() %>%
    fmt_markdown(columns = where(~ is.character(.x))) %>%
    fmt_number(
      columns = where(~ is.numeric(.x)),
      drop_trailing_zeros = TRUE,
      use_seps = FALSE,
      decimals = digits
    ) %>%
    fmt_number(
      columns = where(~ is.numeric(.x) < 0.0001),
      pattern = "<0.0001"
    ) %>%
    fmt_number(
      columns = where(~ is.numeric(.x) < 0.001),
      pattern = "<0.001"
    ) %>%
    fmt_number(
      columns = where(~ is.numeric(.x) < 0.01),
      pattern = "<0.01"
    ) %>%
    fmt_number(
      columns = where(~ is.numeric(.x) < 0.05),
      pattern = "<0.05"
    )

  return(formatted_tbl)
}

#end
##############################################################################
# #start kable Settings                                                             #
##############################################################################

# kable <- function(data, ...){
#   if(!interactive()) {
#   knitr::kable(data) %>%
#   kableExtra::kable_classic(c("striped", "hover", "condensed"), full_width = F, html_font = "Inter")
#   }
#   else {
#      knitr::kable(data)
#   }

# }

kable <- function(data, ...) {
  knitr::kable(data, table.attr = "style = \"color: #C6D0F5;\"")
}

#end

##############################################################################
#  #start Visualisation                                                   #
##############################################################################

#extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

phyla_colours <- c(
  "Actinobacteriota" = '#a6cee3',
  "Bacteroidota" = '#b2df8a',
  "Campylobacterota" = '#33a02c',
  "Bacillota" = '#Fb9a99',
  "Other" = '#E31a1c',
  "Pseudomonadota" = '#Fdbf6f',
  "Unclassified" = '#Ff7f00',
  "Unknown" = '#ffff99',
  "Verrucomicrobiota" = '#Cab2d6'
)

pairedColours <- colorRampPalette(brewer.pal(12, "Paired")) #colour-blind friendly
set2 <- c(
  '#a6cee3',
  '#1f78b4',
  '#b2df8a',
  '#33a02c',
  '#fb9a99',
  '#e31a1c',
  '#fdbf6f',
  '#ff7f00',
  '#cab2d6',
  '#6a3d9a',
  '#ffff99',
  '#b15928'
)

trt_labels <- c('0' = "CPM", '1' = "FeZnPM")
time_labels <- c(
  "1Baseline" = "Baseline",
  "2Midpoint" = "Midpoint",
  "3Endpoint" = "Endpoint"
)

armColours <- c("#FC8D62", "#8DA0CB")
timepoint_colours <- c("#66C2A5", "#E78AC3")

#for use with ggplot tile or ggcorrplot
ggheatmap_palette <- c('#91bfdb', '#ffffbf', '#fc8d59')

okabe_ito <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999" # grey
)
# Function to assign column classes based on a data dictionary
# To view colour-blind friendly palettes run: display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, colorblindFriendly=TRUE)
assign_unique_colors <- function(
  char_vector,
  color_vector = NULL,
  palette = "Set2"
) {
  if (!require("RColorBrewer", quietly = TRUE)) {
    install.packages("RColorBrewer")
  }
  library(RColorBrewer)

  unique_levels <- unique(char_vector)

  # If color_vector is not provided, create an empty one
  if (is.null(color_vector)) {
    color_vector <- vector("character", length = length(unique_levels))
    #names(color_vector) <- unique_levels
  }

  # Find levels without assigned colors
  unassigned_levels <- setdiff(unique_levels, names(color_vector))

  # Generate unique colors for unassigned levels
  if (length(unassigned_levels) < 9) {
    new_colors <- brewer.pal(n = length(unassigned_levels), name = palette)
  } else {
    palette <- colorRampPalette(c(
      '#a6cee3',
      '#1f78b4',
      '#b2df8a',
      '#33a02c',
      '#fb9a99',
      '#e31a1c',
      '#fdbf6f',
      '#ff7f00',
      '#cab2d6',
      '#6a3d9a',
      '#ffff99',
      '#b15928'
    ))
    new_colors <- palette(n = length(unassigned_levels))
  }
  names(new_colors) <- unassigned_levels

  #  Update the existing color vector with new colors for previously unseen levels
  for (level in unassigned_levels) {
    color_vector[level] <- new_colors[level]
  }

  # remove any NAs
  color_vector <- color_vector[!is.na(names(color_vector))]

  # Return the updated color vector
  return(color_vector[unique(char_vector)])
}

patterns <- function(params, boundary_df, aspect_ratio, legend = FALSE) {
  args <- as.list(params)
  pattern_type <- as.factor(as.numeric(params$pattern_type))

  # Define colorblind-friendly colors
  colors <- c(
    '#a6cee3',
    '#1f78b4',
    '#b2df8a',
    '#33a02c',
    '#fb9a99',
    '#e31a1c',
    '#fdbf6f',
    '#ff7f00',
    '#cab2d6',
    '#6a3d9a',
    '#ffff99',
    '#b15928'
  ) #brewer.pal(12, "Paired")

  # Base parameters for the pattern
  base_args <- list(
    x = boundary_df$x,
    y = boundary_df$y,
    id = boundary_df$id,
    prefix = ""
  )

  # Define different patterns for different levels
  if (pattern_type <= 12) {
    # Solid colors for first 12 levels using stripe pattern with 100% density
    args <- c(
      base_args,
      list(
        pattern = "stripe",
        pattern_fill = colors[pattern_type],
        pattern_colour = colors[pattern_type], # Same as fill color
        pattern_spacing = 0.05,
        pattern_density = 1, # Full density makes it appear solid
        pattern_angle = 0
      )
    )
  } else if (pattern_type >= 13) {
    # Striped patterns for remaining levels
    args <- c(
      base_args,
      list(
        pattern = "stripe",
        pattern_fill = colors[pattern_type],
        pattern_colour = colors[round((pattern_type - 12) / 2)],
        pattern_spacing = 0.05,
        pattern_density = 0.7,
      )
    )
  } else {
    # Dot patterns for remaining levels
    args <- c(
      base_args,
      list(
        pattern = "circle",
        pattern_fill = colors[pattern_type],
        pattern_colour = colors[round((pattern_type - 12) / 2)],
        pattern_spacing = 0.05,
        pattern_density = 0.7,
      )
    )
  }

  do.call(gridpattern::patternGrob, args)
}

# Register the pattern
options(ggpattern_geometry_funcs = list(patterns = patterns))

#' Create Correlation Panel with Color-Coded Values
#'
#' @description
#' Creates a correlation panel for use in pairs plot, displaying correlation
#' coefficients with color-coded backgrounds. The color intensity represents
#' the strength and direction of correlation.
#'
#' @param x Numeric vector or matrix for correlation calculation
#' @param y Numeric vector or matrix for correlation calculation
#' @param digits Integer. Number of decimal places to display (default: 2)
#' @param prefix Character string to prepend to correlation value (default: "")
#' @param text_size Numeric. Size of correlation text (default: 3)
#' @param color_palette Character vector of length 3 specifying low, mid, and high colors
#' @param na.rm Logical. Whether to remove NA values (default: TRUE)
#' @param method Character. Correlation method: "pearson", "spearman", or "kendall" (default: "pearson")
#' @param ... Additional arguments passed to text()
#'
#' @return Invisibly returns the correlation coefficient
#'
#' @examples
#' \dontrun{
#' pairs(mtcars,
#'       upper.panel = panel.cor,
#'       lower.panel = panel.smooth,
#'       gap = 0,
#'       row1attop = TRUE)
#' }
#'
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom stats cor
panel.cor <- function(
  x,
  y,
  digits = 2,
  prefix = "",
  text_size = 3,
  color_palette = c('#91bfdb', '#ffffbf', '#fc8d59'),
  na.rm = TRUE,
  method = "pearson",
  ...
) {
  # Input validation
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors or matrices")
  }
  if (!is.numeric(digits) || digits < 0) {
    stop("digits must be a non-negative integer")
  }
  if (!is.character(prefix)) {
    stop("prefix must be a character string")
  }
  if (!method %in% c("pearson", "spearman", "kendall")) {
    stop("method must be one of: 'pearson', 'spearman', 'kendall'")
  }
  if (length(color_palette) != 3) {
    stop("color_palette must be a vector of three colors")
  }

  # Calculate correlation with error handling
  r <- tryCatch(
    {
      cor(
        x,
        y,
        use = if (na.rm) "complete.obs" else "everything",
        method = method
      )
    },
    error = function(e) {
      warning("Could not compute correlation: ", e$message)
      return(NA)
    }
  )

  # Format text
  txt <- if (!is.na(r)) {
    paste0(prefix, format(round(r, digits), nsmall = digits))
  } else {
    "NA"
  }

  # Create color scale
  col <- colorRampPalette(color_palette)(100)

  # Map correlation to color index with NA handling
  col.ind <- if (!is.na(r)) {
    round((r + 1) * 49 + 1)
  } else {
    50 # middle color for NA
  }

  # Ensure color index is within bounds
  col.ind <- pmax(1, pmin(100, col.ind))

  # Fill panel with color using usr coordinates
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = col[col.ind])

  # Add text with appropriate color
  text_col <- if (is.na(r) || abs(r) <= 0.5) "black" else "white"
  text(mean(usr[1:2]), mean(usr[3:4]), txt, col = text_col, cex = text_size)

  # Return correlation coefficient invisibly
  invisible(r)
}

#end
################################################################################
# Plot limma-voom results
################################################################################
plot_limmavoom = function(tt, physeq, alpha = 0.1) {
  require("ggplot2")
  require("phyloseq")
  # Rownames of tt are the OTU IDs
  xdf = cbind(tt, tax_table(physeq)[rownames(tt), ])
  xdf$OTU <- rownames(xdf)
  # The points to highlight
  specialdf = xdf[(xdf$P.Value < alpha), ]
  # The plot
  p = ggplot(
    data = xdf,
    mapping = aes(x = logFC, y = -log10(P.Value), size = -log10(P.Value))
  ) +
    geom_vline(xintercept = 0.0, linetype = 2) +
    # Background points
    geom_point(color = "black", alpha = 0.65) +
    scale_size_continuous(range = c(1, 4)) +
    guides(size = FALSE, colour = guide_legend(nrow = 2, byrow = TRUE)) +
    # axis labels
    labs(y = expression(-log[10](p)), x = expression(log[2](FC)))
  # The interesting OTUs
  if (nrow(specialdf) > 0) {
    p <- p +
      geom_point(data = specialdf, mapping = aes(color = Genus)) +
      geom_text(
        data = specialdf,
        mapping = aes(label = paste("OTU:", OTU)),
        size = 2,
        nudge_x = 0.3
      )
  }
  return(p)
}

##############################################################################
#  #start Functions                                                                #
##############################################################################

#for scaling numeric variables on a 0 to 1 scale
scale_values <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#backtransforms values transformed by scale_values function
backtransform_scaled_values <- function(scaled_value, original_data) {
  original_range <- max(original_data) - min(original_data)
  backtransformed_value <- scaled_value * original_range + min(original_data)
  return(backtransformed_value)
}

#' Calculate Summary Statistics for Numeric Data
#'
#' @description Calculates summary statistics (mean, standard deviation, standard error,
#' median, maximum, and minimum) for numeric data. Works with vectors, matrices, or
#' data frames, with optional grouping variables.
#'
#' @param data A numeric vector, matrix, or data frame
#' @param varname Character string specifying the variable name to summarize.
#' If NULL and data is a vector, uses the entire vector.
#' @param groupnames Character vector of column names to group by.
#' Optional for data frames, ignored for vectors and matrices.
#'
#' @return A data frame containing summary statistics:
#' \itemize{
#'   \item mean
#'   \item sd (standard deviation)
#'   \item SEM (standard error of the mean)
#'   \item median
#'   \item max (maximum value)
#'   \item min (minimum value)
#' }
#'
#' @examples
#' # Vector example
#' vec <- rnorm(100)
#' data_summary(vec)
#'
#' # Data frame example with grouping
#' df <- data.frame(
#'   value = rnorm(100),
#'   group = rep(c("A", "B"), each = 50)
#' )
#' data_summary(df, "value", "group")
#'
#' # Matrix example
#' mat <- matrix(rnorm(100), ncol = 2)
#' data_summary(mat)
#'
#' @importFrom stats sd median
#' @importFrom plyr ddply rename
#' @export
data_summary <- function(data, varname = NULL, groupnames = NULL) {
  # Input validation and conversion
  if (is.vector(data)) {
    if (!is.numeric(data)) {
      stop("Vector must be numeric")
    }
    data <- data.frame(value = data)
    varname <- "value"
  } else if (is.matrix(data)) {
    if (!is.numeric(data)) {
      stop("Matrix must be numeric")
    }
    data <- as.data.frame(data)
    if (is.null(varname)) {
      varname <- colnames(data)[1]
    }
  } else if (!is.data.frame(data)) {
    stop("Input must be a vector, matrix, or data frame")
  }

  # Check if varname exists in data frame
  if (!is.null(varname) && !(varname %in% colnames(data))) {
    stop("Variable name not found in data")
  }

  # Check if groupnames exist in data frame
  if (!is.null(groupnames) && !all(groupnames %in% colnames(data))) {
    stop("One or more grouping variables not found in data")
  }

  # Summary function
  summary_func <- function(x, col) {
    if (all(is.na(x[[col]]))) {
      return(c(
        n = NA,
        mean = NA,
        sd = NA,
        SEM = NA,
        median = NA,
        max = NA,
        min = NA
      ))
    }
    c(
      n = length(x[[col]]),
      mean = mean(x[[col]], na.rm = TRUE),
      sd = sd(x[[col]], na.rm = TRUE),
      SEM = sd(x[[col]], na.rm = TRUE) /
        sqrt(sum(!is.na(x[[col]]))),
      median = median(x[[col]], na.rm = TRUE),
      max = max(x[[col]], na.rm = TRUE),
      min = min(x[[col]], na.rm = TRUE)
    )
  }

  # Calculate summaries
  if (is.null(groupnames)) {
    # No grouping - calculate overall summary
    data_sum <- as.data.frame(t(summary_func(data, varname)))
  } else {
    # With grouping
    data_sum <- ddply(data, groupnames, .fun = summary_func, varname)
  }

  # Rename mean column if varname provided
  if (!is.null(varname)) {
    data_sum <- rename(data_sum, c("mean" = varname))
  }

  return(data_sum)
}

#' Generate Multiple Random Vectors from Various Distributions
#'
#' @description Creates multiple random vectors using specified distribution functions.
#' Can generate vectors of equal or varying lengths using different random distributions.
#'
#' @param n_vec Numeric. Number of vectors to generate
#' @param n_obs Numeric vector or single value. If single value, all vectors will have
#' the same length. If vector, must have length equal to n_vec
#' @param dist Character or function. Distribution to use. Can be "norm" (default), "unif",
#' "exp", "pois", or a custom random generation function
#' @param params List of parameters for the distribution. For normal distribution:
#' list(mean = 0, sd = 1), for uniform: list(min = 0, max = 1), etc.
#' @param names Character vector of names for the vectors. If NULL, default names
#' will be "vec1", "vec2", etc.
#'
#' @return A list of numeric vectors, each containing random values from the specified distribution
#'
#' @examples
#' # Create 3 normal vectors
#' mkvec(n_vec = 3, n_obs = 10, dist = "norm")
#'
#' # Create 3 uniform vectors
#' mkvec(n_vec = 3, n_obs = 10,
#'                      dist = "unif",
#'                      params = list(min = 0, max = 10))
#'
#' # Create 2 exponential vectors with different rates
#' mkvec(n_vec = 2,
#'                      n_obs = c(10, 20),
#'                      dist = "exp",
#'                      params = list(rate = c(1, 2)))
#'
#' # Using a custom distribution function
#' custom_dist <- function(n) rbeta(n, 2, 2)
#' mkvec(n_vec = 2, n_obs = 10, dist = custom_dist)
#'
#' @export
mkvec <- function(
  n_vec,
  n_obs,
  dist = "norm",
  params = list(mean = 0, sd = 1)
) {
  # Input validation
  if (!is.numeric(n_vec) || n_vec < 1) {
    stop("n_vec must be a positive number")
  }

  # Handle n_obs input - now with recycling
  if (length(n_obs) == 1) {
    n_obs <- rep(n_obs, n_vec)
  } else {
    # Recycle n_obs values if needed
    n_obs <- rep_len(n_obs, n_vec)
  }

  # Generate default names
  names <- paste0("vec", 1:n_vec)

  # Handle custom function with multiple parameters
  if (is.function(dist)) {
    # Create list of vectors
    result <- vector("list", n_vec)

    # Ensure all parameters are properly replicated
    params <- lapply(params, function(p) {
      if (length(p) == 1) {
        rep(p, n_vec)
      } else {
        rep_len(p, n_vec) # Recycle parameter values if needed
      }
    })

    for (i in seq_len(n_vec)) {
      # Extract parameters for this iteration
      current_params <- lapply(params, function(p) p[i])

      # Add n parameter
      current_params$n <- n_obs[i]

      # Call custom function with parameters
      tryCatch(
        {
          result[[i]] <- do.call(dist, current_params)
          # Convert to matrix to vector if needed
          if (is.matrix(result[[i]])) {
            result[[i]] <- as.vector(result[[i]])
          }
        },
        error = function(e) {
          stop("Error in custom function for vector ", i, ": ", e$message)
        }
      )
    }
  } else if (is.character(dist)) {
    # Handle built-in distributions (code for built-in distributions...)
    stop("Built-in distributions not yet implemented")
  } else {
    stop("dist must be either a character string or a function")
  }

  # Assign names to the vectors in the list
  names(result) <- names

  return(result)
}

#end
##############################################################################
#  #start LMMs Tools for LMMs                                                       #
##############################################################################

#generate a dataframe with the common goodness of fit metrics from an nlme derived model
model_stats <- function(nlme_model, vif = TRUE, format = "markdown") {
  if (class(nlme_model) != "lme") {
    message("Model must be of class 'lme'.")
    stop()
  }

  x <- summary(nlme_model)

  mod_df <- data.frame(
    row.names = "Value",
    Bhat = x$coefficients$fixed[1],
    intercept_variance = as.numeric(nlme::VarCorr(x)[, 1])[1], #tau
    residual_variance = as.numeric(nlme::VarCorr(x)[, 1])[2], #sigma
    total_varaince = as.numeric(nlme::VarCorr(x)[, 1])[1] +
      as.numeric(nlme::VarCorr(x)[, 1])[2],
    icc = as.numeric(nlme::VarCorr(x)[, 1])[1] /
      (as.numeric(nlme::VarCorr(x)[, 1])[1] +
        as.numeric(nlme::VarCorr(x)[, 1])[2]), #intercalss correlation coefficient
    BIC = x$BIC,
    AIC = x$AIC,
    loglik = x$logLik
  )

  if (vif == TRUE) {
    VIF <- car::vif(nlme_model) #check for multicollineratiy issues
    VIF_Intercept <- NA
    VIF <- c(Intercept = VIF_Intercept, VIF)

    tTable <- as.data.frame(x$tTable)
    tTable$VIF <- VIF
  } else {
    VIF <- car::vif(nlme_model) #check for multicollineratiy issues
    VIF$Intercept <- NA

    tTable <- as.data.frame(x$tTable)
    cbind(tTable, VIF)
  }

  results <- bind_rows(tTable, as.data.frame(t(mod_df)))
  results <- knitr::kable(results, format = format)
  return(results)
}

#plot residuals to test assumption of normal distribution
# ADD Breusch-Pagan test or similar for heteroskedacity test
# use render = FALSE for interactive use
residual_plots <- function(
  model,
  model_name = names(model),
  render = TRUE,
  chol = TRUE
) {
  set.seed(1234)
  library(ggplot2)
  library(dplyr)
  library(mgcv)

  if (chol == TRUE) {
    rawRes <- residuals(model, type = "pearson")
    estCov <- extract.lme.cov(model, model$data)

    # Perform Cholesky decomposition on the residuals' covariance matrix
    message(
      "Cholesky decomposition is being performed on residuals from ",
      deparse(substitute(model))
    )
    residuals <- solve(t(chol(estCov))) %*% rawRes
  } else {
    residuals <- residuals(model, type = "pearson")
  }

  xy_plot <- ggplot(
    data.frame(x = predict(model), y = residuals),
    aes(x = x, y = y)
  ) +
    geom_point(shape = 1, size = 2, color = "black") +
    #geom_smooth(method = "lm") +
    #ggpubr::stat_cor(position = "dodge") +
    geom_smooth(method = "loess") +
    geom_hline(yintercept = 0, color = "red") +
    xlab("Fitted Values") +
    ylab("Model Residuals") +
    labs(title = paste("Residuals vs Fitted Values")) +
    theme(plot.title = element_text(size = 10))

  qqplot <- ggplot(data.frame(resid = residuals), aes(sample = resid)) +
    stat_qq(shape = 1, size = 2, color = "black") +
    stat_qq_line(color = "red") +
    labs(title = paste("QQ Plot")) +
    ylab(NULL) +
    xlab(NULL) +
    theme(plot.title = element_text(size = 10))

  #create histogram and overlay normal curve
  hist_data <- hist(residuals, plot = FALSE)
  data <- data.frame(x = residuals)

  hist <- ggplot(data, aes(x = x)) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = hist_data$breaks[2] - hist_data$breaks[1],
      fill = "lightblue",
      color = "black",
      boundary = 0.5
    ) +
    stat_function(
      fun = dnorm,
      args = list(mean = mean(data$x), sd = sd(data$x)),
      color = "red",
      size = 1
    ) +
    ylab("Frequency") +
    xlab("Model residuals") +
    labs(title = paste("Residuals Histogram")) +
    theme(plot.title = element_text(size = 10))

  norm_test1 <- capture.output(shapiro.test(residuals))
  norm <- ggplot(data.frame(x = 0:1, y = 0:1), aes(x = x, y = y)) +
    geom_blank() +
    xlab(NULL) +
    ylab(NULL) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) +
    annotate(
      geom = "text",
      x = 0.75,
      y = 0.20,
      label = paste(
        "Summary Output:\n",
        paste(norm_test1, collapse = "\n"),
        "P-value <0.05",
        "\n => NOT normally distributed"
      ),
      hjust = 1,
      vjust = 0,
      size = 4
    )

  # Create a list of plots
  if (render == FALSE) {
    X11(display = "", title = model_name)
  }

  plots <- gridExtra::grid.arrange(
    xy_plot,
    qqplot,
    hist,
    norm,
    ncol = 2,
    nrow = 2
  )

  return(plots)
}

residual_norm_test <- function(model_list) {
  norm_test <- lapply(names(model_list), function(model_name) {
    model <- model_list[[model_name]]
    st <- shapiro.test(residuals(model))
    data.frame(Model = model_name, W = st$statistic, P_value = st$p.value)
  })
  norm_df <- do.call(rbind, norm_test)
  return(norm_df)
}

vda <- function(data, variable, by, ...) {
  # Compute VD.A
  A <- effsize::VD.A(as.formula(glue::glue("`{variable}` ~ {by}")), data = data)
  tibble(vda_estimate = A$estimate, vda_magnitude = A$magnitude)
}

sturges_bins <- function(data) {
  n <- length(data)
  k <- ceiling(log2(n) + 1)
  return(k)
}

#' Test for Normality of Data Frame Columns
#'
#' @description Performs Shapiro-Wilk normality test on a specified numeric column
#' of a data frame and returns formatted results suitable for Quarto/RMarkdown documents.
#'
#' @param df A data frame containing the data to test
#' @param j Column name or index to test for normality
#' @param format Output format for the table. One of "auto", "pipe", "simple",
#' "html", "latex", or "rst" (default: "auto")
#' @param digits Number of decimal places for p-value (default: 4)
#' @param alpha Significance level for normality test (default: 0.05)
#'
#' @return A formatted kable object containing test results. The table includes:
#' \itemize{
#'   \item Name: Name of the tested variable
#'   \item W: Shapiro-Wilk test statistic
#'   \item P_value: P-value of the test
#' }
#'
#' @details
#' The function performs a Shapiro-Wilk test for normality on the specified column.
#' Results are formatted using kable() for clean output in RMarkdown/Quarto documents.
#' The null hypothesis is that the data is normally distributed.
#'
#' @note
#' The Shapiro-Wilk test is most appropriate for sample sizes between 3 and 5000.
#'
#' @examples
#' # Create example data
#' df <- data.frame(
#'   normal = rnorm(100),
#'   uniform = runif(100),
#'   categorical = factor(rep(1:4, 25))
#' )
#'
#' # Test normally distributed data
#' norm_test(df, "normal")
#'
#' # Test non-normal data
#' norm_test(df, "uniform")
#'
#' # Test with different format
#' norm_test(df, "normal", format = "html")
#'
#' @importFrom stats shapiro.test
#' @importFrom knitr kable
#' @export
norm_test <- function(df, j, format, digits = 4, alpha = 0.05) {
  # Input validation
  if (!is.data.frame(df)) {
    stop("First argument must be a data frame")
  }

  if (is.character(j) && !(j %in% names(df))) {
    stop("Column '", j, "' not found in data frame")
  }

  if (is.numeric(j) && (j < 1 || j > ncol(df))) {
    stop("Column index out of bounds")
  }

  # Extract column
  y <- df[[j]]
  col_name <- if (is.character(j)) j else names(df)[j]

  # Check if numeric
  if (!is.numeric(y)) {
    warning("Column '", col_name, "' is not numeric. Skipping test.")
    return(NULL)
  }

  # Check sample size
  if (length(y) < 3 || length(y) > 5000) {
    warning(
      "Sample size (",
      length(y),
      ") is outside recommended range (3-5000) for Shapiro-Wilk test"
    )
  }

  # Remove NA values
  y <- na.omit(y)
  if (length(y) == 0) {
    warning("No non-missing values in column '", col_name, "'")
    return(NULL)
  }

  # Perform test
  norm_test <- shapiro.test(y)

  # Create results data frame
  norm_stats <- data.frame(
    "Variable" = col_name,
    "n" = length(y),
    "W" = norm_test$statistic,
    "P_value" = round(norm_test$p.value, digits),
    "Distribution" = ifelse(
      norm_test$p.value > alpha,
      "Likely Normal",
      "Likely Non-normal"
    ),
    check.names = FALSE
  )

  # Create caption
  caption <- sprintf(
    "Shapiro-Wilk Normality Test Results (α = %.2f)",
    alpha
  )

  # Format table
  table <- knitr::kable(
    norm_stats,
    format = format,
    caption = caption,
    align = c('l', 'r', 'r', 'r', 'l'),
    booktabs = TRUE,
    digits = digits
  )

  # Add styling if format is latex and kableExtra is available
  if (format == "latex") {
    if (requireNamespace("kableExtra", quietly = TRUE)) {
      table <- kableExtra::kable_styling(
        table,
        bootstrap_options = c("striped", "hover"),
        full_width = FALSE
      )
    }
  }

  return(table)
}

#' Print method for norm_test output
#' @param x Output from norm_test function
#' @param ... Additional arguments passed to print
#' @export
print.norm_test <- function(x, ...) {
  print(x, ...)
}

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  round(p.mat, 4)
}

#' Summarize Variables for Linear Mixed Model Terms
#'
#' @description
#' Creates a summary table of variables intended for use in linear mixed models,
#' including basic characteristics and summary statistics for both numeric and
#' categorical variables.
#'
#' @param variables A list of variables to be summarized
#' @param variable_names A character vector of names corresponding to the variables
#'
#' @return A kable-formatted table containing:
#' \itemize{
#'   \item Variable: Name of the variable
#'   \item Type: Class of the variable
#'   \item N: Number of observations
#'   \item Missing: Number of missing values
#'   \item Unique values: Number of unique values
#'   \item Levels: Number of levels (for factors)
#'   \item Summary Stats: For numeric variables: mean, SD, and range.
#'                       For non-numeric variables: "NA"
#' }
#'
#' @details
#' The function processes both numeric and categorical variables differently:
#' * Numeric variables include summary statistics (mean, SD, range)
#' * Categorical variables show number of levels but no summary statistics
#' All results are formatted using knitr::kable for nice printing
#'
#' @examples
#' variables <- list(
#'   age = c(25, 30, 35),
#'   group = factor(c("A", "B", "A")),
#'   time = c(1, 2, 3)
#' )
#' var_names <- c("age", "group", "time")
#' summarise_terms(variables, var_names)
#'
#' @aliases summarize_terms
#' @export
#'
#' @import knitr
#'
#' @seealso
#' \code{\link[knitr]{kable}} for the table formatting
#'
#' @note
#' Ensures consistent formatting across variable types and handles missing values
summarise_terms <- function(variables, variable_names) {
  # Input validation
  if (length(variables) != length(variable_names)) {
    stop("Number of variables must match number of variable names")
  }
  if (!is.list(variables)) {
    stop("variables must be a list")
  }
  if (!is.character(variable_names)) {
    stop("variable_names must be a character vector")
  }

  # Initialize empty data frame
  summaries <- data.frame()

  for (i in seq_along(variables)) {
    var <- variables[[i]]
    name <- variable_names[i]

    # Basic info for all variables
    summary <- data.frame(
      Variable = name,
      Type = class(var)[1],
      N = length(var),
      Missing = sum(is.na(var)),
      `Unique values` = length(unique(var)),
      Levels = length(levels(var))
    )

    # Add type-specific summaries
    if (is.numeric(var)) {
      summary$`Summary Stats` <- paste0(
        "Mean=",
        round(mean(var, na.rm = TRUE), 2),
        ", SD=",
        round(sd(var, na.rm = TRUE), 2),
        ", Range=[",
        round(min(var, na.rm = TRUE), 2),
        "-",
        round(max(var, na.rm = TRUE), 2),
        "]"
      )
    } else {
      summary$`Summary Stats` <- "NA"
    }

    # Bind to main results
    summaries <- rbind(summaries, summary)
  }

  return(knitr::kable(summaries, digits = 2))
}

anova_pvalue <- function(x) {
  UseMethod("anova_pvalue", x)
}
anova_pvalue.default <- function(x) {
  terms <- rownames(x)
  # Extract the p-value column from the data frame
  pvals <- x$`p-value`
  # Combine the p-values and terms into a temporary data frame
  data.frame(Terms = terms, Model = model, P.value = pvals)
}
anova_pvalue.list <- function(x) {
  # Initialize an empty data frame to store results
  combined_results <- data.frame()

  pvalues <- lapply(seq_along(x), function(i) {
    # Get p-values and terms for each model in the current list
    terms <- rownames(x[[i]])
    # Extract the p-value column from the data frame
    pvals <- x[[i]]$`p-value`
    model <- names(x)[i]
    # Combine the p-values and terms into a temporary data frame
    data.frame(Terms = terms, Model = model, P.value = pvals)
  })
  # Bind the temporary data frame to the combined results data frame
  combined_results <- dplyr::bind_rows(combined_results, pvalues)
  rownames(combined_results) <- NULL

  # Return the combined results data frame
  return(combined_results)
}
anova_pvalue.listoflists <- function(x) {
  # Initialize an empty data frame to store results
  combined_results <- data.frame()

  for (l in x) {
    # Get p-values and terms for each model in the current list
    pvalues <- lapply(seq_along(l), function(i) {
      terms <- rownames(l[[i]])
      # Extract the p-value column from the data frame
      pvals <- l[[i]]$`p-value`
      model <- names(l)[i]
      # Combine the p-values and terms into a temporary data frame
      data.frame(Terms = terms, Model = model, P.value = pvals)
    })
    # Bind the temporary data frame to the combined results data frame
    combined_results <- dplyr::bind_rows(combined_results, pvalues)
    rownames(combined_results) <- NULL
  }

  # Return the combined results data frame
  return(combined_results)
}

het_breuschpagan_lme <- function(nlme_model, exog_het, robust = TRUE) {
  # Extract residuals and predicted values from the LME model
  resid <- residuals(nlme_model, type = "pearson")
  y <- resid^2
  if (!robust) {
    y <- y / mean(y)
  }

  # Convert exog_het to a matrix
  exog_het <- as.matrix(exog_het)

  # Fit a linear model using OLS
  lm_fit <- lm(y ~ exog_het - 1) # -1 removes intercept

  # Calculate test statistics
  lm <- if (robust) {
    nobs(lm_fit) * summary(lm_fit)$r.squared
  } else {
    summary(lm_fit)$sigma / 2
  }

  # Calculate p-value
  lm_pval <- 1 - pchisq(lm, length(exog_het) - 1)

  # Return the result as a list
  result <- list(
    lm = lm,
    lm_pvalue = lm_pval,
    fvalue = summary(lm_fit)$fstatistic[1],
    f_pvalue = summary(lm_fit)$fstatistic[2]
  )

  return(result)
}

# Function to check and apply transformation
check_and_transform <- function(var) {
  # Example: Check skewness and apply log-transformation if needed
  if (shapiro.test(var)$p.value < 0.05) {
    var <- log(var)
    cat("Log-transformation applied to", response_var, "\n")
  }

  # Return the modified data
  return(data)
}

#' Transform Variable to Achieve Linear Relationships
#'
#' @description
#' Attempts to find a transformation that improves linear relationships between
#' an outcome variable and predictor variables using correlation analysis.
#'
#' @param df Data frame containing variables
#' @param y Character, name of outcome variable
#' @param cor.threshold Numeric, minimum acceptable correlation coefficient
#' @param method Character, correlation method ("pearson", "spearman", "kendall"). Default: "pearson"
#' @param use Character, handling of missing values in correlation
#'
#' @return List containing transformed data frame and transformation details
#'
#' @examples
#' \dontrun{
#' result <- straighten(df, "outcome_variable", cor.threshold = 0.3)
#' print(result$transformation)
#' print(result$correlations)
#' transformed_df <- result$data
#' }
#'
#' @export
straighten <- function(
  df,
  y,
  cor.threshold = 0.2,
  method = "pearson",
  use = "pairwise.complete.obs"
) {
  # Input validation
  if (!is.data.frame(df)) {
    stop("df must be a data frame")
  }
  if (!y %in% names(df)) {
    stop("y must be a column name in df")
  }
  if (!is.numeric(df[[y]])) {
    stop("y must be numeric")
  }
  if (any(df[[y]] <= 0)) {
    warning("y contains non-positive values, some transformations may fail")
  }

  # Define transformations to try
  transforms <- list(
    original = function(x) x,
    log = function(x) log(x),
    sqrt = function(x) sqrt(x),
    square = function(x) x^2,
    cube = function(x) x^3,
    inverse = function(x) 1 / x,
    boxcox = function(x) {
      if (requireNamespace("MASS", quietly = TRUE)) {
        bc <- MASS::boxcox(lm(x ~ 1), plotit = FALSE)
        lambda <- bc$x[which.max(bc$y)]
        if (abs(lambda) < 0.001) log(x) else (x^lambda - 1) / lambda
      } else {
        x
      }
    }
  )

  # Function to check correlations
  check_correlations <- function(transformed_y) {
    # Handle transformation errors
    if (any(is.na(transformed_y)) || any(is.infinite(transformed_y))) {
      return(FALSE)
    }

    # Create temporary df with transformed y
    temp_df <- df
    temp_df[[y]] <- transformed_y

    # Calculate correlations
    cors <- cor(
      temp_df[[y]],
      df[setdiff(names(df), y)],
      method = method,
      use = use
    )

    # Check if any correlation exceeds threshold
    any(abs(cors) >= cor.threshold)
  }

  # Try each transformation
  results <- list()
  for (transform_name in names(transforms)) {
    tryCatch(
      {
        transformed_y <- transforms[[transform_name]](df[[y]])
        if (check_correlations(transformed_y)) {
          # Store successful transformation
          transformed_df <- df
          transformed_df[[y]] <- transformed_y
          results[[transform_name]] <- list(
            data = transformed_df,
            transformation = transform_name,
            correlations = cor(
              transformed_df[[y]],
              df[setdiff(names(df), y)],
              method = method,
              use = use
            )
          )
        }
      },
      error = function(e) {
        warning(paste("Transformation", transform_name, "failed:", e$message))
      }
    )
  }

  # If no transformations worked
  if (length(results) == 0) {
    warning("No transformation achieved the desired correlation threshold")
    return(list(
      data = df,
      transformation = "none",
      correlations = cor(
        df[[y]],
        df[setdiff(names(df), y)],
        method = method,
        use = use
      )
    ))
  }

  # Find best transformation (highest max correlation)
  best_transform <- names(results)[which.max(sapply(results, function(x) {
    max(abs(x$correlations))
  }))]

  return(results[[best_transform]])
}

#' Transform Variables to Achieve Linear Relationships and Test Normality
#'
#' @description
#' Attempts to find transformations that improve linear relationships between
#' variables and tests for normality. Can work with either a single outcome
#' variable or all numeric variables in the dataset.
#'
#' @param df Data frame containing variables
#' @param y Character, name of outcome variable (optional)
#' @param cor.threshold Numeric, minimum acceptable correlation coefficient
#' @param method Character, correlation method ("pearson", "spearman", "kendall")
#' @param use Character, handling of missing values in correlation
#' @param exclude Character vector of column names to exclude
#' @param min_samples Minimum number of samples for normality testing
#' @param alpha Significance level for normality test
#'
#' @return List containing transformed data and analysis results
#'
#' @examples
#' \dontrun{
#' # For single outcome variable
#' result <- straighten(df, y = "outcome_variable")
#'
#' # For all numeric variables
#' result <- straighten(df)
#' }
#'
#' @export
straighten2 <- function(
  df,
  y = NULL,
  cor.threshold = 0.2,
  method = "pearson",
  use = "pairwise.complete.obs",
  exclude = NULL,
  min_samples = 3,
  alpha = 0.05
) {
  # Input validation
  if (!is.data.frame(df)) {
    stop("df must be a data frame")
  }

  # Initialize results list
  results <- list()

  # Define transformations
  transforms <- list(
    original = function(x) x,
    log = function(x) log(x),
    sqrt = function(x) sqrt(x),
    square = function(x) x^2,
    cube = function(x) x^3,
    inverse = function(x) 1 / x,
    boxcox = function(x) {
      if (requireNamespace("MASS", quietly = TRUE)) {
        bc <- MASS::boxcox(lm(x ~ 1), plotit = FALSE)
        lambda <- bc$x[which.max(bc$y)]
        if (abs(lambda) < 0.001) log(x) else (x^lambda - 1) / lambda
      } else {
        x
      }
    }
  )

  if (!is.null(y)) {
    # Single outcome variable mode
    if (!y %in% names(df)) {
      stop("y must be a column name in df")
    }

    if (!is.numeric(df[[y]])) {
      stop("y must be numeric")
    }

    # Check for non-positive values
    if (any(df[[y]] <= 0)) {
      warning("y contains non-positive values, some transformations may fail")
    }

    # Function to check correlations
    check_correlations <- function(transformed_y) {
      if (any(is.na(transformed_y)) || any(is.infinite(transformed_y))) {
        return(FALSE)
      }

      temp_df <- df
      temp_df[[y]] <- transformed_y
      cors <- cor(
        temp_df[[y]],
        df[setdiff(names(df), y)],
        method = method,
        use = use
      )

      any(abs(cors) >= cor.threshold)
    }

    # Try transformations
    transform_results <- list()
    for (transform_name in names(transforms)) {
      tryCatch(
        {
          transformed_y <- transforms[[transform_name]](df[[y]])

          if (check_correlations(transformed_y)) {
            transformed_df <- df
            transformed_df[[y]] <- transformed_y
            transform_results[[transform_name]] <- list(
              data = transformed_df,
              transformation = transform_name,
              correlations = cor(
                transformed_df[[y]],
                df[setdiff(names(df), y)],
                method = method,
                use = use
              )
            )
          }
        },
        error = function(e) {
          warning(paste("Transformation", transform_name, "failed:", e$message))
        }
      )
    }

    if (length(transform_results) == 0) {
      warning("No transformation achieved the desired correlation threshold")

      best_result <- list(
        data = df,
        transformation = "none",
        correlations = cor(
          df[[y]],
          df[setdiff(names(df), y)],
          method = method,
          use = use
        )
      )
    } else {
      best_transform <- names(transform_results)[which.max(sapply(
        transform_results,
        function(x) max(abs(x$correlations))
      ))]
      best_result <- transform_results[[best_transform]]
    }

    # Add normality test results
    best_result$normality <- multi_shapiro(best_result$data)
    return(best_result)
  } else {
    # Multi-variable mode
    numeric_cols <- names(df)[sapply(df, is.numeric)]

    if (!is.null(exclude)) {
      numeric_cols <- setdiff(numeric_cols, exclude)
    }

    if (length(numeric_cols) < 2) {
      stop("Need at least 2 numeric columns for correlation analysis")
    }

    # Initialize results
    transformed_df <- df
    transformation_summary <- data.frame(
      variable = character(),
      transformation = character(),
      max_correlation = numeric(),
      stringsAsFactors = FALSE
    )

    # Process each numeric column
    for (col in numeric_cols) {
      if (any(df[[col]] <= 0)) {
        warning(paste("Column", col, "contains non-positive values"))
      }

      best_transform <- "original"
      best_correlation <- 0
      best_data <- df[[col]]

      for (transform_name in names(transforms)) {
        tryCatch(
          {
            transformed_col <- transforms[[transform_name]](df[[col]])

            if (
              any(is.na(transformed_col)) || any(is.infinite(transformed_col))
            ) {
              next
            }

            other_cols <- setdiff(numeric_cols, col)
            cors <- sapply(other_cols, function(other_col) {
              cor(transformed_col, df[[other_col]], method = method, use = use)
            })

            max_cor <- max(abs(cors))

            if (max_cor > best_correlation) {
              best_correlation <- max_cor
              best_transform <- transform_name
              best_data <- transformed_col
            }
          },
          error = function(e) {
            warning(paste(
              "Transformation",
              transform_name,
              "failed for column",
              col,
              ":",
              e$message
            ))
          }
        )
      }

      transformed_df[[col]] <- best_data
      transformation_summary <- rbind(
        transformation_summary,
        data.frame(
          variable = col,
          transformation = best_transform,
          max_correlation = best_correlation,
          stringsAsFactors = FALSE
        )
      )
    }

    # Calculate final correlation matrix
    final_cors <- cor(transformed_df[numeric_cols], method = method, use = use)

    # Add normality test results
    normality_results <- multi_shapiro(transformed_df)

    return(list(
      data = transformed_df,
      transformations = transformation_summary,
      correlations = final_cors,
      normality = normality_results,
      threshold = cor.threshold,
      method = method
    ))
  }
}

# Print method
print.straighten_results <- function(x, ...) {
  if (!is.null(x$transformations)) {
    # Multi-variable mode
    cat("\nTransformation Summary:\n")
    print(x$transformations)

    cat("\nCorrelation Matrix of Transformed Variables:\n")
    print(round(x$correlations, 3))
  } else {
    # Single variable mode
    cat("\nBest Transformation:", x$transformation, "\n")
    cat("\nCorrelations with predictors:\n")
    print(round(x$correlations, 3))
  }

  cat("\nNormality Test Results:\n")
  print(data.frame(
    Variable = x$normality$variable,
    N = x$normality$n_samples,
    "Shapiro-W" = round(x$normality$statistic, 3),
    "p-value" = sprintf("%.3e", x$normality$p_value),
    "Normal" = x$normality$is_normal,
    "Skewness" = round(x$normality$skewness, 3),
    "Kurtosis" = round(x$normality$kurtosis, 3)
  ))

  cat("\nCorrelation method:", x$method, "\n")
  cat("Correlation threshold:", x$threshold, "\n")
}

#end

pad_columns <- function(df) {
  # Determine the maximum length among all columns
  max_length <- max(sapply(df, length))

  # Pad each column with NA values
  padded_df <- lapply(df, function(x) {
    if (length(x) < max_length) {
      c(x, rep(NA, max_length - length(x)))
    } else {
      x
    }
  })

  # Convert the result back to a data frame
  padded_df <- as.data.frame(padded_df)

  return(padded_df)
}

#determine the coefficient threshold for a particular fold change of interest
coef <- function(desired_fold_change, log_scale = 2) {
  log(desired_fold_change, base = log_scale)
}

##############################################################################
#  #start Wrangling functions for data manipulation                         #
##############################################################################

loadcombine <- function(
  file_dir,
  file_pattern,
  header = TRUE,
  sep,
  skip = 0,
  output
) {
  # Usage

  #loadcombine(file_dir = "C:/Users/Admin/OneDrive - University of Aberdeen (1)/Experiments/Sequencing/Transcriptomics/Nate-transcriptome_analysis/results/exploratory_MoA/GO/REVIGO_tiny/", file_pattern = "*.tsv", sep = "\t", output = "C:/Users/Admin/OneDrive - University of Aberdeen (1)/Experiments/Sequencing/Transcriptomics/Nate-transcriptome_analysis/results/exploratory_MoA/GO/REVIGO_tiny/REVIGO_results_combined.csv", header = T)

  currentDir <- getwd()
  setwd(file_dir)

  #create var with files paths as a character vector
  all_file_paths <- list.files(
    path = file_dir,
    pattern = file_pattern,
    full.names = TRUE
  )

  #load all the files into the environment
  all_files <- lapply(
    all_file_paths,
    read.table,
    header = header,
    sep = sep,
    skip = skip
  )

  all_filenames <-
    all_file_paths %>%
    basename() %>% #basename removes the path from the file name
    as.list() #turns them into a list

  #combines the file content list and file name list i.e. matches the names and the files
  all_lists <- mapply(c, all_files, all_filenames, SIMPLIFY = FALSE)

  all_result <- rbindlist(all_lists, fill = TRUE)
  names(all_result)[length(all_result)] <- "sample"
  #BUG all_result <- all_result[ ,colSums(is.na(all_result))<=ncol(all_result)]
  #export results file
  write.csv(all_result, file = output, row.names = FALSE)
  setwd(currentDir)
}

save_plot <- function(plot, filename) {
  svg(filename = paste0(filename, ".svg"))
  print(plot)
  dev.off()
  png(filename = paste0(filename, ".png"))
  print(plot)
  dev.off()
}

# Function to assign column classes based on a data dictionary
assign_col_classes <- function(
  data,
  data_dict_path,
  variable_name = "Variable.Name",
  classCol = "colClass.R"
) {
  # Read the data dictionary
  data_dict <- read.csv(data_dict_path, stringsAsFactors = FALSE)

  # Read the first row of the data dictionary to get the variable names
  variables <- data_dict[[variable_name]]

  # Read the first row of the data dictionary to get the classes
  classes <- data_dict[[classCol]]

  if (is.object(data) == FALSE) {
    data <- read.csv(data, header = TRUE)
  }

  # Match variable names to column names in the data dictionary
  col_classes <- sapply(colnames(data), function(var_name) {
    class_index <- match(var_name, variables)
    if (!is.na(class_index)) {
      return(classes[class_index])
    }
  })

  # Loop through columns and assign classes
  for (i in seq_along(col_classes)) {
    tryCatch(
      {
        if (is.null(col_classes[[i]])) {
          data[[i]] <- type.convert(data[[i]], as.is = TRUE)
        } else {
          if (col_classes[[i]] == "factor") {
            data[[i]] <- as.factor(data[[i]])
          } else {
            class(data[[i]]) <- col_classes[[i]]
          }
        }
      },
      error = function(e) {
        cat("Error assigning class to column:", names(data)[i], "\n")
        cat("Error message:", conditionMessage(e), "\n")
      }
    )
  }

  return(data)
}
# Function to assign column classes based on a data dictionary
assign_col_classes <- function(
  data,
  data_dict_path,
  variable_name = "Variable.Name",
  classCol = "colClass.R"
) {
  # Read the data dictionary
  data_dict <- read.csv(data_dict_path, stringsAsFactors = FALSE)

  # Read the first row of the data dictionary to get the variable names
  variables <- data_dict[[variable_name]]

  # Read the first row of the data dictionary to get the classes
  classes <- data_dict[[classCol]]

  if (is.object(data) == FALSE) {
    data <- read.csv(data, header = TRUE)
  }

  # Match variable names to column names in the data dictionary
  col_classes <- sapply(colnames(data), function(var_name) {
    class_index <- match(var_name, variables)
    if (!is.na(class_index)) {
      return(classes[class_index])
    }
  })

  # Loop through columns and assign classes
  for (i in seq_along(col_classes)) {
    tryCatch(
      {
        if (is.null(col_classes[[i]])) {
          data[[i]] <- type.convert(data[[i]], as.is = TRUE)
        } else {
          if (col_classes[[i]] == "factor") {
            data[[i]] <- as.factor(data[[i]])
          } else {
            class(data[[i]]) <- col_classes[[i]]
          }
        }
      },
      error = function(e) {
        cat("Error assigning class to column:", names(data)[i], "\n")
        cat("Error message:", conditionMessage(e), "\n")
      }
    )
  }

  return(data)
}

#' @title Utilities For \code{\link{phyloseq-class}} Slots to Data Frames

#' @description Utility to convert phyloseq slots to data frames.

#' @details These functions are a fork of \code{tibble_utilities} from the [\code{microbiome} package](https://github.com/microbiome/microbiome/blob/master/R/tibble_utilities.R). Instead of converting to tibbles, these functions instead convert to data frame format, thus retaining rownames instead of adding the rowname as a new column.
#'
#' Convert different \code{phyloseq} slots into data frames.
#' \code{otu_tibble} gets the otu_table in data frame format.
#' \code{tax_tibble} gets the taxa_table in data frame format.
#' \code{combine_otu_tax} combines otu_table and taxa_table into one data frame.
#' \code{sample_df} is effectively the same as \code{microbime::meta()}, and extracts the sample metadata as a data.frame.

#' @param x \code{\link{phyloseq-class}} object.
#' @param column.id Provide name for the column which will hold the rownames.
#'                  of slot.
#' @return A \code{data frame}
#' @examples
#' library(promicrobialTools)
#' data("dietswap")
#' otu_tib <- otu_tibble(dietswap,column.id="FeatureID")
#' tax_tib <- tax_tibble(dietswap,column.id="FeatureID")
#' sample_tib <- sample_tibble(dietswap,column.id="SampleID")
#' otu_tax <- combine_otu_tax(dietswap,column.id = "FeatureID")
#' head(otu_tax)
#'
#' @name DfUtilites
#' @author Contact: Nathaniel Cole \email{nc564@cornell.edu}
#' @seealso To convert OTU tables to abundacne matrices see \code{microbiome::abundances}. To extract
NULL

#' @rdname DfUtilites
#' @aliases otu_df
#' @export
otu_df <- function(x) {
  if (any(c("phyloseq", "otu_table") %in% is(x))) {
    require(microbiome)
        # Pick OTU matrix
    otu <- abundances(x)
  }
  otu_df <- otu %>%
    as.data.frame()
  return(otu_df)
}


#' @rdname DfUtilites
#' @aliases tax_df
#' @export
tax_df <- function(x, column.id = "FeatureID") {
  if (any(c("phyloseq", "tax_table") %in% is(x))) {
    # Pick OTU matrix
    tax <- phyloseq::tax_table(x)
    tax_df <- tax %>%
      as.matrix() %>%
      as.data.frame()
    return(tax_df)
  } else {
    return(message("tax_table not found"))
  }
}

#' @rdname DfUtilites
#' @aliases sample_df
#' @export
sample_df <- function(x, column.id = "sample_names") {
  if (any(c("phyloseq", "sample_data") %in% is(x))) {
    require(microbiome)
    # Pick OTU matrix
    smd <- meta(x)
    return(smd)
  } else {
    stop("Phyloseq object not supplied.")
  }
}

#' @rdname TibbleUtilites
#' @aliases combine_otu_tax
#' @export
combine_otu_tax <- function(x, column.id = "FeatureID") {
  otu_tb <- tax_tb <- NULL
  otu_tb <- otu_df(x, column.id)
  tax_tb <- tax_df(x, column.id)
  otu_tb <- tax_tb %>%
    dplyr::left_join(otu_tb, by = column.id)
  return(otu_tb)
}

#start check_metadata function

#' Check Metadata Variable Properties
#'
#' @description
#' Checks and reports the class, number of NA values, and length of a metadata variable.
#' This function is designed to validate variables used in microbiome analyses,
#' ensuring they meet expected specifications.
#'
#' @param var The variable to check (vector)
#' @param expected_class Character string specifying the expected class of the variable
#' @param name Character string specifying the name of the variable for reporting
#'
#' @return Prints a summary of the variable's properties:
#'   \itemize{
#'     \item Class (actual and expected)
#'     \item Number of NA values
#'     \item Length of variable
#'   }
#'
#' @examples
#' subject <- c("A1", "A2", "A3")
#' check_metadata(subject, "character", "subject")
#'
#' age <- c(25, 30, NA, 40)
#' check_metadata(age, "numeric", "age")
#'
#' group <- factor(c("Control", "Treatment"))
#' check_metadata(group, "factor", "group")
#'
#' @seealso
#' \code{\link{class}}, \code{\link{is.na}}
#'
#' @export
#'
#' @note
#' This function is particularly useful for checking metadata variables
#' before conducting microbiome analyses where complete cases and
#' correct variable types are crucial.
#'
#' @author Nathaniel Cole
#' @version 1.0.0
check_metadata <- function(var) {
  cat("Checking:", deparse(substitute(var)), "\n")
  cat("Class:", class(var), "\n")
  cat("NAs:", sum(is.na(var)), "\n")
  cat("Length:", length(var), "\n")
  cat("-------------------\n")
}

#end

#' Create a Formatted Table from Multiple Results Objects
#'
#' @description
#' Combines multiple results objects into a single formatted table.
#' Compatible with emmeans, data.frames, matrices, and other objects that can be
#' coerced to data.frames. Provides flexible formatting options and is compatible
#' with knitr/quarto rendering.
#'
#' @param results_list List of objects to combine. Names of list elements used as identifiers.
#' @param id_col Character. Name for the identifier column (default: "Model")
#' @param digits Integer. Number of decimal places for numeric columns (default: 2)
#' @param format Character. Output format: "markdown", "html", or "latex" (default: "markdown")
#' @param caption Character. Optional table caption
#' @param col.names Character vector. Custom column names
#' @param round_cols Character vector. Columns to round (default: all numeric)
#' @param align Character vector. Column alignment (default: auto-detect)
#' @param table_style Character vector. kableExtra styling options
#'
#' @return A knitr_kable object that can be directly used in R Markdown or Quarto
#'
#' @examples
#' \dontrun{
#' # With emmeans objects
#' results1 <- list(
#'   "Model 1" = emmeans(model1, ~treatment),
#'   "Model 2" = emmeans(model2, ~treatment)
#' )
#'
#' # With data frames
#' results2 <- list(
#'   "Analysis 1" = data.frame(x = 1:3, y = 4:6),
#'   "Analysis 2" = data.frame(x = 7:9, y = 10:12)
#' )
#'
#' # Create tables
#' table1 <- table_bind(
#'   results1,
#'   caption = "EMM Comparison",
#'   table_style = c("striped", "hover")
#' )
#'
#' table2 <- table_bind(
#'   results2,
#'   id_col = "Analysis",
#'   align = c('l', 'r', 'r')
#' )
#' }
#'
#' @export
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
table_bind <- function(
  results_list,
  id_col = "Model",
  digits = 2,
  format = c("markdown", "html", "latex"),
  caption = NULL,
  col.names = NULL,
  round_cols = NULL,
  align = NULL,
  table_style = NULL
) {
  # Input validation
  if (!is.list(results_list)) {
    stop("Input must be a list of results objects")
  }

  if (length(results_list) == 0) {
    stop("Empty list provided")
  }

  # Match format argument
  format <- match.arg(format)

  # Convert each object to data frame
  df_list <- lapply(results_list, function(x) {
    tryCatch(
      {
        # Handle different object types
        if (inherits(x, "emmGrid")) {
          as.data.frame(x)
        } else if (inherits(x, "matrix")) {
          as.data.frame(x)
        } else if (is.data.frame(x)) {
          x
        } else {
          as.data.frame(x)
        }
      },
      error = function(e) {
        stop("Error converting object to data frame: ", e$message)
      }
    )
  })

  # Check for consistent columns
  col_names <- lapply(df_list, colnames)
  if (length(unique(lapply(col_names, length))) > 1) {
    stop("Inconsistent number of columns in objects")
  }

  # Create results data frame
  results <- data.frame(
    ID = rep(names(results_list), sapply(df_list, nrow)),
    do.call(rbind, df_list),
    check.names = FALSE
  )
  names(results)[1] <- id_col

  # Round numeric columns
  if (is.null(round_cols)) {
    numeric_cols <- sapply(results, is.numeric)
    results[, numeric_cols] <- round(results[, numeric_cols], digits)
  } else {
    results[, round_cols] <- round(results[, round_cols], digits)
  }

  # Use custom column names if provided
  if (!is.null(col.names)) {
    if (length(col.names) != ncol(results)) {
      stop("Length of col.names must match number of columns")
    }
    colnames(results) <- col.names
  }

  # Auto-detect alignment if not specified
  if (is.null(align)) {
    align <- ifelse(sapply(results, is.numeric), 'r', 'l')
  }

  # Create kable
  table <- knitr::kable(
    results,
    format = format,
    digits = digits,
    caption = caption,
    align = align,
    booktabs = TRUE
  )

  # Apply kableExtra styling if requested
  if (!is.null(table_style) && requireNamespace("kableExtra", quietly = TRUE)) {
    table <- kableExtra::kable_styling(table, bootstrap_options = table_style)
  }

  # Add class for method dispatch
  class(table) <- c("bound_table", class(table))

  return(table)
}

# Print method
#' @export
print.bound_table <- function(x, ...) {
  print(unclass(x))
}

# Knit print method
#' @export
knit_print.bound_table <- function(x, ...) {
  knitr::knit_print(unclass(x), ...)
}
#end

#start EMMs functions for EMMs
#################
#     EMMs      #
#################

#' Bind Multiple Estimated Marginal Means into a Table
#'
#' @description
#' Combines multiple emmeans objects into a single formatted table, with model names
#' as identifiers. The output is compatible with knitr/quarto rendering.
#'
#' @param emm List of emmeans objects. Names of the list elements are used as model identifiers.
#' @param digits Integer. Number of decimal places for numeric columns (default: 2)
#' @param format Character. Output format: "markdown", "html", or "latex" (default: "markdown")
#' @param caption Character. Optional table caption
#' @param col.names Character vector. Custom column names
#'
#' @return A knitr_kable object that can be directly used in R Markdown or Quarto
#'
#' @examples
#' \dontrun{
#' # Create multiple emmeans objects
#' emm1 <- emmeans(model1, ~treatment)
#' emm2 <- emmeans(model2, ~treatment)
#'
#' # Combine into table
#' emm_table <- emm_bind(
#'   list("Model 1" = emm1, "Model 2" = emm2),
#'   digits = 3,
#'   caption = "Estimated Marginal Means Comparison"
#' )
#'
#' # Print table
#' emm_table
#' }
#'
#' @export
#' @importFrom knitr kable
#' @importFrom emmeans emmeans
emm_bind <- function(
  emm,
  digits = 2,
  format = "markdown",
  caption = NULL,
  col.names = NULL,
  kable = TRUE
) {
  # Input validation
  if (!is.list(emm)) {
    stop("'emm' must be a list of emmeans objects")
  }

  if (length(emm) == 0) {
    stop("Empty list provided")
  }

  if (!all(sapply(emm, inherits, "emmGrid"))) {
    stop("All elements must be emmeans objects (class 'emmGrid')")
  }

  # Match format argument
  format <- match.arg(format)

  # Convert each emmeans object to data frame
  df_list <- lapply(emm, function(x) {
    tryCatch(
      {
        as.data.frame(x)
      },
      error = function(e) {
        stop("Error converting emmeans to data frame: ", e$message)
      }
    )
  })

  # Check for consistent columns
  col_names <- lapply(df_list, colnames)
  if (length(unique(lapply(col_names, length))) > 1) {
    warning("Inconsistent number of columns in emmeans objects")
  }

  # Create results data frame
  results <- data.frame(
    Model = rep(names(df_list), sapply(df_list, nrow)),
    data.table::rbindlist(df_list, fill = TRUE),
    check.names = FALSE
  )

  # Round numeric columns to specified digits
  numeric_cols <- sapply(results, is.numeric)
  results[numeric_cols] <- round(results[numeric_cols], digits)

  # Use custom column names if provided
  if (!is.null(col.names)) {
    if (length(col.names) != ncol(results)) {
      stop("Length of col.names must match number of columns")
    }
    colnames(results) <- col.names
  }

  # Create kable with appropriate format
  if(kable){
    table <- knitr::kable(
      results,
      format = format,
      digits = digits,
      caption = caption,
      booktabs = TRUE
    )
  } else {
    table <- results
  }

  # Add class for potential method dispatch
  class(table) <- c("emm_table", class(table))

  return(table)
}

plot.emm_table = function(x, 
                         title = "EMM Contrasts", 
                         subtitle = NULL,
                         sig_levels = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
                         colors = NULL,
                         base_size = 11,comparisons = TRUE,
                         alpha = 0.05,
                         ...) {
    
  # Input validation
      if (!inherits(x, "emm_table")) {
          stop("Input must be of class 'emm_table'")
      }
      
      # Required columns check
      required_cols <- c("Model", "estimate", "SE", "p.value", "Model", "Species")
      if (!all(required_cols %in% colnames(x))) {
          stop("Missing required columns: ", 
              paste(setdiff(required_cols, colnames(x)), collapse = ", "))
      }

      # Create significance levels
  # Create combined labels and add to data
    x <- x %>%
        mutate(
            combined_label = paste0(Model, "; ", Species),
            sig_level = case_when(
                p.value < sig_levels["***"] ~ "***",
                p.value < sig_levels["**"] ~ "**",
                p.value < sig_levels["*"] ~ "*",
                TRUE ~ "ns"
            ),
            # Add y-position for each point
            y_position = match(combined_label, unique(combined_label))
        )

    # Calculate plot dimensions
    max_est <- ceiling(max(abs(x$estimate + x$SE)))
    x_padding <- max_est * 0.15  # 15% padding
 
  # Add comparison arrows if requested
    if(comparisons) {
        # Create comparison data
        comp_data <- x %>%
            filter(p.value < alpha) %>%
            mutate(
                # Calculate arrow positions
                arrow_start = estimate - SE,
                arrow_end = estimate + SE,
                # Y positions for arrows
                y_pos = y_position,
                # Arrow aesthetic parameters based on significance
                arrow_alpha = case_when(
                    p.value < sig_levels["***"] ~ 1,
                    p.value < sig_levels["**"] ~ 0.7,
                    p.value < sig_levels["*"] ~ 0.4,
                    TRUE ~ 0.2
                )
            )
        
        # Add arrows
        p <- p + 
            geom_segment(data = comp_data,
                        aes(x = 0, 
                            xend = estimate,
                            y = y_pos,
                            yend = y_pos,
                            alpha = arrow_alpha),
                        arrow = arrow(length = unit(0.2, "cm"), 
                                    type = "closed"),
                        color = "gray30") +
            scale_alpha_identity()
    }
    
    # Add point ranges
    p <- p + geom_pointrange(
        aes(xmin = estimate - SE, 
            xmax = estimate + SE),
        size = 1,
        position = position_dodge(width = 0.5)
    )
    
    # Add significance indicators
    if (max(nchar(x$sig_level)) > 0) {
        p <- p + geom_text(
            aes(x = max(estimate) + x_padding/2,
                label = sig_level),
            hjust = 0,
            size = 3
        )
    }
    
    # Theming and labels
    p <- p + 
        labs(
            title = title,
            subtitle = subtitle,
            x = "EMM Contrast Estimate",
            y = NULL
        ) +
        theme(
            axis.text.y = element_text(size = base_size * 0.8),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = base_size * 1.2),
            plot.subtitle = element_text(size = base_size * 0.9)
        ) +
        scale_x_continuous(
            limits = c(-max_est - x_padding, max_est + x_padding),
            expand = expansion(mult = 0.1)
        )
    
    # Add custom colors if provided
    if (!is.null(colors)) {
      p <- p + scale_color_manual(values = colors)
    } else {
      p <- p + scale_color_brewer("Set2")
    }

    
    # Add significance legend
    p <- p + 
        annotate(
            "text",
            x = max_est + x_padding/2,
            y = 1,
            label = "Significance:\n*** p<0.001\n** p<0.01\n* p<0.05",
            hjust = 0,
            vjust = 0,
            size = 3
        )
    
    return(p)
}

#' Plot Phylogenetic Estimated Marginal Means
#'
#' @description
#' Creates a forest plot of estimated marginal means (EMMs) for phylogenetic data,
#' with colored labels for species present across multiple timepoints. The plot includes
#' error bars, significance indicators, and phylum-based coloring.
#'
#' @param df A data frame containing the following required columns:
#'   * `contrast`: The comparison being made
#'   * `Model`: Model identifier
#'   * `Species`: Species name
#'   * `Genus`: Genus name
#'   * `Phylum`: Phylum classification
#'   * `estimate`: EMM estimate (or `ratio` if estimate not available)
#'   * `SE`: Standard error
#'   * `p.value`: P-value for statistical significance
#'   * `timepoint`: Timepoint identifier
#'   Optional column:
#'   * `effect.size`: Effect size for magnitude calculation
#'
#' @param comparisons Logical, whether to show comparisons (default: TRUE)
#' @param common_color Character, color for species present in multiple timepoints (default: "red")
#' @param other_color Character, color for species present in single timepoint (default: "black")
#' @param common_weight Character, font weight for common species (default: "bold")
#' @param other_weight Character, font weight for other species (default: "normal")
#' @param text_size Numeric, size of y-axis text (default: 8)
#'
#' @return A ggplot2 object containing the forest plot with:
#'   * Colored species labels (red for common species)
#'   * Point estimates with error bars
#'   * Phylum-based color coding
#'   * Significance indicators
#'   * Reference line at 0
#'   * Magnitude indicators (if effect.size is present)
#'
#' @details
#' The function automatically:
#' * Calculates magnitude categories if effect.size is present:
#'   - negligible: |effect.size| < 0.2
#'   - small: 0.2 ≤ |effect.size| < 0.5
#'   - medium: 0.5 ≤ |effect.size| < 0.8
#'   - large: |effect.size| ≥ 0.8
#' * Adds significance indicators:
#'   - ***: p < 0.001
#'   - **: p < 0.01
#'   - *: p ≤ 0.05
#'   - ns: p > 0.05
#'
#' @examples
#' # Basic usage
#' plotPhylaEMMs(emm_results)
#'
#' # Custom styling
#' plotPhylaEMMs(
#'   emm_results,
#'   common_color = "#E41A1C",
#'   other_color = "#666666",
#'   text_size = 10
#' )
#'
#' # Without effect size
#' plotPhylaEMMs(emm_results_no_effect_size)
#'
#' @import ggplot2
#' @import dplyr
#' @import forcats
#' @import ggtext
#' @importFrom scales pretty_breaks
#'
#' @seealso
#' * [emmeans] for generating estimated marginal means
#' * [ggplot2] for additional plot customization
#'
#' @export
plotPhylaEMMs <- function(df, 
                         comparisons = TRUE,
                         highlight_common = TRUE,
                         common_color = "red",
                         other_color = "black",
                         common_weight = "bold",
                         other_weight = "normal",
                         text_size = 10) {
  
  contrast <- unique(df$contrast)
  
  # Effect size processing
  if(any(colnames(df) == "effect.size")){
    df <- df %>%
      mutate(
        magnitude = case_when(
          abs(effect.size) < 0.2 ~ "negligible",
          abs(effect.size) >= 0.2 & abs(effect.size) < 0.5 ~ "small",
          abs(effect.size) >= 0.5 & abs(effect.size) < 0.8 ~ "medium",
          abs(effect.size) >= 0.8 ~ "large"
        )
      )
    legend = "legend"
  } else {
    df <- df %>%
      mutate(magnitude = "1")
    legend = NULL
  }
  
  # Data preparation
  df <- df %>%
    arrange(desc(Genus)) %>%
    mutate(
      Species = paste0(Model, "; ", 
        case_when(
          Species != "" ~ Species,
          Genus == "" ~ Phylum,
          TRUE ~ Genus
        )),
      sig_level = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value <= 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      estimate = ifelse(is.na(estimate), ratio, estimate)
    )

  if(highlight_common == TRUE){
    # Find common species
    common_species <- df %>%
      group_by(Species) %>%
      summarize(n = n_distinct(timepoint)) %>%
      filter(n > 1) %>%
      pull(Species)
  # Add colored labels
    df <- df %>%
      mutate(
        Species_colored = case_when(
          Species %in% common_species ~ 
            sprintf("<span style='color:%s;font-weight:%s'>%s</span>",
                  common_color, common_weight, Species),
          TRUE ~ 
            sprintf("<span style='color:%s;font-weight:%s'>%s</span>",
                  other_color, other_weight, Species)
        )
      )
    } else {
          df <- df %>% mutate(Species_colored = Species)
  }
  
  # Create plot
  max_estimate <- max(df$estimate, na.rm = TRUE)
  significance_x_pos <- max_estimate + 0.5
  
  p <- ggplot(df, aes(x = estimate, y = forcats::fct_inorder(Species_colored))) +
    # Add reference line at 0
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    # Add horizontal bars for EMM contrasts
    geom_pointrange(
      aes(xmin = estimate - SE, xmax = estimate + SE, 
          color = Phylum, shape = magnitude), 
      size = 1
    ) +
    # Add significance indicators
    geom_text(
      aes(x = significance_x_pos, label = sig_level), 
      hjust = 0
    ) +
    scale_color_brewer(palette = "Paired") +
    scale_shape(guide = legend) +
    theme(
      axis.text.y = element_markdown(size = text_size),
      panel.grid.major.y = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      caption = paste("Contrast:", contrast),
      x = "Estimated Marginal Mean Contrast Estimate",
      y = "",
      color = "Phylum"
    )
  
  # Add significance level annotation
  if(exists("sigArms")) {  # Check if sigArms exists
    p <- p + annotate(
      "text",
      x = max(sigArms$estimate) + 0.5,
      y = 0,
      label = "* p<0.05\n** p<0.01\n*** p<0.001",
      hjust = 0,
      vjust = 1
    )
  }
  
  return(p)
}

#end

#' Color Common Labels
#' 
#' Creates a plot with colored labels for items that appear in multiple groups
#' @param data A data frame
#' @param grp Column name for grouping variable (e.g., "timepoint")
#' @param y Column name for y-axis variable
#' @param x Column name for x-axis variable (default: NULL)
#' @param fill Column name for fill variable (default: same as grp)
#' @param facet Column name for faceting (default: NULL)
#' @param common_color Color for common items (default: "red")
#' @param other_color Color for other items (default: "black")
#' @param common_weight Font weight for common items (default: "normal")
#' @param other_weight Font weight for other items (default: "normal")
#' @param text_size Size of y-axis text (default: 10)
#' @param ... Additional arguments passed to ggplot(aes())
#' @return A ggplot object
color_common_labels <- function(data, 
                              grp, 
                              y, 
                              x = NULL,
                              fill = NULL,
                              facet = NULL,
                              common_color = "red",
                              other_color = "black",
                              common_weight = "normal",
                              other_weight = "normal",
                              text_size = 10,
                              position = "dodge",
                              ...) {
  
  # Input validation
  if (!all(c(grp, y) %in% colnames(data))) {
    stop("Specified columns not found in data")
  }
  
  # Set fill to grp if not specified
  if (is.null(fill)) fill <- grp
  
  # Find common items
  common <- data %>%
    group_by(!!sym(y)) %>%
    summarize(n = n_distinct(!!sym(grp))) %>%
    filter(n > 1) %>%
    pull(!!sym(y))
  
  # Add colored labels to data
  data <- data %>%
    mutate(
      label_colored = case_when(
        !!sym(y) %in% common ~ 
          sprintf("<span style='color:%s;font-weight:%s'>%s</span>",
                 common_color, common_weight, !!sym(y)),
        TRUE ~ 
          sprintf("<span style='color:%s;font-weight:%s'>%s</span>",
                 other_color, other_weight, !!sym(y))
      )
    )
  
  # Create base plot
  if (is.null(x)) {
    # Default to using row numbers if x is not specified
    p <- ggplot(data, aes(y = label_colored, fill = !!sym(fill), ...))
  } else {
    p <- ggplot(data, aes(x = !!sym(x), y = label_colored, fill = !!sym(fill), ...))
  }
  
  # Add geom based on whether x is specified
  if (is.null(x)) {
    p <- p + geom_bar(position = position)
  } else {
    p <- p + geom_col(position = position)
  }
  
  # Add faceting if specified
  if (!is.null(facet)) {
    p <- p + facet_wrap(as.formula(paste("~", facet)), scales = "free_y")
  }
  
  # Add theme elements
  p <- p + 
    theme_minimal() +
    theme(
      axis.text.y = element_markdown(size = text_size),
      panel.grid.major.y = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(y = NULL)
  
  return(p)
}

# #TODO #ADD Example usage

# # Basic usage
# p1 <- color_common_labels(
#   data = data,
#   grp = "timepoint",
#   y = "Model",
#   x = "estimate",
#   text_size = 8
# )

# # With faceting
# p2 <- color_common_labels(
#   data = data,
#   grp = "timepoint",
#   y = "Model",
#   x = "estimate",
#   facet = "Phylum",
#   common_color = "#E41A1C",
#   other_color = "#666666",
#   common_weight = "bold",
#   text_size = 8
# )

# # With custom aesthetics
# p3 <- color_common_labels(
#   data = data,
#   grp = "timepoint",
#   y = "Model",
#   x = "estimate",
#   fill = "change",  # Using different fill variable
#   facet = "Phylum",
#   common_color = "darkred",
#   other_color = "grey30",
#   common_weight = "bold",
#   text_size = 8
# ) +
#   scale_fill_manual(values = c("Increase" = "forestgreen", "Decrease" = "red"))

# # Print plots
# print(p1)
# print(p2)
# print(p3)

# # Example with different position
# p4 <- color_common_labels(
#   data = data,
#   grp = "timepoint",
#   y = "Model",
#   x = "estimate",
#   position = "stack"
# )


calculate_pairwise_vda <- function(data, value_col, group_col) {
  # Get unique groups
  groups <- unique(data[[group_col]])
  n_groups <- length(groups)

  # Create empty results dataframe
  results <- data.frame(
    group1 = character(),
    group2 = character(),
    VDA = numeric(),
    effect_size = character(),
    stringsAsFactors = FALSE
  )

  # Calculate VDA for each pair
  for (i in 1:(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      group1_data <- data[data[[group_col]] == groups[i], ][[value_col]]
      group2_data <- data[data[[group_col]] == groups[j], ][[value_col]]

      # Calculate VDA
      score <- 0
      for (x1 in group1_data) {
        for (x2 in group2_data) {
          if (x1 > x2) {
            score <- score + 1
          }
          if (x1 == x2) score <- score + 0.5
        }
      }

      vda <- score / (length(group1_data) * length(group2_data))

      # Determine effect size
      effect <- case_when(
        vda >= 0.71 ~ "large",
        vda >= 0.64 ~ "medium",
        vda >= 0.56 ~ "small",
        vda <= 0.29 ~ "large negative",
        vda <= 0.36 ~ "medium negative",
        vda <= 0.44 ~ "small negative",
        TRUE ~ "negligible"
      )

      # Add to results
      results <- rbind(
        results,
        data.frame(
          group1 = groups[i],
          group2 = groups[j],
          VDA = vda,
          effect_size = effect
        )
      )
    }
  }

  return(results)
}

#' Create Tertiles from a Continuous Variable
#'
#' @description
#' Divides a continuous variable into three groups (tertiles) while ensuring that
#' equal values are not split across different tertiles. Returns both the categorized
#' data and a summary of the tertile distribution.
#'
#' @param continuous_variable A numeric vector to be divided into tertiles
#'
#' @return A list containing:
#'   \itemize{
#'     \item tertiles: An ordered factor with levels "Low", "Medium", "High"
#'     \item summary: A data frame with columns:
#'       \itemize{
#'         \item Tertile: Character indicating tertile level
#'         \item N: Integer count of observations
#'         \item Percent: Numeric percentage of observations
#'       }
#'     \item breaks: Numeric vector of break points (excluding -Inf and Inf)
#'   }
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' result <- tertiles(x)
#' print(result$summary)
#' print(result$breaks)
#'
#' @export
#'
#' @note
#' The function ensures that equal values are not split across tertiles,
#' which may result in groups that deviate from exact thirds.
#'
#' @seealso
#' \code{\link{cut}}, \code{\link{quantile}}
#'
#' @importFrom stats quantile
#' @importFrom knitr kable
#'
tertiles <- function(continuous_variable) {
  # Input validation
  if (!is.numeric(continuous_variable)) {
    stop("Input must be a numeric vector")
  }

  if (length(continuous_variable) < 3) {
    stop("Input vector must have at least 3 values to create tertiles")
  }

  # Remove NA values
  x <- continuous_variable[!is.na(continuous_variable)]

  if (length(x) == 0) {
    stop("No non-NA values in input vector")
  }

  # Get unique values and sort them
  unique_values <- sort(unique(x))

  # Calculate cumulative percentages for unique values
  cum_pct <- cumsum(table(x)) / length(x)

  # Initialize tertile breaks
  breaks <- c(-Inf, Inf) # Default breaks if logic fails

  # Find the closest values to 33.33% and 66.67% that don't split equal values
  if (length(unique_values) > 2) {
    # Find break points that don't split equal values
    lower_third <- unique_values[which.min(abs(cum_pct - 1 / 3))]
    upper_third <- unique_values[which.min(abs(cum_pct - 2 / 3))]

    breaks <- c(-Inf, lower_third, upper_third, Inf)
  }

  # Create tertiles
  tertiles <- cut(
    continuous_variable,
    breaks = breaks,
    labels = c("Low", "Medium", "High"),
    include.lowest = TRUE,
    ordered = TRUE
  ) # Make it an ordered factor

  # Create summary statistics
  summary_stats <- data.frame(
    Tertile = c("Low", "Medium", "High"),
    N = as.numeric(table(tertiles)),
    Percent = round(as.numeric(prop.table(table(tertiles)) * 100), 1)
  )

  # Print summary
  cat("Tertile Distribution:")
  knitr::kable(summary_stats)
  cat("Break Points:")
  cat(breaks[-c(1, length(breaks))]) # Don't show -Inf and Inf

  return(tertiles)
}

#' Process Heatmap Data with Metadata
#'
#' @description
#' Processes heatmap data by combining it with metadata and calculating summary statistics
#'
#' @param hm A list containing heatmap data with a $data element
#' @param i The name of the grouping variable
#' @param metadata A data frame containing metadata including sample_names and grouping variables
#'
#' @return A data frame with summary statistics grouped by phyla and the specified grouping variable
#'
#' @examples
#' # Example usage:
#' # hm_data(hm = hmaps$b12_tertiles, i = "b12_tertiles", metadata = b12_base_meta)
#'
#' @export
hm_data <- function(hm, i, metadata) {
  # Input validation
  if (!is.list(hm) || is.null(hm$data)) {
    stop("'hm' must be a list containing a 'data' element")
  }
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data frame")
  }
  if (!all(c("sample_names", i) %in% colnames(metadata))) {
    stop(
      "metadata must contain 'sample_names' and the grouping variable '",
      i,
      "'"
    )
  }

  # Process data
  df <- hm$data %>%
    dplyr::select(phyla = Display, Sample, relAbund = Sum) %>%
    # Join with metadata using proper column names
    left_join(
      metadata %>%
        dplyr::select(Sample = sample_names, grouping_var = !!sym(i)),
      by = "Sample"
    ) %>%
    # Calculate summary statistics
    group_by(phyla, grouping_var) %>%
    summarise(
      mean = mean(relAbund, na.rm = TRUE),
      median = median(relAbund, na.rm = TRUE),
      max = max(relAbund, na.rm = TRUE),
      min = min(relAbund, na.rm = TRUE),
      n = n(),
      `Median (IQR; Min-Max)` = paste0(
        round(median(relAbund, na.rm = TRUE), 3),
        "\n(",
        round(IQR(relAbund), 3),
        "; ",
        round(min(relAbund, na.rm = TRUE), 3),
        "-",
        round(max(relAbund, na.rm = TRUE), 3),
        ")"
      ),
      .groups = "drop"
    )

  return(df)
}

#' Create a Phylum-Level Heatmap from ampvis2 Object
#'
#' @description
#' Generates a customized heatmap visualization of microbial community data at the phylum level
#' using the ampvis2 package. The function includes enhanced formatting options and input validation.
#'
#' @param ampvisObject An ampvis2 object containing microbial community data
#' @param group_by Character. The grouping variable name from the metadata
#' @param tax_show Numeric. Number of taxa to show (default: 20)
#' @param plot_values Logical. Whether to show values in cells (default: TRUE)
#' @param value_size Numeric. Size of value labels (default: 8)
#' @param decimals Integer. Number of decimal places for values (default: 3)
#' @param text_size Numeric. Size of axis text (default: 14)
#' @param axis_angle Numeric. Angle of x-axis labels (default: 45)
#' @param color_palette Character vector. Custom color palette (optional)
#'
#' @return A ggplot2 object containing the heatmap visualization
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' phylum_heatmap(amp_data, "Treatment")
#'
#' # With custom parameters
#' phylum_heatmap(amp_data, "Treatment",
#'                tax_show = 15,
#'                plot_values = FALSE,
#'                text_size = 12)
#' }
#'
#' @export
#'
#' @importFrom ampvis2 amp_heatmap
#' @importFrom ggplot2 theme element_text
#' @importFrom assertthat assert_that is.string
#'
#' @seealso
#' \code{\link[ampvis2]{amp_heatmap}}
#'
phylum_heatmap <- function(
  ampvisObject,
  group_by,
  tax_show = 20,
  plot_values = TRUE,
  value_size = 8,
  decimals = 2,
  text_size = 8,
  axis_angle = 45,
  color_palette = NULL,
  order_x_by = NULL
) {
  require(assertthat)

  # Input validation
  assert_that(
    inherits(ampvisObject, "ampvis2"),
    msg = "Input must be an ampvis2 object"
  )

  assert_that(
    is.string(group_by),
    msg = "group_by must be a character string"
  )

  assert_that(
    group_by %in% colnames(ampvisObject$metadata),
    msg = sprintf("'%s' not found in metadata", group_by)
  )

  assert_that(
    is.numeric(tax_show),
    tax_show > 0,
    msg = "tax_show must be a positive number"
  )

  # Create heatmap
  heatmap <- amp_heatmap(
    ampvisObject,
    group_by = group_by,
    tax_aggregate = "Phylum",
    tax_show = tax_show,
    plot_values = plot_values,
    plot_values_size = value_size,
    normalise = TRUE,
    round = decimals,
    tax_empty = "OTU",
    order_x_by = order_x_by,
    order_y_by = "cluster"
  )

  # Add custom theme
  heatmap <- heatmap +
    theme(
      axis.text.x = element_text(
        angle = axis_angle,
        size = text_size,
        vjust = 1,
        hjust = 1
      ),
      axis.text.y = element_text(
        size = text_size
      ),
      legend.position = "right",
      legend.title = element_text(
        color = "black",
        size = text_size * 0.7,
        face = "bold"
      ),
      # Additional theme elements for better visualization
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title = element_text(size = text_size, face = "bold")
    )

  # Apply custom color palette if provided
  if (!is.null(color_palette)) {
    heatmap <- heatmap +
      scale_fill_gradientn(colours = color_palette)
  }

  # Add attributes for Quarto
  attr(heatmap, "caption") <- sprintf(
    "Phylum-level heatmap grouped by %s",
    group_by
  )

  class(heatmap) <- c("phylum_heatmap", class(heatmap))

  return(heatmap)
}

#' Print method for phylum_heatmap objects
#'
#' @param x A phylum_heatmap object
#' @param ... Additional arguments passed to print
#'
#' @export
print.phylum_heatmap <- function(x, ...) {
  print(ggplot2::ggplot(x))
}

#' Save phylum heatmap with optimal dimensions
#'
#' @param heatmap A phylum_heatmap object
#' @param filename Character. Path to save the file
#' @param width Numeric. Width in inches (default: 12)
#' @param height Numeric. Height in inches (default: 8)
#' @param dpi Numeric. Resolution (default: 300)
#'
#' @export
save_phylum_heatmap <- function(
  heatmap,
  filename,
  width = 12,
  height = 8,
  dpi = 300
) {
  ggsave(
    filename = filename,
    plot = heatmap,
    width = width,
    height = height,
    dpi = dpi
  )
}

# phylum_heatmap <- function(ampvisObject, group_by){
#   heatmap <- amp_heatmap(ampvisObject,
#   group_by = group_by,
#   tax_aggregate = "Phylum",
#   tax_show = 20,
#   plot_values = T,
#   normalise = T,
#   plot_values_size = 8,
#   round = 3,
#   tax_empty="OTU",
#   order_x_by = order_x_by,
#   order_y_by = "cluster") +
#   theme(axis.text.x = element_text(angle = 45, size = 14, vjust = 1), axis.text.y = element_text(size = 14),  legend.position = "right", legend.title = element_text(color = "black", size = 10, face = "bold"))

#   return(heatmap)
# }

#default settings plot relative abundances (normalise = TRUE)
phyla_boxplot <- function(ampvis_object, group) {
  amp_boxplot(
    ampvis_object,
    normalise = TRUE,
    group_by = group,
    tax_aggregate = "Phylum",
    tax_show = 5,
    point_size = 1,
    plot_flip = T
  ) +
    scale_color_brewer(palette = "Set2")
}

#' Create a Genus-Level Heatmap from ampvis2 Object
#'
#' @description
#' Generates a customized heatmap visualization of microbial community data at the genus level
#' using the ampvis2 package. The function includes enhanced formatting options and input validation.
#'
#' @param ampvisObject An ampvis2 object containing microbial community data
#' @param group_by Character. The grouping variable name from the metadata
#' @param tax_show Numeric. Number of taxa to show (default: 20)
#' @param plot_values Logical. Whether to show values in cells (default: TRUE)
#' @param value_size Numeric. Size of value labels (default: 8)
#' @param decimals Integer. Number of decimal places for values (default: 3)
#' @param text_size Numeric. Size of axis text (default: 14)
#' @param axis_angle Numeric. Angle of x-axis labels (default: 45)
#' @param color_palette Character vector. Custom color palette (optional)
#'
#' @return A ggplot2 object containing the heatmap visualization
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' genus_heatmap(amp_data, "Treatment")
#'
#' # With custom parameters
#' genus_heatmap(amp_data, "Treatment",
#'                tax_show = 15,
#'                plot_values = FALSE,
#'                text_size = 12)
#' }
#'
#' @export
#'
#' @importFrom ampvis2 amp_heatmap
#' @importFrom ggplot2 theme element_text
#' @importFrom assertthat assert_that is.string
#'
#' @seealso
#' \code{\link[ampvis2]{amp_heatmap}}
#'
genus_heatmap <- function(
  ampvisObject,
  group_by,
  tax_show = 20,
  plot_values = TRUE,
  value_size = 8,
  decimals = 2,
  text_size = 14,
  axis_angle = 45,
  color_palette = NULL,
  order_x_by = NULL
) {
  require(assertthat)

  # Input validation
  assert_that(
    inherits(ampvisObject, "ampvis2"),
    msg = "Input must be an ampvis2 object"
  )

  assert_that(
    is.string(group_by),
    msg = "group_by must be a character string"
  )

  assert_that(
    group_by %in% colnames(ampvisObject$metadata),
    msg = sprintf("'%s' not found in metadata", group_by)
  )

  assert_that(
    is.numeric(tax_show),
    tax_show > 0,
    msg = "tax_show must be a positive number"
  )

  # Create heatmap
  heatmap <- amp_heatmap(
    ampvisObject,
    group_by = group_by,
    tax_aggregate = "Genus",
    tax_add = "Phylum",
    tax_show = tax_show,
    plot_values = plot_values,
    normalise = TRUE,
    plot_values_size = value_size,
    round = decimals,
    tax_empty = "OTU",
    order_x_by = order_x_by,
    order_y_by = "cluster"
  )

  # Add custom theme
  heatmap <- heatmap +
    theme(
      axis.text.x = element_text(
        angle = axis_angle,
        size = text_size,
        vjust = 1,
        hjust = 1
      ),
      axis.text.y = element_text(
        size = text_size
      ),
      legend.position = "right",
      legend.title = element_text(
        color = "black",
        size = text_size * 0.7,
        face = "bold"
      ),
      # Additional theme elements for better visualization
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title = element_text(size = text_size, face = "bold")
    )

  # Apply custom color palette if provided
  if (!is.null(color_palette)) {
    heatmap <- heatmap +
      scale_fill_gradientn(colours = color_palette)
  }

  # Add attributes for Quarto
  attr(heatmap, "caption") <- sprintf(
    "Genus-level heatmap grouped by %s",
    group_by
  )

  class(heatmap) <- c("genus_heatmap", class(heatmap))

  return(heatmap)
}

#' Create a Species-Level Heatmap from ampvis2 Object
#'
#' @description
#' Generates a customized heatmap visualization of microbial community data at the species level
#' using the ampvis2 package. The function includes enhanced formatting options and input validation.
#'
#' @param ampvisObject An ampvis2 object containing microbial community data
#' @param group_by Character. The grouping variable name from the metadata
#' @param tax_show Numeric. Number of taxa to show (default: 20)
#' @param plot_values Logical. Whether to show values in cells (default: TRUE)
#' @param value_size Numeric. Size of value labels (default: 8)
#' @param decimals Integer. Number of decimal places for values (default: 3)
#' @param text_size Numeric. Size of axis text (default: 14)
#' @param axis_angle Numeric. Angle of x-axis labels (default: 45)
#' @param color_palette Character vector. Custom color palette (optional)
#'
#' @return A ggplot2 object containing the heatmap visualization
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' species_heatmap(amp_data, "Treatment")
#'
#' # With custom parameters
#' species_heatmap(amp_data, "Treatment",
#'                tax_show = 15,
#'                plot_values = FALSE,
#'                text_size = 12)
#' }
#'
#' @export
#'
#' @importFrom ampvis2 amp_heatmap
#' @importFrom ggplot2 theme element_text
#' @importFrom assertthat assert_that is.string
#'
#' @seealso
#' \code{\link[ampvis2]{amp_heatmap}}
#'
species_heatmap <- function(
  ampvisObject,
  group_by,
  tax_show = 20,
  plot_values = TRUE,
  value_size = 8,
  decimals = 2,
  text_size = 14,
  axis_angle = 45,
  color_palette = NULL,
  order_x_by = NULL
) {
  require(assertthat)

  # Input validation
  assert_that(
    inherits(ampvisObject, "ampvis2"),
    msg = "Input must be an ampvis2 object"
  )

  assert_that(
    is.string(group_by),
    msg = "group_by must be a character string"
  )

  assert_that(
    group_by %in% colnames(ampvisObject$metadata),
    msg = sprintf("'%s' not found in metadata", group_by)
  )

  assert_that(
    is.numeric(tax_show),
    tax_show > 0,
    msg = "tax_show must be a positive number"
  )

  # Create heatmap
  heatmap <- amp_heatmap(
    ampvisObject,
    group_by = group_by,
    tax_aggregate = "Species",
    tax_add = "Genus",
    tax_show = tax_show,
    plot_values = plot_values,
    normalise = TRUE,
    plot_values_size = value_size,
    round = decimals,
    tax_empty = "OTU",
    order_x_by = order_x_by,
    order_y_by = "cluster"
  )

  # Add custom theme
  heatmap <- heatmap +
    theme(
      axis.text.x = element_text(
        angle = axis_angle,
        size = text_size,
        vjust = 1,
        hjust = 1
      ),
      axis.text.y = element_text(
        size = text_size
      ),
      legend.position = "right",
      legend.title = element_text(
        color = "black",
        size = text_size * 0.7,
        face = "bold"
      ),
      # Additional theme elements for better visualization
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title = element_text(size = text_size, face = "bold")
    )

  # Apply custom color palette if provided
  if (!is.null(color_palette)) {
    heatmap <- heatmap +
      scale_fill_gradientn(colours = color_palette)
  }

  # Add attributes for Quarto
  attr(heatmap, "caption") <- sprintf(
    "Species-level heatmap grouped by %s",
    group_by
  )

  class(heatmap) <- c("species_heatmap", class(heatmap))

  return(heatmap)
}

#' Calculate Omega-Squared Effect Size for PERMANOVA (adonis) Results
#'
#' @description
#' Calculates either partial or non-partial omega-squared effect sizes from the results
#' of a PERMANOVA analysis (from \code{vegan::adonis()} or \code{vegan::adonis2()}).
#' Omega-squared provides a measure of variance explained that is less biased than R-squared,
#' particularly for small sample sizes.
#'
#' @param adonisOutput An object of class "adonis" or "anova.cca" from vegan's adonis/adonis2 function
#' @param partial logical; if TRUE (default) calculates partial omega-squared,
#'                if FALSE calculates non-partial omega-squared
#'
#' @return Returns a modified version of the input object with an additional column for
#'         omega-squared values. For adonis objects, returns the full object with modified
#'         aov.tab. For anova.cca objects, returns just the modified ANOVA table.
#'
#' @details
#' The function calculates omega-squared (ω²) or partial omega-squared (ω²p) effect sizes
#' following the formulas:
#'
#' Partial ω² = (df_effect * (MS_effect - MS_residual)) /
#'              (df_effect * MS_effect + (N - df_effect) * MS_residual)
#'
#' Non-partial ω² = (SS_effect - df_effect * MS_residual) / (SS_total + MS_residual)
#'
#' Where:
#' - df_effect: degrees of freedom for the effect
#' - MS_effect: Mean Square for the effect
#' - MS_residual: Mean Square for the residual
#' - SS_effect: Sum of Squares for the effect
#' - SS_total: Total Sum of Squares
#' - N: Total sample size
#'
#' @note
#' - Omega-squared values will be NA for residual and total rows
#' - Effect sizes can be negative when effect size is very small
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' data(dune)
#' data(dune.env)
#'
#' # Run PERMANOVA
#' adonis_result <- adonis2(dune ~ Management + A1, data = dune.env)
#'
#' # Calculate partial omega-squared
#' result_with_omega <- adonis_OmegaSq(adonis_result)
#'
#' # Calculate non-partial omega-squared
#' result_with_omega <- adonis_OmegaSq(adonis_result, partial = FALSE)
#' }
#'
#' @export
#' @importFrom methods is
#'
#' @author Original code by [original author], modified by [your name]
#' @references
#' Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared statistics:
#' measures of effect size for some common research designs.
#' Psychological Methods, 8(4), 434-447.
adonis_OmegaSq <- function(adonisOutput, partial = TRUE) {
  # Input validation
  if (!requireNamespace("methods", quietly = TRUE)) {
    stop("Package 'methods' is required but not installed")
  }

  if (
    !(methods::is(adonisOutput, "adonis") ||
      methods::is(adonisOutput, "anova.cca"))
  ) {
    stop("Input must be an object of class 'adonis' or 'anova.cca'")
  }

  if (!is.logical(partial)) {
    stop("'partial' must be logical (TRUE/FALSE)")
  }

  if (methods::is(adonisOutput, "anova.cca")) {
    aov_tab <- adonisOutput
    aov_tab$MeanSqs <- aov_tab$SumOfSqs / aov_tab$Df
    aov_tab$MeanSqs[length(aov_tab$Df)] <- NA
  } else {
    aov_tab <- adonisOutput$aov.tab
  }

  heading <- attr(aov_tab, "heading")

  MS_res <- aov_tab[pmatch("Residual", rownames(aov_tab)), "MeanSqs"]
  if (is.na(MS_res)) {
    stop("Cannot find 'Residual' row in ANOVA table")
  }

  SS_tot <- aov_tab[rownames(aov_tab) == "Total", "SumsOfSqs"]

  N <- aov_tab[rownames(aov_tab) == "Total", "Df"] + 1

  if (partial) {
    omega <- apply(aov_tab, 1, function(x) {
      (x["Df"] * (x["MeanSqs"] - MS_res)) /
        (x["Df"] * x["MeanSqs"] + (N - x["Df"]) * MS_res)
    })
    aov_tab$parOmegaSq <- c(omega[1:(length(omega) - 2)], NA, NA)
  } else {
    omega <- apply(aov_tab, 1, function(x) {
      (x["SumsOfSqs"] - x["Df"] * MS_res) / (SS_tot + MS_res)
    })
    aov_tab$OmegaSq <- c(omega[1:(length(omega) - 2)], NA, NA)
  }
  if (is(adonisOutput, "adonis")) {
    cn_order <- c(
      "Df",
      "SumsOfSqs",
      "MeanSqs",
      "F.Model",
      "R2",
      if (partial) "parOmegaSq" else "OmegaSq",
      "Pr(>F)"
    )
  } else {
    cn_order <- c(
      "Df",
      "SumOfSqs",
      "F",
      if (partial) "parOmegaSq" else "OmegaSq",
      "Pr(>F)"
    )
  }
  aov_tab <- aov_tab[, cn_order]
  attr(aov_tab, "names") <- cn_order
  attr(aov_tab, "heading") <- heading
  if (is(adonisOutput, "adonis")) {
    adonisOutput$aov.tab <- aov_tab
  } else {
    adonisOutput <- aov_tab
  }
  return(adonisOutput)
}

gt_strata <- function(
  data,
  strata,
  by,
  continuous_stat = "{median}  \n({IQR}; {min}-{max})",
  categorical_stat = "{n} / {N} ({p}%)",
  digits = 2,
  missing = "no",
  p.adjust = "fdr",
  vda = TRUE
) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (missing(by)) {
    stop("'by' argument must be specified")
  }

  if (!missing %in% c("no", "ifany", "always")) {
    stop("'missing' must be one of: 'no', 'ifany', 'always'")
  }

  # Check for required packages
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' needed. Please install it.")
  } else {
    require(gt)
  }

  if (!requireNamespace("gtsummary", quietly = TRUE)) {
    stop("Package 'gtsummary' needed. Please install it.")
  } else {
    require(gtsummary)
  }

  if (vda == FALSE) {
    # Create summary table
    summary_tbl <- tbl_strata(
      data = data,
      strata = strata,
      .combine_with = c("tbl_stack"),
      .tbl_fun = ~ .x %>%
        tbl_summary(
          by = by,
          missing = missing,
          statistic = list(
            all_continuous() ~ continuous_stat,
            all_categorical() ~ categorical_stat
          )
        ) %>%
        add_p(pvalue_fun = pvalue_format) %>%
        add_q(pvalue_fun = pvalue_format, method = p.adjust) %>%
        bold_p(q = TRUE),
      .header = "**{strata}**, N = {n}"
    )
  } else if (vda == TRUE) {
    summary_tbl <-
      tbl_strata(
        data = data,
        strata = strata,
        .combine_with = c("tbl_stack"),
        .tbl_fun = ~ .x %>%
          tbl_summary(
            by = by,
            digits = everything() ~ 2,
            missing = "no",
            statistic = list(
              all_continuous() ~ "{median}  \n({IQR}; {min}-{max})", # stats and format for continuous columns
              all_categorical() ~ "{n} / {N}  \n({p}%)"
            )
          ) %>% # stats and format for categorical columns
          add_n() %>%
          add_overall() %>%
          add_p(pvalue_fun = function(x) style_pvalue(x, digits = 3)) %>%
          add_q(pvalue_fun = function(x) style_pvalue(x, digits = 3)) %>%
          bold_p(q = TRUE) %>%
          add_stat(fns = all_continuous() ~ vda) %>%
          modify_fmt_fun(list(vda_magnitude ~ as.character)) %>%
          modify_header(list(
            vda_estimate ~ "**Effect Size Estimate**",
            vda_magnitude ~ "**Magnitude**"
          )),
        .header = "**{strata}**, N = {n}"
      )
  } else {
    stop("`vda` must be set to either `TRUE` or `FALSE`")
  }

  # Convert to gt and apply formatting
  formatted_tbl <- summary_tbl %>%
    as_gt() %>%
    fmt_markdown(columns = where(~ is.character(.x))) %>%
    fmt_number(
      columns = where(~ is.numeric(.x)),
      drop_trailing_zeros = TRUE,
      use_seps = FALSE,
      decimals = digits
    ) %>%
    fmt_number(
      columns = where(~ is.numeric(.x) < 0.0001),
      pattern = "<0.0001"
    ) %>%
    fmt_number(
      columns = where(~ is.numeric(.x) < 0.001),
      pattern = "<0.001"
    ) %>%
    fmt_number(
      columns = where(~ is.numeric(.x) < 0.01),
      pattern = "<0.01"
    ) %>%
    fmt_number(
      columns = where(~ is.numeric(.x) < 0.05),
      pattern = "<0.05"
    )

  return(formatted_tbl)
}

#' Create Stratified Summary Tables with Formatting
#'
#' @description
#' Creates formatted stratified summary tables using gtsummary and gt packages.
#' The function supports both continuous and categorical variables, with options
#' for visual display of agreement (VDA) statistics and customizable formatting.
#'
#' @param data A data frame containing the variables to be summarized
#' @param strata Character. Name of the stratification variable
#' @param by Character. Name of the grouping variable for comparisons
#' @param continuous_stat Character. Format string for continuous variable statistics.
#'        Default: "{median} ({IQR}; {min}-{max})"
#' @param categorical_stat Character. Format string for categorical variable statistics.
#'        Default: "{n} / {N} ({p}%)"
#' @param digits Integer. Number of decimal places for numeric formatting. Default: 2
#' @param missing Character. How to handle missing values: "no", "ifany", or "always".
#'        Default: "no"
#' @param p.adjust Character. Method for p-value adjustment. Must be one of the methods
#'        supported by stats::p.adjust. Default: "fdr"
#' @param vda Logical. Whether to include visual display of agreement statistics.
#'        Default: TRUE
#'
#' @return A gt table object with formatted statistics and p-values
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' gt_strata(data = my_data, strata = "group", by = "treatment")
#'
#' # Custom formatting without VDA
#' gt_strata(
#'   data = my_data,
#'   strata = "group",
#'   by = "treatment",
#'   continuous_stat = "{mean} ({sd})",
#'   vda = FALSE
#' )
#' }
#'
#' @import gt
#' @import gtsummary
#' @importFrom dplyr where
#'
#' @export
gt_strata <- function(
  data,
  strata,
  by,
  continuous_stat = "{median} ({IQR}; {min}-{max})",
  categorical_stat = "{n} / {N} ({p}%)",
  digits = 2,
  missing = "no",
  p.adjust = "fdr",
  vda = TRUE
) {
  # Enhanced input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (missing(strata) || !is.character(strata) || !strata %in% names(data)) {
    stop("'strata' must be a valid column name in the data frame")
  }

  if (missing(by) || !is.character(by) || !by %in% names(data)) {
    stop("'by' must be a valid column name in the data frame")
  }

  if (!is.numeric(digits) || digits < 0) {
    stop("'digits' must be a non-negative number")
  }

  if (!missing %in% c("no", "ifany", "always")) {
    stop("'missing' must be one of: 'no', 'ifany', 'always'")
  }

  if (!is.logical(vda)) {
    stop("'vda' must be TRUE or FALSE")
  }

  # Validate p.adjust method
  valid_p_adjust <- c(
    "holm",
    "hochberg",
    "hommel",
    "bonferroni",
    "BH",
    "BY",
    "fdr",
    "none"
  )
  if (!p.adjust %in% valid_p_adjust) {
    stop(sprintf(
      "'p.adjust' must be one of: %s",
      paste(valid_p_adjust, collapse = ", ")
    ))
  }

  # Check for required packages with informative messages
  required_packages <- c("gt", "gtsummary")
  missing_packages <- required_packages[
    !sapply(required_packages, requireNamespace, quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    stop(sprintf(
      "The following required packages are missing: %s\nPlease install them using install.packages()",
      paste(missing_packages, collapse = ", ")
    ))
  }

  # Create the summary table based on VDA setting
  summary_tbl <- if (!vda) {
    create_basic_summary(
      data,
      strata,
      by,
      continuous_stat,
      categorical_stat,
      digits,
      missing,
      p.adjust
    )
  } else {
    create_vda_summary(data, strata, by, digits, missing, p.adjust)
  }

  # Apply consistent formatting
  formatted_tbl <- format_gt_table(summary_tbl, digits)

  # Add metadata for Quarto compatibility
  attr(formatted_tbl, "quarto") <- list(
    format = "html",
    dependencies = c("gt", "gtsummary"),
    options = list(
      digits = digits,
      missing = missing,
      p.adjust = p.adjust,
      vda = vda
    )
  )

  return(formatted_tbl)
}

#' Create Basic Summary Table
#' @noRd
create_basic_summary <- function(
  data,
  strata,
  by,
  continuous_stat,
  categorical_stat,
  digits,
  missing,
  p.adjust
) {
  tbl_strata(
    data = data,
    strata = strata,
    .combine_with = c("tbl_stack"),
    .tbl_fun = ~ .x %>%
      tbl_summary(
        by = by,
        missing = missing,
        statistic = list(
          all_continuous() ~ continuous_stat,
          all_categorical() ~ categorical_stat
        )
      ) %>%
      add_p(pvalue_fun = pvalue_format) %>%
      add_q(pvalue_fun = pvalue_format, method = p.adjust) %>%
      bold_p(q = TRUE),
    .header = "**{strata}**, N = {n}"
  )
}

#' Create VDA Summary Table
#' @noRd
create_vda_summary <- function(data, strata, by, digits, missing, p.adjust) {
  tbl_strata(
    data = data,
    strata = strata,
    .combine_with = c("tbl_stack"),
    .tbl_fun = ~ .x %>%
      tbl_summary(
        by = by,
        digits = everything() ~ digits,
        missing = missing,
        statistic = list(
          all_continuous() ~ "{median}\n({IQR}; {min}-{max})",
          all_categorical() ~ "{n} / {N}\n({p}%)"
        )
      ) %>%
      add_n() %>%
      add_overall() %>%
      add_stat(fns = all_continuous() ~ vda) %>%
      modify_fmt_fun(list(vda_magnitude ~ as.character)) %>%
      modify_header(list(
        vda_estimate ~ "**Effect Size Estimate**",
        vda_magnitude ~ "**Magnitude**"
      )) %>%
      add_p(pvalue_fun = function(x) style_pvalue(x, digits = 3)) %>%
      add_q(pvalue_fun = function(x) style_pvalue(x, digits = 3)) %>%
      bold_p(q = TRUE),
    .header = "**{strata}**, N = {n}"
  )
}

#' Format GT Table
#' @noRd
format_gt_table <- function(tbl, digits) {
  tbl %>%
    as_gt() %>%
    fmt_markdown(columns = where(~ is.character(.x))) %>%
    fmt_number(
      columns = where(~ is.numeric(.x)),
      drop_trailing_zeros = TRUE,
      use_seps = FALSE,
      decimals = digits
    ) %>%
    fmt_number(
      columns = where(~ is.numeric(.x) < 0.0001),
      pattern = "<0.0001"
    ) %>%
    fmt_number(
      columns = where(~ is.numeric(.x) < 0.001),
      pattern = "<0.001"
    ) %>%
    fmt_number(
      columns = where(~ is.numeric(.x) < 0.01),
      pattern = "<0.01"
    ) %>%
    fmt_number(
      columns = where(~ is.numeric(.x) < 0.05),
      pattern = "<0.05"
    )
}


#end
##############################################################################
# #start kable Settings                                                             #
##############################################################################

# kable <- function(data, ...){
#   if(!interactive()) {
#   knitr::kable(data) %>%
#   kableExtra::kable_classic(c("striped", "hover", "condensed"), full_width = F, html_font = "Inter")
#   }
#   else {
#      knitr::kable(data)
#   }

# }

kable <- function(data, ...) {
  knitr::kable(data, table.attr = "style = \"color: #C6D0F5;\"")
}

#end

##############################################################################
#  #start Visualisation                                                   #
##############################################################################

#extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

phyla_colours <- c(
  "Actinobacteriota" = '#a6cee3',
  "Bacteroidota" = '#b2df8a',
  "Campylobacterota" = '#33a02c',
  "Bacillota" = '#Fb9a99',
  "Other" = '#E31a1c',
  "Pseudomonadota" = '#Fdbf6f',
  "Unclassified" = '#Ff7f00',
  "Unknown" = '#ffff99',
  "Verrucomicrobiota" = '#Cab2d6'
)

pairedColours <- colorRampPalette(brewer.pal(12, "Paired")) #colour-blind friendly
set2 <- c(
  '#a6cee3',
  '#1f78b4',
  '#b2df8a',
  '#33a02c',
  '#fb9a99',
  '#e31a1c',
  '#fdbf6f',
  '#ff7f00',
  '#cab2d6',
  '#6a3d9a',
  '#ffff99',
  '#b15928'
)

trt_labels <- c('0' = "CPM", '1' = "FeZnPM")
time_labels <- c(
  "1Baseline" = "Baseline",
  "2Midpoint" = "Midpoint",
  "3Endpoint" = "Endpoint"
)

armColours <- c("#FC8D62", "#8DA0CB")
timepoint_colours <- c("#66C2A5", "#E78AC3")

#for use with ggplot tile or ggcorrplot
ggheatmap_palette <- c('#91bfdb', '#ffffbf', '#fc8d59')

okabe_ito <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999" # grey
)
# Function to assign column classes based on a data dictionary
# To view colour-blind friendly palettes run: display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, colorblindFriendly=TRUE)
assign_unique_colors <- function(
  char_vector,
  color_vector = NULL,
  palette = "Set2"
) {
  if (!require("RColorBrewer", quietly = TRUE)) {
    install.packages("RColorBrewer")
  }
  library(RColorBrewer)

  unique_levels <- unique(char_vector)

  # If color_vector is not provided, create an empty one
  if (is.null(color_vector)) {
    color_vector <- vector("character", length = length(unique_levels))
    #names(color_vector) <- unique_levels
  }

  # Find levels without assigned colors
  unassigned_levels <- setdiff(unique_levels, names(color_vector))

  # Generate unique colors for unassigned levels
  if (length(unassigned_levels) < 9) {
    new_colors <- brewer.pal(n = length(unassigned_levels), name = palette)
  } else {
    palette <- colorRampPalette(c(
      '#a6cee3',
      '#1f78b4',
      '#b2df8a',
      '#33a02c',
      '#fb9a99',
      '#e31a1c',
      '#fdbf6f',
      '#ff7f00',
      '#cab2d6',
      '#6a3d9a',
      '#ffff99',
      '#b15928'
    ))
    new_colors <- palette(n = length(unassigned_levels))
  }
  names(new_colors) <- unassigned_levels

  #  Update the existing color vector with new colors for previously unseen levels
  for (level in unassigned_levels) {
    color_vector[level] <- new_colors[level]
  }

  # remove any NAs
  color_vector <- color_vector[!is.na(names(color_vector))]

  # Return the updated color vector
  return(color_vector[unique(char_vector)])
}

patterns <- function(params, boundary_df, aspect_ratio, legend = FALSE) {
  args <- as.list(params)
  pattern_type <- as.factor(as.numeric(params$pattern_type))

  # Define colorblind-friendly colors
  colors <- c(
    '#a6cee3',
    '#1f78b4',
    '#b2df8a',
    '#33a02c',
    '#fb9a99',
    '#e31a1c',
    '#fdbf6f',
    '#ff7f00',
    '#cab2d6',
    '#6a3d9a',
    '#ffff99',
    '#b15928'
  ) #brewer.pal(12, "Paired")

  # Base parameters for the pattern
  base_args <- list(
    x = boundary_df$x,
    y = boundary_df$y,
    id = boundary_df$id,
    prefix = ""
  )

  # Define different patterns for different levels
  if (pattern_type <= 12) {
    # Solid colors for first 12 levels using stripe pattern with 100% density
    args <- c(
      base_args,
      list(
        pattern = "stripe",
        pattern_fill = colors[pattern_type],
        pattern_colour = colors[pattern_type], # Same as fill color
        pattern_spacing = 0.05,
        pattern_density = 1, # Full density makes it appear solid
        pattern_angle = 0
      )
    )
  } else if (pattern_type >= 13) {
    # Striped patterns for remaining levels
    args <- c(
      base_args,
      list(
        pattern = "stripe",
        pattern_fill = colors[pattern_type],
        pattern_colour = colors[round((pattern_type - 12) / 2)],
        pattern_spacing = 0.05,
        pattern_density = 0.7,
      )
    )
  } else {
    # Dot patterns for remaining levels
    args <- c(
      base_args,
      list(
        pattern = "circle",
        pattern_fill = colors[pattern_type],
        pattern_colour = colors[round((pattern_type - 12) / 2)],
        pattern_spacing = 0.05,
        pattern_density = 0.7,
      )
    )
  }

  do.call(gridpattern::patternGrob, args)
}

# Register the pattern
options(ggpattern_geometry_funcs = list(patterns = patterns))

#' Create Correlation Panel with Color-Coded Values
#'
#' @description
#' Creates a correlation panel for use in pairs plot, displaying correlation
#' coefficients with color-coded backgrounds. The color intensity represents
#' the strength and direction of correlation.
#'
#' @param x Numeric vector or matrix for correlation calculation
#' @param y Numeric vector or matrix for correlation calculation
#' @param digits Integer. Number of decimal places to display (default: 2)
#' @param prefix Character string to prepend to correlation value (default: "")
#' @param text_size Numeric. Size of correlation text (default: 3)
#' @param color_palette Character vector of length 3 specifying low, mid, and high colors
#' @param na.rm Logical. Whether to remove NA values (default: TRUE)
#' @param method Character. Correlation method: "pearson", "spearman", or "kendall" (default: "pearson")
#' @param ... Additional arguments passed to text()
#'
#' @return Invisibly returns the correlation coefficient
#'
#' @examples
#' \dontrun{
#' pairs(mtcars,
#'       upper.panel = panel.cor,
#'       lower.panel = panel.smooth,
#'       gap = 0,
#'       row1attop = TRUE)
#' }
#'
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom stats cor
panel.cor <- function(
  x,
  y,
  digits = 2,
  prefix = "",
  text_size = 3,
  color_palette = c('#91bfdb', '#ffffbf', '#fc8d59'),
  na.rm = TRUE,
  method = "pearson",
  ...
) {
  # Input validation
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors or matrices")
  }
  if (!is.numeric(digits) || digits < 0) {
    stop("digits must be a non-negative integer")
  }
  if (!is.character(prefix)) {
    stop("prefix must be a character string")
  }
  if (!method %in% c("pearson", "spearman", "kendall")) {
    stop("method must be one of: 'pearson', 'spearman', 'kendall'")
  }
  if (length(color_palette) != 3) {
    stop("color_palette must be a vector of three colors")
  }

  # Calculate correlation with error handling
  r <- tryCatch(
    {
      cor(
        x,
        y,
        use = if (na.rm) "complete.obs" else "everything",
        method = method
      )
    },
    error = function(e) {
      warning("Could not compute correlation: ", e$message)
      return(NA)
    }
  )

  # Format text
  txt <- if (!is.na(r)) {
    paste0(prefix, format(round(r, digits), nsmall = digits))
  } else {
    "NA"
  }

  # Create color scale
  col <- colorRampPalette(color_palette)(100)

  # Map correlation to color index with NA handling
  col.ind <- if (!is.na(r)) {
    round((r + 1) * 49 + 1)
  } else {
    50 # middle color for NA
  }

  # Ensure color index is within bounds
  col.ind <- pmax(1, pmin(100, col.ind))

  # Fill panel with color using usr coordinates
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = col[col.ind])

  # Add text with appropriate color
  text_col <- if (is.na(r) || abs(r) <= 0.5) "black" else "white"
  text(mean(usr[1:2]), mean(usr[3:4]), txt, col = text_col, cex = text_size)

  # Return correlation coefficient invisibly
  invisible(r)
}

#end
################################################################################
# Plot limma-voom results
################################################################################
plot_limmavoom = function(tt, physeq, alpha = 0.1) {
  require("ggplot2")
  require("phyloseq")
  # Rownames of tt are the OTU IDs
  xdf = cbind(tt, tax_table(physeq)[rownames(tt), ])
  xdf$OTU <- rownames(xdf)
  # The points to highlight
  specialdf = xdf[(xdf$P.Value < alpha), ]
  # The plot
  p = ggplot(
    data = xdf,
    mapping = aes(x = logFC, y = -log10(P.Value), size = -log10(P.Value))
  ) +
    geom_vline(xintercept = 0.0, linetype = 2) +
    # Background points
    geom_point(color = "black", alpha = 0.65) +
    scale_size_continuous(range = c(1, 4)) +
    guides(size = FALSE, colour = guide_legend(nrow = 2, byrow = TRUE)) +
    # axis labels
    labs(y = expression(-log[10](p)), x = expression(log[2](FC)))
  # The interesting OTUs
  if (nrow(specialdf) > 0) {
    p <- p +
      geom_point(data = specialdf, mapping = aes(color = Genus)) +
      geom_text(
        data = specialdf,
        mapping = aes(label = paste("OTU:", OTU)),
        size = 2,
        nudge_x = 0.3
      )
  }
  return(p)
}

##############################################################################
#  #start Functions                                                                #
##############################################################################

#' Calculate Summary Statistics for Numeric Data
#'
#' @description Calculates summary statistics (mean, standard deviation, standard error,
#' median, maximum, and minimum) for numeric data. Works with vectors, matrices, or
#' data frames, with optional grouping variables.
#'
#' @param data A numeric vector, matrix, or data frame
#' @param varname Character string specifying the variable name to summarize.
#' If NULL and data is a vector, uses the entire vector.
#' @param groupnames Character vector of column names to group by.
#' Optional for data frames, ignored for vectors and matrices.
#'
#' @return A data frame containing summary statistics:
#' \itemize{
#'   \item mean
#'   \item sd (standard deviation)
#'   \item SEM (standard error of the mean)
#'   \item median
#'   \item max (maximum value)
#'   \item min (minimum value)
#' }
#'
#' @examples
#' # Vector example
#' vec <- rnorm(100)
#' data_summary(vec)
#'
#' # Data frame example with grouping
#' df <- data.frame(
#'   value = rnorm(100),
#'   group = rep(c("A", "B"), each = 50)
#' )
#' data_summary(df, "value", "group")
#'
#' # Matrix example
#' mat <- matrix(rnorm(100), ncol = 2)
#' data_summary(mat)
#'
#' @importFrom stats sd median
#' @importFrom plyr ddply rename
#' @export
data_summary <- function(data, varname = NULL, groupnames = NULL) {
  # Input validation and conversion
  if (is.vector(data)) {
    if (!is.numeric(data)) {
      stop("Vector must be numeric")
    }
    data <- data.frame(value = data)
    varname <- "value"
  } else if (is.matrix(data)) {
    if (!is.numeric(data)) {
      stop("Matrix must be numeric")
    }
    data <- as.data.frame(data)
    if (is.null(varname)) {
      varname <- colnames(data)[1]
    }
  } else if (!is.data.frame(data)) {
    stop("Input must be a vector, matrix, or data frame")
  }

  # Check if varname exists in data frame
  if (!is.null(varname) && !(varname %in% colnames(data))) {
    stop("Variable name not found in data")
  }

  # Check if groupnames exist in data frame
  if (!is.null(groupnames) && !all(groupnames %in% colnames(data))) {
    stop("One or more grouping variables not found in data")
  }

  # Summary function
  summary_func <- function(x, col) {
    if (all(is.na(x[[col]]))) {
      return(c(
        n = NA,
        mean = NA,
        sd = NA,
        SEM = NA,
        median = NA,
        max = NA,
        min = NA
      ))
    }
    c(
      n = length(x[[col]]),
      mean = mean(x[[col]], na.rm = TRUE),
      sd = sd(x[[col]], na.rm = TRUE),
      SEM = sd(x[[col]], na.rm = TRUE) /
        sqrt(sum(!is.na(x[[col]]))),
      median = median(x[[col]], na.rm = TRUE),
      max = max(x[[col]], na.rm = TRUE),
      min = min(x[[col]], na.rm = TRUE)
    )
  }

  # Calculate summaries
  if (is.null(groupnames)) {
    # No grouping - calculate overall summary
    data_sum <- as.data.frame(t(summary_func(data, varname)))
  } else {
    # With grouping
    data_sum <- ddply(data, groupnames, .fun = summary_func, varname)
  }

  # Rename mean column if varname provided
  if (!is.null(varname)) {
    data_sum <- rename(data_sum, c("mean" = varname))
  }

  return(data_sum)
}

#' Generate Multiple Random Vectors from Various Distributions
#'
#' @description Creates multiple random vectors using specified distribution functions.
#' Can generate vectors of equal or varying lengths using different random distributions.
#'
#' @param n_vec Numeric. Number of vectors to generate
#' @param n_obs Numeric vector or single value. If single value, all vectors will have
#' the same length. If vector, must have length equal to n_vec
#' @param dist Character or function. Distribution to use. Can be "norm" (default), "unif",
#' "exp", "pois", or a custom random generation function
#' @param params List of parameters for the distribution. For normal distribution:
#' list(mean = 0, sd = 1), for uniform: list(min = 0, max = 1), etc.
#' @param names Character vector of names for the vectors. If NULL, default names
#' will be "vec1", "vec2", etc.
#'
#' @return A list of numeric vectors, each containing random values from the specified distribution
#'
#' @examples
#' # Create 3 normal vectors
#' mkvec(n_vec = 3, n_obs = 10, dist = "norm")
#'
#' # Create 3 uniform vectors
#' mkvec(n_vec = 3, n_obs = 10,
#'                      dist = "unif",
#'                      params = list(min = 0, max = 10))
#'
#' # Create 2 exponential vectors with different rates
#' mkvec(n_vec = 2,
#'                      n_obs = c(10, 20),
#'                      dist = "exp",
#'                      params = list(rate = c(1, 2)))
#'
#' # Using a custom distribution function
#' custom_dist <- function(n) rbeta(n, 2, 2)
#' mkvec(n_vec = 2, n_obs = 10, dist = custom_dist)
#'
#' @export
mkvec <- function(
  n_vec,
  n_obs,
  dist = "norm",
  params = list(mean = 0, sd = 1)
) {
  # Input validation
  if (!is.numeric(n_vec) || n_vec < 1) {
    stop("n_vec must be a positive number")
  }

  # Handle n_obs input - now with recycling
  if (length(n_obs) == 1) {
    n_obs <- rep(n_obs, n_vec)
  } else {
    # Recycle n_obs values if needed
    n_obs <- rep_len(n_obs, n_vec)
  }

  # Generate default names
  names <- paste0("vec", 1:n_vec)

  # Handle custom function with multiple parameters
  if (is.function(dist)) {
    # Create list of vectors
    result <- vector("list", n_vec)

    # Ensure all parameters are properly replicated
    params <- lapply(params, function(p) {
      if (length(p) == 1) {
        rep(p, n_vec)
      } else {
        rep_len(p, n_vec) # Recycle parameter values if needed
      }
    })

    for (i in seq_len(n_vec)) {
      # Extract parameters for this iteration
      current_params <- lapply(params, function(p) p[i])

      # Add n parameter
      current_params$n <- n_obs[i]

      # Call custom function with parameters
      tryCatch(
        {
          result[[i]] <- do.call(dist, current_params)
          # Convert to matrix to vector if needed
          if (is.matrix(result[[i]])) {
            result[[i]] <- as.vector(result[[i]])
          }
        },
        error = function(e) {
          stop("Error in custom function for vector ", i, ": ", e$message)
        }
      )
    }
  } else if (is.character(dist)) {
    # Handle built-in distributions (code for built-in distributions...)
    stop("Built-in distributions not yet implemented")
  } else {
    stop("dist must be either a character string or a function")
  }

  # Assign names to the vectors in the list
  names(result) <- names

  return(result)
}

#end
##############################################################################
#  #start LMMs Tools for LMMs                                                       #
##############################################################################

#generate a dataframe with the common goodness of fit metrics from an nlme derived model
model_stats <- function(nlme_model, vif = TRUE, format = "markdown") {
  if (class(nlme_model) != "lme") {
    message("Model must be of class 'lme'.")
    stop()
  }

  x <- summary(nlme_model)

  mod_df <- data.frame(
    row.names = "Value",
    Bhat = x$coefficients$fixed[1],
    intercept_variance = as.numeric(nlme::VarCorr(x)[, 1])[1], #tau
    residual_variance = as.numeric(nlme::VarCorr(x)[, 1])[2], #sigma
    total_varaince = as.numeric(nlme::VarCorr(x)[, 1])[1] +
      as.numeric(nlme::VarCorr(x)[, 1])[2],
    icc = as.numeric(nlme::VarCorr(x)[, 1])[1] /
      (as.numeric(nlme::VarCorr(x)[, 1])[1] +
        as.numeric(nlme::VarCorr(x)[, 1])[2]), #intercalss correlation coefficient
    BIC = x$BIC,
    AIC = x$AIC,
    loglik = x$logLik
  )

  if (vif == TRUE) {
    VIF <- car::vif(nlme_model) #check for multicollineratiy issues
    VIF_Intercept <- NA
    VIF <- c(Intercept = VIF_Intercept, VIF)

    tTable <- as.data.frame(x$tTable)
    tTable$VIF <- VIF
  } else {
    VIF <- car::vif(nlme_model) #check for multicollineratiy issues
    VIF$Intercept <- NA

    tTable <- as.data.frame(x$tTable)
    cbind(tTable, VIF)
  }

  results <- bind_rows(tTable, as.data.frame(t(mod_df)))
  results <- knitr::kable(results, format = format)
  return(results)
}

#plot residuals to test assumption of normal distribution
# ADD Breusch-Pagan test or similar for heteroskedacity test
# use render = FALSE for interactive use
residual_plots <- function(
  model,
  model_name = names(model),
  render = TRUE,
  chol = TRUE
) {
  set.seed(1234)
  library(ggplot2)
  library(dplyr)
  library(mgcv)

  if (chol == TRUE) {
    rawRes <- residuals(model, type = "pearson")
    estCov <- extract.lme.cov(model, model$data)

    # Perform Cholesky decomposition on the residuals' covariance matrix
    message(
      "Cholesky decomposition is being performed on residuals from ",
      deparse(substitute(model))
    )
    residuals <- solve(t(chol(estCov))) %*% rawRes
  } else {
    residuals <- residuals(model, type = "pearson")
  }

  xy_plot <- ggplot(
    data.frame(x = predict(model), y = residuals),
    aes(x = x, y = y)
  ) +
    geom_point(shape = 1, size = 2, color = "black") +
    #geom_smooth(method = "lm") +
    #ggpubr::stat_cor(position = "dodge") +
    geom_smooth(method = "loess") +
    geom_hline(yintercept = 0, color = "red") +
    xlab("Fitted Values") +
    ylab("Model Residuals") +
    labs(title = paste("Residuals vs Fitted Values")) +
    theme(plot.title = element_text(size = 10))

  qqplot <- ggplot(data.frame(resid = residuals), aes(sample = resid)) +
    stat_qq(shape = 1, size = 2, color = "black") +
    stat_qq_line(color = "red") +
    labs(title = paste("QQ Plot")) +
    ylab(NULL) +
    xlab(NULL) +
    theme(plot.title = element_text(size = 10))

  #create histogram and overlay normal curve
  hist_data <- hist(residuals, plot = FALSE)
  data <- data.frame(x = residuals)

  hist <- ggplot(data, aes(x = x)) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = hist_data$breaks[2] - hist_data$breaks[1],
      fill = "lightblue",
      color = "black",
      boundary = 0.5
    ) +
    stat_function(
      fun = dnorm,
      args = list(mean = mean(data$x), sd = sd(data$x)),
      color = "red",
      size = 1
    ) +
    ylab("Frequency") +
    xlab("Model residuals") +
    labs(title = paste("Residuals Histogram")) +
    theme(plot.title = element_text(size = 10))

  norm_test1 <- capture.output(shapiro.test(residuals))
  norm <- ggplot(data.frame(x = 0:1, y = 0:1), aes(x = x, y = y)) +
    geom_blank() +
    xlab(NULL) +
    ylab(NULL) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) +
    annotate(
      geom = "text",
      x = 0.75,
      y = 0.20,
      label = paste(
        "Summary Output:\n",
        paste(norm_test1, collapse = "\n"),
        "P-value <0.05",
        "\n => NOT normally distributed"
      ),
      hjust = 1,
      vjust = 0,
      size = 4
    )

  # Create a list of plots
  if (render == FALSE) {
    X11(display = "", title = model_name)
  }

  plots <- gridExtra::grid.arrange(
    xy_plot,
    qqplot,
    hist,
    norm,
    ncol = 2,
    nrow = 2
  )

  return(plots)
}

residual_norm_test <- function(model_list) {
  norm_test <- lapply(names(model_list), function(model_name) {
    model <- model_list[[model_name]]
    st <- shapiro.test(residuals(model))
    data.frame(Model = model_name, W = st$statistic, P_value = st$p.value)
  })
  norm_df <- do.call(rbind, norm_test)
  return(norm_df)
}

vda <- function(data, variable, by, ...) {
  # Compute VD.A
  A <- effsize::VD.A(as.formula(glue::glue("`{variable}` ~ {by}")), data = data)
  tibble(vda_estimate = A$estimate, vda_magnitude = A$magnitude)
}

sturges_bins <- function(data) {
  n <- length(data)
  k <- ceiling(log2(n) + 1)
  return(k)
}

#' Test for Normality of Data Frame Columns
#'
#' @description Performs Shapiro-Wilk normality test on a specified numeric column
#' of a data frame and returns formatted results suitable for Quarto/RMarkdown documents.
#'
#' @param df A data frame containing the data to test
#' @param j Column name or index to test for normality
#' @param format Output format for the table. One of "auto", "pipe", "simple",
#' "html", "latex", or "rst" (default: "auto")
#' @param digits Number of decimal places for p-value (default: 4)
#' @param alpha Significance level for normality test (default: 0.05)
#'
#' @return A formatted kable object containing test results. The table includes:
#' \itemize{
#'   \item Name: Name of the tested variable
#'   \item W: Shapiro-Wilk test statistic
#'   \item P_value: P-value of the test
#' }
#'
#' @details
#' The function performs a Shapiro-Wilk test for normality on the specified column.
#' Results are formatted using kable() for clean output in RMarkdown/Quarto documents.
#' The null hypothesis is that the data is normally distributed.
#'
#' @note
#' The Shapiro-Wilk test is most appropriate for sample sizes between 3 and 5000.
#'
#' @examples
#' # Create example data
#' df <- data.frame(
#'   normal = rnorm(100),
#'   uniform = runif(100),
#'   categorical = factor(rep(1:4, 25))
#' )
#'
#' # Test normally distributed data
#' norm_test(df, "normal")
#'
#' # Test non-normal data
#' norm_test(df, "uniform")
#'
#' # Test with different format
#' norm_test(df, "normal", format = "html")
#'
#' @importFrom stats shapiro.test
#' @importFrom knitr kable
#' @export
norm_test <- function(df, j, format, digits = 4, alpha = 0.05) {
  # Input validation
  if (!is.data.frame(df)) {
    stop("First argument must be a data frame")
  }

  if (is.character(j) && !(j %in% names(df))) {
    stop("Column '", j, "' not found in data frame")
  }

  if (is.numeric(j) && (j < 1 || j > ncol(df))) {
    stop("Column index out of bounds")
  }

  # Extract column
  y <- df[[j]]
  col_name <- if (is.character(j)) j else names(df)[j]

  # Check if numeric
  if (!is.numeric(y)) {
    warning("Column '", col_name, "' is not numeric. Skipping test.")
    return(NULL)
  }

  # Check sample size
  if (length(y) < 3 || length(y) > 5000) {
    warning(
      "Sample size (",
      length(y),
      ") is outside recommended range (3-5000) for Shapiro-Wilk test"
    )
  }

  # Remove NA values
  y <- na.omit(y)
  if (length(y) == 0) {
    warning("No non-missing values in column '", col_name, "'")
    return(NULL)
  }

  # Perform test
  norm_test <- shapiro.test(y)

  # Create results data frame
  norm_stats <- data.frame(
    "Variable" = col_name,
    "n" = length(y),
    "W" = norm_test$statistic,
    "P_value" = round(norm_test$p.value, digits),
    "Distribution" = ifelse(
      norm_test$p.value > alpha,
      "Likely Normal",
      "Likely Non-normal"
    ),
    check.names = FALSE
  )

  # Create caption
  caption <- sprintf(
    "Shapiro-Wilk Normality Test Results (α = %.2f)",
    alpha
  )

  # Format table
  table <- knitr::kable(
    norm_stats,
    format = format,
    caption = caption,
    align = c('l', 'r', 'r', 'r', 'l'),
    booktabs = TRUE,
    digits = digits
  )

  # Add styling if format is latex and kableExtra is available
  if (format == "latex") {
    if (requireNamespace("kableExtra", quietly = TRUE)) {
      table <- kableExtra::kable_styling(
        table,
        bootstrap_options = c("striped", "hover"),
        full_width = FALSE
      )
    }
  }

  return(table)
}

#' Print method for norm_test output
#' @param x Output from norm_test function
#' @param ... Additional arguments passed to print
#' @export
print.norm_test <- function(x, ...) {
  print(x, ...)
}

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  round(p.mat, 4)
}

#' Summarize Variables for Linear Mixed Model Terms
#'
#' @description
#' Creates a summary table of variables intended for use in linear mixed models,
#' including basic characteristics and summary statistics for both numeric and
#' categorical variables.
#'
#' @param variables A list of variables to be summarized
#' @param variable_names A character vector of names corresponding to the variables
#'
#' @return A kable-formatted table containing:
#' \itemize{
#'   \item Variable: Name of the variable
#'   \item Type: Class of the variable
#'   \item N: Number of observations
#'   \item Missing: Number of missing values
#'   \item Unique values: Number of unique values
#'   \item Levels: Number of levels (for factors)
#'   \item Summary Stats: For numeric variables: mean, SD, and range.
#'                       For non-numeric variables: "NA"
#' }
#'
#' @details
#' The function processes both numeric and categorical variables differently:
#' * Numeric variables include summary statistics (mean, SD, range)
#' * Categorical variables show number of levels but no summary statistics
#' All results are formatted using knitr::kable for nice printing
#'
#' @examples
#' variables <- list(
#'   age = c(25, 30, 35),
#'   group = factor(c("A", "B", "A")),
#'   time = c(1, 2, 3)
#' )
#' var_names <- c("age", "group", "time")
#' summarise_terms(variables, var_names)
#'
#' @aliases summarize_terms
#' @export
#'
#' @import knitr
#'
#' @seealso
#' \code{\link[knitr]{kable}} for the table formatting
#'
#' @note
#' Ensures consistent formatting across variable types and handles missing values
summarise_terms <- function(variables, variable_names) {
  # Input validation
  if (length(variables) != length(variable_names)) {
    stop("Number of variables must match number of variable names")
  }
  if (!is.list(variables)) {
    stop("variables must be a list")
  }
  if (!is.character(variable_names)) {
    stop("variable_names must be a character vector")
  }

  # Initialize empty data frame
  summaries <- data.frame()

  for (i in seq_along(variables)) {
    var <- variables[[i]]
    name <- variable_names[i]

    # Basic info for all variables
    summary <- data.frame(
      Variable = name,
      Type = class(var)[1],
      N = length(var),
      Missing = sum(is.na(var)),
      `Unique values` = length(unique(var)),
      Levels = length(levels(var))
    )

    # Add type-specific summaries
    if (is.numeric(var)) {
      summary$`Summary Stats` <- paste0(
        "Mean=",
        round(mean(var, na.rm = TRUE), 2),
        ", SD=",
        round(sd(var, na.rm = TRUE), 2),
        ", Range=[",
        round(min(var, na.rm = TRUE), 2),
        "-",
        round(max(var, na.rm = TRUE), 2),
        "]"
      )
    } else {
      summary$`Summary Stats` <- "NA"
    }

    # Bind to main results
    summaries <- rbind(summaries, summary)
  }

  return(knitr::kable(summaries, digits = 2))
}

anova_pvalue <- function(x) {
  UseMethod("anova_pvalue", x)
}
anova_pvalue.default <- function(x) {
  terms <- rownames(x)
  # Extract the p-value column from the data frame
  pvals <- x$`p-value`
  # Combine the p-values and terms into a temporary data frame
  data.frame(Terms = terms, Model = model, P.value = pvals)
}
anova_pvalue.list <- function(x) {
  # Initialize an empty data frame to store results
  combined_results <- data.frame()

  pvalues <- lapply(seq_along(x), function(i) {
    # Get p-values and terms for each model in the current list
    terms <- rownames(x[[i]])
    # Extract the p-value column from the data frame
    pvals <- x[[i]]$`p-value`
    model <- names(x)[i]
    # Combine the p-values and terms into a temporary data frame
    data.frame(Terms = terms, Model = model, P.value = pvals)
  })
  # Bind the temporary data frame to the combined results data frame
  combined_results <- dplyr::bind_rows(combined_results, pvalues)
  rownames(combined_results) <- NULL

  # Return the combined results data frame
  return(combined_results)
}
anova_pvalue.listoflists <- function(x) {
  # Initialize an empty data frame to store results
  combined_results <- data.frame()

  for (l in x) {
    # Get p-values and terms for each model in the current list
    pvalues <- lapply(seq_along(l), function(i) {
      terms <- rownames(l[[i]])
      # Extract the p-value column from the data frame
      pvals <- l[[i]]$`p-value`
      model <- names(l)[i]
      # Combine the p-values and terms into a temporary data frame
      data.frame(Terms = terms, Model = model, P.value = pvals)
    })
    # Bind the temporary data frame to the combined results data frame
    combined_results <- dplyr::bind_rows(combined_results, pvalues)
    rownames(combined_results) <- NULL
  }

  # Return the combined results data frame
  return(combined_results)
}

het_breuschpagan_lme <- function(nlme_model, exog_het, robust = TRUE) {
  # Extract residuals and predicted values from the LME model
  resid <- residuals(nlme_model, type = "pearson")
  y <- resid^2
  if (!robust) {
    y <- y / mean(y)
  }

  # Convert exog_het to a matrix
  exog_het <- as.matrix(exog_het)

  # Fit a linear model using OLS
  lm_fit <- lm(y ~ exog_het - 1) # -1 removes intercept

  # Calculate test statistics
  lm <- if (robust) {
    nobs(lm_fit) * summary(lm_fit)$r.squared
  } else {
    summary(lm_fit)$sigma / 2
  }

  # Calculate p-value
  lm_pval <- 1 - pchisq(lm, length(exog_het) - 1)

  # Return the result as a list
  result <- list(
    lm = lm,
    lm_pvalue = lm_pval,
    fvalue = summary(lm_fit)$fstatistic[1],
    f_pvalue = summary(lm_fit)$fstatistic[2]
  )

  return(result)
}

# Function to check and apply transformation
check_and_transform <- function(var) {
  # Example: Check skewness and apply log-transformation if needed
  if (shapiro.test(var)$p.value < 0.05) {
    var <- log(var)
    cat("Log-transformation applied to", response_var, "\n")
  }

  # Return the modified data
  return(data)
}

#end

pad_columns <- function(df) {
  # Determine the maximum length among all columns
  max_length <- max(sapply(df, length))

  # Pad each column with NA values
  padded_df <- lapply(df, function(x) {
    if (length(x) < max_length) {
      c(x, rep(NA, max_length - length(x)))
    } else {
      x
    }
  })

  # Convert the result back to a data frame
  padded_df <- as.data.frame(padded_df)

  return(padded_df)
}

#determine the coefficient threshold for a particular fold change of interest
#TODO change function name so not to conflict with stats package coef function
coef <- function(desired_fold_change, log_scale = 2) {
  log(desired_fold_change, base = log_scale)
}

##############################################################################
#  #start Wrangling functions for data manipulation                         #
##############################################################################

loadcombine <- function(
  file_dir,
  file_pattern,
  header = TRUE,
  sep,
  skip = 0,
  output
) {
  # Usage

  #loadcombine(file_dir = "C:/Users/Admin/OneDrive - University of Aberdeen (1)/Experiments/Sequencing/Transcriptomics/Nate-transcriptome_analysis/results/exploratory_MoA/GO/REVIGO_tiny/", file_pattern = "*.tsv", sep = "\t", output = "C:/Users/Admin/OneDrive - University of Aberdeen (1)/Experiments/Sequencing/Transcriptomics/Nate-transcriptome_analysis/results/exploratory_MoA/GO/REVIGO_tiny/REVIGO_results_combined.csv", header = T)

  currentDir <- getwd()
  setwd(file_dir)

  #create var with files paths as a character vector
  all_file_paths <- list.files(
    path = file_dir,
    pattern = file_pattern,
    full.names = TRUE
  )

  #load all the files into the environment
  all_files <- lapply(
    all_file_paths,
    read.table,
    header = header,
    sep = sep,
    skip = skip
  )

  all_filenames <-
    all_file_paths %>%
    basename() %>% #basename removes the path from the file name
    as.list() #turns them into a list

  #combines the file content list and file name list i.e. matches the names and the files
  all_lists <- mapply(c, all_files, all_filenames, SIMPLIFY = FALSE)

  all_result <- rbindlist(all_lists, fill = TRUE)
  names(all_result)[length(all_result)] <- "sample"
  #BUG all_result <- all_result[ ,colSums(is.na(all_result))<=ncol(all_result)]
  #export results file
  write.csv(all_result, file = output, row.names = FALSE)
  setwd(currentDir)
}

save_plot <- function(plot, filename) {
  svg(filename = paste0(filename, ".svg"))
  print(plot)
  dev.off()
  png(filename = paste0(filename, ".png"))
  print(plot)
  dev.off()
}

# Function to assign column classes based on a data dictionary
assign_col_classes <- function(
  data,
  data_dict_path,
  variable_name = "Variable.Name",
  classCol = "colClass.R"
) {
  # Read the data dictionary
  data_dict <- read.csv(data_dict_path, stringsAsFactors = FALSE)

  # Read the first row of the data dictionary to get the variable names
  variables <- data_dict[[variable_name]]

  # Read the first row of the data dictionary to get the classes
  classes <- data_dict[[classCol]]

  if (is.object(data) == FALSE) {
    data <- read.csv(data, header = TRUE)
  }

  # Match variable names to column names in the data dictionary
  col_classes <- sapply(colnames(data), function(var_name) {
    class_index <- match(var_name, variables)
    if (!is.na(class_index)) {
      return(classes[class_index])
    }
  })

  # Loop through columns and assign classes
  for (i in seq_along(col_classes)) {
    tryCatch(
      {
        if (is.null(col_classes[[i]])) {
          data[[i]] <- type.convert(data[[i]], as.is = TRUE)
        } else {
          if (col_classes[[i]] == "factor") {
            data[[i]] <- as.factor(data[[i]])
          } else {
            class(data[[i]]) <- col_classes[[i]]
          }
        }
      },
      error = function(e) {
        cat("Error assigning class to column:", names(data)[i], "\n")
        cat("Error message:", conditionMessage(e), "\n")
      }
    )
  }

  return(data)
}
# Function to assign column classes based on a data dictionary
assign_col_classes <- function(
  data,
  data_dict_path,
  variable_name = "Variable.Name",
  classCol = "colClass.R"
) {
  # Read the data dictionary
  data_dict <- read.csv(data_dict_path, stringsAsFactors = FALSE)

  # Read the first row of the data dictionary to get the variable names
  variables <- data_dict[[variable_name]]

  # Read the first row of the data dictionary to get the classes
  classes <- data_dict[[classCol]]

  if (is.object(data) == FALSE) {
    data <- read.csv(data, header = TRUE)
  }

  # Match variable names to column names in the data dictionary
  col_classes <- sapply(colnames(data), function(var_name) {
    class_index <- match(var_name, variables)
    if (!is.na(class_index)) {
      return(classes[class_index])
    }
  })

  # Loop through columns and assign classes
  for (i in seq_along(col_classes)) {
    tryCatch(
      {
        if (is.null(col_classes[[i]])) {
          data[[i]] <- type.convert(data[[i]], as.is = TRUE)
        } else {
          if (col_classes[[i]] == "factor") {
            data[[i]] <- as.factor(data[[i]])
          } else {
            class(data[[i]]) <- col_classes[[i]]
          }
        }
      },
      error = function(e) {
        cat("Error assigning class to column:", names(data)[i], "\n")
        cat("Error message:", conditionMessage(e), "\n")
      }
    )
  }

  return(data)
}

#' @title Utilities For \code{\link{phyloseq-class}} Slots to Data Frames

#' @description Utility to convert phyloseq slots to data frames.

#' @details These functions are a fork of \code{tibble_utilities} from the [\code{microbiome} package](https://github.com/microbiome/microbiome/blob/master/R/tibble_utilities.R). Instead of converting to tibbles, these functions instead convert to data frame format, thus retaining rownames instead of adding the rowname as a new column.
#'
#' Convert different \code{phyloseq} slots into data frames.
#' \code{otu_tibble} gets the otu_table in data frame format.
#' \code{tax_tibble} gets the taxa_table in data frame format.
#' \code{combine_otu_tax} combines otu_table and taxa_table into one data frame.
#' \code{sample_df} is effectively the same as \code{microbime::meta()}, and extracts the sample metadata as a data.frame.

#' @param x \code{\link{phyloseq-class}} object.
#' @param column.id Provide name for the column which will hold the rownames.
#'                  of slot.
#' @return A \code{data frame}
#' @examples
#' library(promicrobialTools)
#' data("dietswap")
#' otu_tib <- otu_tibble(dietswap,column.id="FeatureID")
#' tax_tib <- tax_tibble(dietswap,column.id="FeatureID")
#' sample_tib <- sample_tibble(dietswap,column.id="SampleID")
#' otu_tax <- combine_otu_tax(dietswap,column.id = "FeatureID")
#' head(otu_tax)
#'
#' @name DfUtilites
#' @author Contact: Nathaniel Cole \email{nc564@cornell.edu}
#' @seealso To convert OTU tables to abundacne matrices see \code{microbiome::abundances}. To extract
NULL

#' @rdname DfUtilites
#' @aliases ogu_df
#' @export
otu_df <- function(x) {
  if (any(c("phyloseq", "otu_table") %in% is(x))) {
    # Pick OTU matrix
    otu <- microbiome::abundances(x)
  }
  otu_df <- otu %>%
    as.data.frame()
  return(otu_df)
}


#' @rdname DfUtilites
#' @aliases tax_df
#' @export
tax_df <- function(x, column.id = "FeatureID") {
  if (any(c("phyloseq", "tax_table") %in% is(x))) {
    # Pick OTU matrix
    tax <- phyloseq::tax_table(x)
    tax_df <- tax %>%
      as.matrix() %>%
      as.data.frame()
    return(tax_df)
  } else {
    return(message("tax_table not found"))
  }
}

#' @rdname DfUtilites
#' @aliases sample_df
#' @export
sample_df <- function(x, column.id = "SampleID") {
  if (any(c("phyloseq", "sample_data") %in% is(x))) {
    # Pick OTU matrix
    smd <- microbiome::meta(x)
    return(sample_tibble)
  }
}

#' @rdname TibbleUtilites
#' @aliases combine_otu_tax
#' @export
combine_otu_tax <- function(x, column.id = "FeatureID") {
  otu_tb <- tax_tb <- NULL
  otu_tb <- otu_df(x, column.id)
  tax_tb <- tax_df(x, column.id)
  otu_tb <- tax_tb %>%
    dplyr::left_join(otu_tb, by = column.id)
  return(otu_tb)
}

#start check_metadata function

#' Check Metadata Variable Properties
#'
#' @description
#' Checks and reports the class, number of NA values, and length of a metadata variable.
#' This function is designed to validate variables used in microbiome analyses,
#' ensuring they meet expected specifications.
#'
#' @param var The variable to check (vector)
#' @param expected_class Character string specifying the expected class of the variable
#' @param name Character string specifying the name of the variable for reporting
#'
#' @return Prints a summary of the variable's properties:
#'   \itemize{
#'     \item Class (actual and expected)
#'     \item Number of NA values
#'     \item Length of variable
#'   }
#'
#' @examples
#' subject <- c("A1", "A2", "A3")
#' check_metadata(subject, "character", "subject")
#'
#' age <- c(25, 30, NA, 40)
#' check_metadata(age, "numeric", "age")
#'
#' group <- factor(c("Control", "Treatment"))
#' check_metadata(group, "factor", "group")
#'
#' @seealso
#' \code{\link{class}}, \code{\link{is.na}}
#'
#' @export
#'
#' @note
#' This function is particularly useful for checking metadata variables
#' before conducting microbiome analyses where complete cases and
#' correct variable types are crucial.
#'
#' @author Nathaniel Cole
#' @version 1.0.0
check_metadata <- function(var) {
  cat("Checking:", deparse(substitute(var)), "\n")
  cat("Class:", class(var), "\n")
  cat("NAs:", sum(is.na(var)), "\n")
  cat("Length:", length(var), "\n")
  cat("-------------------\n")
}

#end

#' Create a Formatted Table from Multiple Results Objects
#'
#' @description
#' Combines multiple results objects into a single formatted table.
#' Compatible with emmeans, data.frames, matrices, and other objects that can be
#' coerced to data.frames. Provides flexible formatting options and is compatible
#' with knitr/quarto rendering.
#'
#' @param results_list List of objects to combine. Names of list elements used as identifiers.
#' @param id_col Character. Name for the identifier column (default: "Model")
#' @param digits Integer. Number of decimal places for numeric columns (default: 2)
#' @param format Character. Output format: "markdown", "html", or "latex" (default: "markdown")
#' @param caption Character. Optional table caption
#' @param col.names Character vector. Custom column names
#' @param round_cols Character vector. Columns to round (default: all numeric)
#' @param align Character vector. Column alignment (default: auto-detect)
#' @param table_style Character vector. kableExtra styling options
#'
#' @return A knitr_kable object that can be directly used in R Markdown or Quarto
#'
#' @examples
#' \dontrun{
#' # With emmeans objects
#' results1 <- list(
#'   "Model 1" = emmeans(model1, ~treatment),
#'   "Model 2" = emmeans(model2, ~treatment)
#' )
#'
#' # With data frames
#' results2 <- list(
#'   "Analysis 1" = data.frame(x = 1:3, y = 4:6),
#'   "Analysis 2" = data.frame(x = 7:9, y = 10:12)
#' )
#'
#' # Create tables
#' table1 <- table_bind(
#'   results1,
#'   caption = "EMM Comparison",
#'   table_style = c("striped", "hover")
#' )
#'
#' table2 <- table_bind(
#'   results2,
#'   id_col = "Analysis",
#'   align = c('l', 'r', 'r')
#' )
#' }
#'
#' @export
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
table_bind <- function(
  results_list,
  id_col = "Model",
  digits = 2,
  format = c("markdown", "html", "latex"),
  caption = NULL,
  col.names = NULL,
  round_cols = NULL,
  align = NULL,
  table_style = NULL
) {
  # Input validation
  if (!is.list(results_list)) {
    stop("Input must be a list of results objects")
  }

  if (length(results_list) == 0) {
    stop("Empty list provided")
  }

  # Match format argument
  format <- match.arg(format)

  # Convert each object to data frame
  df_list <- lapply(results_list, function(x) {
    tryCatch(
      {
        # Handle different object types
        if (inherits(x, "emmGrid")) {
          as.data.frame(x)
        } else if (inherits(x, "matrix")) {
          as.data.frame(x)
        } else if (is.data.frame(x)) {
          x
        } else {
          as.data.frame(x)
        }
      },
      error = function(e) {
        stop("Error converting object to data frame: ", e$message)
      }
    )
  })

  # Check for consistent columns
  col_names <- lapply(df_list, colnames)
  if (length(unique(lapply(col_names, length))) > 1) {
    stop("Inconsistent number of columns in objects")
  }

  # Create results data frame
  results <- data.frame(
    ID = rep(names(results_list), sapply(df_list, nrow)),
    do.call(rbind, df_list),
    check.names = FALSE
  )
  names(results)[1] <- id_col

  # Round numeric columns
  if (is.null(round_cols)) {
    numeric_cols <- sapply(results, is.numeric)
    results[, numeric_cols] <- round(results[, numeric_cols], digits)
  } else {
    results[, round_cols] <- round(results[, round_cols], digits)
  }

  # Use custom column names if provided
  if (!is.null(col.names)) {
    if (length(col.names) != ncol(results)) {
      stop("Length of col.names must match number of columns")
    }
    colnames(results) <- col.names
  }

  # Auto-detect alignment if not specified
  if (is.null(align)) {
    align <- ifelse(sapply(results, is.numeric), 'r', 'l')
  }

  # Create kable
  table <- knitr::kable(
    results,
    format = format,
    digits = digits,
    caption = caption,
    align = align,
    booktabs = TRUE
  )

  # Apply kableExtra styling if requested
  if (!is.null(table_style) && requireNamespace("kableExtra", quietly = TRUE)) {
    table <- kableExtra::kable_styling(table, bootstrap_options = table_style)
  }

  # Add class for method dispatch
  class(table) <- c("bound_table", class(table))

  return(table)
}

# Print method
#' @export
print.bound_table <- function(x, ...) {
  print(unclass(x))
}

# Knit print method
#' @export
knit_print.bound_table <- function(x, ...) {
  knitr::knit_print(unclass(x), ...)
}
#end

#start EMMs functions for EMMs
#################
#     EMMs      #
#################

#' Bind Multiple Estimated Marginal Means into a Table
#'
#' @description
#' Combines multiple emmeans objects into a single formatted table, with model names
#' as identifiers. The output is compatible with knitr/quarto rendering.
#'
#' @param emm List of emmeans objects. Names of the list elements are used as model identifiers.
#' @param digits Integer. Number of decimal places for numeric columns (default: 2)
#' @param format Character. Output format: "markdown", "html", or "latex" (default: "markdown")
#' @param caption Character. Optional table caption
#' @param col.names Character vector. Custom column names
#'
#' @return A knitr_kable object that can be directly used in R Markdown or Quarto
#'
#' @examples
#' \dontrun{
#' # Create multiple emmeans objects
#' emm1 <- emmeans(model1, ~treatment)
#' emm2 <- emmeans(model2, ~treatment)
#'
#' # Combine into table
#' emm_table <- emm_bind(
#'   list("Model 1" = emm1, "Model 2" = emm2),
#'   digits = 3,
#'   caption = "Estimated Marginal Means Comparison"
#' )
#'
#' # Print table
#' emm_table
#' }
#'
#' @export
#' @importFrom knitr kable
#' @importFrom emmeans emmeans
# emm_bind <- function(
#   emm,
#   digits = 2,
#   format = "markdown",
#   caption = NULL,
#   col.names = NULL
# ) {
#   # Input validation
#   if (!is.list(emm)) {
#     stop("'emm' must be a list of emmeans objects")
#   }

#   if (length(emm) == 0) {
#     stop("Empty list provided")
#   }

#   if (!all(sapply(emm, inherits, "emmGrid"))) {
#     stop("All elements must be emmeans objects (class 'emmGrid')")
#   }

#   # Match format argument
#   format <- match.arg(format)

#   # Convert each emmeans object to data frame
#   df_list <- lapply(emm, function(x) {
#     tryCatch(
#       {
#         as.data.frame(x)
#       },
#       error = function(e) {
#         stop("Error converting emmeans to data frame: ", e$message)
#       }
#     )
#   })

#   # Check for consistent columns
#   col_names <- lapply(df_list, colnames)
#   if (length(unique(lapply(col_names, length))) > 1) {
#     stop("Inconsistent number of columns in emmeans objects")
#   }

#   # Create results data frame
#   results <- data.frame(
#     Model = rep(names(df_list), sapply(df_list, nrow)),
#     do.call(rbind, df_list),
#     check.names = FALSE
#   )

#   # Round numeric columns to specified digits
#   numeric_cols <- sapply(results, is.numeric)
#   results[numeric_cols] <- round(results[numeric_cols], digits)

#   # Use custom column names if provided
#   if (!is.null(col.names)) {
#     if (length(col.names) != ncol(results)) {
#       stop("Length of col.names must match number of columns")
#     }
#     colnames(results) <- col.names
#   }

#   # Create kable with appropriate format
#   table <- knitr::kable(
#     results,
#     format = format,
#     digits = digits,
#     caption = caption,
#     booktabs = TRUE
#   )

#   # Add class for potential method dispatch
#   class(table) <- c("emm_table", class(table))

#   return(table)
# }

#end

calculate_pairwise_vda <- function(data, value_col, group_col) {
  # Get unique groups
  groups <- unique(data[[group_col]])
  n_groups <- length(groups)

  # Create empty results dataframe
  results <- data.frame(
    group1 = character(),
    group2 = character(),
    VDA = numeric(),
    effect_size = character(),
    stringsAsFactors = FALSE
  )

  # Calculate VDA for each pair
  for (i in 1:(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      group1_data <- data[data[[group_col]] == groups[i], ][[value_col]]
      group2_data <- data[data[[group_col]] == groups[j], ][[value_col]]

      # Calculate VDA
      score <- 0
      for (x1 in group1_data) {
        for (x2 in group2_data) {
          if (x1 > x2) {
            score <- score + 1
          }
          if (x1 == x2) score <- score + 0.5
        }
      }

      vda <- score / (length(group1_data) * length(group2_data))

      # Determine effect size
      effect <- case_when(
        vda >= 0.71 ~ "large",
        vda >= 0.64 ~ "medium",
        vda >= 0.56 ~ "small",
        vda <= 0.29 ~ "large negative",
        vda <= 0.36 ~ "medium negative",
        vda <= 0.44 ~ "small negative",
        TRUE ~ "negligible"
      )

      # Add to results
      results <- rbind(
        results,
        data.frame(
          group1 = groups[i],
          group2 = groups[j],
          VDA = vda,
          effect_size = effect
        )
      )
    }
  }

  return(results)
}

#' Create Tertiles from a Continuous Variable
#'
#' @description
#' Divides a continuous variable into three groups (tertiles) while ensuring that
#' equal values are not split across different tertiles. Returns both the categorized
#' data and a summary of the tertile distribution.
#'
#' @param continuous_variable A numeric vector to be divided into tertiles
#'
#' @return A list containing:
#'   \itemize{
#'     \item tertiles: An ordered factor with levels "Low", "Medium", "High"
#'     \item summary: A data frame with columns:
#'       \itemize{
#'         \item Tertile: Character indicating tertile level
#'         \item N: Integer count of observations
#'         \item Percent: Numeric percentage of observations
#'       }
#'     \item breaks: Numeric vector of break points (excluding -Inf and Inf)
#'   }
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' result <- tertiles(x)
#' print(result$summary)
#' print(result$breaks)
#'
#' @export
#'
#' @note
#' The function ensures that equal values are not split across tertiles,
#' which may result in groups that deviate from exact thirds.
#'
#' @seealso
#' \code{\link{cut}}, \code{\link{quantile}}
#'
#' @importFrom stats quantile
#' @importFrom knitr kable
#'
tertiles <- function(continuous_variable) {
  # Input validation
  if (!is.numeric(continuous_variable)) {
    stop("Input must be a numeric vector")
  }

  if (length(continuous_variable) < 3) {
    stop("Input vector must have at least 3 values to create tertiles")
  }

  # Remove NA values
  x <- continuous_variable[!is.na(continuous_variable)]

  if (length(x) == 0) {
    stop("No non-NA values in input vector")
  }

  # Get unique values and sort them
  unique_values <- sort(unique(x))

  # Calculate cumulative percentages for unique values
  cum_pct <- cumsum(table(x)) / length(x)

  # Initialize tertile breaks
  breaks <- c(-Inf, Inf) # Default breaks if logic fails

  # Find the closest values to 33.33% and 66.67% that don't split equal values
  if (length(unique_values) > 2) {
    # Find break points that don't split equal values
    lower_third <- unique_values[which.min(abs(cum_pct - 1 / 3))]
    upper_third <- unique_values[which.min(abs(cum_pct - 2 / 3))]

    breaks <- c(-Inf, lower_third, upper_third, Inf)
  }

  # Create tertiles
  tertiles <- cut(
    continuous_variable,
    breaks = breaks,
    labels = c("Low", "Medium", "High"),
    include.lowest = TRUE,
    ordered = TRUE
  ) # Make it an ordered factor

  # Create summary statistics
  summary_stats <- data.frame(
    Tertile = c("Low", "Medium", "High"),
    N = as.numeric(table(tertiles)),
    Percent = round(as.numeric(prop.table(table(tertiles)) * 100), 1)
  )

  # Print summary
  print("Tertile Distribution:")
  print(summary_stats)
  print("Break Points:")
  print(breaks[-c(1, length(breaks))]) # Don't show -Inf and Inf

  return(tertiles)
}

#' Process Heatmap Data with Metadata
#'
#' @description
#' Processes heatmap data by combining it with metadata and calculating summary statistics
#'
#' @param hm A list containing heatmap data with a $data element
#' @param i The name of the grouping variable
#' @param metadata A data frame containing metadata including sample_names and grouping variables
#'
#' @return A data frame with summary statistics grouped by phyla and the specified grouping variable
#'
#' @examples
#' # Example usage:
#' # hm_data(hm = hmaps$b12_tertiles, i = "b12_tertiles", metadata = b12_base_meta)
#'
#' @export
hm_data <- function(hm, i, metadata) {
  # Input validation
  if (!is.list(hm) || is.null(hm$data)) {
    stop("'hm' must be a list containing a 'data' element")
  }
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data frame")
  }
  if (!all(c("sample_names", i) %in% colnames(metadata))) {
    stop(
      "metadata must contain 'sample_names' and the grouping variable '",
      i,
      "'"
    )
  }

  # Process data
  df <- hm$data %>%
    dplyr::select(phyla = Display, Sample, relAbund = Sum) %>%
    # Join with metadata using proper column names
    left_join(
      metadata %>%
        dplyr::select(Sample = sample_names, grouping_var = !!sym(i)),
      by = "Sample"
    ) %>%
    # Calculate summary statistics
    group_by(phyla, grouping_var) %>%
    summarise(
      mean = mean(relAbund, na.rm = TRUE),
      median = median(relAbund, na.rm = TRUE),
      max = max(relAbund, na.rm = TRUE),
      min = min(relAbund, na.rm = TRUE),
      n = n(),
      `Median (IQR; Min-Max)` = paste0(
        round(median(relAbund, na.rm = TRUE), 3),
        "\n(",
        round(IQR(relAbund), 3),
        "; ",
        round(min(relAbund, na.rm = TRUE), 3),
        "-",
        round(max(relAbund, na.rm = TRUE), 3),
        ")"
      ),
      .groups = "drop"
    )

  return(df)
}

#' Create a Phylum-Level Heatmap from ampvis2 Object
#'
#' @description
#' Generates a customized heatmap visualization of microbial community data at the phylum level
#' using the ampvis2 package. The function includes enhanced formatting options and input validation.
#'
#' @param ampvisObject An ampvis2 object containing microbial community data
#' @param group_by Character. The grouping variable name from the metadata
#' @param tax_show Numeric. Number of taxa to show (default: 20)
#' @param plot_values Logical. Whether to show values in cells (default: TRUE)
#' @param value_size Numeric. Size of value labels (default: 8)
#' @param decimals Integer. Number of decimal places for values (default: 3)
#' @param text_size Numeric. Size of axis text (default: 14)
#' @param axis_angle Numeric. Angle of x-axis labels (default: 45)
#' @param color_palette Character vector. Custom color palette (optional)
#'
#' @return A ggplot2 object containing the heatmap visualization
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' phylum_heatmap(amp_data, "Treatment")
#'
#' # With custom parameters
#' phylum_heatmap(amp_data, "Treatment",
#'                tax_show = 15,
#'                plot_values = FALSE,
#'                text_size = 12)
#' }
#'
#' @export
#'
#' @importFrom ampvis2 amp_heatmap
#' @importFrom ggplot2 theme element_text
#' @importFrom assertthat assert_that is.string
#'
#' @seealso
#' \code{\link[ampvis2]{amp_heatmap}}
#'
phylum_heatmap <- function(
  ampvisObject,
  group_by,
  tax_show = 20,
  plot_values = TRUE,
  value_size = 8,
  decimals = 2,
  text_size = 14,
  axis_angle = 45,
  color_palette = NULL,
  order_x_by = NULL
) {
  require(assertthat)

  # Input validation
  assert_that(
    inherits(ampvisObject, "ampvis2"),
    msg = "Input must be an ampvis2 object"
  )

  assert_that(
    is.string(group_by),
    msg = "group_by must be a character string"
  )

  assert_that(
    group_by %in% colnames(ampvisObject$metadata),
    msg = sprintf("'%s' not found in metadata", group_by)
  )

  assert_that(
    is.numeric(tax_show),
    tax_show > 0,
    msg = "tax_show must be a positive number"
  )

  # Create heatmap
  heatmap <- amp_heatmap(
    ampvisObject,
    group_by = group_by,
    tax_aggregate = "Phylum",
    tax_show = tax_show,
    plot_values = plot_values,
    normalise = TRUE,
    plot_values_size = value_size,
    round = decimals,
    tax_empty = "OTU",
    order_x_by = order_x_by,
    order_y_by = "cluster"
  )

  # Add custom theme
  heatmap <- heatmap +
    theme(
      axis.text.x = element_text(
        angle = axis_angle,
        size = text_size,
        vjust = 1,
        hjust = 1
      ),
      axis.text.y = element_text(
        size = text_size
      ),
      legend.position = "right",
      legend.title = element_text(
        color = "black",
        size = text_size * 0.7,
        face = "bold"
      ),
      # Additional theme elements for better visualization
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title = element_text(size = text_size, face = "bold")
    )

  # Apply custom color palette if provided
  if (!is.null(color_palette)) {
    heatmap <- heatmap +
      scale_fill_gradientn(colours = color_palette)
  }

  # Add attributes for Quarto
  attr(heatmap, "caption") <- sprintf(
    "Phylum-level heatmap grouped by %s",
    group_by
  )

  class(heatmap) <- c("phylum_heatmap", class(heatmap))

  return(heatmap)
}

#' Print method for phylum_heatmap objects
#'
#' @param x A phylum_heatmap object
#' @param ... Additional arguments passed to print
#'
#' @export
print.phylum_heatmap <- function(x, ...) {
  print(ggplot2::ggplot(x))
}

#' Save phylum heatmap with optimal dimensions
#'
#' @param heatmap A phylum_heatmap object
#' @param filename Character. Path to save the file
#' @param width Numeric. Width in inches (default: 12)
#' @param height Numeric. Height in inches (default: 8)
#' @param dpi Numeric. Resolution (default: 300)
#'
#' @export
save_phylum_heatmap <- function(
  heatmap,
  filename,
  width = 12,
  height = 8,
  dpi = 300
) {
  ggsave(
    filename = filename,
    plot = heatmap,
    width = width,
    height = height,
    dpi = dpi
  )
}

# phylum_heatmap <- function(ampvisObject, group_by){
#   heatmap <- amp_heatmap(ampvisObject,
#   group_by = group_by,
#   tax_aggregate = "Phylum",
#   tax_show = 20,
#   plot_values = T,
#   normalise = T,
#   plot_values_size = 8,
#   round = 3,
#   tax_empty="OTU",
#   order_x_by = order_x_by,
#   order_y_by = "cluster") +
#   theme(axis.text.x = element_text(angle = 45, size = 14, vjust = 1), axis.text.y = element_text(size = 14),  legend.position = "right", legend.title = element_text(color = "black", size = 10, face = "bold"))

#   return(heatmap)
# }

#default settings plot relative abundances (normalise = TRUE)
phyla_boxplot <- function(ampvis_object, group) {
  amp_boxplot(
    ampvis_object,
    normalise = TRUE,
    group_by = group,
    tax_aggregate = "Phylum",
    tax_show = 5,
    point_size = 1,
    plot_flip = T
  ) +
    scale_color_brewer(palette = "Set2")
}

#' Calculate Omega-Squared Effect Size for PERMANOVA (adonis) Results
#'
#' @description
#' Calculates either partial or non-partial omega-squared effect sizes from the results
#' of a PERMANOVA analysis (from \code{vegan::adonis()} or \code{vegan::adonis2()}).
#' Omega-squared provides a measure of variance explained that is less biased than R-squared,
#' particularly for small sample sizes.
#'
#' @param adonisOutput An object of class "adonis" or "anova.cca" from vegan's adonis/adonis2 function
#' @param partial logical; if TRUE (default) calculates partial omega-squared,
#'                if FALSE calculates non-partial omega-squared
#'
#' @return Returns a modified version of the input object with an additional column for
#'         omega-squared values. For adonis objects, returns the full object with modified
#'         aov.tab. For anova.cca objects, returns just the modified ANOVA table.
#'
#' @details
#' The function calculates omega-squared (ω²) or partial omega-squared (ω²p) effect sizes
#' following the formulas:
#'
#' Partial ω² = (df_effect * (MS_effect - MS_residual)) /
#'              (df_effect * MS_effect + (N - df_effect) * MS_residual)
#'
#' Non-partial ω² = (SS_effect - df_effect * MS_residual) / (SS_total + MS_residual)
#'
#' Where:
#' - df_effect: degrees of freedom for the effect
#' - MS_effect: Mean Square for the effect
#' - MS_residual: Mean Square for the residual
#' - SS_effect: Sum of Squares for the effect
#' - SS_total: Total Sum of Squares
#' - N: Total sample size
#'
#' @note
#' - Omega-squared values will be NA for residual and total rows
#' - Effect sizes can be negative when effect size is very small
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' data(dune)
#' data(dune.env)
#'
#' # Run PERMANOVA
#' adonis_result <- adonis2(dune ~ Management + A1, data = dune.env)
#'
#' # Calculate partial omega-squared
#' result_with_omega <- adonis_OmegaSq(adonis_result)
#'
#' # Calculate non-partial omega-squared
#' result_with_omega <- adonis_OmegaSq(adonis_result, partial = FALSE)
#' }
#'
#' @export
#' @importFrom methods is
#'
#' @author Original code by [original author], modified by [your name]
#' @references
#' Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared statistics:
#' measures of effect size for some common research designs.
#' Psychological Methods, 8(4), 434-447.
adonis_OmegaSq <- function(adonisOutput, partial = TRUE) {
  # Input validation
  if (!requireNamespace("methods", quietly = TRUE)) {
    stop("Package 'methods' is required but not installed")
  }

  if (
    !(methods::is(adonisOutput, "adonis") ||
      methods::is(adonisOutput, "anova.cca"))
  ) {
    stop("Input must be an object of class 'adonis' or 'anova.cca'")
  }

  if (!is.logical(partial)) {
    stop("'partial' must be logical (TRUE/FALSE)")
  }

  if (methods::is(adonisOutput, "anova.cca")) {
    aov_tab <- adonisOutput
    aov_tab$MeanSqs <- aov_tab$SumOfSqs / aov_tab$Df
    aov_tab$MeanSqs[length(aov_tab$Df)] <- NA
  } else {
    aov_tab <- adonisOutput$aov.tab
  }

  heading <- attr(aov_tab, "heading")

  MS_res <- aov_tab[pmatch("Residual", rownames(aov_tab)), "MeanSqs"]
  if (is.na(MS_res)) {
    stop("Cannot find 'Residual' row in ANOVA table")
  }

  SS_tot <- aov_tab[rownames(aov_tab) == "Total", "SumsOfSqs"]

  N <- aov_tab[rownames(aov_tab) == "Total", "Df"] + 1

  if (partial) {
    omega <- apply(aov_tab, 1, function(x) {
      (x["Df"] * (x["MeanSqs"] - MS_res)) /
        (x["Df"] * x["MeanSqs"] + (N - x["Df"]) * MS_res)
    })
    aov_tab$parOmegaSq <- c(omega[1:(length(omega) - 2)], NA, NA)
  } else {
    omega <- apply(aov_tab, 1, function(x) {
      (x["SumsOfSqs"] - x["Df"] * MS_res) / (SS_tot + MS_res)
    })
    aov_tab$OmegaSq <- c(omega[1:(length(omega) - 2)], NA, NA)
  }
  if (is(adonisOutput, "adonis")) {
    cn_order <- c(
      "Df",
      "SumsOfSqs",
      "MeanSqs",
      "F.Model",
      "R2",
      if (partial) "parOmegaSq" else "OmegaSq",
      "Pr(>F)"
    )
  } else {
    cn_order <- c(
      "Df",
      "SumOfSqs",
      "F",
      if (partial) "parOmegaSq" else "OmegaSq",
      "Pr(>F)"
    )
  }
  aov_tab <- aov_tab[, cn_order]
  attr(aov_tab, "names") <- cn_order
  attr(aov_tab, "heading") <- heading
  if (is(adonisOutput, "adonis")) {
    adonisOutput$aov.tab <- aov_tab
  } else {
    adonisOutput <- aov_tab
  }
  return(adonisOutput)
}

hm_data_alpha <- function(hm, i, metadata) {
  # Input validation
  if (!is.list(hm) || is.null(hm$data)) {
    stop("'hm' must be a list containing a 'data' element")
  }
  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data frame")
  }
  if (!all(c("sample_names", i) %in% colnames(metadata))) {
    stop(
      "metadata must contain 'sample_names' and the grouping variable '",
      i,
      "'"
    )
  }

  # Process data
  df <- hm$data %>%
    dplyr::select(
      phyla = Display,
      Sample,
      `Chao 1 Index` = chao1,
      `Shannon Diversity Index` = diversity_shannon,
      Evenness = evenness_pielou,
      `Faith's Phylogenetic Diversity` = pd
    ) %>%
    # Join with metadata using proper column names
    left_join(
      metadata %>%
        dplyr::select(Sample = sample_names, grouping_var = !!sym(i)),
      by = "Sample"
    ) %>%
    # Calculate summary statistics
    group_by(phyla, grouping_var) %>%
    summarise(
      # mean = mean(relAbund, na.rm = TRUE),
      # median = median(relAbund, na.rm = TRUE),
      # max = max(relAbund, na.rm = TRUE),
      # min = min(relAbund, na.rm = TRUE),
      # n = n(),
      `Median (IQR; Min-Max)` = paste0(
        round(median(relAbund, na.rm = TRUE), 3),
        "\n(",
        round(IQR(relAbund), 3),
        "; ",
        round(min(relAbund, na.rm = TRUE), 3),
        "-",
        round(max(relAbund, na.rm = TRUE), 3),
        ")"
      ),
      .groups = "drop"
    )

  return(df)
}

#' Perform Shapiro-Wilk normality test on all numeric columns in a dataframe
#'
#' @param df A dataframe containing the data to be tested
#' @param min_samples Minimum number of samples required (default = 3, but note that shapiro.test requires n >= 3)
#' @param exclude_cols Character vector of column names to exclude (optional)
#' @param alpha Significance level for the test (default = 0.05)
#' @return A dataframe containing test results including statistics, p-values, and normality conclusion
#'
#' @examples
#' # Basic usage
#' results <- multi_shapiro(my_dataframe)
#'
#' # Exclude specific columns and set minimum samples
#' results <- multi_shapiro(my_dataframe,
#'                              min_samples = 10,
#'                              exclude_cols = c("ID", "Group"))
multi_shapiro <- function(
  df,
  min_samples = 3,
  exclude_cols = NULL,
  alpha = 0.05
) {
  # Input validation
  if (!is.data.frame(df)) {
    stop("Input must be a dataframe")
  }

  if (min_samples < 3) {
    warning(
      "min_samples set below 3. Shapiro-Wilk test requires at least 3 samples."
    )
    min_samples <- 3
  }

  # Initialize results dataframe
  results <- data.frame(
    variable = character(),
    n_samples = numeric(),
    statistic = numeric(),
    p_value = numeric(),
    is_normal = logical(),
    stringsAsFactors = FALSE
  )

  # Get columns to process
  cols_to_process <- names(df)
  if (!is.null(exclude_cols)) {
    cols_to_process <- setdiff(cols_to_process, exclude_cols)
  }

  # Process each column
  for (col in cols_to_process) {
    # Skip if not numeric
    if (!is.numeric(df[[col]])) {
      next
    }

    # Remove NA values
    valid_data <- na.omit(df[[col]])
    n_samples <- length(valid_data)

    # Skip if not enough samples
    if (n_samples < min_samples) {
      next
    }

    # Perform Shapiro-Wilk test
    test_result <- tryCatch(
      shapiro.test(valid_data),
      error = function(e) NULL,
      warning = function(w) NULL
    )

    # Add results if test was successful
    if (!is.null(test_result)) {
      results <- rbind(
        results,
        data.frame(
          variable = col,
          n_samples = n_samples,
          statistic = test_result$statistic,
          p_value = test_result$p.value,
          is_normal = test_result$p.value >= alpha
        )
      )
    }
  }

  # Add descriptive statistics if results exist
  if (nrow(results) > 0) {
    # Calculate skewness and kurtosis
    results$skewness <- sapply(results$variable, function(col) {
      data <- na.omit(df[[col]])
      mean((data - mean(data))^3) / (sd(data)^3)
    })

    results$kurtosis <- sapply(results$variable, function(col) {
      data <- na.omit(df[[col]])
      mean((data - mean(data))^4) / (sd(data)^4) - 3
    })

    # Format p-values for printing
    results$p_value_fmt <- sprintf("%.3e", results$p_value)

    # Add interpretation
    results$interpretation <- ifelse(
      results$is_normal,
      "Normal distribution",
      "Non-normal distribution"
    )
  }

  # Sort by p-value
  results <- results[order(results$p_value), ]

  # Return results
  return(results)
}

# Print method for nicer output
print.shapiro_results <- function(x, ...) {
  print(data.frame(
    Variable = x$variable,
    N = x$n_samples,
    "Shapiro-W" = round(x$statistic, 3),
    "p-value" = x$p_value_fmt,
    "Normality" = x$interpretation
  ))
}

#' Test if Data Follows a Gamma Distribution
#'
#' @description
#' Tests whether data follows a gamma distribution using multiple methods:
#' 1. Visual QQ plot comparison
#' 2. Kolmogorov-Smirnov test
#' 3. Anderson-Darling test
#' 4. Shape and scale parameter estimation
#' 5. Skewness and kurtosis analysis
#'
#' @param x Numeric vector of data
#' @param plot Logical, whether to create diagnostic plots
#' @param alpha Significance level for tests
#' @return List containing test results and diagnostics
#'
#' @import fitdistrplus
#' @import ggplot2
test_gamma_distribution <- function(x, plot = TRUE, alpha = 0.05) {
  if (!is.numeric(x)) {
    stop("Data must be numeric")
  }
  if (any(x <= 0)) {
    stop("Gamma distribution requires positive values")
  }

  # Initialize results list
  results <- list()

  # Fit gamma distribution
  tryCatch(
    {
      require(fitdistrplus)
      fit <- fitdistrplus::fitdist(x, "gamma")

      # Store parameter estimates
      results$parameters <- list(
        shape = fit$estimate["shape"],
        rate = fit$estimate["rate"],
        scale = 1 / fit$estimate["rate"]
      )

      # Calculate theoretical quantiles
      theoretical_quantiles <- qgamma(
        ppoints(length(x)),
        shape = results$parameters$shape,
        rate = results$parameters$rate
      )

      # Perform Kolmogorov-Smirnov test
      ks_test <- ks.test(
        x,
        "pgamma",
        shape = results$parameters$shape,
        rate = results$parameters$rate
      )
      results$ks_test <- list(
        statistic = ks_test$statistic,
        p_value = ks_test$p.value,
        is_gamma = ks_test$p.value > alpha
      )

      # Calculate Anderson-Darling test
      ad_test <- tryCatch(
        {
          require(ADGofTest)
          ADGofTest::ad.test(
            x,
            null = "pgamma",
            shape = results$parameters$shape,
            rate = results$parameters$rate
          )
        },
        error = function(e) NULL
      )

      if (!is.null(ad_test)) {
        results$ad_test <- list(
          statistic = ad_test$statistic,
          p_value = ad_test$p.value,
          is_gamma = ad_test$p.value > alpha
        )
      }

      # Calculate summary statistics
      results$summary <- list(
        n = length(x),
        mean = mean(x),
        sd = sd(x),
        skewness = mean((x - mean(x))^3) / sd(x)^3,
        kurtosis = mean((x - mean(x))^4) / sd(x)^4 - 3,
        cv = sd(x) / mean(x) # Coefficient of variation
      )

      # Compare theoretical vs observed moments
      theoretical_cv <- 1 / sqrt(results$parameters$shape)
      results$cv_comparison <- list(
        observed = results$summary$cv,
        theoretical = theoretical_cv,
        difference = abs(results$summary$cv - theoretical_cv)
      )

      # Create diagnostic plots if requested
      if (plot) {
        require(ggplot2)

        # QQ plot
        qq_data <- data.frame(
          theoretical = theoretical_quantiles,
          observed = sort(x)
        )

        qq_plot <- ggplot(qq_data, aes(theoretical, observed)) +
          geom_point(alpha = 0.5) +
          geom_abline(
            intercept = 0,
            slope = 1,
            color = "red",
            linetype = "dashed"
          ) +
          labs(
            title = "Q-Q Plot against Gamma Distribution",
            x = "Theoretical Quantiles",
            y = "Sample Quantiles"
          ) +
          theme_minimal()

        # Density plot
        x_seq <- seq(min(x), max(x), length.out = 100)
        density_data <- data.frame(
          x = x_seq,
          theoretical = dgamma(
            x_seq,
            shape = results$parameters$shape,
            rate = results$parameters$rate
          )
        )

        density_plot <- ggplot() +
          geom_density(aes(x), data = data.frame(x = x)) +
          geom_line(
            data = density_data,
            aes(x = x, y = theoretical),
            color = "red",
            linetype = "dashed"
          ) +
          labs(title = "Density Plot Comparison", x = "Value", y = "Density") +
          theme_minimal()

        results$plots <- list(
          qq_plot = qq_plot,
          density_plot = density_plot
        )
      }

      # Overall assessment
      results$conclusion <- list(
        is_gamma = all(c(
          results$ks_test$is_gamma,
          if (!is.null(ad_test)) results$ad_test$is_gamma else TRUE,
          results$cv_comparison$difference < 0.1
        )),
        confidence = "medium",
        reasons = character(0)
      )

      # Add reasons for conclusion
      if (results$ks_test$p_value < alpha) {
        results$conclusion$reasons <- c(
          results$conclusion$reasons,
          "Failed Kolmogorov-Smirnov test"
        )
      }
      if (!is.null(ad_test) && ad_test$p.value < alpha) {
        results$conclusion$reasons <- c(
          results$conclusion$reasons,
          "Failed Anderson-Darling test"
        )
      }
      if (results$cv_comparison$difference >= 0.1) {
        results$conclusion$reasons <- c(
          results$conclusion$reasons,
          "Coefficient of variation differs from theoretical"
        )
      }
    },
    error = function(e) {
      results$error <- e$message
      results$conclusion <- list(
        is_gamma = FALSE,
        confidence = "low",
        reasons = "Error in fitting gamma distribution"
      )
    }
  )

  class(results) <- "gamma_test"
  return(results)
}

# Print method for gamma_test results
print.gamma_test <- function(x, ...) {
  cat("\nGamma Distribution Test Results\n")
  cat("==============================\n\n")

  if (!is.null(x$error)) {
    cat("Error:", x$error, "\n")
    return(invisible(x))
  }

  cat("Parameter Estimates:\n")
  cat("Shape:", round(x$parameters$shape, 4), "\n")
  cat("Scale:", round(x$parameters$scale, 4), "\n")

  cat("\nTest Statistics:\n")
  cat(
    "Kolmogorov-Smirnov test p-value:",
    format.pval(x$ks_test$p_value, digits = 4),
    "\n"
  )
  if (!is.null(x$ad_test)) {
    cat(
      "Anderson-Darling test p-value:",
      format.pval(x$ad_test$p_value, digits = 4),
      "\n"
    )
  }

  cat("\nSummary Statistics:\n")
  cat("Mean:", round(x$summary$mean, 4), "\n")
  cat("SD:", round(x$summary$sd, 4), "\n")
  cat("Skewness:", round(x$summary$skewness, 4), "\n")
  cat("Kurtosis:", round(x$summary$kurtosis, 4), "\n")
  cat(
    "Coefficient of Variation:",
    round(x$summary$cv, 4),
    "(Expected:",
    round(x$cv_comparison$theoretical, 4),
    ")\n"
  )

  cat("\nConclusion:\n")
  cat(
    "Data appears to",
    if (x$conclusion$is_gamma) "" else "not",
    "follow a gamma distribution\n"
  )
  if (length(x$conclusion$reasons) > 0) {
    cat("Reasons:\n")
    cat(paste("-", x$conclusion$reasons), sep = "\n")
  }

  if (!is.null(x$plots)) {
    cat("\nPlots available in results$plots\n")
  }
}

# Usage example:
# result <- test_gamma_distribution(your_data)
# print(result)
# if (require(gridExtra) && !is.null(result$plots)) {
#     grid.arrange(result$plots$qq_plot, result$plots$density_plot, ncol = 2)
# }

#' Test if Data Follows a Negative Binomial Distribution
#'
#' @description
#' Tests whether count data follows a negative binomial distribution using:
#' 1. Chi-square goodness of fit test
#' 2. Visual comparison plots
#' 3. Mean-variance relationship check
#' 4. Parameter estimation
#' 5. Comparison with Poisson distribution
#'
#' @param x Numeric vector of count data
#' @param plot Logical, whether to create diagnostic plots
#' @param alpha Significance level for tests
#' @return List containing test results and diagnostics
#'
#' @import MASS
#' @import ggplot2
test_negative_binomial <- function(x, plot = TRUE, alpha = 0.05) {
  if (!is.numeric(x)) {
    stop("Data must be numeric")
  }
  if (any(x %% 1 != 0)) {
    stop("Data must be counts (integers)")
  }
  if (any(x < 0)) {
    stop("Negative binomial distribution requires non-negative values")
  }

  results <- list()

  # Basic statistics
  results$summary <- list(
    n = length(x),
    mean = mean(x),
    variance = var(x),
    overdispersion = var(x) / mean(x)
  )

  # Fit negative binomial distribution
  tryCatch(
    {
      require(MASS)
      nb_fit <- MASS::fitdistr(x, "negative binomial")
      results$parameters <- list(
        size = nb_fit$estimate["size"],
        mu = nb_fit$estimate["mu"]
      )

      # Fit Poisson for comparison
      pois_fit <- MASS::fitdistr(x, "poisson")
      results$poisson_lambda <- pois_fit$estimate

      # Calculate theoretical probabilities
      max_x <- max(x)
      x_range <- 0:max_x

      nb_probs <- dnbinom(
        x_range,
        size = results$parameters$size,
        mu = results$parameters$mu
      )
      pois_probs <- dpois(x_range, lambda = results$poisson_lambda)

      # Observed frequencies
      obs_freq <- table(factor(x, levels = x_range))

      # Expected frequencies
      n <- length(x)
      nb_exp <- n * nb_probs
      pois_exp <- n * pois_probs

      # Chi-square test
      # Combine cells with expected frequency < 5
      combine_cells <- function(obs, exp) {
        while (any(exp < 5) && length(exp) > 1) {
          idx <- which(exp < 5)[1]
          if (idx == length(exp)) {
            exp[idx - 1] <- exp[idx - 1] + exp[idx]
            obs[idx - 1] <- obs[idx - 1] + obs[idx]
            exp <- exp[-idx]
            obs <- obs[-idx]
          } else {
            exp[idx + 1] <- exp[idx + 1] + exp[idx]
            obs[idx + 1] <- obs[idx + 1] + obs[idx]
            exp <- exp[-(idx)]
            obs <- obs[-(idx)]
          }
        }
        return(list(obs = obs, exp = exp))
      }

      combined <- combine_cells(obs_freq, nb_exp)

      chi_sq_stat <- sum((combined$obs - combined$exp)^2 / combined$exp)
      df <- length(combined$obs) - length(nb_fit$estimate) - 1
      p_value <- pchisq(chi_sq_stat, df, lower.tail = FALSE)

      results$goodness_of_fit <- list(
        method = "Chi-square",
        statistic = chi_sq_stat,
        df = df,
        p_value = p_value,
        is_negbin = p_value > alpha
      )

      # AIC comparison with Poisson
      results$model_comparison <- list(
        negbin_aic = 2 * nb_fit$loglik + 2 * length(nb_fit$estimate),
        poisson_aic = 2 * pois_fit$loglik + 2 * length(pois_fit$estimate),
        preferred_model = if (nb_fit$loglik > pois_fit$loglik) {
          "Negative Binomial"
        } else {
          "Poisson"
        }
      )

      # Create diagnostic plots if requested
      if (plot) {
        require(ggplot2)

        # Probability mass function comparison
        pmf_data <- data.frame(
          x = rep(x_range, 3),
          probability = c(nb_probs, pois_probs, (as.vector(obs_freq) / n)),
          type = factor(rep(
            c("Negative Binomial", "Poisson", "Observed"),
            each = length(x_range)
          ))
        )

        pmf_plot <- ggplot(
          pmf_data,
          aes(x = x, y = probability, color = type)
        ) +
          geom_line() +
          geom_point() +
          labs(
            title = "Probability Mass Function Comparison",
            x = "Count",
            y = "Probability"
          ) +
          theme_minimal()

        # Q-Q plot
        theoretical_quantiles <- qnbinom(
          ppoints(length(x)),
          size = results$parameters$size,
          mu = results$parameters$mu
        )

        qq_data <- data.frame(
          theoretical = sort(theoretical_quantiles),
          observed = sort(x)
        )

        qq_plot <- ggplot(qq_data, aes(theoretical, observed)) +
          geom_point(alpha = 0.5) +
          geom_abline(
            intercept = 0,
            slope = 1,
            color = "red",
            linetype = "dashed"
          ) +
          labs(
            title = "Q-Q Plot against Negative Binomial Distribution",
            x = "Theoretical Quantiles",
            y = "Sample Quantiles"
          ) +
          theme_minimal()

        results$plots <- list(
          pmf_plot = pmf_plot,
          qq_plot = qq_plot
        )
      }

      # Overall assessment
      results$conclusion <- list(
        is_negbin = results$goodness_of_fit$is_negbin &&
          results$model_comparison$preferred_model == "Negative Binomial",
        confidence = "medium",
        reasons = character(0)
      )

      if (results$goodness_of_fit$p_value < alpha) {
        results$conclusion$reasons <- c(
          results$conclusion$reasons,
          "Failed chi-square goodness of fit test"
        )
      }

      if (results$model_comparison$preferred_model != "Negative Binomial") {
        results$conclusion$reasons <- c(
          results$conclusion$reasons,
          "Poisson distribution provides better fit based on AIC"
        )
      }
    },
    error = function(e) {
      results$error <- e$message
      results$conclusion <- list(
        is_negbin = FALSE,
        confidence = "low",
        reasons = "Error in fitting negative binomial distribution"
      )
    }
  )

  class(results) <- "negbin_test"
  return(results)
}

# Print method for negbin_test results
print.negbin_test <- function(x, ...) {
  cat("\nNegative Binomial Distribution Test Results\n")
  cat("=========================================\n\n")

  if (!is.null(x$error)) {
    cat("Error:", x$error, "\n")
    return(invisible(x))
  }

  cat("Summary Statistics:\n")
  cat("Sample size:", x$summary$n, "\n")
  cat("Mean:", round(x$summary$mean, 4), "\n")
  cat("Variance:", round(x$summary$variance, 4), "\n")
  cat("Overdispersion ratio:", round(x$summary$overdispersion, 4), "\n")

  cat("\nParameter Estimates:\n")
  cat("Size (r):", round(x$parameters$size, 4), "\n")
  cat("Mean (μ):", round(x$parameters$mu, 4), "\n")

  cat("\nGoodness of Fit Test:\n")
  cat("Chi-square statistic:", round(x$goodness_of_fit$statistic, 4), "\n")
  cat("Degrees of freedom:", x$goodness_of_fit$df, "\n")
  cat("P-value:", format.pval(x$goodness_of_fit$p_value, digits = 4), "\n")

  cat("\nModel Comparison (AIC):\n")
  cat("Negative Binomial:", round(x$model_comparison$negbin_aic, 2), "\n")
  cat("Poisson:", round(x$model_comparison$poisson_aic, 2), "\n")
  cat("Preferred model:", x$model_comparison$preferred_model, "\n")

  cat("\nConclusion:\n")
  cat(
    "Data appears to",
    if (x$conclusion$is_negbin) "" else "not",
    "follow a negative binomial distribution\n"
  )
  if (length(x$conclusion$reasons) > 0) {
    cat("Reasons:\n")
    cat(paste("-", x$conclusion$reasons), sep = "\n")
  }

  if (!is.null(x$plots)) {
    cat("\nPlots available in results$plots\n")
  }
}

# Usage example:
# result <- test_negative_binomial(your_count_data)
# print(result)
# if (require(gridExtra) && !is.null(result$plots)) {
#     grid.arrange(result$plots$pmf_plot, result$plots$qq_plot, ncol = 2)
# }

#' Find Best Fitting Distribution for Data
#'
#' @description
#' Tests data against multiple theoretical distributions and ranks them by goodness of fit.
#' Uses multiple criteria including AIC, BIC, and Kolmogorov-Smirnov test.
#'
#' @param x Numeric vector of data
#' @param plot Logical, whether to create diagnostic plots
#' @param distributions Character vector of distributions to test
#' @return List containing ranked distributions and diagnostic plots
#'
#' @import fitdistrplus
#' @import ggplot2
find_dist <- function(
  x,
  plot = TRUE,
  distributions = c(
    "norm",
    "lnorm",
    "gamma",
    "weibull",
    "exp",
    "beta",
    "logis",
    "cauchy",
    "nbinom",
    "pois"
  )
) {
  # Load required packages
  require(fitdistrplus)
  require(ggplot2)

  # Initialize results list
  results <- list()

  # Store original data in results
  results$data <- x

  # Store basic statistics
  results$summary <- list(
    n = length(x),
    mean = mean(x),
    sd = sd(x),
    min = min(x),
    max = max(x),
    skewness = mean((x - mean(x))^3) / sd(x)^3,
    kurtosis = mean((x - mean(x))^4) / sd(x)^4 - 3
  )

  # Initialize data frame for storing fitting results
  fits <- data.frame(
    distribution = character(),
    AIC = numeric(),
    BIC = numeric(),
    KS_stat = numeric(),
    KS_p = numeric(),
    AD_stat = numeric(),
    AD_p = numeric(),
    convergence = logical(),
    stringsAsFactors = FALSE
  )

  # Function to handle special cases
  prepare_data <- function(x, dist) {
    if (dist == "beta") {
      # Scale to (0,1) for beta
      return((x - min(x)) / (max(x) - min(x)))
    }
    if (dist %in% c("gamma", "lnorm", "weibull", "exp")) {
      # Handle non-positive values
      if (any(x <= 0)) {
        warning(paste("Non-positive values found. Adding offset for", dist))
        return(x - min(x) + 0.001)
      }
    }
    return(x)
  }

  # Improved error messages
  distribution_requirements <- list(
    norm = "none",
    lnorm = "positive values",
    gamma = "positive values",
    weibull = "positive values",
    exp = "positive values",
    beta = "values between 0 and 1",
    logis = "none",
    cauchy = "none",
    nbinom = "non-negative integers",
    pois = "non-negative integers"
  )

  # Function to check if data meets distribution requirements
  check_data_requirements <- function(x, dist) {
    switch(
      dist,
      "norm" = TRUE,
      "lnorm" = all(x > 0),
      "gamma" = all(x > 0),
      "weibull" = all(x > 0),
      "exp" = all(x > 0),
      "beta" = all(x >= 0) && all(x <= 1),
      "logis" = TRUE,
      "cauchy" = TRUE,
      "nbinom" = all(x >= 0) && all(x == floor(x)),
      "pois" = all(x >= 0) && all(x == floor(x)),
      TRUE
    )
  }

  # Store fitting messages
  fit_messages <- character()

  # Try fitting each distribution
  for (dist in distributions) {
    tryCatch(
      {
        # Check basic requirements first
        if (!check_data_requirements(x, dist)) {
          fits <- rbind(
            fits,
            data.frame(
              distribution = dist,
              AIC = NA,
              BIC = NA,
              KS_stat = NA,
              KS_p = NA,
              convergence = FALSE,
              message = paste(
                "Data does not meet",
                dist,
                "distribution requirements:",
                distribution_requirements[[dist]]
              ),
              stringsAsFactors = FALSE
            )
          )
          next
        }

        # Prepare data if needed
        x_prepared <- prepare_data(x, dist)

        # Check for ties before K-S test
        has_ties <- any(duplicated(x_prepared))

        # Fit distribution
        fit <- switch(
          dist,
          "norm" = fitdistrplus::fitdist(x_prepared, "norm"),
          "lnorm" = fitdistrplus::fitdist(x_prepared, "lnorm"),
          "gamma" = fitdistrplus::fitdist(x_prepared, "gamma"),
          "weibull" = fitdistrplus::fitdist(x_prepared, "weibull"),
          "exp" = fitdistrplus::fitdist(x_prepared, "exp"),
          "beta" = fitdistrplus::fitdist(x_prepared, "beta"),
          "logis" = fitdistrplus::fitdist(x_prepared, "logis"),
          "cauchy" = fitdistrplus::fitdist(x_prepared, "cauchy"),
          "nbinom" = {
            if (all(x_prepared == floor(x_prepared))) {
              fitdistrplus::fitdist(x_prepared, "nbinom")
            } else {
              NULL
            }
          },
          "pois" = {
            if (all(x_prepared == floor(x_prepared))) {
              fitdistrplus::fitdist(x_prepared, "pois")
            } else {
              NULL
            }
          }
        )

        if (!is.null(fit)) {
          # Perform goodness of fit tests
          if (!has_ties) {
            ks_test <- suppressWarnings(
              ks.test(x_prepared, paste0("p", dist), ...$fit$estimate)
            )
            ks_stat <- ks_test$statistic
            ks_p <- ks_test$p.value
          } else {
            # Use alternative test or note ties
            ks_stat <- NA
            ks_p <- NA
          }

          # Add results to fits dataframe
          fits <- rbind(
            fits,
            data.frame(
              distribution = dist,
              AIC = fit$aic,
              BIC = fit$bic,
              KS_stat = ks_stat,
              KS_p = ks_p,
              convergence = TRUE,
              message = "Successfully fitted",
              stringsAsFactors = FALSE
            )
          )

          # Store fitted parameters
          results[[paste0(dist, "_fit")]] <- fit
        }
      },
      error = function(e) {
        msg <- paste(
          "Failed to fit",
          dist,
          "distribution:",
          gsub("Error.*?: ", "", e$message)
        )
        fit_messages <- c(fit_messages, msg)

        fits <- rbind(
          fits,
          data.frame(
            distribution = dist,
            AIC = NA,
            BIC = NA,
            KS_stat = NA,
            KS_p = NA,
            convergence = FALSE,
            message = msg,
            stringsAsFactors = FALSE
          )
        )
      }
    )
  }

  # Modified print method to include helpful messages
  results$fit_messages <- fit_messages
  results$fits <- fits

  class(results) <- "distfit"
  return(results)
}

print.distfit <- function(x, ...) {
  cat("\nDistribution Fitting Results\n")
  cat("==========================\n\n")

  cat("Data Summary:\n")
  cat("n:", x$summary$n, "\n")
  cat("Mean:", round(x$summary$mean, 4), "\n")
  cat("SD:", round(x$summary$sd, 4), "\n")
  cat("Skewness:", round(x$summary$skewness, 4), "\n")
  cat("Kurtosis:", round(x$summary$kurtosis, 4), "\n\n")

  # Show successful fits - CORRECTED THIS SECTION
  if (!is.null(x$fits) && nrow(x$fits) > 0) {
    successful_fits <- x$fits[x$fits$convergence == TRUE, ]
    if (nrow(successful_fits) > 0) {
      cat("Successful fits:\n")
      print(head(successful_fits[, c("distribution", "AIC", "BIC", "KS_p")], 5))
    } else {
      cat("No distributions were successfully fitted.\n")
    }

    # Print messages from failed fits
    failed_fits <- x$fits[!x$fits$convergence, ]
    if (nrow(failed_fits) > 0) {
      cat("\nFailed fits:\n")
      for (i in 1:nrow(failed_fits)) {
        cat(
          "- ",
          failed_fits$distribution[i],
          ": ",
          failed_fits$message[i],
          "\n",
          sep = ""
        )
      }
    }
  } else {
    cat("No distribution fitting results available.\n")
  }

  # Data characteristics
  cat("\nData Characteristics:\n")
  if (any(duplicated(x$data))) {
    cat("- Contains duplicate values (ties)\n")
  }
  if (any(x$data <= 0)) {
    cat("- Contains zero or negative values\n")
  }
  if (all(x$data == floor(x$data))) {
    cat("- All values are integers\n")
  }
  if (x$summary$skewness > 2) {
    cat("- Highly skewed (skewness > 2)\n")
  }
  if (x$summary$kurtosis > 7) {
    cat("- Heavy-tailed (high kurtosis)\n")
  }

  # Distribution-specific messages
  cat("\nFitting Notes:\n")
  if (any(x$data <= 0)) {
    cat(
      "- Gamma, lognormal, Weibull, and exponential distributions require positive values\n"
    )
    cat("- Consider adding a constant to make all values positive\n")
  }
  if (x$summary$skewness > 2) {
    cat("- High skewness may indicate need for transformation\n")
  }
  if (x$summary$kurtosis > 7) {
    cat("- High kurtosis suggests heavy tails\n")
  }

  # Suggestions
  cat("\nSuggestions:\n")
  if (x$summary$skewness > 2) {
    cat("- Consider log or power transformation\n")
  }
  if (all(x$data >= 0) && all(x$data == floor(x$data))) {
    cat("- Data appears to be counts, consider Poisson or Negative Binomial\n")
  }
  if (all(x$data > 0) && x$summary$skewness > 2) {
    cat(
      "- For positive, skewed data, consider gamma or lognormal distributions\n"
    )
  }
  if (x$summary$kurtosis > 7) {
    cat("- For heavy-tailed data, consider stable or t distributions\n")
  }

  # KS test warning if ties present
  if (any(duplicated(x$data))) {
    cat(
      "\nNote: Kolmogorov-Smirnov p-values may be unreliable due to ties in the data\n"
    )
  }

  if (!is.null(x$plots)) {
    cat("\nPlots available in results$plots\n")
  }
}

#' Create a Volcano Plot for Differential Analysis
#'
#' This function creates a volcano plot to visualize the results of differential
#' analysis, plotting log fold change against -log10 FDR p-values with
#' customizable significance thresholds.
#'
#' @param df A data frame containing differential analysis results.
#' @param title Character string. Optional title for the plot. Default is NULL.
#' @param p.fdr Column name (unquoted) containing FDR-adjusted p-values.
#'   Default is `p.adj.fdr`.
#' @param lfc Column name (unquoted) containing log fold change values.
#'   Default is `lfc`.
#' @param p.threshold Numeric. FDR p-value threshold for significance.
#'   Default is 0.05.
#' @param lfc.threshold Numeric. Log fold change threshold for significance.
#'   Default is log(2).
#' @param point.size Numeric. Size of points in the plot. Default is 2.
#' @param sig.color Character. Color for significant points. Default is "red".
#' @param nonsig.color Character. Color for non-significant points. Default is "blue".
#' @param width Numeric. Plot width in inches for output. Default is 8.
#' @param height Numeric. Plot height in inches for output. Default is 6.
#' @param dpi Numeric. Resolution for raster output. Default is 300.
#'
#' @return A ggplot2 object representing a volcano plot with proper theming
#'   for both HTML and PDF output.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' sample_data <- data.frame(
#'   p.adj.fdr = runif(100, 0, 1),
#'   lfc = rnorm(100, 0, 2)
#' )
#' plotDAvol(sample_data, title = "Differential Expression Analysis")
#'
#' # With custom thresholds
#' plotDAvol(sample_data, 
#'           title = "Custom Analysis", 
#'           p.threshold = 0.01, 
#'           lfc.threshold = log(1.5))
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang enquo quo_name
#' @export
plotDAvol <- function(df, 
                      title = NULL, 
                      p.fdr = p.adj.fdr, 
                      lfc = lfc,
                      p.threshold = 0.05,
                      lfc.threshold = log(2),
                      point.size = 2,
                      sig.color = "#2166AC",
                      nonsig.color = "#B2182B",
                      width = 8,
                      height = 6,
                      dpi = 300) {
  
  # Input validation
  if (!is.data.frame(df)) {
    stop("df must be a data frame")
  }
  
  if (nrow(df) == 0) {
    stop("df cannot be empty")
  }
  
  # Capture column names using non-standard evaluation
  p.fdr_col <- enquo(p.fdr)
  lfc_col <- enquo(lfc)
  
  p.fdr_name <- quo_name(p.fdr_col)
  lfc_name <- quo_name(lfc_col)
  
  # Check if required columns exist
  if (!p.fdr_name %in% names(df)) {
    stop(paste("Column", p.fdr_name, "not found in df"))
  }
  
  if (!lfc_name %in% names(df)) {
    stop(paste("Column", lfc_name, "not found in df"))
  }
  
  # Validate numeric columns
  if (!is.numeric(df[[p.fdr_name]])) {
    stop(paste("Column", p.fdr_name, "must be numeric"))
  }
  
  if (!is.numeric(df[[lfc_name]])) {
    stop(paste("Column", lfc_name, "must be numeric"))
  }
  
  # Check for valid p-values
  if (any(df[[p.fdr_name]] < 0 | df[[p.fdr_name]] > 1, na.rm = TRUE)) {
    warning("Some p-values are outside the range [0, 1]")
  }
  
  # Validate parameters
  if (!is.null(title) && !is.character(title)) {
    stop("title must be a character string or NULL")
  }
  
  if (!is.numeric(p.threshold) || p.threshold <= 0 || p.threshold > 1) {
    stop("p.threshold must be a number between 0 and 1")
  }
  
  if (!is.numeric(lfc.threshold) || lfc.threshold < 0) {
    stop("lfc.threshold must be a positive number")
  }
  
  # Create the plot with improved styling for Quarto output
  p <- df %>%
    mutate(
      logP = -log10(!!p.fdr_col),
      significance = case_when(
        logP > -log10(p.threshold) & abs(!!lfc_col) > lfc.threshold ~ "Significant",
        TRUE ~ "Not Significant"
      )
    ) %>%
    ggplot(aes(x = !!lfc_col, y = logP, color = significance)) +
    geom_point(size = point.size, alpha = 0.7) +
    scale_color_manual(
      values = c("Significant" = sig.color, "Not Significant" = nonsig.color),
      name = "Significance"
    ) +
    labs(
      x = bquote("Log Fold Change" ~ (log[e](FC))),
      y = bquote("FDR" ~ (-log[10] ~ "P-value")),
      title = title
    ) +
    geom_hline(yintercept = -log10(p.threshold), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = c(-lfc.threshold, lfc.threshold), linetype = "dashed", alpha = 0.5) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      # Ensure compatibility with both HTML and PDF output
      text = element_text(family = "sans"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # Add attributes for Quarto compatibility
  attr(p, "width") <- width
  attr(p, "height") <- height
  attr(p, "dpi") <- dpi
  
  return(p)
}