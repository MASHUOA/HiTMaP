.hitmap_load_packages <- function(..., packages = NULL) {
  if (is.null(packages)) {
    packages <- as.character(substitute(list(...)))[-1L]
  }

  packages <- unique(as.character(packages))
  packages <- packages[nzchar(packages)]
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing)) {
    stop(
      "Missing required package(s): ",
      paste(missing, collapse = ", "),
      ". Install dependencies before running HiTMaP.",
      call. = FALSE
    )
  }

  invisible(lapply(packages, function(package) {
    suppressPackageStartupMessages(
      library(package, character.only = TRUE)
    )
  }))
}

.hitmap_content_bin <- function(x, nbins, method = "content") {
  if (!identical(method, "content")) {
    stop("HiTMaP currently supports content binning only.", call. = FALSE)
  }

  nbins <- as.integer(nbins)[1L]
  if (is.na(nbins) || nbins < 1L) {
    stop("nbins must be positive.", call. = FALSE)
  }

  valid <- !is.na(x)
  if (!any(valid)) {
    return(factor(rep(NA_character_, length(x))))
  }

  if (nbins == 1L) {
    result <- rep(NA_character_, length(x))
    result[valid] <- "1"
    return(factor(result))
  }

  if (length(unique(x[valid])) <= nbins) {
    result <- rep(NA_character_, length(x))
    result[valid] <- as.character(factor(x[valid]))
    return(factor(result))
  }

  probabilities <- seq(1 / nbins, 1 - 1 / nbins, length.out = nbins - 1L)
  breaks <- unique(c(
    -Inf,
    stats::quantile(x[valid], probs = probabilities, names = FALSE),
    Inf
  ))

  cut(x, breaks = breaks, include.lowest = TRUE)
}

.hitmap_rolling_mean <- function(x, width = 2L) {
  width <- as.integer(width)[1L]
  if (is.na(width) || width < 1L) {
    stop("width must be a positive integer.", call. = FALSE)
  }
  if (length(x) < width) {
    return(numeric(0))
  }

  starts <- seq_len(length(x) - width + 1L)
  vapply(
    starts,
    function(start) mean(x[start:(start + width - 1L)]),
    numeric(1)
  )
}

.hitmap_cardinal_peak_pick_methods <- function() {
  c("diff", "sd", "mad", "quantile", "filter", "cwt")
}

.hitmap_normalize_disable_method <- function(method) {
  if (is.null(method)) {
    return(NULL)
  }

  method <- as.character(method[1L])
  if (tolower(method) == "disable") {
    return("Disable")
  }
  method
}

.hitmap_normalize_peak_pick_method <- function(method, fallback = "mad", warn = TRUE) {
  method <- .hitmap_normalize_disable_method(method)
  if (is.null(method)) {
    return(fallback)
  }

  lowered <- tolower(method)
  if (method == "Disable") {
    return(method)
  }
  if (lowered %in% c("default", "default|default")) {
    return("Default")
  }

  valid_methods <- .hitmap_cardinal_peak_pick_methods()
  if (lowered %in% c("adaptive", "simple")) {
    if (isTRUE(warn)) {
      warning(
        "Cardinal no longer supports peakPick method '", method,
        "'. Using '", fallback, "' instead.",
        call. = FALSE
      )
    }
    return(fallback)
  }
  if (!lowered %in% valid_methods) {
    if (isTRUE(warn)) {
      warning(
        "Invalid peakPick method '", method, "'. Using '", fallback,
        "'. Valid Cardinal methods are: ",
        paste(valid_methods, collapse = ", "),
        call. = FALSE
      )
    }
    return(fallback)
  }

  lowered
}

.hitmap_normalize_preprocess_methods <- function(preprocess, warn = TRUE) {
  if (is.null(preprocess)) {
    return(preprocess)
  }

  if (is.null(preprocess$smoothSignal)) {
    preprocess$smoothSignal <- list(method = "Disable")
  } else {
    preprocess$smoothSignal$method <- .hitmap_normalize_disable_method(preprocess$smoothSignal$method)
    if (is.null(preprocess$smoothSignal$method)) {
      preprocess$smoothSignal$method <- "Disable"
    }
  }

  if (is.null(preprocess$reduceBaseline)) {
    preprocess$reduceBaseline <- list(method = "locmin")
  } else {
    preprocess$reduceBaseline$method <- .hitmap_normalize_disable_method(preprocess$reduceBaseline$method)
    if (is.null(preprocess$reduceBaseline$method)) {
      preprocess$reduceBaseline$method <- "locmin"
    }
  }

  if (!is.null(preprocess$normalize)) {
    preprocess$normalize$method <- .hitmap_normalize_disable_method(preprocess$normalize$method)
  }

  if (is.null(preprocess$peakPick)) {
    preprocess$peakPick <- list(method = "mad")
  } else {
    preprocess$peakPick$method <- .hitmap_normalize_peak_pick_method(
      preprocess$peakPick$method,
      warn = warn
    )
  }

  if (!is.null(preprocess$peakAlign)) {
    preprocess$peakAlign$method <- .hitmap_normalize_disable_method(preprocess$peakAlign$method)
    if (is.null(preprocess$peakAlign$method)) {
      preprocess$peakAlign$method <- "Enable"
    }
  }

  preprocess
}

.hitmap_cart2pol <- function(x, y, degrees = FALSE) {
  theta <- atan2(y, x)
  if (isTRUE(degrees)) {
    theta <- theta * 180 / pi
  }
  data.frame(r = sqrt(x^2 + y^2), theta = theta)
}

.hitmap_percent <- function(x, digits = 0) {
  paste0(format(round(x * 100, digits), nsmall = digits, trim = TRUE), "%")
}

.hitmap_formula_atoms <- function(formula) {
  if (is.null(formula) || length(formula) == 0L) {
    return(NULL)
  }

  formula <- as.character(formula[1L])
  if (!nzchar(formula) || identical(toupper(formula), "FALSE")) {
    return(NULL)
  }

  tokens <- regmatches(formula, gregexpr("[A-Z][a-z]*[0-9.]*", formula, perl = TRUE))[[1]]
  if (!length(tokens)) {
    return(NULL)
  }

  elements <- sub("[0-9.]+$", "", tokens)
  counts <- sub("^[A-Z][a-z]*", "", tokens)
  counts <- ifelse(nzchar(counts), as.numeric(counts), 1)
  totals <- tapply(counts, elements, sum)
  totals <- totals[totals != 0]
  as.list(totals)
}

.hitmap_formula_string <- function(formula) {
  atoms <- .hitmap_formula_atoms(formula)
  if (is.null(atoms)) {
    return("")
  }
  paste0(names(atoms), unlist(atoms, use.names = FALSE), collapse = "")
}

.hitmap_formula_mass <- function(formula, charge = 0) {
  atoms <- .hitmap_formula_atoms(formula)
  if (is.null(atoms)) {
    return(0)
  }
  OrgMassSpecR::MonoisotopicMass(atoms, charge = charge)
}

.hitmap_theme_article <- function() {
  ggplot2::theme_classic(base_size = 11)
}
