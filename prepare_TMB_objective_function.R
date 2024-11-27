#' Prepare a TMB Objective Function for Optimizing Model Parameters
#'
#' This function returns an objective function that can be automatically
#' differentiated, created with `TMB::MakeADFun`, that, given a set of parameters,
#' calculates the negative log-likelihood of the catch and adds it to a penalty
#' term that penalizes deviations from the observed yield.
#'
#' The reason that we have a function that returns the objective function
#' instead of simply defining the objective function directly is that we want
#' to do some pre-processing once and then use the pre-processed data in the
#' objective function. This is more efficient than doing the pre-processing
#' every time the objective function is called.
#'
#' The main preprocessing makes sure that we have a comprehensive set of bins
#' that cover the entire size range, even though there will not be observations
#' at all sizes. Missing observations should be interpreted as a 0 count.
#'
#' @param params A MizerParams object
#' @param df A data frame containing the observed binned count data. It must contain
#'   the following columns:
#'   * `length`: The start of each bin.
#'   * `dl`: The width of each bin.
#'   * `count`: The observed count for each bin.
#' @param pars A named list of starting values for the parameters that will be
#'   optimized.
#'
#' @return The objective function
#'
prepare_TMB_objective_function <- function(params, df, pars) {

    # Validate MizerParams object
    params <- validParams(params)
    sp <- params@species_params

    # Validate data frame
    if (!all(c('length', 'dl', 'count') %in% names(df))) {
        stop("Data frame 'df' must contain columns 'length', 'dl', and 'count'.")
    }

    # Extract observed bin starts, ends, and counts
    observed_bins <- data.frame(
        bin_start = df$length,
        bin_end = df$length + df$dl,
        count = df$count
    )

    # Create a comprehensive set of bin edges covering all observed bins
    bin_edges <- sort(unique(c(observed_bins$bin_start, observed_bins$bin_end)))

    # Define full bins covering the observed range
    full_bins <- data.frame(
        bin_start = bin_edges[-length(bin_edges)],
        bin_end = bin_edges[-1],
        count = 0  # Initialize counts to zero
    )

    # Map observed counts to the corresponding bins in full_bins
    bin_key <- paste0(full_bins$bin_start, "_", full_bins$bin_end)
    observed_bin_key <- paste0(observed_bins$bin_start, "_", observed_bins$bin_end)
    matched_indices <- match(observed_bin_key, bin_key)
    full_bins$count[matched_indices] <- observed_bins$count

    # Extract counts and compute total observations
    counts <- full_bins$count

    # When we calculate the likelihood, we will need to have values of the
    # modelled catch density at all bin boundaries. We will use interpolation
    # for that purpose
    # Collect all unique bin boundaries for interpolation
    l_bin_boundaries <- unique(c(full_bins$bin_start, full_bins$bin_end))
    # Give this in terms of weights
    w_bin_boundaries <- sp$a * l_bin_boundaries^sp$b

    # Precompute bin widths
    w_bin_widths <- diff(w_bin_boundaries)

    # Interpolate EReproAndGrowth and repro_prop to all bin boundaries
    EReproAndGrowth <-
        approx(w(params), getEReproAndGrowth(params), xout = w_bin_boundaries)$y
    repro_prop <-
        approx(w(params), repro_prop(params), xout = w_bin_boundaries)$y

    # Prepare data
    data_list <- list(
        counts = counts,
        bin_widths = w_bin_widths,
        bin_boundaries = w_bin_boundaries,
        bin_boundary_lengths = l_bin_boundaries,
        yield = params@gear_params$yield_observed,
        biomass = sp$biomass_observed,
        EReproAndGrowth = EReproAndGrowth,
        repro_prop = repro_prop,
        w_mat = sp$w_mat,
        d = sp$d
    )

    MakeADFun(data = data_list,
              parameters = pars,
              DLL = "objective_function")
}
