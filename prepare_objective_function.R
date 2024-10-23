#' Prepare an Objective Function for Optimizing Model Parameters
#'
#' This function returns an objective function that, given a set of parameters,
#' calculates he negative log-likelihood of the catch and adds it to a penalty
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
#' @param yield_lambda A parameter that controls the strength of the penalty for
#'  deviation from the observed yield.
#'
#' @return The objective function
#'
prepare_objective_function <- function(params, df, yield_lambda) {

    # Validate MizerParams object
    params <- validParams(params)
    sp <- params@species_params

    # We use lengths instead of weights in the catch data
    lengths <- (params@w / sp$a)^(1/sp$b)

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
    N <- sum(counts)

    # Some terms in the log-likelihood formula for the multinomial
    # distribution are independent of the probabilities and can be precomputed
    data_log_likelihood_constant <- lgamma(N + 1) - sum(lgamma(counts + 1))

    # When we calculate the likelihood, we will need to have values of the
    # modelled catch density at all bin boundaries. We will use interpolation
    # for that purpose
    # Collect all unique bin boundaries for interpolation
    all_bin_boundaries <- unique(c(full_bins$bin_start, full_bins$bin_end))

    # Precompute bin widths
    bin_widths <- full_bins$bin_end - full_bins$bin_start

    # Function that computes the log-likelihood given a probability density function
    log_likelihood_fn <- function(lengths, pdf_values) {
        # lengths: Numeric vector of lengths at which the density function is
        #          calculated by mizer
        # pdf_values: Numeric vector of PDF values at those lengths

        # Validate inputs
        if (length(lengths) != length(pdf_values)) {
            stop("lengths and pdf_values must be vectors of the same length.")
        }

        # Interpolate the PDF at all bin boundaries at once
        interpolated_pdf <- approx(lengths, pdf_values,
                                   xout = all_bin_boundaries, rule = 2)$y

        # Create a named vector for easy lookup
        names(interpolated_pdf) <- all_bin_boundaries

        # Get the densities at bin starts and ends
        density_start <- interpolated_pdf[as.character(full_bins$bin_start)]
        density_end <- interpolated_pdf[as.character(full_bins$bin_end)]

        # To get the probability in each bin, we need to integrate the density.
        # Approximate the integral over the bin using the trapezoidal rule
        bin_probs <- ((density_start + density_end) / 2) * bin_widths

        # Ensure probabilities are positive and normalize them
        if (any(bin_probs <= 0)) {
            stop("All bin probabilities must be positive.")
        }
        bin_probs <- bin_probs / sum(bin_probs)

        # Compute the log-likelihood using the multinomial distribution formula
        LL <- data_log_likelihood_constant + sum(counts * log(bin_probs))

        return(LL)
    }

    # Return the objective function
    function(pars) {

        # Update the model using the new values of the parameters
        # We use `pars` for the parameters that we are optimizing
        # and `params` for the full MizerParams object
        params <- update_params(params, pars)

        # Get the catch density as function of weight from the model
        pdf_values <- params@initial_n * getFMort(params)
        # Convert to density as a function of length
        pdf_values <- pdf_values * sp$b * params@w / lengths

        # Compute the negative log-likelihood of the observed catch
        neg_log_likelihood <- -log_likelihood_fn(lengths, pdf_values)

        # Handle cases where the log-likelihood is infinite or NaN
        if (is.infinite(neg_log_likelihood) || is.nan(neg_log_likelihood)) {
            neg_log_likelihood <- .Machine$double.xmax
        }

        # Add a penalty for deviation from observed yield
        penalty <- neg_log_likelihood +
            yield_lambda * log(getYield(params) / params@gear_params$yield_observed)^2

        return(penalty)
    }
}
