catch_pdf <- function(pars) {
    # Update the model using the new values of the parameters ----
    gp$l50 <- pars["l50"]
    gp$l25 <- pars["ratio"] * gp$l50
    sp$M <- pars["M"]
    ext_mort(p)[] <- sp$M * w(p)^sp$d

    # We'll test for the existence of the following parameters so that we do no
    # need to include them in our optimisation unless we want to.
    if (hasName(pars, "w_mat")) {
        sp$w_mat <- pars["w_mat"]
    }
    if (hasName(pars, "U")) {
        sp$w_mat25 <- sp$w_mat / 3^(1 / pars["U"])
    }

    p@species_params <- sp
    p <- setReproduction(p)
    gear_params(p) <- gp

    # Calculate the PDF of the catch in the new steady state ----
    p <- steadySingleSpecies(p)
    catch_pdf <- initialN(p) * getFMort(p) * dw(p)

    return(catch_pdf)
}

objective_function <- function(params) {
    # Get the catch PDF from the model
    pdf_values <- catch_pdf(params)

    # Compute the negative log-likelihood
    neg_log_likelihood <- -log_likelihood_fn(lengths, pdf_values)

    # Handle cases where the log-likelihood is infinite or NaN
    if (is.infinite(neg_log_likelihood) || is.nan(neg_log_likelihood)) {
        neg_log_likelihood <- .Machine$double.xmax
    }

    return(neg_log_likelihood)
}
