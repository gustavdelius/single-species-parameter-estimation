# Function that updates params object with new parameter values
update_params <- function(params, pars) {
    sp <- params@species_params
    gp <- params@gear_params

    gp$l50 <- pars["l50"]
    gp$l25 <- pars["ratio"] * gp$l50
    gp$catchability <- pars["catchability"]
    sp$M <- pars["M"]
    ext_mort(params)[] <- sp$M * params@w^sp$d

    # We'll test for the existence of the following parameters so that we do no
    # need to include them in our optimisation unless we want to.
    if (hasName(pars, "w_mat")) {
        sp$w_mat <- pars["w_mat"]
    }
    if (hasName(pars, "U")) {
        sp$w_mat25 <- sp$w_mat / 3^(1 / pars["U"])
    }

    params@species_params <- sp
    params <- setReproduction(params)
    gear_params(params) <- gp

    # Calculate the new steady state ----
    params <- steadySingleSpecies(params)
    # Rescale it to get the observed biomass
    total <- sum(params@initial_n * params@w * params@dw)
    factor <- sp$biomass_observed / total
    params@initial_n <- params@initial_n * factor

    return(params)
}
