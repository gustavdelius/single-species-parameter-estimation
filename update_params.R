#' Function that updates params object with new parameter values
#'
#' This function currently updates the gear selectivity parameters `l50` and
#' `l25, the catchability, the steepness `U` of the maturity ogive and the
#' coefficient `M` of the power-law mortality rate. It then recalculates the
#' steady state of the model and rescales it to match the observed biomass.
#'
#' @param params The MizerParams object to update
#' @param pars A named numeric vector of parameter values
#' @return The updated MizerParams object
update_params <- function(params, pars) {
    sp <- params@species_params
    gp <- params@gear_params

    # Update the gear parameters
    gp$l50 <- pars["l50"]
    gp$l25 <- pars["ratio"] * gp$l50
    gp$catchability <- pars["catchability"]
    gear_params(params) <- gp

    # recalculate the power-law mortality rate
    sp$M <- pars["M"]
    ext_mort(params)[] <- sp$M * params@w^sp$d

    # Update the steepness of the maturity ogive
    sp$w_mat25 <- sp$w_mat / 3^(1 / pars["U"])
    params@species_params <- sp
    params <- setReproduction(params)

    # Calculate the new steady state ----
    params <- steadySingleSpecies(params)
    # Rescale it to get the observed biomass
    total <- sum(params@initial_n * params@w * params@dw)
    factor <- sp$biomass_observed / total
    params@initial_n <- params@initial_n * factor

    return(params)
}
