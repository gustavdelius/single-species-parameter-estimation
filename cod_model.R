# Prepare the model
#
# This is just a simple model for cod, using some numbers from the Celtic Sea.
cod_model <- function() {
    p <- newSingleSpeciesParams(
        species_name = "cod",
        w_max = 42500,
        w_mat = 4000
    )
    sp <- species_params(p)

    # Set predation kernel
    # While we change the predation kernel we want to not change the consumption
    # rate Q, so we'll calculate it before the change and after the change and then
    # adjust gamma to get back to the original value of Q.
    q <- sweep(getEncounter(p) * (1 - getFeedingLevel(p)) *
                   initialN(p), 2, dw(p), "*")
    Q <- rowSums(q)
    sp$pred_kernel_type <- "power_law"
    sp$kernel_exp <- -0.87
    sp$kernel_l_l <- 4.6
    sp$kernel_u_l <- 3
    sp$kernel_l_r <- 12.5
    sp$kernel_u_r <- 4.3
    species_params(p) <- sp
    q_new <- sweep(getEncounter(p) * (1 - getFeedingLevel(p)) *
                       initialN(p), 2, dw(p), "*")
    Q_new <- rowSums(q_new)
    sp$gamma <- sp$gamma * Q / Q_new

    # Set weight-length relationship
    sp$a <- 0.0079
    sp$b <- 3.05

    # Observed biomass
    sp$biomass_observed <- 0.713

    # age at maturity (for calibrating growth)
    sp$age_mat <- 6.3

    # exponent and coefficient of mortality power law
    sp$d <- sp$n - 1
    M <- getExtMort(p) / w(p)^sp$d
    sp$M <- mean(M)
    species_params(p) <- sp

    # Set gear parameters
    gp <- data.frame(
        species = "cod",
        gear = "total",
        sel_func = "sigmoid_length",
        l50 = 40, # made-up number
        l25 = 35, # made-up number
        catchability = 0.2, # made-up number
        yield_observed = 0.415
    )
    gear_params(p) <- gp
    initial_effort(p) <- 1

    p <- steadySingleSpecies(p)
    p <- matchGrowth(p)
    p <- matchBiomasses(p)
    p <- matchYield(p)
    p <- steadySingleSpecies(p)
    return(p)
}
