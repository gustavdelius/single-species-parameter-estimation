# Project population over time using the dynamic model
gpt_project <- function(N0, G, Z, R, n_steps_per_year) {
    # Number of years
    n_years <- length(R)
    # Number of size bins
    n_bins <- length(N0)
    # Time step size
    delta_t <- 1 / n_steps_per_year
    # Initialize the population matrix
    N <- matrix(0, nrow = n_years * n_steps_per_year + 1, ncol = n_bins)
    N[1, ] <- N0
    # Project the population over time
    for (y in 1:n_years) {
        for (t in 1:n_steps_per_year) {
            i <- (y - 1) * n_steps_per_year + t
            for (j in 1:n_bins) {
                growth_contrib <- (if (j == 1) R[y] else (G[y, j - 1] * N[i + 1, j - 1])) * delta_t / bin_width[j]
                denominator <- 1 + (G[y, j] * delta_t / bin_width[j]) + (Z[y, j] * delta_t)
                N[i + 1, j] <- (N[i, j] + growth_contrib) / denominator
            }
        }
    }
    return(N)
}

gpt_steady <- function(N00, g, z, bin_width) {
    n_bins <- length(g)
    N <- rep(0, n_bins)
    N[1] <- N00
    for (j in 2:n_bins) {
        N[j] <- g[j - 1] * N[j - 1] / (g[j] + z[j] * bin_width[j])
    }
    return(N)
}
