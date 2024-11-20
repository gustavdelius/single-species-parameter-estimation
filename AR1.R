generate_separable_AR1 <- function(n_years, n_bins, rho_t, rho_s, sigma_t, sigma_s) {
    # Temporal AR(1) correlation matrix
    temporal_corr <- outer(1:n_years, 1:n_years, function(i, j) rho_t^abs(i - j))
    temporal_cov <- sigma_t^2 * temporal_corr # Temporal covariance

    # Spatial AR(1) correlation matrix
    spatial_corr <- outer(1:n_bins, 1:n_bins, function(i, j) rho_s^abs(i - j))
    spatial_cov <- sigma_s^2 * spatial_corr # Spatial covariance

    # Full spatiotemporal covariance as Kronecker product
    full_cov <- kronecker(temporal_cov, spatial_cov)

    # Generate a sample from a multivariate normal distribution
    library(MASS)
    sample_vector <- mvrnorm(n = 1, mu = rep(0, n_years * n_bins), Sigma = full_cov)

    # Reshape the sample into a matrix of dimensions (n_years x n_bins)
    matrix(sample_vector, nrow = n_years, ncol = n_bins)
}
