# This code was created with ChatGPT, see
# https://chatgpt.com/share/673b5cfa-c0dc-8007-a835-4bbd70e794c0
# I added the code to use the cod model in the creation of synthetic data

# Prepare data for TMB
library(TMB)
library(mizerExperimental)
source("cod_model.R")
p <- cod_model()

# Create simulated datasets using the dynamic model
set.seed(123)  # Set seed for reproducibility

# Simulate years and size bins
n_years <- 20
years <- seq(2000, 2000 + n_years - 1)
bin_start <- seq(10, 10000, by = 40)
n_bins <- length(bin_start)
bin_width <- rep(40, n_bins)

# Define given growth and mortality rates
g_given <- approx(w(p), getEGrowth(p), xout = bin_start)$y
m_given <- approx(w(p), getExtMort(p), xout = bin_start)$y
f <- approx(w(p), getFMort(p), xout = bin_start)$y
plot(bin_start, f, type = "l", xlab = "Weight", ylab = "Fishing Mortality Rate")

# Use random walk for recruitment
R0 <- getRDD(p)
sigma_R <- R0  # Standard deviation for recruitment
log_R <- numeric(n_years)
log_R[1] <- log(100)
for (y in 2:n_years) {
    log_R[y] <- rho_R * log_R[y - 1] + rnorm(1, mean = 0, sd = sigma_R)
}
R <- exp(log_R)

# Define fishing mortality scaling factor
log_F0 <- rep(log(0.1), n_years)
F0 <- exp(log_F0)

# Time step size
delta_t <- 0.1
n_steps_per_year <- 1 / delta_t

# Initialize population matrix
N <- matrix(0, nrow = n_years * n_steps_per_year + 1, ncol = n_bins)
N[1, ] <- N0

# Project population over time using the dynamic model
for (y in 1:n_years) {
    for (t in 1:n_steps_per_year) {
        i <- (y - 1) * n_steps_per_year + t
        N_current <- N[i, ]
        F_i <- f * F0[y]
        Z_i <- F_i + m_given

        # Update population using the new dynamic equation
        N_next <- numeric(n_bins)
        for (j in 1:n_bins) {
            growth_contrib <- if (j == 1) R[y] else (g_given[j - 1] * delta_t / bin_width[j - 1]) * N[i + 1, j - 1]
            denominator <- 1 + (g_given[j] * delta_t / bin_width[j]) + (Z_i[j] * delta_t)
            N_next[j] <- (N_current[j] + growth_contrib) / denominator
        }
        N[i + 1, ] <- N_next
    }
}

# Calculate total catches for each year and size bin
C_y <- matrix(0, nrow = n_years, ncol = n_bins)
for (y in 1:n_years) {
    for (t in 1:n_steps_per_year) {
        i <- (y - 1) * n_steps_per_year + t
        N_current <- N[i, ]
        F_i <- f * F0[y]
        Z_i <- F_i + m_given

        # Calculate catches
        C_i <- N_current * (F_i / Z_i) * (1 - exp(-Z_i * delta_t))
        C_y[y, ] <- C_y[y, ] + C_i
    }
}

# Simulate sample sizes for each year
n_sample <- rpois(n_years, lambda = 1000)  # Average sample size of 1000 fish per year

# Simulate counts for each year and size bin based on predicted proportions
simulated_counts <- matrix(NA, nrow = n_years, ncol = n_bins)
for (y in 1:n_years) {
    probs <- C_y[y, ] / sum(C_y[y, ])
    simulated_counts[y, ] <- rmultinom(1, size = n_sample[y], prob = probs)
}

# Simulate yield for each year
yield <- rowSums(C_y) + rnorm(n_years, mean = 0, sd = exp(log_sigma_yield))  # Adding noise with flexible standard deviation

# Prepare simulated landings data
landings <- data.frame(
    year = rep(years, each = n_bins),
    bin_start = rep(seq(1, n_bins), times = n_years),
    count = as.vector(t(simulated_counts))
)

# Prepare data for TMB
# Organize counts into a list or matrix
counts_list <- split(landings$count, landings$year)
# Or create a matrix of counts [n_years x n_bins]
counts_matrix <- simulated_counts  # Use simulated counts matrix

# Prepare data for TMB
data_list <- list(
    years = unique(landings$year),
    bin_start = unique(landings$bin_start),
    bin_width = diff(c(unique(landings$bin_start), max(landings$bin_start) + 1)),  # Assuming equal widths
    obs_counts = counts_matrix,   # Matrix of observed counts [n_years x n_bins]
    n_sample = n_sample,           # Vector of sample sizes [n_years]
    yield = yield,                 # Simulated total yield
    f = f,                         # Selectivity vector
    given_g = g_given,             # Given growth rates for each weight class
    given_m = m_given,             # Given natural mortality rates for each weight class
    delta_t = delta_t              # Time step size
)

# Define parameters and initial values for TMB
parameters <- list(
    log_N0 = rep(log(1000), n_bins),
    log_F0 = rep(log(0.1), length(unique(landings$year))),
    log_R = rep(log(100), length(unique(landings$year))),
    epsilon_g = matrix(0, length(unique(landings$year)), n_bins),
    epsilon_m = matrix(0, length(unique(landings$year)), n_bins),
    log_sigma_g = log(0.2),
    log_sigma_m = log(0.2),
    log_sigma_yield = log(5000),  # Initial guess for standard deviation of yield
    rho_t_g = 0.5,  # Initial guess for temporal correlation for growth rate errors
    rho_s_g = 0.5,  # Initial guess for spatial correlation for growth rate errors
    rho_t_m = 0.5,  # Initial guess for temporal correlation for mortality rate errors
    rho_s_m = 0.5   # Initial guess for spatial correlation for mortality rate errors
)

# Define random effects
random_effects <- c("log_R", "epsilon_g", "epsilon_m")

# Compile and load the TMB model
compile("gpt.cpp")
dyn.load(dynlib("gpt"))

# Create the TMB objective function
obj <- MakeADFun(
    data = data_list,
    parameters = parameters,
    random = random_effects,
    DLL = "gpt"
)

# Optimize
opt <- nlminb(obj$par, obj$fn, obj$gr)

# Report results
rep <- sdreport(obj)
summary(rep)

# Plot results if necessary
plot(unique(landings$year), exp(parameters$log_F0), type = "b", ylab = "Fishing Mortality Scaling Factor (F0)", xlab = "Year")
