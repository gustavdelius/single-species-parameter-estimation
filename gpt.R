# This code was created with ChatGPT, see
# https://chatgpt.com/share/673b5cfa-c0dc-8007-a835-4bbd70e794c0
# I added the code to use the cod model in the creation of synthetic data

library(TMB)
library(mizerExperimental)
source("gpt_animate_N.R")
source("AR1.R")
source("gpt_project.R")

# Create synthetic data ----

# This creates an example mizer model
source("cod_model.R")
p <- cod_model()

set.seed(123)  # Set seed for reproducibility

# Set years and size bins
n_years <- 20
years <- seq(2000, 2000 + n_years - 1)
bin_start <- exp(seq(log(10), log(p@species_params$w_max), length.out = 25))
n_bins <- length(bin_start)
bin_width <- diff(bin_start)
bin_width <- c(bin_width, bin_width[length(bin_width)])  # Add the last bin width
# Time step size
delta_t <- 0.1
n_steps_per_year <- 1 / delta_t

# Set growth and mortality rates
g_given <- approx(w(p), getEGrowth(p), xout = bin_start)$y
plot(bin_start, g_given, type = "l", log = "x", xlab = "Weight", ylab = "Growth rate")
m_given <- approx(w(p), getExtMort(p), xout = bin_start)$y
plot(bin_start, m_given, type = "l", log = "xy", xlab = "Weight", ylab = "Natural mortality rate")
f <- approx(w(p), getFMort(p), xout = bin_start)$y
plot(bin_start, f, type = "l", log = "x", xlab = "Weight", ylab = "Fishing mortality")

# Set initial population size
Ns <- approx(w(p), initialN(p)[1, ], xout = bin_start)$y
plot(bin_start, Ns, type = "l", log = "xy", xlab = "Weight", ylab = "Abundance")
Ns <- gpt_steady(Ns[1], g_given, m_given + f, bin_width)
lines(bin_start, Ns, col = "red")

# Set recruitment
Rs <- Ns[1] * (g_given[1] + (m_given[1] + f[1]) * bin_width[1])
R <- rep(Rs, n_years)

# Check that this is at steady state
G <- matrix(rep(g_given, each = n_years), nrow = n_years, byrow = FALSE)
Z <- matrix(rep(m_given + f, each = n_years), nrow = n_years, byrow = FALSE)
N <- gpt_project(Ns, G, Z, R, n_steps_per_year)
all.equal(N[1, ], N[n_years * n_steps_per_year + 1, ])

# Random fishing mortality scaling factor
sigma_F0 <- 0.1
F0 <- exp(rnorm(n_years, sd = sigma_F0) - 0.5 * sigma_F0^2)
plot(years, F0, type = "b", xlab = "Year", ylab = "Fishing mortality scaling factor")

# Define separate parameters for AR(1) processes in time and bins
rho_t_g <- 0.5   # Temporal correlation for growth rate errors
rho_s_g <- 0.3   # Spatial correlation for growth rate errors
sigma_t_g <- 0.4 # Standard deviation for temporal errors in growth
sigma_s_g <- 0.4 # Standard deviation for spatial errors in growth

rho_t_m <- 0.8   # Temporal correlation for mortality rate errors
rho_s_m <- 0.2   # Spatial correlation for mortality rate errors
sigma_t_m <- 0.4 # Standard deviation for temporal errors in mortality
sigma_s_m <- 0.2 # Standard deviation for spatial errors in mortality

# Generate AR(1) errors for growth and mortality rates
epsilon_g <- generate_separable_AR1(n_years, n_bins, rho_t_g, rho_s_g, sigma_t_g, sigma_s_g)
epsilon_m <- generate_separable_AR1(n_years, n_bins, rho_t_m, rho_s_m, sigma_t_m, sigma_s_m)

# Rate matrices
G <- sweep(exp(epsilon_g), 2, g_given, "*")
plot(bin_start, G[1, ], type = "l", log = "x", xlab = "Weight", ylab = "Growth rate",
     main = "Growth rate for the first year")
lines(bin_start, G[2, ], col = "blue")
lines(bin_start, g_given, col = "red")
FMort <- outer(F0, f)
Z <- sweep(exp(epsilon_m), 2, m_given, "*") + FMort
plot(bin_start, Z[1, ], type = "l", log = "xy", xlab = "Weight", ylab = "Total mortality rate",
     main = "Total mortality rate for the first year")
lines(bin_start, Z[2, ], col = "blue")
lines(bin_start, m_given + FMort[1, ], col = "red")

# Random recruitment
sigma_R <- 0.2 # Standard deviation for recruitment
R <- Rs * exp(rnorm(n_years, sd = sigma_R) - 0.5 * sigma_R^2)
# Plot recruitment over time
plot(years, R, type = "b", xlab = "Year", ylab = "Recruitment")

# Calculate population matrix
N <- gpt_project(Ns, G, Z, R, n_steps_per_year)

# Animate N over time
animation <- gpt_animate_N(N, years, bin_start)
animation
# anim_save("animation.gif", animation)
# Let's animate the deviation from steady state instead
N_factor <- sweep(N, 2, Ns, "/")
animation <- gpt_animate_N(N_factor, years, bin_start, logy = FALSE)
animation

# Calculate total catches for each year and size bin
C_y <- matrix(0, nrow = n_years, ncol = n_bins)
for (y in 1:n_years) {
    for (t in 1:n_steps_per_year) {
        i <- (y - 1) * n_steps_per_year + t
        # Calculate catches
        C_i <- N[i, ] * FMort[y, ] * delta_t
        C_y[y, ] <- C_y[y, ] + C_i
    }
}

# Simulate sample sizes for each year
n_sample <- rpois(n_years, lambda = 10000)  # Average sample size of 1000 fish per year

# Simulate counts for each year and size bin based on predicted proportions
simulated_counts <- matrix(NA, nrow = n_years, ncol = n_bins)
for (y in 1:n_years) {
    probs <- C_y[y, ] / sum(C_y[y, ])
    simulated_counts[y, ] <- rmultinom(1, size = n_sample[y], prob = probs)
}
# We now put the simulated counts into the data frame that real data will come in
landings <- data.frame(
    year = rep(years, each = n_bins),
    bin_start = rep(bin_start, times = n_years),
    count = as.vector(t(simulated_counts))
)
y <- 5
plot(bin_start, simulated_counts[y, ], type = "b", log = "xy",
     xlab = "Weight", ylab = "Simulated Counts",
     main = paste("Year", years[y]))

# Simulate yield for each year
sigma_yield <- 0.1
log_sigma_yield <- log(sigma_yield)
yield <- rowSums(C_y) * exp(rnorm(n_years, mean = 0, sd = sigma_yield))
# Plot yield against years
plot(years, yield, type = "b", xlab = "Year", ylab = "Yield", main = "Simulated Yield Over Time")
lines(years, rowSums(C_y), col = "blue")

# Prepare data for TMB
data_list <- list(
    years = unique(landings$year),
    bin_start = unique(landings$bin_start),
    bin_width = bin_width,
    obs_counts = simulated_counts, # Matrix of observed counts [n_years x n_bins]
    n_sample = n_sample,           # Vector of sample sizes [n_years]
    yield = yield,                 # Simulated total yield
    f = f,                         # Selectivity vector
    given_g = g_given,             # Given growth rates for each weight class
    given_m = m_given,             # Given natural mortality rates for each weight class
    n_steps_per_year = n_steps_per_year,
    Ns = Ns,                       # Steady state population density
    Rs = Rs                        # Steady state recruitment
)

# Define parameters and initial values for TMB
parameters <- list(
    log_N0_factor = rep(0, n_bins),
    log_F0 = rep(log(0.1), length(unique(landings$year))),
    log_R_factor = rep(0, length(unique(landings$year))),
    epsilon_g = array(0, dim = c(length(unique(landings$year)), n_bins)),
    epsilon_m = array(0, dim = c(length(unique(landings$year)), n_bins)),
    log_sigma_yield =  log(0.2),  # Initial guess for standard deviation of yield
    logit_rho_t_g = 0,  # Initial guess for temporal correlation for growth rate errors
    logit_rho_s_g = 0,  # Initial guess for spatial correlation for growth rate errors
    logit_rho_t_m = 0,  # Initial guess for temporal correlation for mortality rate errors
    logit_rho_s_m = 0,  # Initial guess for spatial correlation for mortality rate errors
    log_sigma_t_g = log(sigma_t_g),
    log_sigma_s_g = log(sigma_s_g),
    log_sigma_t_m = log(sigma_t_m),
    log_sigma_s_m = log(sigma_s_m)
)

# Define random effects
random_effects <- c("epsilon_g", "epsilon_m")

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
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 3))

# Report results
rep <- sdreport(obj)
summary(rep)

best_par <- obj$env$last.par.best

# Plot the estimates for N0 against size
estimated_N0_factor <- exp(best_par[names(best_par) == "log_N0_factor"])
plot(bin_start, estimated_N0_factor, type = "b", log = "x", xlab = "Weight", ylab = "Estimated Initial Abundance (N0)", main = "Estimated N0 vs Size")

# Plot the estimates for F0 against years
estimated_F0 <- exp(best_par[names(best_par) == "log_F0"])
plot(years, estimated_F0, type = "b", ylab = "Estimated Fishing Mortality Scaling Factor (F0)", xlab = "Year")
lines(years, F0, col = "blue")

# Plot the estimates for R against years
estimated_R <- exp(best_par[names(best_par) == "log_R_factor"]) * Rs
plot(years, estimated_R, type = "b", ylab = "Estimated Recruitment (R)", xlab = "Year")
lines(years, R, col = "blue")

