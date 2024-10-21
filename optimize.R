library(mizer)
source("prepareLogLikelihoodFn.R")
source("objective_function.R")

# Prepare the model
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
# We use lengths instead of weights in the catch data
lengths <- (w(p) / sp$a)^(1/sp$b)

# Observed biomass
sp$biomass_observed <- 0.713

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
    yield_observed = 0.465
)
gear_params(p) <- gp
initial_effort(p) <- 1

p <- steadySingleSpecies(p)

# Load catch size distribution
catch <- readRDS("catch.rds")
# Prepare the log-likelihood function
log_likelihood_fn <- prepare_log_likelihood(catch)

# Initial parameter estimates
initial_params <- c(l50 = gp$l50, ratio = gp$l25 / gp$l50, M = 0, U = 10)

# Set parameter bounds (if necessary)
lower_bounds <- c(l50 = 5, ratio = 0.1, M = 0, U = 1)
upper_bounds <- c(l50 = Inf, ratio = 0.99, M = Inf, U = 20)

# Perform the optimization
optim_result <- optim(
    par = initial_params,
    fn = objective_function,
    method = "L-BFGS-B",  # Allows for parameter bounds
    lower = lower_bounds,
    upper = upper_bounds,
    control = list(fnscale = 1, maxit = 100)
)

# Extract the optimal parameters
optimal_params <- optim_result$par
print(optimal_params)

# Optimal negative log-likelihood value
optimal_neg_log_likelihood <- optim_result$value
print(optimal_neg_log_likelihood)

# Plot observed counts
observed_data <- catch
hist_data <- data.frame(
    length = observed_data$length + observed_data$dl / 2,
    count = observed_data$count
)

plot(hist_data$length, hist_data$count, type = 'h', lwd = 2, col = 'blue',
     xlab = 'Length', ylab = 'Count', main = 'Observed Data and Fitted PDF')

# Generate fitted PDF with optimal parameters
pdf_values <- catch_pdf(optimal_params)

# Scale the PDF for visualization
pdf_values_scaled <- pdf_values * max(hist_data$count) / max(pdf_values)

# Add the fitted PDF to the plot
lines(lengths, pdf_values_scaled, col = 'red', lwd = 2)
legend('topright', legend = c('Observed Counts', 'Fitted PDF'),
       col = c('blue', 'red'), lwd = 2)

