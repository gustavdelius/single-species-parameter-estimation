install_github("sizespectrum/mizerExperimental")
library(mizerExperimental)
source("prepare_objective_function.R")
source("update_params.R")
source("cod_model.R")

# Prepare the model
p <- cod_model()

# Load cod catch size distribution
catch <- readRDS("cod_catch.rds")
# Prepare the objective function
objective_function <- prepare_objective_function(p, catch)

# Initial parameter estimates
gp <- gear_params(p)
initial_params <- c(l50 = gp$l50, ratio = gp$l25 / gp$l50, M = 0, U = 10,
                    catchability = gp$catchability)

# Set parameter bounds (if necessary)
lower_bounds <- c(l50 = 5, ratio = 0.1, M = 0, U = 1, catchability = 0)
upper_bounds <- c(l50 = Inf, ratio = 0.99, M = Inf, U = 20, catchability = Inf)

# Perform the optimization
optim_result <- optim(
    par = initial_params,
    fn = objective_function,
    method = "L-BFGS-B",  # Allows for parameter bounds
    lower = lower_bounds,
    upper = upper_bounds,
    control = list(fnscale = 1, maxit = 100,
                   trace = 1, REPORT = 1)
)

# Plot observed catches
hist_data <- data.frame(
    length = catch$length + catch$dl / 2,
    count = catch$count
)
plot(hist_data$length, hist_data$count, type = 'h', lwd = 2, col = 'blue',
     xlab = 'Length', ylab = 'Count', main = 'Observed Data and Fitted PDF')

# Generate fitted catch with optimal parameters
optimal_params <- update_params(p, optim_result$par)
model_catch <- optimal_params@initial_n * getFMort(optimal_params) * p@dw
# Scale the catch for visualization
model_catch <- model_catch * max(hist_data$count) / max(model_catch)
# We use lengths instead of weights in the catch data
lengths <- (p@w / p@species_params$a)^(1/p@species_params$b)

# Add the fitted PDF to the plot
lines(lengths, model_catch, col = 'red', lwd = 2)
legend('topright', legend = c('Observed Counts', 'Fitted PDF'),
       col = c('blue', 'red'), lwd = 2)

gp$yield_observed
getYield(optimal_params)


