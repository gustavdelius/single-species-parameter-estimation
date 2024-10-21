install_github("sizespectrum/mizerExperimental")
library(mizerExperimental)
source("prepare_objective_function.R")
source("update_params.R")
source("plot_catch.R")
source("cod_model.R")

# Prepare the model
p <- cod_model()

# Load cod catch size distribution
catch <- readRDS("cod_catch.rds")

# Model does not fit the observed catch yet:
plot_catch(p, catch)
# and it has the wrong total yield:
p@gear_params$yield_observed
getYield(p)

# Prepare the objective function
objective_function <- prepare_objective_function(p, catch)

# Initial parameter estimates
gp <- gear_params(p)
initial_params <- c(l50 = gp$l50, ratio = gp$l25 / gp$l50, M = 0, U = 10,
                    catchability = gp$catchability)

# Set parameter bounds
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

# Plot fit of calibrated model to catch
optimal_params <- update_params(p, optim_result$par)
plot_catch(optimal_params, catch)

# Also the yield is matched:
gp$yield_observed
getYield(optimal_params)


