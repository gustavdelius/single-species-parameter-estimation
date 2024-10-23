# In this script we calibrate the cod model to the observed catch size distribution
# and total yield. We use the L-BFGS-B optimization algorithm to find the parameter
# values that minimize the negative log-likelihood of the observed catch data.
# The objective function also includes a penalty for deviation from the observed
# total yield.

install_github("sizespectrum/mizerExperimental")
library(mizerExperimental)
source("plot_catch.R")
source("cod_model.R")

use_TMB <- FALSE # Set to TRUE to use the TMB version of the optimization
use_TMB <- TRUE
if (use_TMB) {
    # This is the TMB version of the optimization. It is faster than the optim
    # version below, and can handle more complex models.
    library(TMB)

    # Compile the model
    compile("objective_function.cpp")
    dyn.load(dynlib("objective_function"))

    source("prepare_TMB_objective_function.R")
} else {
    source("prepare_objective_function.R")
    source("update_params.R")
}

# We demonstrate this with a single-species model for cod
p <- cod_model()

# Load cod catch size distribution
catch <- readRDS("cod_catch.rds")

# Model does not fit the observed catch yet:
plot_catch(p, catch)
# and it has the wrong total yield:
p@gear_params$yield_observed
getYield(p)

# Initial parameter estimates
gp <- gear_params(p)
initial_params <- c(l50 = gp$l50, ratio = gp$l25 / gp$l50, M = 0, U = 10,
                    catchability = gp$catchability)

# Set parameter bounds
lower_bounds <- c(l50 = 5, ratio = 0.1, M = 0, U = 1, catchability = 0)
upper_bounds <- c(l50 = Inf, ratio = 0.99, M = Inf, U = 20, catchability = Inf)

if (use_TMB) {
    # Prepare the objective function. See prepare_TMB_objective_function.R for details.
    obj <- prepare_TMB_objective_function(p, catch, yield_lambda = 1e7,
                                          pars = initial_params)
    # Perform the optimization. This starts with the initial parameter estimates and
    # iteratively updates them to minimize the objective function.
    optim_result <- nlminb(obj$par, obj$fn, obj$gr,
                           lower = lower_bounds, upper = upper_bounds,
                           control = list(trace = 1))
} else {
    # This version uses numerical differentiation. It is slower than the TMB
    # version above, but is easier to use and understand.
    # Prepare the objective function. See prepare_objective_function.R for details.
    objective_function <- prepare_objective_function(p, catch, yield_lambda = 1e7)

    # Perform the optimization. This starts with the initial parameter estimates and
    # iteratively updates them to minimize the objective function.
    optim_result <- optim(
        par = initial_params,
        fn = objective_function,
        method = "L-BFGS-B",  # Allows for parameter bounds
        lower = lower_bounds,
        upper = upper_bounds,
        control = list(fnscale = 1, maxit = 100,
                       # We print the value of the objective function at each iteration
                       trace = 1, REPORT = 1)
    )
    # After the last iteration there is a pause. That is normal. Be patient.
}

# Set model to use the optimal parameters
optimal_params <- update_params(p, optim_result$par)
# and plot the model catch again against the observed catch
plot_catch(optimal_params, catch)
# The fit is quite good.

# Also the yield is approximately matched:
gp$yield_observed
getYield(optimal_params)
# If you want a better match you can increase the `yield_lambda` parameter


