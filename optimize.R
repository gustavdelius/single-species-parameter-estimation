# In this script we calibrate the cod model to the observed catch size distribution
# and total yield. We use the L-BFGS-B optimization algorithm to find the parameter
# values that minimize the negative log-likelihood of the observed catch data.
# The objective function also includes a penalty for deviation from the observed
# total yield.

install_github("sizespectrum/mizerExperimental")
library(mizerExperimental)
library(TMB)
library(dplyr)
source("plot_catch.R")
source("cod_model.R")
source("update_params.R")
source("prepare_TMB_objective_function.R")

# Compile the model
compile("objective_function.cpp")
dyn.load(dynlib("objective_function"))

p <- readParams("celtic_params.rds")
sp <- species_params(p)
gp <- gear_params(p)

# Load cod catch size distribution
catch <- readRDS("celtic_catch.rds")
# Remove 4 very large fish
catch <- catch |>
    filter(length < 161)
# Fix incorrect maximum sizes
max_lengths <- catch |>
    filter(species != "Pilchard") |>
    group_by(species) |>
    summarize(max_length = max(length + dl))
sp <- left_join(sp, max_lengths, by = "species")
sp$max_weight <- l2w(sp$max_length, p)
sp |> select(species, w_max, max_weight)
wrong <- sp$w_repro_max < sp$max_weight
species_params(p)$w_max[wrong] <- sp$max_weight[wrong]
species_params(p)$w_repro_max[wrong] <- sp$max_weight[wrong]
sp <- species_params(p)
p <- steadySingleSpecies(p) |> matchGrowth() |> matchBiomasses()
plotlySpectra(p)

species <- valid_species_arg(p, 3)

sp_select <- sp$species == species
sps <- sp[sp_select, ]
gps <- gp[gp$species == species, ]

# Model does not fit the observed catch yet:
plot_catch(p, species, catch)
# and it has the wrong total yield:
gps$yield_observed
getYield(p)[sp_select]

# Set parameter bounds
lower_bounds <- c(l50 = 5, ratio = 0.1, catchability = 0, M = 0,
                  U = 1, w_repro_max = sps$w_mat)
upper_bounds <- c(l50 = Inf, ratio = 0.99, catchability = Inf, M = Inf,
                  U = 20, w_repro_max = Inf)

# Initial parameter estimates
initial_params <- c(l50 = gps$l50, ratio = gps$l25 / gps$l50,
                    catchability = gps$catchability, M = 0,
                    U = 10, w_repro_max = sps$w_repro_max)

# Prepare the objective function. See prepare_TMB_objective_function.R for details.
obj <- prepare_TMB_objective_function(p, species, catch, yield_lambda = 1e7,
                                      pars = initial_params)
# Perform the optimization. This starts with the initial parameter estimates and
# iteratively updates them to minimize the objective function.
optim_result <- nlminb(obj$par, obj$fn, obj$gr,
                       lower = lower_bounds, upper = upper_bounds,
                       control = list(trace = 1))

optim_result$par
report <- obj$report()

# Set model to use the optimal parameters
optimal_params <- update_params(p, species, optim_result$par)
# and plot the model catch again against the observed catch
plot_catch(optimal_params, species, catch)

# Also the yield is approximately matched:
gps$yield_observed
report$model_yield
getYield(optimal_params)[sp_select]
# If you want a better match you can increase the `yield_lambda` parameter

# Biomass is matched perfectly, by design
sps$biomass_observed
getBiomass(optimal_params)[sp_select]
