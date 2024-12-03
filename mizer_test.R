# Load necessary libraries
library(TMB)
library(mizer)

params <- NS_params

# Prepare data for TMB
data <- list(
    encounter = getEncounter(params),
    intake_max = intake_max(params),
    feeding_level = getFeedingLevel(params),
    alpha = species_params(params)$alpha,
    metab = metab(params),
    psi = params@psi,
    selectivity = selectivity(params),
    catchability = catchability(params),
    effort = initial_effort(params),
    mu_b = ext_mort(params),
    pred_mort = getPredMort(params),
    R_max = species_params(params)$R_max,
    rdi = getRDI(params)
)

# Parameters (not of interest here but needed for TMB)
parameters <- list(dummy = rnorm(2))

# Compile the C++ code
compile("mizer.cpp")
dyn.load(dynlib("mizer"))

# Create the TMB object
obj <- MakeADFun(data = data, parameters = parameters, DLL = "mizer")

opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(trace = 3,
                             rel.tol = 1e-3))

# Run the TMB model to get the outputs
rep <- obj$report()


# Compare the outputs
all.equal(getFeedingLevel(params), rep$feeding_level,
          tolerance = 1e-6, check.attributes = FALSE)
all.equal(getEReproAndGrowth(params),rep$e,
          tolerance = 1e-6, check.attributes = FALSE)
all.equal(getERepro(params), rep$e_repro,
          tolerance = 1e-6, check.attributes = FALSE)
all.equal(getEGrowth(params), rep$e_growth,
          tolerance = 1e-6, check.attributes = FALSE)
all.equal(getFMort(params), rep$f_mort,
          tolerance = 1e-6, check.attributes = FALSE)
all.equal(getMort(params), rep$mort, tolerance = 1e-6, check.attributes = FALSE)
all.equal(getRDD(params), rep$rdd, tolerance = 1e-6, check.attributes = FALSE)

