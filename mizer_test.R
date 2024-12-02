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

# Parameters (empty because we're not estimating any parameters here)
parameters <- list()

# Compile the C++ code
compile("mizer.cpp")
dyn.load(dynlib("mizer"))

# Create the TMB object
obj <- MakeADFun(data = data, parameters = parameters, DLL = "mizer", silent = TRUE)

# Run the TMB model to get the outputs
result <- obj$report()

# Extract the outputs from TMB
feeding_level_Cpp <- result$feeding_level_Cpp
e_Cpp <- result$e_Cpp
e_repro_Cpp <- result$e_repro_Cpp
e_growth_Cpp <- result$e_growth_Cpp
f_mort_Cpp <- result$f_mort_Cpp
mort_Cpp <- result$mort_Cpp
rdd_Cpp <- result$rdd_Cpp

# Compare the outputs
all.equal(feeding_level_R, feeding_level_Cpp, tolerance = 1e-6)
all.equal(e_R, e_Cpp, tolerance = 1e-6)
all.equal(e_repro_R, e_repro_Cpp, tolerance = 1e-6)
all.equal(e_growth_R, e_growth_Cpp, tolerance = 1e-6)
all.equal(f_mort_R, f_mort_Cpp, tolerance = 1e-6)
all.equal(mort_R, mort_Cpp, tolerance = 1e-6)
all.equal(rdd_R, as.vector(rdd_Cpp), tolerance = 1e-6)

# Clean up temporary files
dyn.unload(dynlib(cpp_file))
file.remove(cpp_file)
file.remove(paste0(tools::file_path_sans_ext(cpp_file), ".o"))
file.remove(paste0(tools::file_path_sans_ext(cpp_file), ".dll"))  # On Windows
file.remove(paste0(tools::file_path_sans_ext(cpp_file), ".so"))   # On Unix
