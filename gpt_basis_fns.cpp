// Include necessary headers
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
    // Data inputs
    DATA_VECTOR(years);
    DATA_VECTOR(bin_start);
    DATA_VECTOR(bin_width);
    DATA_MATRIX(obs_counts);    // Observed counts [n_years x n_bins]
    DATA_VECTOR(n_sample);      // Sample sizes [n_years]
    DATA_VECTOR(yield);         // Yield data [n_years]
    DATA_VECTOR(f);             // Selectivity vector [n_bins]
    DATA_VECTOR(given_g);
    DATA_VECTOR(given_m);
    DATA_INTEGER(n_steps_per_year);
    DATA_VECTOR(Ns);            // Steady-state abundances [n_bins]
    DATA_SCALAR(Rs);            // Steady-state recruitment [n_years]

    // Basis functions for time and size
    DATA_MATRIX(B_time);        // Time basis matrix [n_years x K_time]
    DATA_MATRIX(B_size);        // Size basis matrix [n_bins x K_size]

    // Parameters and random effects
    PARAMETER_VECTOR(log_N0_factor);    // Log of deviation factor from steady state abundance
    PARAMETER_VECTOR(log_F0);           // Log of fishing mortality scaling factor for each year
    PARAMETER_VECTOR(log_R_factor);     // Log of deviation factor from steady state recruitment
    PARAMETER(log_sigma_R_factor);      // Add parameter for standard deviation
    PARAMETER(log_sigma_yield);         // Log of standard deviation for yield deviations

    // Coefficients for basis functions
    PARAMETER_MATRIX(coeff_g);          // Growth coefficients [K_time x K_size]
    PARAMETER_MATRIX(coeff_m);          // Mortality coefficients [K_time x K_size]
    // Standard deviations of coefficients
    PARAMETER(log_sigma_coeff_g);       // Log standard deviation for growth coefficients
    PARAMETER(log_sigma_coeff_m);       // Log standard deviation for mortality coefficients


    // Transform parameters
    vector<Type> N0 = exp(log_N0_factor) * Ns;  // Initial population abundances
    vector<Type> F0 = exp(log_F0);              // Fishing mortality scaling factor for each year
    Type sigma_R_factor = exp(log_sigma_R_factor);
    vector<Type> R = exp(log_R_factor) * Rs;    // Recruitment for each year
    Type sigma_yield = exp(log_sigma_yield);
    Type delta_t = Type(1.0) / n_steps_per_year;
    Type sigma_coeff_g = exp(log_sigma_coeff_g);
    Type sigma_coeff_m = exp(log_sigma_coeff_m);

    // Initialize negative log-likelihood
    Type nll = Type(0.0);

    // Number of years and bins
    int n_years = years.size();
    int n_bins = bin_start.size();

    // Check dimensions of basis matrices
    int K_time = B_time.cols();  // Number of time basis functions
    int K_size = B_size.cols();  // Number of size basis functions

    // Initialize population array N
    array<Type> N(n_years * n_steps_per_year + 1, n_bins);
    for (int j = 0; j < n_bins; ++j) {
        N(0, j) = N0(j);
    }

    // Calculate rate matrices
    array<Type> G(n_years, n_bins);
    array<Type> F(n_years, n_bins);
    array<Type> Z(n_years, n_bins);

    // Compute epsilon_g and epsilon_m using tensor product of basis functions
    // Initialize epsilon_g and epsilon_m
    array<Type> epsilon_g(n_years, n_bins);
    array<Type> epsilon_m(n_years, n_bins);

    for(int y = 0; y < n_years; ++y){
        for(int j = 0; j < n_bins; ++j){
            Type epsilon_g = Type(0.0);
            Type epsilon_m = Type(0.0);
            for(int k_time = 0; k_time < K_time; ++k_time){
                for(int k_size = 0; k_size < K_size; ++k_size){
                    epsilon_g += coeff_g(k_time, k_size) * B_time(y, k_time) * B_size(j, k_size);
                    epsilon_m += coeff_m(k_time, k_size) * B_time(y, k_time) * B_size(j, k_size);
                }
            }
            G(y, j) = exp(epsilon_g) * given_g(j);
            F(y, j) = f(j) * F0(y);
            Z(y, j) = exp(epsilon_m) * given_m(j) + F(y, j);
        }
    }

    // Project population over time using the dynamic model equation
    for(int y = 0; y < n_years; ++y) {
        int year_offset = y * n_steps_per_year;
        for(int t = 0; t < n_steps_per_year; ++t) {
            int i = year_offset + t;
            for(int j = 0; j < n_bins; ++j) {
                Type growth_contrib = (j == 0 ? R(y) : G(y, j - 1) * N(i, j - 1)) * delta_t / bin_width(j);
                Type denominator = Type(1.0) + (G(y, j) * delta_t / bin_width(j)) + (Z(y, j) * delta_t);
                N(i + 1, j) = (N(i, j) + growth_contrib) / denominator;
            }
        }
    }

    // Observation likelihoods
    for(int y = 0; y < n_years; ++y) {
        // Predicted catches in numbers
        vector<Type> C_y = vector<Type>::Zero(n_bins);
        for (int t = 0; t < n_steps_per_year; ++t) {
            int i = y * n_steps_per_year + t;
            for (int j = 0; j < n_bins; ++j) {
                C_y(j) += N(i, j) * F(y, j) * delta_t;
            }
        }
        // Calculate expected proportions
        Type total_C_y = C_y.sum();
        if (total_C_y > Type(0.0)) {
            vector<Type> p_y = C_y / total_C_y;
            // Observation likelihood: Multinomial distribution
            vector<Type> obs_counts_y = obs_counts.row(y).transpose();
            nll -= dmultinom(obs_counts_y, p_y, true);
        }

        // Yield likelihood
        Type predicted_yield_y = total_C_y;
        nll -= dnorm(log(yield(y)), log(predicted_yield_y), sigma_yield, true);
    }

    // Random effects likelihoods
    using namespace density;

    // Priors for coefficients_g and coefficients_m
    for(int i = 0; i < coeff_g.rows(); ++i){
        for(int j = 0; j < coeff_g.cols(); ++j){
            nll -= dnorm(coeff_g(i, j), Type(0.0), sigma_coeff_g, true);
            nll -= dnorm(coeff_m(i, j), Type(0.0), sigma_coeff_m, true);
        }
    }

    // Adjust the mean in the likelihood for log_R_factor
    nll -= sum(dnorm(log_R_factor, -0.5 * sigma_R_factor * sigma_R_factor, sigma_R_factor, true));

    // Report variables for diagnostics
    REPORT(N); // Predicted population abundances

    // Return negative log-likelihood
    return nll;
}
