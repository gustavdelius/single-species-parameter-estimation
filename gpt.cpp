// This code was created with ChatGPT, see
//  https://chatgpt.com/share/673b5cfa-c0dc-8007-a835-4bbd70e794c0

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
    DATA_VECTOR(yield);
    DATA_VECTOR(f);             // Selectivity vector [n_bins]
    DATA_VECTOR(given_g);
    DATA_VECTOR(given_m);
    DATA_SCALAR(delta_t);

    // Parameters and random effects
    PARAMETER_VECTOR(log_N0);           // Log of initial population size in each bin
    PARAMETER_VECTOR(log_F0);           // Log of fishing mortality scaling factor for each year
    PARAMETER_VECTOR(log_R);            // Log of recruitment for each year
    PARAMETER_MATRIX(epsilon_g);        // Growth rate errors [n_years x n_bins]
    PARAMETER_MATRIX(epsilon_m);        // Mortality rate errors [n_years x n_bins]
    PARAMETER(log_sigma_g);             // Log of standard deviation for growth rate errors
    PARAMETER(log_sigma_m);             // Log of standard deviation for mortality rate errors
    PARAMETER(log_sigma_yield);         // Log of standard deviation for yield deviations
    PARAMETER(rho_t_g);                 // Temporal correlation for growth rate errors
    PARAMETER(rho_s_g);                 // Spatial correlation for growth rate errors
    PARAMETER(rho_t_m);                 // Temporal correlation for mortality rate errors
    PARAMETER(rho_s_m);                 // Spatial correlation for mortality rate errors

    // Transform parameters
    vector<Type> N0 = exp(log_N0);      // Initial population size in each bin
    vector<Type> F0 = exp(log_F0);      // Fishing mortality scaling factor for each year
    vector<Type> R = exp(log_R);        // Recruitment for each year
    Type sigma_g = exp(log_sigma_g);    // Standard deviation for growth rate errors
    Type sigma_m = exp(log_sigma_m);    // Standard deviation for mortality rate errors
    Type sigma_yield = exp(log_sigma_yield); // Standard deviation for yield deviations

    // Initialize negative log-likelihood
    Type nll = Type(0.0);

    // Process model
    int n_years = years.size();
    int n_bins = bin_start.size();
    int n_steps_per_year = 1 / delta_t;

    // Initialize population matrix
    matrix<Type> N(n_years * n_steps_per_year + 1, n_bins);
    N.row(0) = N0;

    // Calculate rate matrices
    matrix<Type> G(n_years, n_bins);
    matrix<Type> F(n_years, n_bins);
    matrix<Type> Z(n_years, n_bins);
    for(int y = 0; y < n_years; ++y) {
        G.row(y) = given_g + epsilon_g.row(y);
        F.row(y) = f * F0(y);
        Z.row(y) = F.row(y) + given_m + epsilon_m.row(y);
    }

    // Project population over time using the new dynamic model equation
    for(int y = 0; y < n_years; ++y) {
        int year_offset = y * n_steps_per_year;
        for(int t = 0; t < n_steps_per_year; ++t) {
            int i = year_offset + t;
            for(int j = 0; j < n_bins; ++j) {
                Type growth_contrib = ((j == 0) ? R(y) : (G(y, j - 1) * N(i + 1, j - 1))) * delta_t / bin_width(j);
                Type denominator = Type(1.0) + (G(y, j) * delta_t / bin_width(j)) + (Z(y, j) * delta_t);
                N(i + 1, j) = (N(i, j) + growth_contrib) / denominator;
            }
        }
    }

    // Observation likelihoods
    for(int y = 0; y < n_years; ++y) {
        // Predicted catches in numbers
        vector<Type> C_y(n_bins, Type(0.0));

        for(int t = 0; t < n_steps_per_year; ++t) {
            int i = y * n_steps_per_year + t;
            // Accumulate catches
            C_y += N.row(i) * F.row(y);
        }

        // Calculate expected proportions
        Type total_C_y = C_y.sum();
        vector<Type> p_y = C_y / total_C_y;

        // Observation likelihood: Multinomial distribution
        vector<Type> obs_counts_y = obs_counts.row(y);
        nll -= dmultinom(obs_counts_y, p_y, true);

        // Yield likelihood
        Type predicted_yield_y = C_y.sum();
        nll -= dnorm(yield(y), predicted_yield_y, sigma_yield, true);
    }

    // Random effects likelihoods
    using namespace density;
    AR1_t<Type> ar1_t_g(rho_t_g);
    AR1_t<Type> ar1_t_m(rho_t_m);
    nll += SCALE(ar1_t_g, sigma_g)(epsilon_g);
    nll += SCALE(ar1_t_m, sigma_m)(epsilon_m);

    // Return negative log-likelihood
    return nll;
}
