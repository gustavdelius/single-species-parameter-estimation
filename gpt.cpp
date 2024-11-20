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
    DATA_INTEGER(n_steps_per_year);

    // Parameters and random effects
    PARAMETER_VECTOR(log_N0);           // Log of initial population size in each bin
    PARAMETER_VECTOR(log_F0);           // Log of fishing mortality scaling factor for each year
    PARAMETER_VECTOR(log_R);            // Log of recruitment for each year
    PARAMETER_ARRAY(epsilon_g);         // Growth rate errors [n_years x n_bins]
    PARAMETER_ARRAY(epsilon_m);         // Mortality rate errors [n_years x n_bins]
    PARAMETER(log_sigma_g);             // Log of standard deviation for growth rate errors
    PARAMETER(log_sigma_m);             // Log of standard deviation for mortality rate errors
    PARAMETER(log_sigma_yield);         // Log of standard deviation for yield deviations
    PARAMETER(rho_t_g);                 // Temporal correlation for growth rate errors
    PARAMETER(rho_s_g);                 // Spatial correlation for growth rate errors
    PARAMETER(rho_t_m);                 // Temporal correlation for mortality rate errors
    PARAMETER(rho_s_m);                 // Spatial correlation for mortality rate errors
    PARAMETER(log_sigma_t_g);           // Log of temporal std for growth rate errors
    PARAMETER(log_sigma_s_g);           // Log of spatial std for growth rate errors
    PARAMETER(log_sigma_t_m);           // Log of temporal std for mortality rate errors
    PARAMETER(log_sigma_s_m);           // Log of spatial std for mortality rate errors

    // Transform parameters
    vector<Type> N0 = exp(log_N0);          // Initial population size in each bin
    vector<Type> F0 = exp(log_F0);          // Fishing mortality scaling factor for each year
    vector<Type> R = exp(log_R);            // Recruitment for each year
    Type sigma_yield = exp(log_sigma_yield);
    Type sigma_t_g = exp(log_sigma_t_g);
    Type sigma_s_g = exp(log_sigma_s_g);
    Type sigma_t_m = exp(log_sigma_t_m);
    Type sigma_s_m = exp(log_sigma_s_m);
    Type delta_t = Type(1.0) / n_steps_per_year;

    // Initialize negative log-likelihood
    Type nll = Type(0.0);

    // Number of years and bins
    int n_years = years.size();
    int n_bins = bin_start.size();

    // Initialize population array N
    array<Type> N(n_years * n_steps_per_year + 1, n_bins);
    for (int j = 0; j < n_bins; ++j) {
        N(0, j) = N0(j);
    }

    // Calculate rate matrices
    array<Type> G(n_years, n_bins);
    array<Type> F(n_years, n_bins);
    array<Type> Z(n_years, n_bins);
    for (int y = 0; y < n_years; ++y) {
        for (int j = 0; j < n_bins; ++j) {
            G(y, j) = exp(epsilon_g(y, j)) * given_g(j);
            F(y, j) = f(j) * F0(y);
            Z(y, j) = exp(epsilon_m(y, j)) * given_m(j) + F(y, j);
        }
    }

    // Project population over time using the dynamic model equation
    for(int y = 0; y < n_years; ++y) {
        int year_offset = y * n_steps_per_year;
        for(int t = 0; t < n_steps_per_year; ++t) {
            int i = year_offset + t;
            for(int j = 0; j < n_bins; ++j) {
                Type growth_contrib = (j == 0 ? R(y) : G(y, j - 1) * N(i + 1, j - 1)) * delta_t / bin_width(j);
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

    // Growth rate errors likelihood
    nll += SEPARABLE(
        SCALE(AR1(rho_t_g), sigma_t_g), // Time dimension
        SCALE(AR1(rho_s_g), sigma_s_g)  // Space dimension
        )(epsilon_g);

    // Mortality rate errors likelihood
    nll += SEPARABLE(
        SCALE(AR1(rho_t_m), sigma_t_m), // Time dimension
        SCALE(AR1(rho_s_m), sigma_s_m)  // Space dimension
        )(epsilon_m);

    // Return negative log-likelihood
    return nll;
}