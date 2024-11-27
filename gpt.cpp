// This code was created with ChatGPT, see
//  https://chatgpt.com/share/673b5cfa-c0dc-8007-a835-4bbd70e794c0

// Include necessary headers
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
    // Data inputs
    DATA_VECTOR(years);       // Years
    DATA_VECTOR(bin_start);   // Bin start points
    DATA_VECTOR(bin_width);   // Bin widths
    DATA_MATRIX(obs_counts);  // Observed counts [n_years x n_bins]
    DATA_VECTOR(n_sample);    // Sample sizes [n_years]
    DATA_VECTOR(yield);       // Yield data [n_years]
    DATA_VECTOR(f);           // Selectivity vector [n_bins]
    DATA_VECTOR(g);           // Growth rate vector [n_bins]
    DATA_VECTOR(m);           // Natural mortality rate vector [n_bins]
    DATA_INTEGER(n_steps_per_year);
    DATA_VECTOR(Ns);          // Steady-state abundances [n_bins]
    DATA_SCALAR(Rs);          // Steady-state recruitment [n_years]

    // Parameters and random effects
    PARAMETER_VECTOR(log_N0_factor);    // Log of deviation factor from steady state abundance
    PARAMETER_VECTOR(log_F0);           // Log of fishing mortality scaling factor for each year
    PARAMETER_VECTOR(log_R_factor);     // Log of deviation factor from steady state recruitment
    PARAMETER(log_sigma_R_factor);      // Add parameter for standard deviation
    PARAMETER_ARRAY(log_epsilon_N);     // Abundance errors [n_years x (n_bins - 1)]
    PARAMETER(logit_rho_N);             // Correlation for errors in log_epsilon_N
    PARAMETER(log_sigma_N);             // Log of sd for errors in log_epsilon_N
    PARAMETER(log_sigma_yield);         // Log of standard deviation for yield deviations

    // Transform parameters
    vector<Type> N0 = exp(log_N0_factor) * Ns;  // Initial population abundances
    vector<Type> F0 = exp(log_F0);              // Fishing mortality scaling factor for each year
    Type sigma_R_factor = exp(log_sigma_R_factor);
    vector<Type> R = exp(log_R_factor) * Rs;    // Recruitment for each year
    Type sigma_yield = exp(log_sigma_yield);
    Type rho_N = invlogit(logit_rho_N);
    Type sigma_N = exp(log_sigma_N);
    Type mu_N = - sigma_N * sigma_N / (1 + rho_N) / Type(2.0);
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

    // Calculate fishing mortality rate matrix
    array<Type> F(n_years, n_bins);
    for (int y = 0; y < n_years; ++y) {
        for (int j = 0; j < n_bins; ++j) {
            F(y, j) = f(j) * F0(y);
        }
    }

    // Project population over time using the dynamic model equation
    int year_start = 0;
    for(int y = 0; y < n_years; ++y) {
        for(int t = 0; t < n_steps_per_year; ++t) {
            int i = year_start + t;
            for(int j = 0; j < n_bins; ++j) {
                Type growth_contrib = (j == 0 ? R(y) : g(j - 1) * N(i + 1, j - 1)) * delta_t / bin_width(j);
                Type denominator = Type(1.0) + (g(j) * delta_t / bin_width(j)) + ((m(j) + F(y, j)) * delta_t);
                N(i + 1, j) = (N(i, j) + growth_contrib) / denominator;
            }
        }
        year_start += n_steps_per_year;
        for(int j = 1; j < n_bins; ++j) {
            N(year_start, j) = N(year_start, j) * exp(log_epsilon_N(y, j - 1));
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

    // likelihood for log_epsilon_N
    for (int y = 0; y < n_years; ++y) {
        nll -= dnorm(log_epsilon_N(y, 0), mu_N, sigma_N, true);
        for (int j = 1; j < n_bins; ++j) {
            nll -= dnorm(log_epsilon_N(y, j) - rho_N * log_epsilon_N(y, j - 1),
                         mu_N, sigma_N, true);
        }
    }

    // Likelihood for log_R_factor
    nll -= sum(dnorm(log_R_factor, -0.5 * sigma_R_factor * sigma_R_factor,
                     sigma_R_factor, true));

    // Report variables for diagnostics
    REPORT(N); // Predicted population abundances

    // Return negative log-likelihood
    return nll;
}
