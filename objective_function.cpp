#include <TMB.hpp>

// Declare the functions//

// Calculate fishing mortality rate
template<class Type>
vector<Type> calculate_F_mort(Type l50, Type ratio, Type catchability,
                              vector<Type> bin_boundaries);

// Calculate total mortality rate
template<class Type>
vector<Type> calculate_mort(vector<Type> F_mort, Type M,
                            vector<Type> bin_boundaries);

// Calculate growth rate
template<class Type>
vector<Type> calculate_growth(vector<Type>EReproAndGrowth, Type w_mat, Type U,
                              vector<Type> bin_boundaries);

// Calculate steady-state number density
template<class Type>
vector<Type> calculate_N(vector<Type> mort, vector<Type> growth, Type biomass,
                         vector<Type> bin_boundaries);

// Calculate model yield
template<class Type>
vector<Type> calculate_yield(vector<Type> N, vector<Type> F_mort,
                             vector<Type> bin_widths,
                             vector<Type> bin_boundaries);

// Calculate catch probabilities
template<class Type>
vector<Type> calculate_catch_probabilities(vector<Type> N, vector<Type> F_mort,
                                           vector<Type> bin_widths,
                                           vector<Type> bin_boundaries);


template<class Type>
Type objective_function<Type>::operator() ()
{
    // **Data Section**
    DATA_VECTOR(counts);           // Count observations in bins
    DATA_VECTOR(bin_widths);
    DATA_VECTOR(bin_boundaries);
    DATA_SCALAR(yield);            // Observed yield
    DATA_SCALAR(biomass);          // Observed biomass
    DATA_VECTOR(EReproAndGrowth);  // The rate at which energy is available for growth
                                   // and reproduction
    DATA_SCALAR(w_mat);
    DATA_SCALAR(yield_lambda);     // controls the strength of the penalty for
                                   // deviation from the observed yield.

    // **Parameter Section**
    PARAMETER(l50);
    PARAMETER(ratio);
    PARAMETER(M);
    PARAMETER(U);
    PARAMETER(catchability);

    // **Calculate fishing mortality rate**
    vector<Type> F_mort = calculate_F_mort(l50, ratio, catchability,
                                           bin_boundaries);

    // **Calculate total mortality rate**
    vector<Type> mort = calculate_mort(F_mort, M, bin_boundaries);

    // **Calculate growth rate**
    vector<Type> growth = calculate_growth(EReproAndGrowth, w_mat, U,
                                           bin_boundaries);

    // **Calculate steady-state number density**
    vector<Type> N = calculate_N(mort, growth, bin_boundaries);

    // **Calculate model yield**
    vector<Type> model_yield = calculate_yield(N, F_mort, bin_widths,
                                               bin_boundaries);

    // **Calculate catch probabilities**
    vector<Type> probs = calculate_catch_probabilities(N, F_mort, bin_widths,
                                                       bin_boundaries);

    // **Negative Log-Likelihood Calculation**
    // Compute the negative log-likelihood using the multinomial distribution
    Type nll = -dmultinom(counts, probs, true);

    // **Add penalty for deviation from observed yield**
    nll += yield_lambda * log(model_yield / yield).square()

    return nll;
}
