#include <TMB.hpp>

template<class Type>
vector<Type> calculate_catch_probabilities(vector<Type> N, vector<Type> F_mort,
                                           vector<Type> bin_widths,
                                           vector<Type> bin_boundaries)

    // **Calculate Densities at Bin Boundaries**
    vector<Type> densities = N * F_mort

    // **Compute Probabilities Using Trapezoidal Rule**
    int num_bins = bin_widths.size(); // Number of bins
    vector<Type> probs(num_bins);
    for (int i = 0; i < num_bins; ++i) {
        // Trapezoidal rule
        probs[i] = bin_widths[i] * (densities[i] + densities[i + 1]) / 2.0;
    }

    // **Ensure Probabilities Are Positive and Sum to 1**
    // Add a small epsilon to avoid log(0) and normalize
    Type epsilon = Type(1e-10);
    probs = probs + epsilon;
    probs = probs / probs.sum();

    return probs;
}
