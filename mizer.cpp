#include <TMB.hpp>

// mizerFeedingLevel function
template<class Type>
matrix<Type> mizerFeedingLevel(matrix<Type> &encounter,
                               matrix<Type> &intake_max) {
    matrix<Type> feeding_level = encounter.array() / (encounter.array() + intake_max.array());
    return feeding_level;
}

// mizerEReproAndGrowth function
template<class Type>
matrix<Type> mizerEReproAndGrowth(matrix<Type> &encounter,
                                  matrix<Type> &feeding_level,
                                  const vector<Type> &alpha,
                                  matrix<Type> &metab) {
    matrix<Type> e = ((Type(1.0) - feeding_level.array()) * encounter.array()).matrix();

    // Multiply each species (row) by alpha
    e.array().colwise() *= alpha.array();

    // Subtract metab
    e = e - metab;
    return e;
}

// mizerERepro function
template<class Type>
matrix<Type> mizerERepro(matrix<Type> &e,
                         matrix<Type> &psi) {
    // Set negative e values to zero
    matrix<Type> e_non_neg = e.cwiseMax(Type(0));
    // Multiply element-wise by psi
    matrix<Type> e_repro = e_non_neg.array() * psi.array();
    return e_repro;
}

// mizerEGrowth function
template<class Type>
matrix<Type> mizerEGrowth(matrix<Type> &e,
                          matrix<Type> &e_repro) {
    // Set negative e values to zero
    matrix<Type> e_non_neg = e.cwiseMax(Type(0));
    // energy for growth is intake - energy for reproduction
    matrix<Type> e_growth = e_non_neg - e_repro;
    return e_growth;
}

// mizerFMort function
template<class Type>
matrix<Type> mizerFMort(array<Type> &selectivity,
                        matrix<Type> &catchability,
                        const vector<Type> &effort) {
    int ngear = selectivity.dim[0];
    int nspecies = selectivity.dim[1];
    int nsize = selectivity.dim[2];

    matrix<Type> fmort(nspecies, nsize);
    fmort.setZero();

    // Sum over gears
    for (int g = 0; g < ngear; ++g) {
        Type eff = effort(g);
        for (int i = 0; i < nspecies; ++i) {
            Type catch_eff = eff * catchability(g, i);
            for (int w = 0; w < nsize; ++w) {
                fmort(i, w) += catch_eff * selectivity(g, i, w);
            }
        }
    }
    return fmort;
}

// mizerMort function
template<class Type>
matrix<Type> mizerMort(matrix<Type> &pred_mort,
                       matrix<Type> &mu_b,
                       matrix<Type> &f_mort) {
    matrix<Type> mort = pred_mort + mu_b + f_mort;
    return mort;
}

// BevertonHoltRDD function (unchanged)
template<class Type>
vector<Type> BevertonHoltRDD(const vector<Type> &rdi,
                             const vector<Type> &R_max) {
    vector<Type> rdd = rdi / (Type(1.0) + rdi / R_max);
    return rdd;
}

// Objective function required by TMB
template<class Type>
Type objective_function<Type>::operator() () {
    // Data input
    DATA_MATRIX(encounter);      // dimensions: species x size
    DATA_MATRIX(intake_max);     // dimensions: species x size
    DATA_VECTOR(alpha);          // length: species
    DATA_MATRIX(metab);          // dimensions: species x size
    DATA_MATRIX(psi);            // dimensions: species x size
    DATA_ARRAY(selectivity);     // dimensions: gear x species x size
    DATA_MATRIX(catchability);   // dimensions: gear x species
    DATA_VECTOR(effort);         // length: gear
    DATA_MATRIX(mu_b);           // dimensions: species x size
    DATA_MATRIX(pred_mort);      // dimensions: species x size
    DATA_VECTOR(R_max);          // length: species
    DATA_VECTOR(rdi);            // length: species

    PARAMETER_VECTOR(dummy);

    // Variables to store results
    matrix<Type> feeding_level = mizerFeedingLevel(encounter, intake_max);
    matrix<Type> e = mizerEReproAndGrowth(encounter, feeding_level, alpha, metab);
    matrix<Type> e_repro = mizerERepro(e, psi);
    matrix<Type> e_growth = mizerEGrowth(e, e_repro);
    matrix<Type> f_mort = mizerFMort(selectivity, catchability, effort);
    matrix<Type> mort = mizerMort(pred_mort, mu_b, f_mort);
    vector<Type> rdd = BevertonHoltRDD(rdi, R_max);

    // Report the results so they can be accessed from R
    REPORT(feeding_level);
    REPORT(e);
    REPORT(e_repro);
    REPORT(e_growth);
    REPORT(f_mort);
    REPORT(mort);
    REPORT(rdd);

    // No likelihood as we're just testing
    return dummy(0)*dummy(0)+dummy(1)*dummy(1);
}
