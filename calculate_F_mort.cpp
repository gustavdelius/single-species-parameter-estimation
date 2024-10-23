template<class Type>
vector<Type> calculate_F_mort(Type l50, Type ratio, Type catchability,
                              vector<Type> bin_boundaries)
{
    int size = bin_boundaries.size();
    vector<Type> F_mort(size);

    // **Calculate F_mort at each bin boundary**
    for (int i = 0; i < size; ++i) {
        // Example density calculation (replace with your actual model)
        Type x = bin_boundaries[i];
        densities[i] = 0; // TODO
    }

    return F_mort;
}
