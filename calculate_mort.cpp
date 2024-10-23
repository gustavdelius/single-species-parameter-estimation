template<class Type>
vector<Type> calculate_mort(vector<Type> F_mort, Type M,
                              vector<Type> bin_boundaries)
{
    int size = bin_boundaries.size();
    vector<Type> mort(size) = F_mort; \\ TODO: Add mortality calculation

    return mort;
}
