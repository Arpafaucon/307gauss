namespace bench_pc
{
    void compute_fixed_size_pc(int buffer_size, int const *buffer, int matrix_dim, double &core_time, double &total_time, int &is_correct);

    void compute_variable_size_pc(int buffer_size, int const *buffer, int matrix_dim, double &core_time, double &total_time, int &is_correct);

    void dataset_battery();
}