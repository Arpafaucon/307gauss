#include "gauss.h"

namespace bench_pc
{

    void fill_matrices_fixed(int const* buffer, idx_t matrix_dim, float_t *exmat, float_t* inmat);
    void compute_fixed_size_pc(int buffer_size, int const *buffer, idx_t matrix_dim, double &core_time, double &total_time, int &is_correct);

    void compute_variable_size_pc(int buffer_size, int const *buffer, idx_t matrix_dim, double &core_time, double &total_time, int &is_correct);

    void dataset_battery();
}