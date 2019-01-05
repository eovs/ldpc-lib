#ifndef __BASE_MATRIX_H
#define __BASE_MATRIX_H

#include <vector>
#include "data_structures.h"
#include "settings.h"

struct base_matrix {
    matrix< int > data;
    int score[4];
    bool solution_exists;
};

base_matrix random_base_matrix(
    int num_check_symbols,
    int code_max_length,
    std::vector< int > const &information_column_weights,
    int num_trials_base,
    int num_trials_column,
    settings const &checking_column_settings,
    settings const &row_constraints
);

#endif
