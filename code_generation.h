#ifndef __CODE_GENERATION_H
#define __CODE_GENERATION_H

#include <vector>

#include "data_structures.h"
#include "equations.h"

bool generate_code(
    matrix< int > const &code_matrix,
    int target_girth,
    int max_module,
    int max_codes_to_test,
    matrix< int > &result_matrix
);

#endif
