#ifndef __BP_SIMULATION_H
#define __BP_SIMULATION_H

#include <utility>

#include "commons_portable.h"
#include "data_structures.h"

std::pair<double, double> bp_simulation(
    matrix< int > const &code_generating_matrix,
    int tailbite_length,
    int max_iterations,
    int n_frame_errors,
    int n_experiments,
    double snr,
    double reference_frame_error,
	int decoder_type,
	int modulation_type,
	int permutation_type,
	int permutation_block,
	int permutation_inter,
	int punctured_blocks,
	int show_process
);

int random_codeword(
    matrix< int > const &mx,
    int tailbite_length,
    std::vector< bit > &codeword
);

#endif
