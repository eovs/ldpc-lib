#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <algorithm>
#include <exception>
#include <string>
#include <utility>
#include <vector>

#include "main_modules.h"

#include "base_matrix.h"
#include "bp_simulation.h"
#include "code_generation.h"
#include "commons_portable.h"
#include "combinatorics.h"
#include "data_structures.h"
#include "decoders.h"
#include "settings.h"

using std::string;
using std::max_element;
using std::vector;
using std::pair;
using std::make_pair;
using std::sort;
using std::find;

class bad_config_interruption_exception : public std::exception {};
class bad_matrix_interruption_exception : public std::exception {};

pair<string, string> generate_paths(string file) {
    if (file.rfind(".jsonx") == file.length() - 6) {
        file = file.substr(0, file.length() - 6);
    }
    file += "_codes.jsonx";
    string absolute_file = file;
    size_t last_slash = file.rfind('/');
    if (last_slash != string::npos) {
        file = file.substr(last_slash + 1);
    }
    last_slash = file.rfind('\\');
    if (last_slash != string::npos) {
        file = file.substr(last_slash + 1);
    }
    return make_pair(absolute_file, file);
}

int main_good_code_search(int argc, char *argv[]) {
    int num_check_symbols;      // former r
    int code_max_length;        // former n
    int target_girth;           // former g
    int tailbite_length;        // former M
    vector< double > snrs;      // former SNRS
    int num_trials_base;        // former trials[0]
    int num_trials_column;      // former trials[1]
    int num_random_codes;       // former N_codes
    int num_bp_iterations;      // former maxiter
    int num_experiments;        // former Nexp
    int num_frame_errors;       // former Nerr
    int num_start_config;       // former Jstart
    vector< int > deg_indices;  // former degs

	int decoder_type;

	int modulation_type;     	// Modulation type:  0 - skip, 1 - QAM4, 2 - QAM16, 3 - QAM64, 4 - QAM256 
	int permutation_type;		// Permutation type:    
	int permutation_block;
	int permutation_inter;
	int punctured_blocks;

	int codes_to_test;          // former constant CODES_TO_TEST
    int matrices_to_test;       // how much matrix to test for each configuration

    int error_index;
    string error_name;
    double error_rate_threshold;
    double good_code_multiple;
    vector< pair< int, double > > bad_code_multiples;

    vector< vector< int > > column_weight_constraints; // triples <what; min; max> for configuration filtering

    vector< string > error_names;
    error_names.push_back("BER");
    error_names.push_back("FER");

    // 0. Helper info about portability choices from commons_portable.h
    //commons_portable_print_info();

    if (argc != 3) {
        die("Expected two arguments: the scenario file name and the result file name");
    }

    // 1. Reading scenario
    settings f_scenario = settings::from_file(argv[1]);
    settings checking_columns;
    settings row_constraints;
    {
        // 1.0. These were compile-time constants, but now they are configurable
        f_scenario.select("random_seed").cast_to(initial_random_seed);
        f_scenario.select("codes_to_test").cast_to(codes_to_test);
        f_scenario.select("matrices_to_test").cast_to(matrices_to_test);

        // 1.1. Reading the necessary fields
        f_scenario.select("rows").cast_to(num_check_symbols);
        f_scenario.select("columns").cast_to(code_max_length);
        f_scenario.select("target_girth").cast_to(target_girth);
        f_scenario.select("target_tailbite_length").cast_to(tailbite_length);
        f_scenario.select("snrs").cast_to(snrs);
        f_scenario.select("trials/base_matrix").cast_to(num_trials_base);
        f_scenario.select("trials/column").cast_to(num_trials_column);
        f_scenario.select("num_random_codes").cast_to(num_random_codes);
        f_scenario.select("max_bp_iterations").cast_to(num_bp_iterations);
        f_scenario.select("num_codewords").cast_to(num_experiments);
        f_scenario.select("error_blocks").cast_to(num_frame_errors);
        f_scenario.select("start_config_number").cast_to(num_start_config);
        f_scenario.select("decoder_type").cast_to(decoder_type);

        f_scenario.select("error_minimization/name").cast_to(error_name);
        f_scenario.select("error_minimization/threshold").cast_to(error_rate_threshold);
        f_scenario.select("error_minimization/good_code_multiple").cast_to(good_code_multiple);

        vector<settings> bad_code_settings;
        f_scenario.select("error_minimization/bad_code_tuning").cast_to(bad_code_settings);
        for (unsigned i = 0; i < bad_code_settings.size(); ++i) {
            int codes;
            double multiple;
            bad_code_settings[i].select("codes").cast_to(codes);
            bad_code_settings[i].select("multiple").cast_to(multiple);
            bad_code_multiples.push_back(make_pair(codes, multiple));
        }
        sort(bad_code_multiples.begin(), bad_code_multiples.end());

        vector<string>::const_iterator name_itr = std::find(error_names.begin(), error_names.end(), error_name);
        if (name_itr == error_names.end()) {
            die("Unknown error name: '%s'", error_name.c_str());
        }
        error_index = (int)(name_itr - error_names.begin());

        settings const &info_config = f_scenario.select("information_columns");
        info_config.select("possible_weights").cast_to(deg_indices);
        if (info_config.can_select("weight_constraints")) {
            info_config.select("weight_constraints").cast_to(column_weight_constraints);
        }

        checking_columns = f_scenario.select("checking_columns");
        row_constraints = f_scenario.select("row_constraints");

		f_scenario.select("modulation_type").cast_to(modulation_type);
		f_scenario.select("permutation_type").cast_to(permutation_type);
		f_scenario.select("permutation_block").cast_to(permutation_block);
		f_scenario.select("permutation_inter").cast_to(permutation_inter);

		f_scenario.select("punctured_blocks").cast_to(punctured_blocks);


        // 1.2. Making sure the correct random seed is written into the config.
        ensure_random_is_initialized();

        // 1.3. Augmenting the config with some useful details (actual random seed, decoder name).
        f_scenario.open("random_seed").set(initial_random_seed);
        f_scenario.open("decoder_name").set(DEC_FULL_NAME[decoder_type]);
    }
    int max_nonzero_degree = *max_element(deg_indices.begin(), deg_indices.end());

    // 2. Writing the main output file. All other output is appended to another file, which is referenced from here.
    pair<string, string> file_with_codes = generate_paths(string(argv[2]));
    settings f_log;
    f_log.open("settings").set(f_scenario);
    f_log.open("results").set_array_ref(file_with_codes.second);
    f_log.to_file(argv[2]);

    // 3. Outer configuration loop
    addend_splitter degree_config_gen(code_max_length - num_check_symbols, (int)deg_indices.size());

    // Best frame error rates for each SNR
    vector< double > best_errors(snrs.size(), error_rate_threshold);

    for (long long config_index = degree_config_gen.n_splits() - 1 - num_start_config; config_index >= 0; --config_index) {
        // 4. Initializing column weights
        vector< int > column_weights(max_nonzero_degree + 1);
        {
            vector< int > degree_config;
            degree_config_gen.split(config_index, degree_config);
            for (unsigned i = 0, i_max = (unsigned int)deg_indices.size(); i < i_max; ++i) {
                column_weights[deg_indices[i]] += degree_config[i];
            }
            for (int i = 0; i <= max_nonzero_degree; ++i) {
                printf("%d ", column_weights[i]);
            }

            bool constraints_violated = false;
            for (unsigned w_index = 0, w_index_max = (unsigned int)deg_indices.size(); w_index < w_index_max; ++w_index) {
                for (unsigned c_index = 0, c_index_max = (unsigned int)column_weight_constraints.size(); c_index < c_index_max; ++c_index) {
                    vector< int > const &cic = column_weight_constraints[c_index];
                    if (cic[0] == deg_indices[w_index]) {
                        if (degree_config[w_index] < cic[1] || degree_config[w_index] > cic[2]) {
                            constraints_violated = true;
                        }
                    }
                }
            }
            if (constraints_violated) {
                printf("=> violates constraints, skipping\n");
                continue;
            }
            printf("\n");
        }

        try {
            console_exception_hook n_hook('n', "skips current configuration", bad_config_interruption_exception());
            for (int matrix_index = 0; matrix_index < matrices_to_test; ++matrix_index) {
                try {
                    console_exception_hook m_hook('m', "skips current matrix", bad_matrix_interruption_exception());
                    // 5. Building a random base matrix
                    base_matrix bm = random_base_matrix(
                            num_check_symbols,
                            code_max_length,
                            column_weights,
                            num_trials_base,
                            num_trials_column,
                            checking_columns,
                            row_constraints
                    );
                    if (!bm.solution_exists) {
                        printf("Base matrix was not found, continuing...\n");
                        continue;
                    }

                    // 6. Generating various codes and test then with different SNRs
                    for (int code_to_test = 0; code_to_test < codes_to_test; ++code_to_test) {
                        double bad_code_multiple = bad_code_multiples[0].second;
                        for (unsigned i = 0; i < bad_code_multiples.size(); ++i) {
                            if (bad_code_multiples[i].first <= code_to_test) {
                                bad_code_multiple = std::min(bad_code_multiple, bad_code_multiples[i].second);
                            }
                        }
                        // 6.0. Printing where we are.
                        printf("---------------------\nCurrent configuration:");
                        for (int i = 0; i <= max_nonzero_degree; ++i) {
                            printf(" %d", column_weights[i]);
                        }
                        printf("\nCurrent matrix: #%d\nCurrent code: #%d\n---------------------\n", matrix_index, code_to_test);

                        matrix< int > result;
                        int min_module = 0;

                        // 6.1. Testing various target girths to find the one which works.
                        int init_girth = bm.score[3] > 300 ? 6 : target_girth; // MB: why?
                        int gs;
                        for (gs = init_girth; gs > 4; gs -= 2) {
                            min_module = generate_code(bm.data, gs, tailbite_length, num_random_codes, result);
                            if (min_module > 0) {
                                break;
                            }
                            min_module = generate_code(bm.data, gs, 399, num_random_codes, result);
                            if (min_module > 0) {
                                // MB: not sure why this works, but the original has it.
                                for (int r = 0; r < result.n_rows(); ++r) {
                                    for (int c = 0; c < result.n_cols(); ++c) {
                                        result(r, c) %= tailbite_length;
                                    }
                                }
                                break;
                            }
                        }
                        if (min_module == 0) {
                            break;
                        }

                        // 6.2. Modeling the code we constructed.
                        bool good_code = false;
                        double tt1 = clock();
						reset_random();
                        matrix< double > pr_err((int)snrs.size(), 2);
                        for (int s = 0, s_max = (int) snrs.size(); s < s_max; ++s) {
                            printf("Modeling started for SNR = %lf\n", snrs[s]);
                            pair<double, double> bp_result = bp_simulation(
                                    result,
                                    tailbite_length,
                                    num_bp_iterations,
                                    num_frame_errors,
                                    num_experiments,
                                    snrs[s],
                                    best_errors[s],
									decoder_type,
									modulation_type,
									permutation_type,
									permutation_block,
									permutation_inter,
									punctured_blocks,
									1
                            );
                            printf("Modeling finished for SNR = %lf\n", snrs[s]);
                            if (bp_result.first < 0) {
                                pr_err(s, 0) = 1;
                                pr_err(s, 1) = 1;
                            } else if (bp_result.first > 1) {
                                f_log.open("BAD_ENCODING_CERTIFICATE").set(result);
                                f_log.to_file(argv[2]);
                                return 0;
                            } else {
                                pr_err(s, 0) = bp_result.first;
                                pr_err(s, 1) = bp_result.second;
                            }
                            printf("BEST_%s=%5.3le, code=%d, pr_err[0] = %lf, pr_err[1] = %lf\n",
                                error_name.c_str(), best_errors[s], code_to_test, pr_err(s, 0), pr_err(s, 1));
                            for (int i = 0; i <= max_nonzero_degree; ++i) {
                                printf("%d ", column_weights[i]);
                            }
                            printf("\n");

                            if (pr_err(s, error_index) > best_errors[s] * bad_code_multiple) {
                                good_code = false;
                                break;
                            }
                            if (pr_err(s, error_index) < best_errors[s] * good_code_multiple) {
                                good_code = true;
                            }
                            if (pr_err(s, error_index) < best_errors[s]) {
                                best_errors[s] = pr_err(s, error_index);
                            }
                        }
                        if (good_code) {
                            double tt2 = clock();
                            double time_sec = (tt2 - tt1) / CLOCKS_PER_SEC;

							int rows = result.n_rows();
							int columns = result.n_cols();

                            vector<settings> snr_results;
                            for (int s = 0, s_max = (int) snrs.size(); s < s_max; ++s) {
								settings snr_curr;
								snr_curr.open("SNR_per_bit___").set(snrs[s]);
								double bitrate = (double)(columns - rows) / (columns - punctured_blocks);
								double EsNo = snrs[s] + 10.0 * log10( 2.0 * bitrate );
								snr_curr.open("SNR_per_symbol").set( EsNo );
								snr_curr.open("BER").set(pr_err(s, 0));
								snr_curr.open("FER").set(pr_err(s, 1));
								snr_curr.open("punctured_blocks").set(punctured_blocks);
								snr_curr.open("bitrate").set(bitrate);
								snr_results.push_back(snr_curr);
                            }

                            settings descriptor;
                            descriptor.open("code").set(result);
                            descriptor.open("column_weights").set(column_weights);
                            descriptor.open("girth").set(gs);
                            descriptor.open("config_index").set(config_index);
                            descriptor.open("matrix_index").set(matrix_index);
                            descriptor.open("code_index").set(code_to_test);
                            descriptor.open("simulation_logs").set(snr_results);
                            descriptor.open("time").set(time_sec);
                            descriptor.to_file(file_with_codes.first, true);

                            printf("Good code written!\n");
                            printf("========>>  TIME = %lf\n", time_sec);
                        }
                    }
                } catch (bad_matrix_interruption_exception &) {
                    printf("Current matrix has been skipped\n");
                }
            }
        } catch (bad_config_interruption_exception &) {
            printf("Current configuration has been skipped\n");
        }
    }

    return 0;
}
