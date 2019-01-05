
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <algorithm>
#include <exception>
#include <string>
#include <utility>
#include <vector>
#include <iostream>

#include "main_modules.h"

#include "base_matrix.h"
#include "bp_simulation.h"
#include "code_generation.h"
#include "commons_portable.h"
#include "combinatorics.h"
#include "data_structures.h"
#include "decoders.h"
#include "settings.h"


typedef struct
{
	int ncol;
	int nrow;
	int size;
	int min_modulo;
	int *data;
} MARK_PARAM;

MARK_PARAM create_mark_param( int nrow, int ncol )
{
	MARK_PARAM mp;

	mp.ncol = ncol;
	mp.nrow = nrow;
	mp.min_modulo = 0;
	mp.size = 0;
	mp.data = (int*) malloc(ncol * nrow * sizeof(mp.data[0]));
	return mp;
}

void delete_mark_param( MARK_PARAM *mp )
{
	free(mp->data);
}

int get_mark_param( MARK_PARAM *mp, FILE *fp )
{
	int c;
	char line[128];
	int i;

	while( (c = fgetc( fp )) != EOF )
	{
		if( c == '{' )
			break;
	}

	if( c == EOF )
		return -1;

	fscanf( fp, "%120s", line );	// read data
	fscanf( fp, "%120s", line );	// read =
	fscanf( fp, "%120s", line );	// read array
	fscanf( fp, "%120s", line );	// read {

	i = 0;
	while( fscanf( fp, "%d", &c ) != 0 )
	{
		mp->data[i++] = c;
	}
	mp->size = i;

	fscanf( fp, "%120s", line );	// read }
	fscanf( fp, "%120s", line );	// read min_modulo
	fscanf( fp, "%120s", line );	// read =
	fscanf( fp, "%d", &mp->min_modulo );	// read data

	fscanf( fp, "%120s", line );	// read }

	return 0;
}

using std::string;
using std::max_element;
using std::vector;
using std::pair;
using std::make_pair;
using std::sort;
using std::find;

class bad_matrix_interruption_exception : public std::exception {};

static pair<string, string> generate_paths(string file) {
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

int main_simulation(int argc, char *argv[]) {
    int tailbite_length;        // former M
    vector< double > snrs;      // former SNRS
	vector< double > snrs_per_symbol;
    int num_bp_iterations;      // former maxiter
    int num_experiments;        // former Nexp
    int num_frame_errors;       // former Nerr
    vector< int > deg_indices;  // former degs

	int decoder_type;

	int modulation_type;     	// Modulation type:  0 - skip, 1 - QAM4, 2 - QAM16, 3 - QAM64, 4 - QAM256
	int permutation_type;		// Permutation type:
	int permutation_block;
	int permutation_inter;
	int punctured_blocks;

    int error_index;
    string error_name;
    double error_rate_threshold;
    double good_code_multiple;
	bool good_code = false;

    vector< vector< int > > column_weight_constraints; // triples <what; min; max> for configuration filtering

    vector< string > error_names;
    error_names.push_back("BER");
    error_names.push_back("FER");

    if (argc != 3 && argc != 4) {
        die("Expected arguments: <scenario file name> <result file name> [marking file]");
    }

    // 1. Reading scenario
    settings f_scenario = settings::from_file(argv[1]);
	int mark_flag = argc > 3 ? 1 : 0;

    {
        // 1.0. These were compile-time constants, but now they are configurable
        f_scenario.select("settings/random_seed").cast_to(initial_random_seed);

        // 1.1. Reading the necessary fields
        f_scenario.select("settings/target_tailbite_length").cast_to(tailbite_length);
        f_scenario.select("settings/snrs").cast_to(snrs);
        f_scenario.select("settings/max_bp_iterations").cast_to(num_bp_iterations);
        f_scenario.select("settings/num_codewords").cast_to(num_experiments);
        f_scenario.select("settings/error_blocks").cast_to(num_frame_errors);
        f_scenario.select("settings/decoder_type").cast_to(decoder_type);

        f_scenario.select("settings/error_minimization/name").cast_to(error_name);
        f_scenario.select("settings/error_minimization/threshold").cast_to(error_rate_threshold);
        f_scenario.select("settings/error_minimization/good_code_multiple").cast_to(good_code_multiple);

        vector<string>::const_iterator name_itr = std::find(error_names.begin(), error_names.end(), error_name);
        if (name_itr == error_names.end()) {
            die("Unknown error name: '%s'", error_name.c_str());
        }
        error_index = (int)(name_itr - error_names.begin());

		f_scenario.select("settings/modulation_type").cast_to(modulation_type);
		f_scenario.select("settings/permutation_type").cast_to(permutation_type);
		f_scenario.select("settings/permutation_block").cast_to(permutation_block);
		f_scenario.select("settings/permutation_inter").cast_to(permutation_inter);
		f_scenario.select("settings/punctured_blocks").cast_to(punctured_blocks);

        // 1.2. Making sure the correct random seed is written into the config.
        ensure_random_is_initialized();

        // 1.3. Augmenting the config with some useful details (actual random seed, decoder name).
        f_scenario.open("decoder_name").set(DEC_FULL_NAME[decoder_type]);
    }

    // 2. Writing the main output file. All other output is appended to another file, which is referenced from here.

    pair<string, string> file_with_codes = generate_paths(string(argv[1]));

	vector<settings> input_codes;
	settings f_code = settings::from_file(argv[1]);
	f_code.select("results").cast_to(input_codes);

    settings f_log;
    f_log.open("settings").set(f_scenario);

    // Best frame error rates for each SNR
    vector< double > best_errors(snrs.size(), error_rate_threshold);

    try
    {
        console_exception_hook m_hook('m', "skips current matrix", bad_matrix_interruption_exception());
        matrix< double > pr_err((int)snrs.size(), 2);
        matrix< int > result;
        int rows = -1, columns = -1;
        matrix<int> best_HM;

        int girth;
        int code_index;
        int matrix_index;
        int config_index;
        vector< int > column_weights;
        matrix<int> org_HM;
        matrix<int> current_HM;

        unsigned i_max = mark_flag ? 1 : (unsigned)input_codes.size();

        for (unsigned i = 0; i < i_max; ++i)
        {
            vector<settings> logs;
            input_codes[i].select("code").cast_to(org_HM);

            input_codes[i].select("column_weights").cast_to(column_weights);
            input_codes[i].select("config_index").cast_to(config_index);
            input_codes[i].select("simulation_logs").cast_to(logs);
            input_codes[i].select("code_index").cast_to(code_index);
            input_codes[i].select("matrix_index").cast_to(matrix_index);
            input_codes[i].select("girth").cast_to(girth);
            if (rows == -1)
			{
                rows = org_HM.n_rows();
                columns = org_HM.n_cols();
            }
            else if (rows != org_HM.n_rows() || columns != org_HM.n_cols())
            {
                std::cout << "Warning: unequal matrices in the input!" << std::endl;
                continue;
            }

            int s_max = mark_flag ? 1 : (int) snrs.size();
            for (int s = 0; s < s_max; ++s)
            {
//								settings mark = settings::from_file(argv[3]);
//								vector< int > marking_data;
//								int marking_min_modulo = 0;
                double tt1 = clock();
                int mark_num  = -1;
                int best_mark = -1;
                pair<double, double> bp_result_org, bp_result_best;

                FILE *fp = NULL;
                MARK_PARAM mark_param = create_mark_param(rows, columns);

                bp_result_best.first  = 1.0;
                bp_result_best.second = 1.0;

                if (mark_flag)
                {
                    if ((fp = fopen(argv[3], "rt")) == NULL)
                        mark_flag = 0;

					//				marking[0] = settings::from_file(argv[3]);
					//				marking[0].select("data").cast_to(marking_data);
					//				marking[0].select("min_modulo").cast_to(marking_min_modulo);

					//				marking[1].select("data").cast_to(marking_data);
					//				marking[1].select("min_modulo").cast_to(marking_min_modulo);
                }

                current_HM = org_HM;

                do
                {
                    reset_random(); // all codes are tested with same noise

                    if (mark_num == -1)
                        printf( "original matrix is being processed\r" );
                    else
                        printf( "marked matrix #%d is being processed\r", mark_num );

                    pair<double, double> bp_result = bp_simulation(
                            current_HM,
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
                            0
                    );
                    if (mark_num == -1)
                        bp_result_org = bp_result;

                    if (bp_result.first < 0)
                    {
                        pr_err(s, 0) = 1;
                        pr_err(s, 1) = 1;
                    }
                    else if (bp_result.first > 1)
                    {
                        f_log.open("BAD_ENCODING_CERTIFICATE").set(result);
                        f_log.to_file(argv[2]);
                        return 0;
                    }
                    else
                    {
                        pr_err(s, 0) = bp_result.first;
                        pr_err(s, 1) = bp_result.second;
                    }

                    good_code = false;
                    if (bp_result.second < bp_result_best.second * good_code_multiple)
                        good_code = true;

                    if (bp_result.second < bp_result_best.second)
                    {
                        best_mark = mark_num;
                        bp_result_best.second = bp_result.second;
						bp_result_best.first  = bp_result.first;
                    }

                    if (good_code)
                    {
                        double tt2 = clock();
                        double time_sec = (tt2 - tt1) / CLOCKS_PER_SEC;

                        vector<settings> snr_results;
                        for (int s = 0, s_max = (int) snrs.size(); s < s_max; ++s)
                        {
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
                        descriptor.open("code").set(current_HM);
                        descriptor.open("column_weights").set(column_weights);
                        descriptor.open("decoder_type").set(decoder_type);
                        descriptor.open("decoder_name").set(DEC_FULL_NAME[decoder_type]);
                        descriptor.open("girth").set(girth);
                        descriptor.open("code_index").set(code_index);
                        descriptor.open("config_index").set(config_index);
                        descriptor.open("matrix_index").set(matrix_index);
                        descriptor.open("simulation_logs").set(snr_results);
                        descriptor.open("time").set(time_sec);
                        descriptor.to_file(argv[2], true);
                    }

                    if (good_code)
                    {
                        printf("\n");
						printf("FER best:   %e, BER best:   %e, best_mark: %d\n", bp_result_best.second, bp_result_best.first, best_mark);
                        printf("FER org:    %e, BER org:    %e,\n",               bp_result_org.second,  bp_result_org.first);
                        printf("FER marked: %e, BER marked: %e,\n\n",             bp_result.second,      bp_result.first);
                    }

                    mark_num++;
                    if (mark_flag)
                    {
                        int res = get_mark_param(&mark_param, fp);
                        if (res < 0)
                            mark_flag = 0;
                        else
                        {
                            for (int i = 0, k = 0; i < mark_param.ncol && k < mark_param.size; i++)
                            {
                                for (int j = 0; j < mark_param.nrow && k < mark_param.size; j++)
                                {
                                    if (current_HM(j, i) != -1)
                                        current_HM(j, i) = mark_param.data[k++];
                                }
                            }
                        }
                    }
                }
                while (mark_flag);

				if (mark_flag)
				    fclose(fp);
            }
        }
    }
    catch (bad_matrix_interruption_exception &)
    {
        printf("Current matrix has been skipped\n");
    }

    return 0;
}
