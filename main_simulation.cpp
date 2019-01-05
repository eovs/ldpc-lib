
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

#define DEFAULT_MARKING "skip"
#define TARGET_ERR 0.01

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
    int num_iterations;			// former maxiter
    int num_experiments;        // former Nexp
    int num_frame_errors;       // former Nerr

	int decoder_type;

	int modulation_type;     	// Modulation type:  0 - skip, 1 - QAM4, 2 - QAM16, 3 - QAM64, 4 - QAM256
	int permutation_type;		// Permutation type:
	int permutation_block;
	int permutation_inter;
	int punctured_blocks;
	char mark_file[512];

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
	int mark_flag;

    {
        // 1.0. These were compile-time constants, but now they are configurable
        f_scenario.select("settings/random_seed").cast_to(initial_random_seed);

        // 1.1. Reading the necessary fields
        f_scenario.select("settings/num_codewords").cast_to(num_experiments);
        f_scenario.select("settings/error_blocks").cast_to(num_frame_errors);

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

        // 1.2. Making sure the correct random seed is written into the config.
        ensure_random_is_initialized();

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
        matrix< int > result;
        int rows = -1, columns = -1;

        int girth;
        int code_index;
        int matrix_index;
        int config_index;
        vector< int > column_weights;
        matrix<int> org_HM;
        matrix<int> current_HM;

        int max_code_idx = (unsigned)input_codes.size();
    
		for( int code_idx = 0; code_idx < max_code_idx; code_idx++ )
        {
            vector<settings> logs;
			vector< double > snrs;      // former SNRS
			int best_snr_idx;
			double best_err = 1.0;
            int best_mark = -1;

			input_codes[code_idx].select("_SNRs").cast_to(snrs);
			int s_max = (int) snrs.size();
			vector< double > best_errors(snrs.size(), error_rate_threshold);
			vector< double > EsN0(s_max);
			vector< double > BER(s_max);
			vector< double > FER(s_max); 
			vector< double > best_BER(s_max);
			vector< double > best_FER(s_max); 
			vector< double > org_BER(s_max);
			vector< double > org_FER(s_max); 

			best_snr_idx = s_max;

			input_codes[code_idx].select("_decoder_type").cast_to(decoder_type);
			input_codes[code_idx].select("_lifting").cast_to(tailbite_length);
			input_codes[code_idx].select("_punctured_blocks").cast_to(punctured_blocks);
			input_codes[code_idx].select("_iterations").cast_to(num_iterations);
			input_codes[code_idx].select("_marking").cast_to(mark_file);
	        
			
			input_codes[code_idx].select("code").cast_to(org_HM);

            input_codes[code_idx].select("column_weights").cast_to(column_weights);
            input_codes[code_idx].select("config_index").cast_to(config_index);
            input_codes[code_idx].select("simulation_logs").cast_to(logs);
            input_codes[code_idx].select("code_index").cast_to(code_index);
            input_codes[code_idx].select("matrix_index").cast_to(matrix_index);
            input_codes[code_idx].select("girth").cast_to(girth);
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

			double bitrate = (double)(columns - rows) / (columns - punctured_blocks);


//          double tt1 = clock();
            int mark_num  = -1;
			double curr_BER, curr_FER;

            FILE *fp = NULL;
            MARK_PARAM mark_param = create_mark_param(rows, columns);

            best_err = 1.0;
			mark_flag = !(strcmp(mark_file, DEFAULT_MARKING) == 0);

            if (mark_flag)
            {
                if ((fp = fopen(mark_file, "rt")) == NULL)
                    mark_flag = 0;
            }

#if 01
			for( int i = 0; i < rows; i++ ) 
			{
				for( int j = 0; j < columns; j++ ) 
				{
					if( org_HM(i, j) > 0 ) 
					{
						int t = org_HM(i, j) % tailbite_length;
						if( j == rows - 1 ) 
						{ // special last column int bidiagonal part of matrix
							t = t == 0 ? 1 : t;
						}
						org_HM(i, j) = t;
					}
				}
			}
#endif

            current_HM = org_HM;


			// Counting row weights - to get compability with ggp
			vector< int > row_weights;
			for (int i = 0; i < rows; ++i) {
				unsigned count = 0;
				for (int j = 0; j < columns; ++j) {
					count += org_HM(i, j) != -1;
				}
				while (row_weights.size() <= count) {
					row_weights.push_back(0);
				}
				++row_weights[count];
			}

			vector< vector< int > > row_weights_write;
			for (unsigned i = 0, i_max = row_weights.size(); i < i_max; ++i) {
				if (row_weights[i] != 0) {
					row_weights_write.push_back(vector< int >());
					row_weights_write.back().push_back(i);
					row_weights_write.back().push_back(row_weights[i]);
				}
			}
			

            do
            {
				int s;

				for( s = 0; s < s_max; s++ )
				{
					pair<double, double> result;

					if (mark_num == -1)
					{
						if( s == 0 )
							printf("====================================================\n");
						printf( "code #%d, original matrix is being processed, SNR = %6.3f\r", code_idx, snrs[s] );
					}
					else
						printf( "code #%d, marked matrix #%d is being processed, SNR = %6.3f\r", code_idx, mark_num, snrs[s] );

					EsN0[s] = snrs[s] + 10.0 * log10( 2.0 * bitrate );

					reset_random(); // all codes are tested with same noise


	                result = bp_simulation(
                            current_HM,
                            tailbite_length,
                            num_iterations,
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

					if( result.first > 1.0 || result.second > 1.0 )
					{
						f_log.open("BAD_ENCODING_CERTIFICATE").set(current_HM);
						f_log.to_file(argv[2]);
						return 0;
					}

					if( result.first < 0 || result.second < 0 )
					{
						curr_BER = 1.0;
						curr_FER = 1.0;
					}
					else
					{
						curr_BER = result.first;
						curr_FER = result.second;
					}

					BER[s] = curr_BER;
					FER[s] = curr_FER;
				}


				for( s = 0; s < s_max; s++ )
				{
					if( FER[s] < TARGET_ERR )					
						break;
				}

				if( s == s_max )
					s--;


				good_code = false;
				if( (s < best_snr_idx) || (s == best_snr_idx && FER[s] < best_err * good_code_multiple ) )
					good_code = true;


				if( (s < best_snr_idx) || (s == best_snr_idx && FER[s] < best_err) )
				{
					best_mark = mark_num;
					best_snr_idx = s;
					best_err = FER[s];
					best_FER = FER;
					best_BER = BER;
				}

				if( mark_num == -1 )
				{
					org_FER = FER;
					org_BER = BER;
				}

				if (good_code )
				{
//				double tt2 = clock();
//					double time_sec = (tt2 - tt1) / CLOCKS_PER_SEC;

					vector<settings> snr_results;

					settings snr_curr;
					snr_curr.open("SNR_per_bit___").set(snrs);
					snr_curr.open("SNR_per_symbol").set( EsN0 );
					snr_curr.open("BER").set(BER);
					snr_curr.open("FER").set(FER);
					
					snr_results.push_back(snr_curr);

					settings descriptor;
					descriptor.open("_decoder_name").set(DEC_FULL_NAME[decoder_type]);
					descriptor.open("_decoder_type").set(decoder_type);
					descriptor.open("_lifting").set(tailbite_length);
					descriptor.open("_SNRs").set(snrs);
					descriptor.open("_punctured_blocks").set(punctured_blocks);
					descriptor.open("_iterations").set(num_iterations);
					descriptor.open("_marking").set(DEFAULT_MARKING);
					descriptor.open("code_bitrate").set(bitrate);
					descriptor.open("code").set(current_HM);
					descriptor.open("column_weights").set(column_weights);
					descriptor.open("row_weights").set(row_weights_write);
					descriptor.open("girth").set(girth);
					descriptor.open("config_index").set(config_index);
					descriptor.open("matrix_index").set(matrix_index);
					descriptor.open("code_index").set(code_index);
					descriptor.open("simulation_logs").set(snr_results);
					//						descriptor.open("time").set(time_sec);
					descriptor.to_file(argv[2], true);
				}


                if (good_code)
                {
                        printf("\nbest mark: %d\n", best_mark);
						printf("        |        FER        |       BER \n"  );
						printf("  SNR   |    orig     best  |   orig     best\n"  );
						for( int s = 0; s < s_max; s++ )
						{
							printf("%7.3f | ", snrs[s] );
							printf("%8.6f %8.6f | ", org_FER[s], best_FER[s]);
							printf("%8.6f %8.6f", org_BER[s], best_BER[s]);
							printf("\n");
						}
						printf("\n");
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

			if( fp != NULL )
			    fclose(fp);
        }
		printf("\n");
    }
    catch (bad_matrix_interruption_exception &)
    {
        printf("Current matrix has been skipped\n");
    }

    return 0;
}
