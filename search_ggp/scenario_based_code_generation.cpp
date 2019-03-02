// This macro is defined by the includer, main_search_ggp.cpp.
// This is done to prevent a possibly dumb Visual Studio from compiling this source,
// which it would definitely fail to do.
#ifdef INTERNAL_SOURCE_INCLUDE

#ifdef WIN32
#define NOMINMAX
#include <Windows.h>
#endif	// WIN32

#include "trace_pm.h"

typedef struct
{
	int ncol;
	int nrow;
	int size;
	int min_modulo;
	int *data;
} MARK_PARAM;

typedef struct
{
	vector< int > a;
	int min_ok_modulo;
	double fer;
} MARK_CANDIDATE;

typedef struct
{
	vector< int > a;
	vector< int > coef;
	int min_ok_modulo;
}MARK_RECORD;


int random_codeword( matrix< int > const &mx, int tailbite_length, vector< bit > &codeword );
static int trace_matrix( matrix<int> current_HM, vector<int> min_ACE, vector<int> max_ACE_spec, int M, int gtarget, vector<int> &ACE, vector<int> &girth_spectrum );
static void show_matrix_property( vector<int> ACE, vector<int> girth_spectrum, int girth, int num );

class scenario_based_code_generation : public ggp_function {
private:
    string scenario_file;
    string output_file_mask;
    string extra_output_file;
    int min_number_width;
    int output_mask_extra_width;

    string output_file_name(int number) const {
        // should have been something like snprintf, but we have only C++03.
        int number_width = 0;
        for (int nn = number == 0 ? 1 : number; nn > 0; ++number_width, nn /= 10);
        int string_size = output_mask_extra_width + std::max(number_width, min_number_width);
        char *space = new char[string_size + 1];
        sprintf(space, output_file_mask.c_str(), number);
        string rv = space;
        delete[] space;
        return rv;
    }

private: // state moved from run
    int rows,               // r
        columns,            // n
        subcode_columns,    // ns
        target_girth,       // g
        subcode_girth,      // gs
        tailbite_length,    // M
        exact_tailbite,     // exact
        starting_length,    // nb
        min_codes,          // min_codes
        max_codes,          // max_codes
        num_cycles;         // N
	
	    vector< int > min_ACE;
		vector< int > max_SPEC;
		int use_ACE;
		int use_SPEC;
		int q_mod;
		matrix< int > HCM;        // HC modified
		vector< int > coef_mask;  // the mask for HC matrix 
		int all_col_2;            // all columns have weight equal to 2
		int bidiagonalFlag;		  // check matrix is bidiagonal

		matrix< int > HM;
		matrix< int > HM_orig;
		double snr_for_selection;
		double min_fer, max_fer;
		int candidate_cnt;
		vector< int > wcol, cumwcol, mask;

		int decoder_type;
		int number_of_iterations;
		int number_of_codewords;
		int number_of_err_blocks;
		vector< int > tested_code_lengths;
		double start_snr_for_estimation;
		double snr_step;
		double snr_bound;
		int number_of_selected_codes;


private:
    bool set_output_file_mask(string const &output_file_mask) {
        this->output_file_mask = output_file_mask;
        string::size_type percent = output_file_mask.find('%', 0);
        string::size_type next_percent = output_file_mask.find('%', percent + 1);
        string::size_type nextD = output_file_mask.find('d', percent + 1);
        if (percent == string::npos || nextD == string::npos || next_percent != string::npos) {
            cout << "Wrong format string in output file mask: '" << output_file_mask << "'" << endl;
            return false;
        }
        min_number_width = -1;
        stringstream(output_file_mask.substr(percent + 1, nextD - percent - 1)) >> min_number_width;
        if (min_number_width == -1) {
            cout << "Wrong format string in output file mask: '" << output_file_mask << "'" << endl;
            return false;
        }
        output_mask_extra_width = (int)(percent + output_file_mask.length() - nextD - 1);
        return true;
    }

	bool check_code_length( int len )
	{
		int size = (int)tested_code_lengths.size();
		for( int i = 0; i < size; i++ )
		{
			if( len == tested_code_lengths[i] )
				return true;
		}
		return false;
	}

	void improve_coefs( vector< int > &coefs, int bidiagonalFlag, int all_w2, int rows, int q_mod )
	{
		if( bidiagonalFlag )
		{
			for( int i = 0; i < rows-1; i++ )
				coefs[i*2+1] = coefs[i*2];
		}

		if( all_w2 )
		{
			if( coefs[rows*2] == coefs[rows*2+1] )
				coefs[rows*2+1] = (coefs[rows*2+1] + 1) % q_mod;
		}
	}


	int check_bidiagonal( matrix< int > matr )
	{
		int rows = matr.n_rows();
		int flag = 1;

		for( int i = 0; i < rows-1; i++ ) 
		{
			int w = 0;
			for( int j = 0; j < rows; j++ )
				w += matr(i, j) != -1;

			if( w != 2 )
			{
				bidiagonalFlag = 0;
				break;
			}

			if( matr(i,i) != 0 )
			{
				bidiagonalFlag = 0;
				break;
			}

			if( matr(i+1,i) != 0 )
			{
				bidiagonalFlag = 0;
				break;
			}
		}

		return flag;
	}

public:
    scenario_based_code_generation(string const &scenario_file) : scenario_file(scenario_file) {}

    virtual bool can_run(ggp_state const &) const {
        return true;
    }

    virtual void write_info(ggp_state const &, ostream &out) const {
        out << "Creating a code according to '" << get_native_path(scenario_file)
            << "'.";
    }

    bool single_iteration(int prev_length,
                          char const *interrupt_info_message,
                          char const *final_message_prefix,
                          char const *final_message_suffix,
                          int codes_tested_log_period,
                          int codes_selected_log_period,
                          int code_expansion_max_trials,
                          int my_min_codes) {

#ifdef WIN32
		HANDLE hStdout= GetStdHandle (STD_OUTPUT_HANDLE) ;
		SetConsoleTextAttribute (hStdout,FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED );
#endif
		int punctured_blocks = 0;
		int tarGirth = 0;
		int tarACEbyACE   = 0;
		int tarSPECbyACE  = 1000000;
		int tarACEbySPEC  = 0;
		int tarSPECbySPEC = 1000000;
        string prev_file_name = output_file_name(prev_length);
        ifstream prev_data(prev_file_name.c_str());
        if (prev_length != 0 && !prev_data) {
            cout << "Cannot open file '" << prev_file_name << "' for reading" << endl;
            return false;
        }
        string curr_file_name = output_file_name(starting_length);
        ofstream curr_data(curr_file_name.c_str());
        if (!curr_data) {
            cout << "Cannot open file '" << curr_file_name << "' for writing" << endl;
            return false;
        }

		if( prev_length == 0 )
		{
			char ans;
			cout << " file " << curr_file_name << " exist. Rewrite it ? ( 'y' or 'n' )\n";
			cin >> ans;
			if( ans == 'n' || ans == 'N' )
				return false;
		}


        if (prev_length) {
            cout << endl << "Extending " << (starting_length - 1) << " columns to " << starting_length << " columns" << endl;
        } else {
            cout << endl << "Creating random " << starting_length << " columns" << endl;
        }
		//========================================================
		vector < MARK_RECORD > vmark_record;
		vector < MARK_CANDIDATE > vmark_candidate;
		MARK_RECORD mark_record;
		MARK_CANDIDATE mark_candidate;
		//========================================================
        vector< int > nonzero_mapping;
        for (int i = starting_length - subcode_columns + 1; i <= starting_length; ++i) {
            for (int j = cumwcol[i - 1]; j < cumwcol[i]; ++j) {
                nonzero_mapping.push_back(j);
            }
        }

        // 1. Constructing "first equations"
        cout << endl
             << " n=" << starting_length
             << ", r=" << rows
             << ", g=" << target_girth
             << ", M=" << tailbite_length << endl << endl;

        // 1.1. Base matrix
        graph< vector_bag > first_equations_base;
        int first_equations_base_edges = 0;
        for (int i = 0; i < starting_length; ++i) {
            for (int j = 0; j < rows; ++j) {
                if (HM(j, i) > 0) {
                    vector_bag fw(first_equations_base_edges++);
                    first_equations_base.add_bidi_edge(i, j + starting_length, fw, -fw);
                }
            }
        }
        cout << " Tanner graph size: nt=" << first_equations_base.n_edges() / 2
             << ", rt=" << first_equations_base.n_vertices() << endl;
        cout << "Generating first equations" << endl;
        // 1.2. Equations
        equation_builder first_equations(first_equations_base, target_girth);
        if (first_equations.tree_girth() < target_girth) {
            cout << "problem: required girth " << target_girth
                 << ", found only " << first_equations.tree_girth() << endl;
            return false;
        }
        cout << "ok: g=" << first_equations.tree_girth()
             << ", eqs=" << first_equations.n_equations()
             << ", nodes=" << first_equations.n_vertices() << endl;

        // 2. Constructing "second equations"
        cout << endl
             << " n=" << starting_length
             << ", r=" << rows
             << ", g=" << subcode_girth
             << ", M=" << tailbite_length << endl << endl;
        // 2.1. Base matrix
        graph< vector_bag > second_equations_base;
        int second_equations_base_edges = 0;
        for (int i = starting_length - subcode_columns; i < starting_length; ++i) {
            for (int j = 0; j < rows; ++j) {
                if (HM(j, i) > 0) {
                    vector_bag fw(second_equations_base_edges++);
                    second_equations_base.add_bidi_edge(
                        i - starting_length + subcode_columns,
                        j + subcode_columns,
                        fw, -fw
                    );
                }
            }
        }
        cout << " Tanner graph size: nt=" << second_equations_base.n_edges() / 2
             << ", rt=" << second_equations_base.n_vertices() << endl;
        cout << "Generating second equations" << endl;
        // 2.2. Equations
        equation_builder second_equations(second_equations_base, subcode_girth);
        if (second_equations.tree_girth() < subcode_girth) {
            cout << "problem: required girth " << subcode_girth
                 << ", found only " << second_equations.tree_girth() << endl;
            return false;
        }
        cout << "ok: g=" << second_equations.tree_girth()
             << ", eqs=" << second_equations.n_equations()
             << ", nodes=" << second_equations.n_vertices() << endl;

        vector< int > a(first_equations_base_edges);
        vector< int > as(second_equations_base_edges);

		vector< int > coef(first_equations_base_edges);    // used for non-binary codes only 

        int m1 = tailbite_length;

        cout << endl << " Press 'x' to stop, 'i' for info" << endl;

        int codes_selected = 0;
        int iteration = 0;
        try {
            console_exception_hook x_hook('x', interrupt_info_message, local_interruption_exception());
            vector< int > previous_modulo;
            vector< vector<int> > previous_codes;
			vector< vector<int> > previous_coefs;
            if (prev_length) {
                while (true) {
                    settings current = settings::from_stream(prev_data);
                    if (!current.can_select("min_modulo")) {
                        break;
                    }
                    int m0;
                    current.select("min_modulo").cast_to(m0);
                    if (m0 > tailbite_length) {
                        continue;
                    }
                    vector<int> current_data;
                    current.select("data").cast_to(current_data);
                    previous_modulo.push_back(m0);
                    previous_codes.push_back(current_data);

					if( q_mod > 2 )
					{
						vector<int> current_coefs;
						current.select("coef").cast_to(current_coefs);
						previous_coefs.push_back(current_coefs);
					}
                }
            }
            unsigned previous_ptr = 0;
            for ( ; iteration < num_cycles && codes_selected < max_codes; ++iteration) {
                if (codes_tested_log_period > 0 && (iteration + 1) % codes_tested_log_period == 0) {
                    cout << iteration + 1 << " tested" << ", " << codes_selected << " selected" << endl;
                }
                int prev_elements = prev_length ? cumwcol[prev_length] : 0;

                int m0 = 1000000000;
                if (prev_length) 
				{
                    m0 = previous_modulo[previous_ptr];
                    vector< int > &previous_v = previous_codes[previous_ptr];
                    for (int i = 0; i < prev_elements; ++i) {
                        a[i] = previous_v[i];
                    }

					if( q_mod > 2 )
					{
						vector< int > &previous_c = previous_coefs[previous_ptr];
						for (int i = 0; i < prev_elements; ++i) {
							coef[i] = previous_c[i];
						}
					}

                    previous_ptr = (previous_ptr + 1) % previous_modulo.size();
                }

                for (int internal = 0; internal < std::min(code_expansion_max_trials, m0); ++internal) 
				{
                    for (int i = prev_elements, i_max = (int) a.size(); i < i_max; ++i) {
                        if (mask[i] == 1) {
                            a[i] = next_random_int(0, tailbite_length);
                        } else {
                            a[i] = 0;
                        }
                    }
                    for (int i = 0, i_max = (int) as.size(); i < i_max; ++i) {
                        as[i] = a[nonzero_mapping[i]];
                    }

                    if( q_mod > 2 )
					{
						for (int i = prev_elements, i_max = (int) a.size(); i < i_max; ++i) {
							if (coef_mask[i] == 1) {
								coef[i] = next_random_int(0, q_mod);
							} else {
								coef[i] = 0;
							}
						}

						improve_coefs( coef, bidiagonalFlag, all_col_2, rows, q_mod );
					}

                    voltage_check_result check;
					
					if( q_mod == 2 )
						check = first_equations.check_voltages(a, tailbite_length);
					else
						check = first_equations.check_voltages_and_coefs(a, tailbite_length, coef, q_mod);


                    if (subcode_girth > target_girth &&
                        check.min_ok_modulo <= tailbite_length &&
                        (exact_tailbite == 0 || check.check_modulo(exact_tailbite))
                    ) 
					{
                        check = second_equations.check_voltages(as, tailbite_length);
                    }
                    if (check.min_ok_modulo <= tailbite_length &&
                        (exact_tailbite == 0 || check.check_modulo(exact_tailbite))
                    ) 
					{
//------------------------------  check matrix --------------------------------------------
						bool bad_matrix = false;
						{
							int min_ok_modulo = check.min_ok_modulo;
							MARK_PARAM mark_param;
							mark_param.ncol = starting_length;
							mark_param.nrow = rows;
							mark_param.size = (int)a.size();

							matrix<int> current_HM(rows,starting_length);
							matrix<int> current_HC(rows,starting_length);

							// get part of original matrix
							for( int i = 0; i < mark_param.nrow; i++ )
							{
								for( int j = 0; j < mark_param.ncol; j++ )
								{
									current_HM(i,j) = HM_orig(i,j);
								}
							}
							// use current mark
							for (int i = 0, k = 0; i < mark_param.ncol && k < mark_param.size; i++)
							{
								for (int j = 0; j < mark_param.nrow && k < mark_param.size; j++)
								{
									if (current_HM(j, i) != -1)
										current_HM(j, i) = a[k++];
								}
							}

                            if( q_mod > 2 )
							{
								// get part of original matrix
								for( int i = 0; i < mark_param.nrow; i++ )
								{
									for( int j = 0; j < mark_param.ncol; j++ )
									{
										current_HC(i,j) = HM_orig(i,j);     // ????????????????????????????????????????
									}
								}
								// use current mark
								for (int i = 0, k = 0; i < mark_param.ncol && k < mark_param.size; i++)
								{
									for (int j = 0; j < mark_param.nrow && k < mark_param.size; j++)
									{
										if (current_HC(j, i) != -1)
											current_HC(j, i) = coef[k++];
									}
								}
							}

							vector< bit > codeword;
							int codeword_exitcode;
							if( q_mod > 2 )
							{
								codeword_exitcode = 0;
							}
							else
							{
								codeword_exitcode = random_codeword(current_HM, tailbite_length, codeword);
							}
							if (codeword_exitcode < 0) {
#if 0
								printf("Bad matrix\n");
								
								for( int i = 0; i < mark_param.nrow; i++ )
								{
									for( int j = 0; j < mark_param.ncol; j++ )
									{
										printf("%3d ", current_HM(i,j) );	// = HM_orig(i,j);
									}
									printf("\n");
								}
#endif
								bad_matrix = true;
							} else if (codeword_exitcode > 0) {
//								printf("Bad encoding\n");

								bad_matrix = true;
							}

							if( !bad_matrix )
							{
								if( use_ACE | use_SPEC )
								{
									vector< int > ACE(GTARGET);
									vector< int > SPEC(GTARGET);

									int curr_girth = trace_matrix( current_HM, min_ACE, max_SPEC, tailbite_length, GTARGET, ACE, SPEC ); 

									if( tarGirth > curr_girth ) 
										bad_matrix = true;
									else
									{
										bool bad_by_ACE  = false;
										bool bad_by_SPEC = false;

										if( use_ACE )
										{
											int n = min_ACE.size();
											int bad_by_ACE = tarACEbyACE > ACE[0];

											for( int k = 0; k < n; k++ )
												bad_by_ACE |= ACE[k] < min_ACE[k];

											if( !bad_by_ACE ) 
											{
												if( tarSPECbyACE * 1.3 < SPEC[0] ) 
													bad_by_ACE = true;
											}

											if( !bad_by_ACE )
											{
												show_matrix_property( ACE, SPEC, curr_girth, min_ACE.size() );

												if( tarACEbyACE < ACE[0] ) 
												{	
													tarGirth = curr_girth;
													tarACEbyACE = ACE[0];
													tarSPECbyACE = SPEC[0];

													cout << "better ACE found!!! tarGirth: " << tarGirth << ", tarACE: " << tarACEbyACE << endl;
												}
											}
										}

										if( use_SPEC )
										{
											//cout << "ACE[0] " << ACE[0] << " SPEC[0] " << SPEC[0] << endl;

											int n = max_SPEC.size();
											int bad_by_SPEC = tarSPECbySPEC < SPEC[0];

											for( int k = 0; k < n; k++ )
												bad_by_SPEC |= SPEC[k] > max_SPEC[k];

											if( !bad_by_SPEC ) 
											{
												//if( tarACEbySPEC * 0.15 > ACE[0] ) 	bad_by_SPEC = true;
											}

											if( !bad_by_SPEC )
											{
												show_matrix_property( ACE, SPEC, curr_girth, min_ACE.size() );


												if( tarSPECbySPEC > SPEC[0] ) 
												{	
													tarGirth = curr_girth;
													tarACEbySPEC = ACE[0];
													tarSPECbySPEC = SPEC[0];

													cout << "better SPEC found!!! tarGirth: " << tarGirth << ", tarSPEC: " << tarSPECbySPEC << endl;
												}
											}
										}
										bad_matrix = bad_by_ACE & bad_by_SPEC;
									}
								}
							}
						}
//-----------------------------------------------------------------------------------------
						if( !bad_matrix )
						{

							settings result;
							//===================================================================
							mark_record.a = a;
							mark_record.min_ok_modulo = check.min_ok_modulo;
							if( q_mod > 2 )
								mark_record.coef = coef;
							vmark_record.push_back(mark_record);
							//===================================================================
							result.open("min_modulo").set(check.min_ok_modulo);
							result.open("data").set(a);
						
							if( q_mod > 2 )
								result.open("coef").set(coef);

							result.to_stream(curr_data);

							++codes_selected;

							if (codes_selected_log_period && codes_selected % codes_selected_log_period == 0) {
								cout << codes_selected << " codes written to '" << curr_file_name << "'" << endl;
							}
							if (check.min_ok_modulo < m1) {
								m1 = check.min_ok_modulo;
								ofstream extra(extra_output_file.c_str(), std::ios_base::out | std::ios_base::app);

								ogroup(extra, cout) << " g=" << target_girth << " for M=" << m1;
	                            if (prev_length) {
		                            ogroup(extra, cout) << " from M=" << m0 << "\n";
			                    } else {
				                    ogroup(extra, cout) << "\n";
					            }

						        for (int i = 0, i_max = (int) a.size(); i < i_max; ++i) {
							        if (i > 0) {
								        extra << " ";
									}
									extra << a[i];
								}
								extra << endl;
							}
						}
                    }
                    console_message_hook i_hook('i', "prints how many codes were tested",
                        format_to_string(" nb=%d, %d tested, %d found", starting_length, iteration + 1, codes_selected));
                    console_check_hooks();
                }
            }
        } 
		catch (local_interruption_exception &) 
		{ 
			//cout << "---------------->local_exception\n"; 
		}
        cout << endl << "  " << codes_selected << final_message_prefix << iteration << final_message_suffix << endl;

		//std::cout << "vmark_record contains " << vmark_record.size() << " elements.\n";

        if (codes_selected < my_min_codes) {
            cout << "Too few codes of length " << starting_length << endl;
            return false;
        }

		bool checked_length = check_code_length( starting_length );
		if( !checked_length )
		{
			cout << "\n----> Simulation skiped\n\n";
			return true;
		}
#ifdef WIN32
		SetConsoleTextAttribute (hStdout, FOREGROUND_GREEN | FOREGROUND_INTENSITY );
#endif
		std::pair<double, double> bp_result;
		MARK_CANDIDATE m_cand;

		int vmarksize = (int)vmark_record.size();
		cout << endl << "Press 's' to stop simulation" << endl;

		int number_of_simulations = 1;
		bool cr_flag = false;
		try
		{
			bool snr_est_done = false;
			for( int i = 0; i < vmarksize; i++ )
			{
				if( !snr_est_done )	// i == 0 )
				{
					vector< int > a = vmark_record[0].a;
					int min_ok_modulo = vmark_record[0].min_ok_modulo;
					MARK_PARAM mark_param;
					mark_param.ncol = columns;	//starting_length;
					mark_param.nrow = rows;
					mark_param.size = (unsigned int)a.size();

					matrix<int> current_HM(rows,columns);
					matrix<int> current_HC(rows,columns);

					// get part of original matrix
					for( int i = 0; i < mark_param.nrow; i++ )
					{
						for( int j = 0; j < mark_param.ncol; j++ )
						{
							current_HM(i,j) = HM_orig(i,j);
						}
					}
					// use current mark
					for (int i = 0, k = 0; i < mark_param.ncol && k < mark_param.size; i++)
					{
						for (int j = 0; j < mark_param.nrow && k < mark_param.size; j++)
						{
							if (current_HM(j, i) != -1)
								current_HM(j, i) = a[k++];
						}
					}

					if( q_mod > 2 )
					{
						// get part of original matrix
						for( int i = 0; i < mark_param.nrow; i++ )
						{
							for( int j = 0; j < mark_param.ncol; j++ )
							{
								current_HC(i,j) = HM_orig(i,j);  // ?????????????????????????
							}
						}
						// use current mark
						for (int i = 0, k = 0; i < mark_param.ncol && k < mark_param.size; i++)
						{
							for (int j = 0; j < mark_param.nrow && k < mark_param.size; j++)
							{
								if (current_HC(j, i) != -1)
									current_HC(j, i) = coef[k++];
							}
						}
					}
					
					double best_fer = 10.0;
					double snr = start_snr_for_estimation;
					cout << "snr estimation for codes of length " << starting_length << endl;
					while( best_fer > snr_bound )
					{
						//reset_random(); // all codes are tested with same noise
						bp_result = bp_simulation(
							            q_mod,
										current_HM,
										current_HC,
										0, // ????????????????????????????????????????????????????
										tailbite_length,
										number_of_iterations,	//num_bp_iterations,
										number_of_err_blocks,	//num_frame_errors,
										number_of_codewords,	//num_experiments,
										snr,					//snrs[s],
										1.0,					//best_errors[s],
										decoder_type,			//decoder_type,
										0,						//modulation_type,
										0,						//permutation_type,
										tailbite_length,		//permutation_block,
										tailbite_length,		//permutation_inter,
										punctured_blocks,
										0
									);
						if( bp_result.second < 0.0 || bp_result.second == 10.0 )
							goto next_mark;	//continue;
						snr_est_done = true;
						best_fer = bp_result.second;
						cout << "snr = " << snr << "  FER = "<< best_fer<< endl;
						snr += snr_step;
					}
					snr_for_selection = snr-snr_step;
					min_fer  = best_fer;
					max_fer = min_fer * 1.5;
					m_cand.a = a;
					m_cand.fer = best_fer;
					m_cand.min_ok_modulo = min_ok_modulo;
					vmark_candidate.push_back( m_cand );
					candidate_cnt = 1;

					//std::cout << "mark_candidate contains " << mark_candidate.size() << " elements.\n";
					cout << "\n\nsimulation for codes of length " << starting_length << endl<<endl;
					cout << "min_fer = " << min_fer <<  " max_fer =  " << max_fer << "  number_of_selected_codes = " << candidate_cnt << "  vmark_candidate contains " <<  vmark_candidate.size() << endl;

				}
				else
				{
					console_exception_hook x_hook('s', interrupt_info_message, local_interruption_exception());			
					vector< int > a = vmark_record[i].a;
					int min_ok_modulo = vmark_record[i].min_ok_modulo;

					MARK_PARAM mark_param;
					mark_param.ncol = columns;	//starting_length;
					mark_param.nrow = rows;
					mark_param.size = (unsigned int)a.size();

					matrix<int> current_HM(rows,columns);
					matrix<int> current_HC(rows,columns);

					// get part of original matrix
					for( int i = 0; i < mark_param.nrow; i++ )
					{
						for( int j = 0; j < mark_param.ncol; j++ )
						{
							current_HM(i,j) = HM_orig(i,j);
						}
					}
					// use current mark
					for (int i = 0, k = 0; i < mark_param.ncol && k < mark_param.size; i++)
					{
						for (int j = 0; j < mark_param.nrow && k < mark_param.size; j++)
						{
							if (current_HM(j, i) != -1)
								current_HM(j, i) = a[k++];
						}
					}

					if( q_mod > 2 )
					{
						// get part of original matrix
						for( int i = 0; i < mark_param.nrow; i++ )
						{
							for( int j = 0; j < mark_param.ncol; j++ )
							{
								current_HC(i,j) = HM_orig(i,j);///// ?????????????????????????????
							}
						}
						// use current mark
						for (int i = 0, k = 0; i < mark_param.ncol && k < mark_param.size; i++)
						{
							for (int j = 0; j < mark_param.nrow && k < mark_param.size; j++)
							{
								if (current_HC(j, i) != -1)
									current_HC(j, i) = coef[k++];
							}
						}
					}

					reset_random(); // all codes are tested with same noise
					bp_result = bp_simulation(
						                q_mod,
										current_HM,
										current_HC,
										0, //?????????????????????????????????????????????????
										tailbite_length,
										number_of_iterations,		//num_bp_iterations,
										number_of_err_blocks,		//num_frame_errors,
										number_of_codewords,		//num_experiments,
										snr_for_selection,			//snrs[s],
										1.0,						//best_errors[s],
										decoder_type,				//decoder_type,
										0,							//modulation_type,
										0,							//permutation_type,
										tailbite_length,			//permutation_block,
										tailbite_length,			//permutation_inter,
										punctured_blocks,			
										0
									);
					if( bp_result.second == 10.0 || bp_result.second == -10.0 )
						break;
					number_of_simulations++; 
					if( bp_result.second < 0.0085 )
					{
						while( bp_result.second < 0.0085 )
						{
							if( cr_flag ) { cout << endl; cr_flag = false; }
							snr_for_selection -= snr_step/2.;
							cout << " fer = " << bp_result.second << "  ";
							cout << " SNR reestimaition :  " << snr_for_selection << endl;
							//reset_random(); // all codes are tested with same noise
							bp_result = bp_simulation(
								        q_mod,
										current_HM,
										current_HC,
										0, //??????????????????????????????????????????????????????????
										tailbite_length,
										number_of_iterations,		//num_bp_iterations,
										number_of_err_blocks,		//num_frame_errors,
										number_of_codewords,		//num_experiments,
										snr_for_selection,			//snrs[s],
										1.0,						//best_errors[s],
										decoder_type,				//decoder_type,
										0,							//modulation_type,
										0,							//permutation_type,
										tailbite_length,			//permutation_block,
										tailbite_length,			//permutation_inter,
										punctured_blocks,
										0
									);
						


						}
						min_fer = bp_result.second;
						max_fer = min_fer * 1.5;
						vmark_candidate.clear();

					}
					if( bp_result.second <= min_fer && bp_result.second != -1.0 ) //|| (bp_result.second > min_fer && bp_result.second < max_fer ) )
					{
						min_fer = bp_result.second;
						max_fer = min_fer * 1.5;
						m_cand.a = a;
						m_cand.fer = min_fer;
						m_cand.min_ok_modulo = min_ok_modulo;
						//std::cout << "1 - mark_candidate contains " << mark_candidate.size() << " elements.\n";
						vmark_candidate.push_back( m_cand );
						//std::cout << "2 - mark_candidate contains " << mark_candidate.size() << " elements.\n";

						//candidate_cnt =  recalc_candidates( mark_candidate );
						{
							vector< MARK_CANDIDATE >::iterator iter;
							int sum = 0;
							for( iter = vmark_candidate.begin(); iter != vmark_candidate.end(); iter++ )
							{
								double fer = (*iter).fer;
								if( fer <= max_fer )
									sum++;
							}
							candidate_cnt = sum;
							if( cr_flag ) { cout << endl; cr_flag = false; }
							cout << "min_fer = " << min_fer <<  " max_fer =  " << max_fer << "  number_of_selected_codes = " << candidate_cnt << "  vmark_candidate contains " <<  vmark_candidate.size() << endl;

						}

					}
					else
					{	
						if( bp_result.second > min_fer && bp_result.second < max_fer )
						{
							m_cand.a = a;
							m_cand.fer = bp_result.second;
							m_cand.min_ok_modulo = min_ok_modulo;
							//std::cout << "3 - mark_candidate contains " << mark_candidate.size() << " elements.\n";
							vmark_candidate.push_back( m_cand );
							//std::cout << "4 - mark_candidate contains " << mark_candidate.size() << " elements.\n";
							candidate_cnt++;	// =  recalc_candidates( mark_candidate );
						}
					}
//					if( vmark_candidate.size() % 2 == 0 )
//						cout << "min_fer = " << min_fer <<  " max_fer =  " << max_fer << "  number_of_selected_codes = " << candidate_cnt << "  vmark_candidate contains " <<  vmark_candidate.size() << endl;
//					printf("# of simulations = %4d ( fer = %f  # of selected = %4d)          \r", number_of_simulations, bp_result.second, candidate_cnt );
					cout << "# of simulations = " << number_of_simulations << "( fer = " << bp_result.second << " # of selected = " << candidate_cnt << ")             \r";
					cr_flag = true;

				}
				if( number_of_selected_codes <= candidate_cnt )
				{
					cout << "min_fer = " << min_fer <<  " max_fer =  " << max_fer << "  number_of_selected_codes = " << candidate_cnt << "  vmark_candidate contains " <<  vmark_candidate.size() << endl;
					break;
				}
				next_mark: ;
			}// cycle over marks
		}	// try block
		catch (local_interruption_exception &) 
		{ 
			//cout << " s -> local_exception\n"; 
		}

		if( candidate_cnt )
		{
			int i;
			string sym_curr_file_name = "";
			sym_curr_file_name += "sim_";
			sym_curr_file_name += curr_file_name;
			//cout << " sym_mark_name = " << sym_curr_file_name << endl;
			
			{
				settings result;
				ofstream curr_data(sym_curr_file_name.c_str());
		        if (!curr_data) 
				{
				    cout << "Cannot open file '" << sym_curr_file_name << "' for writing" << endl;
					return false;
				}
				int imax = (int)vmark_candidate.size();
				for( i = 0; i < imax; i++ )
				{
					if( vmark_candidate[i].fer < max_fer )
				{
					result.open("min_modulo").set( vmark_candidate[i].min_ok_modulo);
                    result.open("data").set(vmark_candidate[i].a);
                    result.to_stream(curr_data);
				}
			}
		}

		}
		else
		{
			cout << "Too few codes of length " << starting_length << endl;
		}

#ifdef WIN32
		SetConsoleTextAttribute (hStdout,FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED );
#endif


        return true;
    }

    virtual void run(ggp_state &) 
	{
        bool continue_existing;
        vector<settings> input_codes;
        string output_file_mask;
        string error_name;
//      string modeling_config_file;
//      string modeling_values_file;
//      vector<string> modeling_additional_lines;

        settings scenario = settings::from_file(scenario_file);
        scenario.select("subcode_columns").cast_to(subcode_columns);
        scenario.select("target_girth").cast_to(target_girth);
        scenario.select("target_girth_subcode").cast_to(subcode_girth);
        scenario.select("tailbite_length").cast_to(tailbite_length);
        scenario.select("exact_tailbite").cast_to(exact_tailbite);
        scenario.select("starting_columns").cast_to(starting_length);
        scenario.select("trials/min_candidates").cast_to(min_codes);
        scenario.select("trials/enough_candidates").cast_to(max_codes);
        scenario.select("trials/max_cycles_per_iteration").cast_to(num_cycles);
        scenario.select("continue_existing").cast_to(continue_existing);
        scenario.select("input").cast_to(input_codes);
        scenario.select("output_file_mask").cast_to(output_file_mask);
        scenario.select("extra_output_file").cast_to(extra_output_file);
        scenario.select("error_name").cast_to(error_name);

		scenario.select("min_ace").cast_to(min_ACE);
		scenario.select("max_ace_spec").cast_to(max_SPEC);

		scenario.select("q_mod").cast_to(q_mod);

//        scenario.select("modeling/config_file").cast_to(modeling_config_file);
//        scenario.select("modeling/values_file").cast_to(modeling_values_file);
//        scenario.select("modeling/additional_lines").cast_to(modeling_additional_lines);
/*
simulation = {
  decoder_type = 2
  number_of_iterations = 15
  number_of_codewords = 10000000
  number_of_err_blocks = 25
  tested_code_lengths = { 21 24 27 30 32}
  start_snr_for_estimation = 2.1
  snr_step = 0.1
  snr_bound = 0.02
  number_of_selected_codes = 100
  }*/

		scenario.select("simulation/decoder_type").cast_to(decoder_type);
		scenario.select("simulation/number_of_iterations").cast_to(number_of_iterations);
		scenario.select("simulation/number_of_codewords").cast_to(number_of_codewords);
		scenario.select("simulation/number_of_err_blocks").cast_to(number_of_err_blocks);
		scenario.select("simulation/tested_code_lengths").cast_to(tested_code_lengths);
		scenario.select("simulation/start_snr_for_estimation").cast_to(start_snr_for_estimation);
		scenario.select("simulation/snr_step").cast_to(snr_step);
		scenario.select("simulation/snr_bound").cast_to(snr_bound);
		scenario.select("simulation/number_of_selected_codes").cast_to(number_of_selected_codes);

//        modeling_config_file = get_native_path(modeling_config_file);
//        modeling_values_file = get_native_path(modeling_values_file);

        set_output_file_mask(get_native_path(output_file_mask));
        extra_output_file = get_native_path(extra_output_file);

        exact_tailbite *= tailbite_length;
        rows = columns = -1;

		if( q_mod <= 2 )
		{
			use_ACE  = min_ACE[0]  < 0 ? 0 : 1;
			use_SPEC = max_SPEC[0] < 0 ? 0 : 1;
		}
		else
		{
			use_ACE  = 0;
			use_SPEC = 0;
		}
		


        double best_metric = 1e9;
        for (unsigned i = 0, i_max = (unsigned int)input_codes.size(); i < i_max; ++i) {
            matrix<int> current_HM;
            vector<settings> logs;
            input_codes[i].select("code").cast_to(current_HM);
            input_codes[i].select("simulation_logs").cast_to(logs);
            if (rows == -1) {
                rows = current_HM.n_rows();
                columns = current_HM.n_cols();
            } else {
                if (rows != current_HM.n_rows() || columns != current_HM.n_cols()) {
                    cout << "Warning: unequal matrices in the input!" << endl;
                    continue;
                }
            }
            double maximum_error = 0;
            for (unsigned j = 0, j_max = (unsigned int)logs.size(); j < j_max; ++j) {
                double current_error;
                logs[j].select(error_name).cast_to(current_error);
                maximum_error = std::max(maximum_error, current_error);
            }
            if (best_metric > maximum_error) {
                best_metric = maximum_error;
                HM = current_HM;
				//================================================
				HM_orig = HM;
				//================================================
            }
        }

        // Raw matrices contain -1's, prepared matrices do not.
        bool is_raw_matrix = false;
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < columns; ++c) {
                is_raw_matrix |= HM(r, c) < 0;
            }
        }
        vector< int > original_values;
        if (is_raw_matrix) 
		{
			for (int c = 0; c < columns; ++c) 
			{
                for (int r = 0; r < rows; ++r) 
				{
                    int cur = HM(r, c);
                    if (cur < 0) 
					{
                        HM(r, c) = 0;
                    } 
					else 
					{
                        original_values.push_back(cur);
                        if (cur == 0) 
						{
                            HM(r, c) = 2;
                        } 
						else 
						{
                            HM(r, c) = 1;
                        }
                    }
                }
            }

			if( q_mod > 2 )
			{
				HCM = HM_orig;
				for (int c = 0; c < columns; ++c) 
				{
					for (int r = 0; r < rows; ++r) 
					{
						if( HM(r, c) < 0 )		// exactly HM !!! 
							HCM(r, c) = 0;
						else 
							HCM(r, c) = 1;
					}
				}
			}
		}
		else
			cout << "is not raw matrix" << endl;

        wcol = vector< int >(columns);
        cumwcol = vector< int >(columns + 1);
        for (int i = 0; i < columns; ++i) {
            for (int j = 0; j < rows; ++j) {
                wcol[i] += HM(j, i) > 0 ? 1 : 0;
            }
            cumwcol[i + 1] = cumwcol[i] + wcol[i];
        }

		all_col_2 = 1;
		for( int i = 0; i < columns; i++ ) 
		{
			if( wcol[i] != 2 )
				all_col_2 = 0;
		}

        mask = vector< int >();
        for (int i = 0; i < columns; ++i) {
            for (int j = 0; j < rows; ++j) {
                if (HM(j, i) > 0) {
                    mask.push_back(HM(j, i));
                }
            }
        }

		if( q_mod > 2 )
		{
			coef_mask = vector< int >();
			for (int i = 0; i < columns; ++i) 
			{
				for (int j = 0; j < rows; ++j) {
					if (HCM(j, i) > 0) {
						coef_mask.push_back(HCM(j, i));
					}
				}
			}
		}

		bidiagonalFlag = check_bidiagonal( HM_orig );

		if (!continue_existing) {
            if (!single_iteration(0, "interrupts initial code generation", " codes selected, ",
                                  " tested", 100000, 1000, 1, 0)) {
                return;
            }
        }


// Check if files exist
		{
			string curr_file_name = output_file_name(starting_length);
			string sym_curr_file_name = "";
			sym_curr_file_name += "sim_";
			sym_curr_file_name += curr_file_name;

			FILE *fold = fopen( curr_file_name.c_str(), "r" );
			FILE *fnew = fopen( sym_curr_file_name.c_str(), "r" );

			if( fnew != NULL && fold != NULL )
			{
				fclose(fold);
				fclose(fnew);
				int ret = remove( curr_file_name.c_str() );
				if( ret )
				{
					cout << "Can't delete file " << curr_file_name << endl;
					return;
				}
				ret = rename( sym_curr_file_name.c_str(), curr_file_name.c_str() );
				if( ret )
				{
					cout << "can't rename file " << sym_curr_file_name << " to " << curr_file_name << endl;
					return;
				}
			}

			if( fold == NULL )
			{
				cout << " Can't open " << curr_file_name << endl;
				return;
			}

			if( fnew == NULL )
			{
				// do nothing
			}
		}

        while (starting_length < columns) {
            ++starting_length;
            if (!single_iteration(starting_length - 1, "interrupts current code expansion",
                                  " codes are found from ", "", 0, 0, 5, min_codes)) {
                return;
            }
			{
				string curr_file_name = output_file_name(starting_length);
				string sym_curr_file_name = "";
				sym_curr_file_name += "sim_";
				sym_curr_file_name += curr_file_name;
				FILE *fold = fopen( curr_file_name.c_str(), "r" );
				FILE *fnew = fopen( sym_curr_file_name.c_str(), "r" );

				if( fnew != NULL && fold != NULL )
				{
					fclose(fold);
					fclose(fnew);
					int ret = remove( curr_file_name.c_str() );
					if( ret )
					{
						cout << "Can't delete file " << curr_file_name << endl;
						return;
					}
					ret = rename( sym_curr_file_name.c_str(), curr_file_name.c_str() );
					if( ret )
					{
						cout << "can't rename file " << sym_curr_file_name << " to " << curr_file_name << endl;
						return;
					}
				}

				if( fold == NULL )
				{
					cout << " Can't open " << curr_file_name << endl;
					return;
				}

				if( fnew == NULL )
				{
					// do nothing
				}
			}

        }

#if 0
        cout << "Writing modeling configuration to " << modeling_config_file << "... ";
        {
            ofstream modeling_config_stream(modeling_config_file.c_str());
            modeling_config_stream << HM.n_rows() << " " << HM.n_cols() << " //  size of base matrix\n";
            for (unsigned i = 0; i < modeling_additional_lines.size(); ++i) {
                modeling_config_stream << modeling_additional_lines[i] << "\n";
            }
            for (int r = 0, r_max = HM.n_rows(); r < r_max; ++r) {
                for (int c = 0, c_max = HM.n_cols(); c < c_max; ++c) {
                    modeling_config_stream << " " << HM(r, c);
                }
                modeling_config_stream << "\n";
            }
        }
        cout << "done!" << endl;
#endif  //0
#if 0
        cout << "Writing modeling values to " << modeling_values_file << "... ";
        {
            string last_file_name = output_file_name(columns);
            ifstream last_data(last_file_name.c_str());
            if (!last_data) {
                cout << "    Error: cannot read the last data file!" << endl;
            } else {
                ofstream modeling_values_stream(modeling_values_file.c_str());
                if (is_raw_matrix) {
                    for (unsigned i = 0, i_max = (unsigned int)original_values.size(); i < i_max; ++i) {
                        modeling_values_stream << original_values[i];
                        if (i + 1 != i_max) {
                            modeling_values_stream << " ";
                        } else {
                            modeling_values_stream << "\n";
                        }
                    }
                }
                while (true) {
                    settings current = settings::from_stream(last_data);
                    if (!current.can_select("data")) {
                        break;
                    }
                    vector<int> data;
                    current.select("data").cast_to(data);
                    for (unsigned i = 0, i_max = (unsigned int)data.size(); i < i_max; ++i) {
                        modeling_values_stream << data[i];
                        if (i + 1 != i_max) {
                            modeling_values_stream << " ";
                        } else {
                            modeling_values_stream << "\n";
                        }
                    }
                }
                modeling_values_stream << "-1\n";
            }
        }
        cout << "done!" << endl;
#endif
    }

    virtual ~scenario_based_code_generation() {}
};

static int trace_matrix( matrix<int> current_HM, vector<int> min_ACE, vector<int> max_ACE_spec, int M, int gtarget, vector<int> &ACE, vector<int> &girth_spectrum )
{
	int i, j;
	int girth;
	ARRAY matr;
	int S[GMAX];
	int SA[GMAX];
	int minACE[GMAX];
	int maxACEspec[GMAX];
	int HDrow = current_HM.n_rows();
	int HDcol = current_HM.n_cols();
	int gmax = GMAX;
	int flag;
	int t0 = min_ACE.size();
	int t1 = max_ACE_spec.size();
	int *minACE_ptr;
	int *maxACEspec_ptr;

	for( i = 0; i < GMAX; i++ ) S[i]  = 0;	
	for( i = 0; i < GMAX; i++ ) SA[i] = 0;	
	for( i = 0; i < GMAX; i++ ) minACE[i] = 0;	
	for( i = 0; i < GMAX; i++ ) maxACEspec[i] = 100000000;	

	if( t0 > GMAX ) 
		t0 = GMAX;
	for( i = 0; i < t0; i++ ) minACE[i] = min_ACE[i];	

	if( t1 > GMAX ) 
		t1 = GMAX;
	for( i = 0; i < t1; i++ ) maxACEspec[i] = max_ACE_spec[i];	

	matr.ndim = 2;

	put_nrow( &matr, HDrow );
	put_ncol( &matr, HDcol );
	put_addr( &matr, Alloc2d_int( HDrow, HDcol ) );

	for( i = 0; i < HDrow; i++ )
		for( j = 0; j < HDcol; j++ )
			matr.addr[i][j] = current_HM(i,j);

	minACE_ptr     =     minACE[0] < 0 ? NULL : minACE;
	maxACEspec_ptr = maxACEspec[0] < 0 ? NULL : maxACEspec;

	flag = trace_bound_pol_mon_pm( matr, M, gmax, gtarget, S, SA, minACE_ptr, maxACEspec_ptr );

	for( girth = 1; girth < GMAX+1; girth++ )
	{
		if( S[girth-1] )
			break;
	}

	free( matr.addr );

	for( int j = 0, i = girth-1; i < gmax; i++ )
	{
		if( SA[i] )
		{
			ACE[j++] = SA[i];
			if( j == GTARGET ) break;
		}
	}

	for( int j = 0, i = girth-1; i < gmax; i++ )
	{
		if( S[i] )
		{
			girth_spectrum[j++] = S[i];
			if( j == GTARGET ) break;
		}
	}


	return flag ? girth : -girth;
}

static void show_matrix_property( vector<int> ACE, vector<int> girth_spectrum, int girth, int num )
{
	printf("girth: %3d, ", girth);

	printf("ACE: "); 

	for( int i = 0; i < num; i++ )
		printf("%3d ", ACE[i] ); 

	printf(",  Girth_spectrum: "); 

	for( int i = 0; i < num; i++ )
		printf("%5d ", girth_spectrum[i] ); 
	printf("\n");
}

#endif
