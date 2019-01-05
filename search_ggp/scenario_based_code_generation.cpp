// This macro is defined by the includer, main_search_ggp.cpp.
// This is done to prevent a possibly dumb Visual Studio from compiling this source,
// which it would definitely fail to do.
#ifdef INTERNAL_SOURCE_INCLUDE

#ifdef WIN32
#define NOMINMAX
#include <Windows.h>
#endif	// WIN32

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
	int min_ok_modulo;
}MARK_RECORD;


int random_codeword( matrix< int > const &mx, int tailbite_length, vector< bit > &codeword );

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
        int m1 = tailbite_length;

        cout << endl << " Press 'x' to stop, 'i' for info" << endl;

        int codes_selected = 0;
        int iteration = 0;
        try {
            console_exception_hook x_hook('x', interrupt_info_message, local_interruption_exception());
            vector< int > previous_modulo;
            vector< vector<int> > previous_codes;
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
                }
            }
            unsigned previous_ptr = 0;
            for ( ; iteration < num_cycles && codes_selected < max_codes; ++iteration) {
                if (codes_tested_log_period > 0 && (iteration + 1) % codes_tested_log_period == 0) {
                    cout << iteration + 1 << " tested" << endl;
                }
                int prev_elements = prev_length ? cumwcol[prev_length] : 0;

                int m0 = 1000000000;
                if (prev_length) {
                    m0 = previous_modulo[previous_ptr];
                    vector< int > &previous_v = previous_codes[previous_ptr];
                    for (int i = 0; i < prev_elements; ++i) {
                        a[i] = previous_v[i];
                    }
                    previous_ptr = (previous_ptr + 1) % previous_modulo.size();
                }

                for (int internal = 0; internal < std::min(code_expansion_max_trials, m0); ++internal) {
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
                    voltage_check_result check = first_equations.check_voltages(a, tailbite_length);
                    if (subcode_girth > target_girth &&
                        check.min_ok_modulo <= tailbite_length &&
                        (exact_tailbite == 0 || check.check_modulo(exact_tailbite))
                    ) {
                        check = second_equations.check_voltages(as, tailbite_length);
                    }
                    if (check.min_ok_modulo <= tailbite_length &&
                        (exact_tailbite == 0 || check.check_modulo(exact_tailbite))
                    ) {
//------------------------------  check matrix --------------------------------------------
						bool bad_matrix = false;
						{
							int min_ok_modulo = check.min_ok_modulo;
							MARK_PARAM mark_param;
							mark_param.ncol = starting_length;
							mark_param.nrow = rows;
							mark_param.size = (int)a.size();

							matrix<int> current_HM(rows,starting_length);
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

							vector< bit > codeword;
							int codeword_exitcode = random_codeword(current_HM, tailbite_length, codeword);

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

						}
//-----------------------------------------------------------------------------------------
						if( !bad_matrix )
						{
                        settings result;
						//===================================================================
						mark_record.a = a;
						mark_record.min_ok_modulo = check.min_ok_modulo;
						vmark_record.push_back(mark_record);
						//===================================================================
                        result.open("min_modulo").set(check.min_ok_modulo);
                        result.open("data").set(a);
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
					
					double best_fer = 10.0;
					double snr = start_snr_for_estimation;
					cout << "snr estimation for codes of length " << starting_length << endl;
					while( best_fer > snr_bound )
					{
						//reset_random(); // all codes are tested with same noise
						bp_result = bp_simulation(
										current_HM,
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

					reset_random(); // all codes are tested with same noise
					bp_result = bp_simulation(
										current_HM,
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
										current_HM,
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

    virtual void run(ggp_state &) {
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
        if (is_raw_matrix) {
            for (int c = 0; c < columns; ++c) {
                for (int r = 0; r < rows; ++r) {
                    int cur = HM(r, c);
                    if (cur < 0) {
                        HM(r, c) = 0;
                    } else {
                        original_values.push_back(cur);
                        if (cur == 0) {
                            HM(r, c) = 2;
                        } else {
                            HM(r, c) = 1;
                        }
                    }
                }
            }
        }

        wcol = vector< int >(columns);
        cumwcol = vector< int >(columns + 1);
        for (int i = 0; i < columns; ++i) {
            for (int j = 0; j < rows; ++j) {
                wcol[i] += HM(j, i) > 0 ? 1 : 0;
            }
            cumwcol[i + 1] = cumwcol[i] + wcol[i];
        }

        mask = vector< int >();
        for (int i = 0; i < columns; ++i) {
            for (int j = 0; j < rows; ++j) {
                if (HM(j, i) > 0) {
                    mask.push_back(HM(j, i));
                }
            }
        }

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

#endif
