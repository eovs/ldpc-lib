// This macro is defined by the includer, main_search_ggp.cpp.
// This is done to prevent a possibly dumb Visual Studio from compiling this source,
// which it would definitely fail to do.
#ifdef INTERNAL_SOURCE_INCLUDE

class random_extension : public ggp_function {
private:
    string input_file;
    string output_file;
    string extra_output_file;

public:
    random_extension(
        string const &input_file,
        string const &output_file,
        string const &extra_output_file
    )
        : input_file(input_file)
        , output_file(output_file)
        , extra_output_file(extra_output_file)
    {}

    virtual bool can_run(ggp_state const &state) const {
        return state.has_matrix_HM();
    }

    virtual void write_info(ggp_state const &, ostream &out) const {
        out << "Random extending '" << input_file
            << "' by one or more columns to '" << output_file
            << "', extra output to '" << extra_output_file
            << "'.";
    }

    virtual void run(ggp_state &state) {
        for (int i = 0; i < state.base_code_length; ++i) {
            if ((i == 0 && state.p_column[i] != 1) || (i > 0 && state.p_column[i - 1] + 1 != state.p_column[i])) {
                cerr << " Wrong base matrix. It should consist of nb first columns. Edit configuration. " << endl;
                return;
            }
        }
        cout << "Extending " << state.num_input_columns << " to " << state.base_code_length << " columns" << endl;
        if (state.num_input_columns >= state.base_code_length) {
            cerr << "Nothing to extend! Check your config!" << endl;
            return;
        }

        vector< int > am;
        for (int i = 0; i < state.base_code_length; ++i) {
            for (int j = 0; j < state.rows; ++j) {
                if (state.HM(j, i) > 0) {
                    am.push_back(state.HM(j, i));
                }
            }
        }

        ifstream input(input_file.c_str());
        if (!input) {
            cerr << "Problem opening '" << input_file << "' for reading" << endl;
            return;
        }
        ofstream out(output_file.c_str());
        if (!out) {
            cerr << "Problem opening '" << output_file << "' for writing" << endl;
            return;
        }

        cout << "Press 'x' to stop" << endl;

        vector< int > a(am.size());
        int m1 = state.tailbite_length;
        try {
            console_exception_hook x_hook('x', "interrupts current random code extension", local_interruption_exception());
            while (true) {
                int m0 = -1;
                input >> m0;
                if (m0 < 0) {
                    input.seekg(0);
                    input >> m0;
                }
                for (int i = 0; i < state.cumwcol[state.num_input_columns - 1]; ++i) {
                    input >> a[i];
                }
                for (int round = 0; round < std::min(m0, 10); ++round) {
                    for (int i = state.cumwcol[state.num_input_columns - 1], i_max = (int) (a.size()); i < i_max; ++i) {
                        a[i] = am[i] == 1 ? next_random_int(0, state.tailbite_length) : 0;
                    }
                    voltage_check_result res = state.equations->check_voltages(a);
                    if (res.first_failed_equation == -1
                            && res.min_ok_modulo <= state.tailbite_length
                            && (state.exact_tailbite == 0 || res.check_modulo(state.exact_tailbite))) {
                        out << res.min_ok_modulo;
                        for (int i = 0, i_max = (int) (a.size()); i < i_max; ++i) {
                            out << " " << a[i];
                        }
                        out << endl;
                        if (m1 > res.min_ok_modulo) {
                            m1 = res.min_ok_modulo;
                            cout << "g=" << state.equations->tree_girth() << " for M=" << m1 << " from M=" << m0 << "\n";
                            ofstream extra(extra_output_file.c_str(), std::ios_base::out | std::ios_base::app);
                            if (extra) {
                                extra << "g=" << state.equations->tree_girth() << " for M=" << m1 << "\n";
                                // in the original code, the upper limit is not "nr" but "n". Why? Not doing it here.
                                for (int i = 0, i_max = (int) (a.size()); i < i_max; ++i) {
                                    if (i > 0) {
                                        extra << " ";
                                    }
                                    extra << a[i];
                                }
                                extra << endl;
                            } else {
                                cerr << "Problem with opening '" << extra_output_file << "' for writing" << endl;
                            }
                        }
                    }
                }
                console_check_hooks();
            }
        } catch (local_interruption_exception &) {
            out << "-1" << endl;
        }
    }

    virtual ~random_extension() {}
};

#endif
