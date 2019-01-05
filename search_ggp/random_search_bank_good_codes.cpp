// This macro is defined by the includer, main_search_ggp.cpp.
// This is done to prevent a possibly dumb Visual Studio from compiling this source,
// which it would definitely fail to do.
#ifdef INTERNAL_SOURCE_INCLUDE

class random_search_bank_good_codes : public ggp_function {
private:
    string output_file;
    string extra_output_file;

public:
    random_search_bank_good_codes(string const &output_file, string const &extra_output_file)
        : output_file(output_file)
        , extra_output_file(extra_output_file)
    {}

    virtual bool can_run(ggp_state const &state) const {
        return state.has_matrix_HM();
    }

    virtual void write_info(ggp_state const &, ostream &out) const {
        out << "Random search with creating a bank of good codes (files '"
            << output_file << "', '" << extra_output_file << "').";
    }

    virtual void run(ggp_state &state) {
        for (int i = 0; i < state.base_code_length; ++i) {
            if ((i == 0 && state.p_column[i] != 1) || (i > 0 && state.p_column[i - 1] + 1 != state.p_column[i])) {
                cerr << " Wrong base matrix. It should consist of nb first columns. Edit configuration. " << endl;
                return;
            }
        }
        cout << "New search. File '" << output_file << "' will be rewritten. Press any key to continue" << endl;
        get_character();

        ofstream out(output_file.c_str());
        if (!out) {
            cerr << "Problem opening '" << output_file << "' for writing" << endl;
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

        cout << "Press 'x' to stop" << endl;

        vector< int > a(am.size());
        try {
            console_exception_hook x_hook('x', "interrupts current code search", local_interruption_exception());
            for (int round = 0, ok_codes = 0, m1 = state.tailbite_length; ; ++round) {
                if ((round + 1) % 1000000000 == 0) {
                    cout << round << " tested" << endl;
                }
                for (int i = 0, i_max = (int) (a.size()); i < i_max; ++i) {
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
                    if (++ok_codes % 1000 == 0) {
                        cout << ok_codes << " codes written to '" << output_file << "'" << endl;
                    }
                    if (m1 > res.min_ok_modulo) {
                        m1 = res.min_ok_modulo;
                        ofstream extra(extra_output_file.c_str(), std::ios_base::out | std::ios_base::app);
                        if (extra) {
                            ogroup(cout, extra) << "g=" << state.equations->tree_girth() << " for M=" << m1 << "\n";
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
                console_check_hooks();
            }
        } catch (local_interruption_exception &) {
            out << "-1" << endl;
        }
    }

    virtual ~random_search_bank_good_codes() {}
};

#endif
