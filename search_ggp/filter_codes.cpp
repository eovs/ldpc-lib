// This macro is defined by the includer, main_search_ggp.cpp.
// This is done to prevent a possibly dumb Visual Studio from compiling this source,
// which it would definitely fail to do.
#ifdef INTERNAL_SOURCE_INCLUDE

class filter_codes : public ggp_function {
private:
    string input_file;
    string output_file;

public:
    filter_codes(string const &input_file, string const &output_file)
        : input_file(input_file)
        , output_file(output_file)
    {}

    virtual bool can_run(ggp_state const &state) const {
        return state.has_matrix_HM();
    }

    virtual void write_info(ggp_state const &, ostream &out) const {
        out << "Filter the set in '" << input_file
            << "' to '" << output_file
            << "'.";
    }

    virtual void run(ggp_state &state) {
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

        int codes_tested = 0;
        int codes_found = 0;
        try {
            console_exception_hook x_hook('x', "interrupts code filtering", local_interruption_exception());
            while (true) {
                int m0 = -1;
                input >> m0;
                if (m0 < 0) {
                    break;
                }
                vector< int > unmapped_a(state.cumwcol[state.num_input_columns - 1]);
                for (int i = 0, i_max = (int) unmapped_a.size(); i < i_max; ++i) {
                    input >> unmapped_a[i];
                }
                vector< int > mapped_a(state.base_code_length);
                for (int i = 0, i_max = (int) mapped_a.size(); i < i_max; ++i) {
                    mapped_a[i] = unmapped_a[state.nonzero_mapping[i]];
                }
                ++codes_tested;
                voltage_check_result res = state.equations->check_voltages(mapped_a);
                if (res.first_failed_equation == -1
                        && res.min_ok_modulo <= state.tailbite_length
                        && (state.exact_tailbite == 0 || res.check_modulo(state.exact_tailbite))) {
                    out << res.min_ok_modulo;
                    for (int i = 0, i_max = (int) (unmapped_a.size()); i < i_max; ++i) {
                        out << " " << unmapped_a[i];
                    }
                    out << endl;
                    ++codes_found;
                }
                console_check_hooks();
            }
        } catch (local_interruption_exception &) {}
        out << "-1" << endl;
        cout << codes_found << " codes found from " << codes_tested << endl;
    }

    virtual ~filter_codes() {}
};

#endif
