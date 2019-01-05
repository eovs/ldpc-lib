// This macro is defined by the includer, main_search_ggp.cpp.
// This is done to prevent a possibly dumb Visual Studio from compiling this source,
// which it would definitely fail to do.
#ifdef INTERNAL_SOURCE_INCLUDE

class check_voltages_from_input : public ggp_function {
private:
    string input_file;
    string output_file;

    void problem_reading_file() {
        cerr << "Problem reading file '" << input_file << "'." << endl;
    }

    bool check_equations_on_voltages(equation_builder const &equations, vector< int > const &voltages) const {
        // Opening the output file for appending
        ofstream out(output_file.c_str(), std::ios_base::out | std::ios_base::app);
        // Grouper writes to both output file and stdout
        ogroup xout(cout, out);

        voltage_check_result res = equations.check_voltages(voltages);
        if (res.first_failed_equation != -1) {
            xout << " gfree < " << equations.tree_girth() << ", wrong eq_number = " << (res.first_failed_equation + 1) << "\n";
            return false;
        } else {
            if (res.failed_modulos.size() == 0) {
                // A very unlikely case, but...
                xout << " g >= " << equations.tree_girth() << " for all modulos\n";
            } else {
                xout << " g >= " << equations.tree_girth() << " for M = " << res.min_ok_modulo << " and all M > Mmax = " << res.max_nonok_modulo << "\n";
                for (int i = res.min_ok_modulo, cnt = 0; i < res.max_nonok_modulo; ++i) {
                    if (res.failed_modulos.count(i) == 0) {
                        xout << i << " ";
                        if (++cnt % 20 == 0) {
                            xout << "\n";
                        }
                    }
                }
                xout << "\n...\n";
            }
            return true;
        }
    }

public:
    check_voltages_from_input(string const &input_file, string const &output_file)
        : input_file(input_file)
        , output_file(output_file)
    {}

    virtual bool can_run(ggp_state const &state) const {
        return state.has_matrix_HM();
    }

    virtual void write_info(ggp_state const &, ostream &out) const {
        out << "Check voltages from '" << input_file << "' (degrees by columns).";
    }

    virtual void run(ggp_state &state) {
        if (state.p_column[state.base_code_length] > state.num_input_columns) {
            cerr << "Number of input columns ninp too small, aborting." << endl;
        } else {
            ifstream in(input_file.c_str());
            if (in) {
                // Raw input_a.txt contents
                vector< int > unmapped_a(state.cumwcol[state.num_input_columns - 1]);
                for (int i = 0; i < (int) unmapped_a.size(); ++i) {
                    in >> unmapped_a[i];
                }
                if (!in) {
                    problem_reading_file();
                    return;
                }
                // Subsequence for analysis
                vector< int > mapped_a(state.equations->n_vertices());  // seems to be an equivalent of original "nt"
                for (int i = 0; i < (int) mapped_a.size(); ++i) {
                    mapped_a[i] = unmapped_a[state.nonzero_mapping[i]];
                }
                int curr_girth = state.equations->tree_girth();
                if (check_equations_on_voltages(*state.equations, mapped_a)) {
                    while (check_equations_on_voltages(equation_builder(state.equations_graph, curr_girth + 2), mapped_a)) {
                        curr_girth += 2;
                    }
                } else {
                    while (!check_equations_on_voltages(equation_builder(state.equations_graph, curr_girth - 2), mapped_a)) {
                        curr_girth -= 2;
                    }
                }
            } else {
                problem_reading_file();
            }
        }
    }

    virtual ~check_voltages_from_input() {}
};

#endif
