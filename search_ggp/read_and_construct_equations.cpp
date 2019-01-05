// This macro is defined by the includer, main_search_ggp.cpp.
// This is done to prevent a possibly dumb Visual Studio from compiling this source,
// which it would definitely fail to do.
#ifdef INTERNAL_SOURCE_INCLUDE

class read_and_construct_equations : public ggp_function {
private:
    string input_file;

    void problem_reading_file(ggp_state const &state) {
        cerr << "Problem reading file '" << input_file << "'.";
        if (state.has_matrix_HM()) {
            cerr << " The existing H and equations are left intact.";
        }
        cerr << endl;
    }

public:
    read_and_construct_equations(string const &input_file)
        : input_file(input_file)
    {}

    virtual bool can_run(ggp_state const &) const {
        return true;
    }

    virtual void write_info(ggp_state const &state, ostream &out) const {
        if (state.has_matrix_HM()) {
            out << "Read the matrix H from file '" << input_file << "' (the existing matrix would be destroyed) and construct the equations.";
        } else {
            out << "Read the matrix H from file '" << input_file << "' and construct the equations.";
        }
    }

    virtual void run(ggp_state &state) {
        ifstream in(input_file.c_str()); // what a badass is C++'s library before C++11...
        if (in) {
            int rows, columns, target_girth, tailbite_length, exact_tailbite, base_code_length;

            in >> rows >> line_eater >> columns >> line_eater >> target_girth >> line_eater;
            in >> tailbite_length >> line_eater >> exact_tailbite >> line_eater;
            exact_tailbite *= tailbite_length;
            in >> base_code_length >> line_eater;

            // A good place before the first error state check.
            if (!in) {
                problem_reading_file(state);
                return;
            }

            // Now may allocate resources.

            vector< int > p_column(base_code_length);
            for (int i = 0; i < base_code_length; ++i) {
                in >> p_column[i];
            }

            int num_input_columns;
            in >> line_eater >> num_input_columns >> line_eater;

            matrix< int > input_matrix(rows, columns);
            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < columns; ++c) {
                    int curr;
                    in >> curr;
                    input_matrix(r, c) = curr;
                }
            }

            // Final error state check.
            if (!in) {
                problem_reading_file(state);
                return;
            }

            cout << "    [debug] [read_and_construct_equations] File contents scanned OK" << endl;

            // column weights (for all, not only used columns)
            vector< int > wcol(columns), cumwcol(columns + 1);
            for (int i = 0; i < columns; ++i) {
                for (int j = 0; j < rows; ++j) {
                    wcol[i] += input_matrix(j, i) > 0;
                    cumwcol[i + 1] = cumwcol[i] + wcol[i];
                }
            }

            // obtain base matrix of proper size for constructing equations
            matrix< int > proper_base_matrix(rows, num_input_columns);
            for (int c = 0; c < num_input_columns; ++c) {
                for (int r = 0; r < rows; ++r) {
                    proper_base_matrix(r, c) = input_matrix(r, p_column[c] - 1);
                }
            }

            // mapping for nonzero positions
            vector< int > nonzero_mapping;
            for (int i = 0; i < base_code_length; ++i) {
                for (int j = cumwcol[p_column[i] - 1]; j < cumwcol[p_column[i]]; ++j) {
                    nonzero_mapping.push_back(j);
                }
            }

            cout << "n=" << columns << ", r=" << rows << ", g=" << target_girth << ", M=" << tailbite_length << endl;

            // create Tanner graph for the base matrix
            graph< vector_bag > base_matrix_graph;
            int tg_edges = 0;
            for (int c = 0; c < num_input_columns; ++c) {
                for (int r = 0; r < rows; ++r) {
                    if (proper_base_matrix(r, c) != 0) {
                        vector_bag fw(tg_edges++);
                        base_matrix_graph.add_bidi_edge(c, r + num_input_columns, fw, -fw);
                    }
                }
            }

            cout << "Tanner graph size: nt=" << tg_edges << ", rt=" << num_input_columns + rows << endl;

            // Cleaning up the previous state.
            if (state.equations) {
                delete state.equations;
                state.equations = NULL;
                cout << "    [debug] [read_and_construct_equations] Old equations destroyed" << endl;
            }

            state.rows = rows;
            state.columns = columns;
            state.target_girth = target_girth;
            state.tailbite_length = tailbite_length;
            state.exact_tailbite = exact_tailbite;
            state.base_code_length = base_code_length;
            state.num_input_columns = num_input_columns;

            state.HM = input_matrix;
            state.wcol = wcol;
            state.cumwcol = cumwcol;
            state.nonzero_mapping = nonzero_mapping;
            state.p_column = p_column;
            state.equations_graph = base_matrix_graph;

            cout << "Generating equations" << endl;
            state.equations = new equation_builder(base_matrix_graph, target_girth);
            cout << "ok: g=" << state.equations->tree_girth()
                 << ", eqs=" << state.equations->n_equations()
                 << ", nodes=" << state.equations->n_vertices()
                 << endl;
        } else {
            cerr << "Problem opening file '" << input_file << "'.";
            if (state.has_matrix_HM()) {
                cerr << " The existing H and equations are left intact.";
            }
            cerr << endl;
        }
    }

    virtual ~read_and_construct_equations() {}
};

#endif
