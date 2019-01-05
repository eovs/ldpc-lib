#define _CRT_SECURE_NO_WARNINGS // to prohibit VS offering their fscanf_s and friends

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "base_matrix.h"
#include "bp_simulation.h"
#include "code_generation.h"
#include "commons_portable.h"
#include "combinatorics.h"
#include "data_structures.h"
#include "settings.h"

#include "main_modules.h"

using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ifstream;
using std::istringstream;
using std::ostream;
using std::ofstream;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

// The state of the search_ggp.

class ggp_state {
public:
    int rows;
    int columns;
    int target_girth;
    int tailbite_length;
    int exact_tailbite;
    int base_code_length;
    int num_input_columns;

    matrix< int >        HM;               // zero size means not read in.
    vector< int >        wcol, cumwcol;
    vector< int >        nonzero_mapping;
    vector< int >        p_column;
    graph< vector_bag >  equations_graph;
    equation_builder    *equations;        // pointer so far, as no default/copy constructors yet.

    ggp_state()
        : rows(-1)
        , columns(-1)
        , target_girth(-1)
        , tailbite_length(-1)
        , exact_tailbite(-1)
        , base_code_length(-1)
        , num_input_columns(-1)
        , equations(NULL)
    {}

    bool has_matrix_HM() const { return HM.n_rows() > 0; }

    ~ggp_state() {
        if (equations)  { delete equations; equations = NULL; }
    }
};

// A common exception for local interruption.

class local_interruption_exception : public std::exception {};

// A base class for all functions provided.

class ggp_function {
public:
    virtual bool can_run(ggp_state const &state) const = 0;
    virtual void write_info(ggp_state const &state, ostream &out) const = 0;
    virtual void run(ggp_state &state) = 0;
    virtual ~ggp_function() {}
};

// An output-grouping helper
struct ogroup {
    ostream &out1;
    ostream &out2;

    ogroup(ostream &out1, ostream &out2) : out1(out1), out2(out2) {}

    template< typename T > ogroup &operator << (T const &value) {
        out1 << value;
        out2 << value;
        return *this;
    }
};

// We keep the submodules in separate files and include them as follows.
#define INTERNAL_SOURCE_INCLUDE
#include "search_ggp/read_and_construct_equations.cpp"
#include "search_ggp/check_voltages_from_input.cpp"
#include "search_ggp/random_search_bank_good_codes.cpp"
#include "search_ggp/random_extension.cpp"
#include "search_ggp/filter_codes.cpp"
#include "search_ggp/scenario_based_code_generation.cpp"
#undef INTERNAL_SOURCE_INCLUDE

bool read_option_from_args(int &argp, int argc, char *argv[], unsigned &option) {
    if (argp < argc) {
        istringstream iss(argv[argp++]);
        return !!(iss >> option);
    } else {
        return false;
    }
}


int main_search_ggp(int argc, char *argv[]) {
    ggp_state state;
    vector< ggp_function* > functions;

    // All functions are to be added there
    functions.push_back(new read_and_construct_equations(
        get_native_path("input_h.txt")
    ));
    functions.push_back(new check_voltages_from_input(
        get_native_path("input_a.txt"),
        get_native_path("output_a.txt")
    ));
    functions.push_back(new random_search_bank_good_codes(
        get_native_path("data.txt"),
        get_native_path("in_out.txt")
    ));
    functions.push_back(new random_extension(
        get_native_path("data.txt"),
        get_native_path("data_ext.txt"),
        get_native_path("in_out.txt")
    ));
    functions.push_back(new filter_codes(
        get_native_path("data.txt"),
        get_native_path("data_f.txt")
    ));
    functions.push_back(new scenario_based_code_generation("scenario.jsonx"));

    int argp = 1;
 	char *pargv[3] = { "ggp","6","0" };
	int pargc = 3;
    while (true) 
	{
#if 0

        cout << endl;
        cout << "Choose option:" << endl;
        cout << " 0: Quit program." << endl;
        for (unsigned i = 1; i <= functions.size(); ++i) {
            if (functions[i - 1]->can_run(state)) {
                cout << " " << i << ": ";
                functions[i - 1]->write_info(state, cout);
                cout << endl;
            } else {
                cout << " " << i << ": [not available now] ";
                functions[i - 1]->write_info(state, cout);
                cout << endl;
            }
        }
#endif	//0
        unsigned option;
        if (read_option_from_args(argp, pargc, pargv, option) || (cin >> option)) {
            if (option > functions.size()) {
                cout << "Error: input too large" << endl;
            } else if (option == 0) {
                break;
            } else if (functions[option - 1]->can_run(state)) {
                functions[option - 1]->run(state);
            } else {
                cout << "Error: option " << option << " is not available now" << endl;
            }
        } else {
            // what a badass is C++'s library ever: "clear" of a stream does not clear anything!
            cin.clear(cin.rdstate() & ~cin.failbit);
            string tmp;
            cin >> tmp;
            cout << "Error: '" << tmp << "' is not a number" << endl;
        }
    }

    while (functions.size() > 0) {
        delete functions.back();
        functions.pop_back();
    }

    return 0;
}
