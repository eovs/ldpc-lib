#include <algorithm>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <iterator>

#include "main_modules.h"

#include "commons_portable.h"
#include "combinatorics.h"
#include "equations.h"
#include "settings.h"

using std::endl;
using std::set;
using std::istringstream;
using std::make_pair;
using std::ostringstream;
using std::pair;
using std::string;
using std::vector;

template<typename T>
class vector_builder {
    vector< T > data;
public:
    vector_builder(T first) {
        data.push_back(first);
    }
    vector_builder<T> &operator () (T next) {
        data.push_back(next);
        return *this;
    }
    vector< T > operator ()() const {
        return data;
    }
    operator vector< T > () const {
        return data;
    }
};

template<typename T>
vector_builder<T> new_vector(T const &value) {
    return vector_builder<T>(value);
}

template<typename T>
class vector_builder< pair< T, T > > {
    vector< pair< T, T > > data;
public:
    vector_builder(T first, T second) {
        data.push_back(make_pair(first, second));
    }
    vector_builder< pair< T, T > > &operator () (T first, T second) {
        data.push_back(make_pair(first, second));
        return *this;
    }
    vector< pair< T, T > > operator ()() const {
        return data;
    }
    operator vector< pair< T, T > > () const {
        return data;
    }
};

template<typename T>
vector_builder< pair< T, T> > new_vector(T const &first, T const &second) {
    return vector_builder< pair< T, T > >(first, second);
}

void println_1(char const *msg, char const *arg) {
    printf(msg, arg);
    printf("\n");
    fflush(stdout);
}

void println(char const *msg) {
    printf("%s\n", msg);
    fflush(stdout);
}

template<typename T>
std::ostream &operator << (std::ostream &lhs, matrix<T> const &rhs) {
    lhs << "{";
    for (int i = 0, i_max = rhs.n_rows(); i < i_max; ++i) {
        lhs << "{";
        for (int j = 0, j_max = rhs.n_cols(); j < j_max; ++j) {
            lhs << rhs(i, j);
            if (j + 1 != j_max) {
                lhs << ",";
            }
        }
        lhs << "}";
    }
    lhs << "}";
    return lhs;
}

template<typename T>
std::ostream &operator << (std::ostream &lhs, vector<T> const &rhs) {
    lhs << "{";
    for (int i = 0, i_max = (int) rhs.size(); i < i_max; ++i) {
        lhs << rhs[i];
        if (i + 1 != i_max) {
            lhs << ",";
        }
    }
    lhs << "}";
    return lhs;
}

/************************** Combinatorics unit tests ***************************/

void test_binom() {
    println("  test_binom()...");
    if (binom(5, 2) != 10) {
        die("binom(5, 2) = %d != 10", binom(5, 2));
    }
    println("    binom(5, 2) passed");
    if (binom(10, 2) != 45) {
        die("binom(100, 2) = %d != 45", binom(100, 2));
    }
    println("    binom(10, 2) passed");
    if (binom(14, 4) != 1001) {
        die("binom(14, 4) = %d != 1001", binom(14, 4));
    }
    println("    binom(14, 4) passed");
}

void test_divisors() {
    println("  test_divisors()...");
    vector< int > const &div12_found = divisors(12);
    vector< int > div12_expected = new_vector(1)(2)(3)(4)(6)(12);
    if (div12_expected != div12_found) {
        for (int i = 0; i < (int) div12_found.size(); ++i) {
            fprintf(stderr, "%d ", div12_found[i]);
        }
        die("Divisors for 12 are wrong");
    }
    println("    divisors(12) passed");

    vector< int > const &div13_found = divisors(13);
    vector< int > div13_expected = new_vector(1)(13);
    if (div13_expected != div13_found) {
        for (int i = 0; i < (int) div13_found.size(); ++i) {
            fprintf(stderr, "%d ", div13_found[i]);
        }
        die("Divisors for 13 are wrong");
    }
    println("    divisors(13) passed");

    vector< int > const &div74_found = divisors(74);
    vector< int > div74_expected = new_vector(1)(2)(37)(74);
    if (div74_expected != div74_found) {
        for (int i = 0; i < (int) div74_found.size(); ++i) {
            fprintf(stderr, "%d ", div74_found[i]);
        }
        die("Divisors for 74 are wrong");
    }
    println("    divisors(74) passed");
}

void test_addend_splitter() {
    println("  test_addend_splitter()...");
    addend_splitter splitter_13_1(13, 1);
    vector< int > expected_13_1 = new_vector(13);
    vector< int > found_13_1;
    if (splitter_13_1.n_splits() != 1) {
        die("Number of (13, 1) splits = %d != 1", splitter_13_1.n_splits());
    }
    splitter_13_1.split(0, found_13_1);
    if (expected_13_1 != found_13_1) {
        for (int i = 0; i < (int) found_13_1.size(); ++i) {
            fprintf(stderr, "%d ", found_13_1[i]);
        }
        die("0th split of (13, 1) is wrong");
    }
    println("    splits(13, 1) passed");

    addend_splitter splitter_13_2(13, 2);
    vector< int > expected_13_2 = new_vector(6)(7);
    vector< int > found_13_2;
    if (splitter_13_2.n_splits() != 14) {
        die("Number of (13, 2) splits = %d != 14", splitter_13_2.n_splits());
    }
    splitter_13_2.split(6, found_13_2);
    if (expected_13_2 != found_13_2) {
        for (int i = 0; i < (int) found_13_2.size(); ++i) {
            fprintf(stderr, "%d ", found_13_2[i]);
        }
        die("6th split of (13, 2) is wrong");
    }
    println("    splits(13, 2) passed");

    addend_splitter splitter_13_3(13, 3);
    vector< int > expected_13_3 = new_vector(2)(3)(8);
    vector< int > found_13_3;
    if (splitter_13_3.n_splits() != 105) {
        die("Number of (13, 3) splits = %d != 105", splitter_13_3.n_splits());
    }
    splitter_13_3.split(30, found_13_3);
    if (expected_13_3 != found_13_3) {
        for (int i = 0; i < (int) found_13_3.size(); ++i) {
            fprintf(stderr, "%d ", found_13_3[i]);
        }
        die("30th split of (13, 3) is wrong");
    }
    println("    splits(13, 3) passed");
}

/************************** Equation unit tests ***************************/

void test_one_equation(char const *name, char const *test, int vertices, int equations) {
    int R, N, G, M;
    istringstream inputh(test);
    inputh >> R >> N >> G >> M;
    vector< vector< int > > incident_to_edges(N);
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < N; ++j) {
            int cell;
            inputh >> cell;
            if (cell == 1) {
                incident_to_edges[j].push_back(i);
            }
        }
    }
    graph< vector_bag > g;
    int edge_count = 0;
    for (int edge = 0; edge < N; ++edge) {
        vector< int > const &nonzero = incident_to_edges[edge];
        for (int i = 0; i < (int) nonzero.size(); ++i) {
            vector_bag fw(edge_count++);
            g.add_bidi_edge(edge, N + nonzero[i], fw, -fw);
        }
    }
    equation_builder b(g, G);
    if (b.n_vertices() != vertices) {
        die("Equations for test %s are wrong: vertices = %d != %d", name, b.n_vertices(), vertices);
    }
    if (b.n_equations() != equations) {
        die("Equations for test %s are wrong: equations = %d != %d", name, b.n_equations(), equations);
    }
    println_1("    %s passed", name);
}

void test_equations() {
    println("  test_equations()...");

    test_one_equation("graph-9", "9 12 14 10000\n"
                                 "1 0 0 1 0 0 1 0 0 1 0 0\n"
                                 "1 0 0 0 1 0 0 1 0 0 1 0\n"
                                 "1 0 0 0 0 1 0 0 1 0 0 1\n"
                                 "0 1 0 1 0 0 0 0 1 0 1 0\n"
                                 "0 1 0 0 1 0 1 0 0 0 0 1\n"
                                 "0 1 0 0 0 1 0 1 0 1 0 0\n"
                                 "0 0 1 1 0 0 0 1 0 0 0 1\n"
                                 "0 0 1 0 1 0 0 0 1 1 0 0\n"
                                 "0 0 1 0 0 1 1 0 0 0 1 0\n",
                      1736, 5058);
    test_one_equation("graph-12", "12 24 10 10000\n"
                                  "0 1 1 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0\n"
                                  "0 1 0 0 0 1 1 1 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 0\n"
                                  "0 0 0 1 1 1 0 1 0 0 0 1 0 0 1 1 0 0 0 0 0 0 0 0\n"
                                  "1 0 1 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 0 0\n"
                                  "0 0 1 0 0 0 1 0 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 0\n"
                                  "0 0 0 0 1 1 0 1 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 0\n"
                                  "0 0 1 1 0 0 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 0 0 0\n"
                                  "0 1 1 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0\n"
                                  "1 0 0 0 1 1 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0\n"
                                  "0 0 0 0 0 1 0 1 0 0 1 1 0 0 0 0 0 0 0 0 0 1 1 0\n"
                                  "0 0 1 1 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1\n"
                                  "1 0 0 0 0 1 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1\n",
                      2491, 7234);
}

/************************** vector_bag unit tests ***************************/

void check_fail(vector< pair< int, int > > const &expected, vector< pair< int, int > > const &found) {
    if (expected != found) {
        printf("Expected:");
        for (int i = 0; i < (int) expected.size(); ++i) {
            printf(" (%d; %d)", expected[i].first, expected[i].second);
        }
        printf("\nFound:");
        for (int i = 0; i < (int) found.size(); ++i) {
            printf(" (%d; %d)", found[i].first, found[i].second);
        }
        printf("\n");
        die("Unequal vectors");
    }
}

void test_vector_bag() {
    println("  test_vector_bag()...");

    vector< pair< int, int > > control;

    vector_bag b(4);

    control.clear();
    control.push_back(make_pair(4, 1));
    check_fail(control, b.data());

    b += vector_bag(3);

    control.clear();
    control.push_back(make_pair(3, 1));
    control.push_back(make_pair(4, 1));
    check_fail(control, b.data());

    b += vector_bag(4);

    control.clear();
    control.push_back(make_pair(3, 1));
    control.push_back(make_pair(4, 2));
    check_fail(control, b.data());

    b -= vector_bag(3);

    control.clear();
    control.push_back(make_pair(4, 2));
    check_fail(control, b.data());

    b -= vector_bag(3);

    control.clear();
    control.push_back(make_pair(3, -1));
    control.push_back(make_pair(4, 2));
    check_fail(control, b.data());

    vector_bag another(3);
    another -= vector_bag(4);
    b += another;
    control.clear();
    control.push_back(make_pair(4, 1));
    check_fail(control, b.data());
}

template<typename T>
void check_setting(settings const &s, T const &expected)
{
    T found;
    s.cast_to(found);
    if (found != expected) {
        ostringstream oss;
        oss << "For setting " << s.get_full_path() << " expected is '" << expected << "' but found '" << found << "'" << endl;
        println(oss.str().c_str());
        die("Settings are parsed wrong");
    }
}

void check_setting(settings const &s, char const *expected)
{
    check_setting(s, string(expected));
}

void test_settings() {
    println("  test_settings()...");

    settings s = settings::from_file("files/tests/main.jsonx");

    check_setting(s.select("integer_number"), 3);
    check_setting(s.select("integer_number"), 3.0); // should also parse like this
    check_setting(s.select("integer_number"), "3"); // ... and like this

    check_setting(s.select("real_number"), 3.1415);
    check_setting(s.select("real_number"), "3.1415");

    check_setting(s.select("example_string"), "test test test");
    check_setting(s.select("example_integer_string"), "42");
    check_setting(s.select("example_integer_string"), 42);

    check_setting(s.select("simple_array"), new_vector(123)(234)(345)(456)());
    check_setting(s.select("simple_array"), new_vector(123.0)(234.0)(345.0)(456.0)());

    check_setting(s.select("external_array"), new_vector(10)(12)(14)(16)(18)(20)());
    check_setting(s.select("external_array"), new_vector(string("10"))("12")("14")("16")("18")("20")());

    matrix<int> simple_matrix(3, 4);
    simple_matrix(0, 0) = 1;
    simple_matrix(0, 1) = 2;
    simple_matrix(0, 2) = 3;
    simple_matrix(0, 3) = 4;
    simple_matrix(1, 0) = 2;
    simple_matrix(1, 1) = 3;
    simple_matrix(1, 2) = 4;
    simple_matrix(1, 3) = 5;
    simple_matrix(2, 0) = 3;
    simple_matrix(2, 1) = 4;
    simple_matrix(2, 2) = 5;
    simple_matrix(2, 3) = 6;
    check_setting(s.select("simple_matrix"), simple_matrix);

    matrix<int> sparse_matrix(4, 3);
    sparse_matrix(0, 0) = 2;
    sparse_matrix(1, 1) = 3;
    sparse_matrix(2, 2) = 4;
    sparse_matrix(2, 0) = -5;
    check_setting(s.select("sparse_matrix"), sparse_matrix);

    matrix<string> sparse_matrix_2(4, 3);
    sparse_matrix_2(0, 0) = "2";
    sparse_matrix_2(1, 1) = "3";
    sparse_matrix_2(2, 2) = "4";
    sparse_matrix_2(2, 0) = "-5";
    check_setting(s.select("sparse_matrix"), sparse_matrix_2);

    matrix<string> string_matrix(4, 3);
    string_matrix(1, 0) = "seventeen";
    string_matrix(2, 1) = "qwerty";
    string_matrix(3, 2) = "white space";
    check_setting(s.select("string_matrix"), string_matrix);

    check_setting(s.select("objref"), 42);

    check_setting(s.select("struct/a"), 1);
    check_setting(s.select("struct/b"), 47.0);
    check_setting(s.select("struct/c"), new_vector(2)(3)());
    // these come from "defaults"
    check_setting(s.select("struct/speed_of_light"), 299792458);
    check_setting(s.select("struct/proton_mass"), 1.6726219e-27);
    check_setting(s.select("struct/one"), "1");
}

void check_interval(string const &str, vector< int > const &expected) {
    vector< int > found;
    try {
        read_intervals(str, std::back_inserter(found));
    } catch (interval_exception &ex) {
        die("Test '%s' failed with exception: '%s'", str.c_str(), ex.what());
    }
    if (found != expected) {
        std::ostringstream oss;
        oss << "Test '" << str << "' failed: expected [";
        for (unsigned i = 0, i_max = (unsigned int)expected.size(); i < i_max; ++i) {
            oss << expected[i];
            if (i + 1 != i_max) {
                oss << " ";
            }
        }
        oss << "] found [";
        for (unsigned i = 0, i_max = (unsigned int)found.size(); i < i_max; ++i) {
            oss << found[i];
            if (i + 1 != i_max) {
                oss << " ";
            }
        }
        oss << "]";
        std::string str = oss.str();
        die(str.c_str());
    }
}

void check_interval_fails(string const &str) {
    try {
        vector< int > found;
        read_intervals(str, std::back_inserter(found));
        die("Test '%s' expected to throw exception but did not", str.c_str());
    } catch (interval_exception &) {}
}

void test_intervals() {
    println("  test_intervals()...");
    check_interval("1", new_vector(1)());
    check_interval("1000", new_vector(1000)());
    check_interval("1-1", new_vector(1)());
    check_interval("1-2", new_vector(1)(2)());
    check_interval("1-3", new_vector(1)(2)(3)());
    check_interval("9-11", new_vector(9)(10)(11)());
    check_interval("1 3", new_vector(1)(3)());
    check_interval("3 1", new_vector(3)(1)());
    check_interval("1 1-3", new_vector(1)(1)(2)(3)());
    check_interval("1-3 6-7", new_vector(1)(2)(3)(6)(7)());
    check_interval("6-8 1-3", new_vector(6)(7)(8)(1)(2)(3)());
    check_interval("1    -     3", new_vector(1)(2)(3)());
    check_interval(" 2   4-6 8 10- 12  ", new_vector(2)(4)(5)(6)(8)(10)(11)(12)());
    check_interval("", vector< int >());
    check_interval("   ", vector< int >());
    check_interval_fails("-");
    check_interval_fails("-1");
    check_interval_fails("1-");
    check_interval_fails("8-1");
    check_interval_fails("1-2-3");
    check_interval_fails("1-2 -");
}

void check_knapsack(double count_solutions, int num_elements, int total_weight, vector< pair< int, int > > const &constraints) {
    // Prepare the failure message in advance
    std::ostringstream oss;
    oss << "Test failed:\n";
    oss << "    num_elements = " << num_elements << "\n";
    oss << "    total_weight = " << total_weight << "\n";
    oss << "    constraints = {\n";
    for (unsigned i = 0; i < constraints.size(); ++i) {
        oss << "        " << i << " => [" << constraints[i].first << ";" << constraints[i].second << "]\n";
    }
    oss << "    }\n";

    // Check the number of solutions
    random_knapsack_solution sol(constraints, total_weight, num_elements);

    if (count_solutions != sol.count_solutions()) {
        oss << "Expected " << count_solutions << " found " << sol.count_solutions() << "\n";
        std::string str = oss.str();
        die(str.c_str());
    }

    if (count_solutions == 0) {
        return;
    }

    set< vector< int > > different_weights;

    // Okay, start checking random solutions
    for (int time = 0; time < 10000; ++time) {
        vector< int > weights(num_elements);
        sol.sample_solution(weights);
        different_weights.insert(weights);
        int my_total_weight = 0;
        for (int i = 0; i < num_elements; ++i) {
            my_total_weight += weights[i];
        }
        if (my_total_weight != total_weight) {
            oss << "Sample #" << time + 1 << ": total weight is wrong.\n";
            oss << "    Expected " << total_weight << " found " << my_total_weight << "\n";
            oss << "    Weights: [";
            for (int i = 0; i < num_elements; ++i) {
                oss << weights[i];
                if (i + 1 != num_elements) {
                    oss << " ";
                }
            }
            oss << "]\n";
            std::string str = oss.str();
            die(str.c_str());
        }

        for (int w = 0, w_max = (int) (constraints.size()); w < w_max; ++w) {
            int countw = (int)(std::count(weights.begin(), weights.end(), w));
            if (countw < constraints[w].first || countw > constraints[w].second) {
                oss << "Sample #" << time + 1 << ": constraint for w = " << w << " is violated\n";
                oss << "    Expected [" << constraints[w].first << "; " << constraints[w].second << "] found " << countw << "\n";
                oss << "    Weights: [";
                for (int i = 0; i < num_elements; ++i) {
                    oss << weights[i];
                    if (i + 1 != num_elements) {
                        oss << " ";
                    }
                }
                oss << "]\n";
                std::string str = oss.str();
                die(str.c_str());
            }
        }
    }

    if (count_solutions > 1 && different_weights.size() == 1) {
        oss << "Random solution sampling is not random\n";
        std::string str = oss.str();
        die(str.c_str());
    }
    printf("    count_solutions = %lf => %u different weights\n", count_solutions, (unsigned)different_weights.size());
}

void test_random_knapsack() {
    println("  test_random_knapsack()...");
    check_knapsack(1, 10, 10, new_vector(0, 0)(0, 20)());
    check_knapsack(1, 10, 20, new_vector(0, 0)(0, 0)(0, 20)());
    check_knapsack(0, 10, 20, new_vector(0, 0)(0, 20)());
    check_knapsack(1, 10, 11, new_vector(0, 0)(0, 20)(0, 1)());
    check_knapsack(1, 10, 11, new_vector(0, 0)(0, 20)(1, 1)());
    check_knapsack(0, 10, 11, new_vector(0, 0)(0, 20)(0, 0)());
    check_knapsack(2, 2, 9, new_vector(0, 0)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)());
    check_knapsack(3, 2, 9, new_vector(0, 0)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)());
    check_knapsack(0, 16,  90, new_vector(0, 0)(0, 0)(0, 0)(0, 0)(0, 0)(0, 0)(0, 20)(0, 20)(0, 20)(0, 20)(0, 0)(0, 0)(0, 0)(0, 20)());
    // I have not checked the numbers of solutions here, but they seem legitimate
    check_knapsack(14,     16, 105, new_vector(0, 0)(0, 0)(0, 0)(0, 0)(0, 0)(0, 0)(0, 20)(0, 20)(0, 20)(0, 20)(0, 0)(0, 0)(0, 0)(0, 20)());
    check_knapsack(117,    16, 130, new_vector(0, 0)(0, 0)(0, 0)(0, 0)(0, 0)(0, 0)(0, 20)(0, 20)(0, 20)(0, 20)(0, 0)(0, 0)(0, 0)(0, 20)());
    check_knapsack(398811, 16, 130, new_vector(0, 0)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)(0, 20)());
}

int main_unit_tests(int, char *[]) {
    println("Starting unit tests...");
    test_binom();
    test_divisors();
    test_addend_splitter();
    test_equations();
    test_vector_bag();
    test_settings();
    test_intervals();
    test_random_knapsack();
    println("All tests OK");

    return 0;
}
