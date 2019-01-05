#include <algorithm>
#include <exception>
#include <limits>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <utility>
#include <queue>
#include <iterator>

#include "commons_portable.h"
#include "data_structures.h"
#include "combinatorics.h"
#include "base_matrix.h"
#include "graph_spider.h"
#include "settings.h"

#undef min
#undef max

const int HIGHLY_IMPROBABLE_INVERSE_PROBABILITY = 1000000;

using std::set;
using std::vector;
using std::string;
using std::pair;
using std::map;
using std::make_pair;
using std::min;
using std::max;

struct cp_random {
    int operator ()(int arg) const {
        return next_random_int(0, arg);
    }
} cp_random_instance;

pair< int, int > intersect_pairs(pair< int, int > const &l, pair< int, int > const &r) {
    int lb = max(l.first, r.first);
    int rb = min(l.second, r.second);
    if (lb > rb) {
        die("Intersection of [%d; %d] and [%d; %d] is empty", l.first, l.second, r.first, r.second);
    }
    return make_pair(lb, rb);
}

pair< bool, vector< int > > can_assign_information_columns(vector< int > info_weights, int current_column, vector< int > row_weights, vector< int > const &selected) {
    for (unsigned i = 0, i_max = (unsigned int)selected.size(); i < i_max; ++i) {
        row_weights[selected[i]]--;
    }

    int sum_rows = 0;
    for (unsigned i = 0, i_max = (unsigned int)row_weights.size(); i < i_max; ++i) {
        if (row_weights[i] < 0) {
            die("row_weights[%d] = %d", i, row_weights[i]);
        }
        sum_rows += row_weights[i];
    }
    int sum_cols = 0;
    for (int i = (int) info_weights.size() - 1; i > current_column; --i) {
        if (info_weights[i] < 0) {
            die("info_weights[%d] = %d", i, info_weights[i]);
        }
        sum_cols += info_weights[i];
    }
//    if (sum_rows != sum_cols) {
//        die("wrong arguments to can_assign: sum_rows = %d, sum_cols = %d", sum_rows, sum_cols);
//    }

    std::sort(info_weights.begin(), info_weights.end());
    std::priority_queue< pair< int, int > > queue;
    for (unsigned i = 0, i_max = (unsigned int)row_weights.size(); i < i_max; ++i) {
        if (row_weights[i] > 0) {
            queue.push(make_pair(row_weights[i], i));
        }
    }

    for (int i = (int) info_weights.size() - 1; i > current_column; --i) {
        vector< pair< int, int > > new_weights;
        if ((int) queue.size() < info_weights[i]) {
            return make_pair(false, std::vector< int >());
        }
        for (int j = 0; j < info_weights[i]; ++j) {
            pair< int, int > curr = queue.top();
            queue.pop();
            if (curr.first > 1) {
                new_weights.push_back(make_pair(curr.first - 1, curr.second));
            }
        }
        for (int j = 0; j < (int) new_weights.size(); ++j) {
            queue.push(new_weights[j]);
        }
    }
    set< int > used_columns;
    vector< int > rv;
    while (queue.size() > 0) {
        pair< int, int > top = queue.top();
        queue.pop();
        if (used_columns.count(top.second) == 0) {
            used_columns.insert(top.second);
            rv.push_back(top.second);
        }
    }
//    if (queue.size() > 0) {
//        int qs = queue.size();
//        printf("[DEBUG] queue:");
//        while (queue.size()) {
//            int top = queue.top();
//            queue.pop();
//            printf(" %d", top);
//        }
//        printf("\n");
//        die("can_assign is broken: queue.size() is %d", qs);
//    }
    return make_pair(true, rv);
}

struct ace_evaluator : graph_spider< int > {
    int girth;
    int min_ace;
    int num_cycles;

    ace_evaluator(graph< int > g)
        : graph_spider< int >(g)
        , girth(std::numeric_limits< int >::max())
        , min_ace(std::numeric_limits< int >::max())
    {
        run_spider();
    }

    virtual void new_cycle_hook(int smaller_end, int larger_end, int cycle_length) {
        if (girth > cycle_length) {
            girth = cycle_length;
            min_ace = std::numeric_limits< int >::max();
        }
        int current_ace = nodes[smaller_end].path_info + nodes[larger_end].path_info;
        if (min_ace > current_ace) {
            min_ace = current_ace;
            num_cycles = 0;
        }
        if (min_ace == current_ace) {
            ++num_cycles;
        }
    }

    virtual bool makes_sense_processing(int cycle_length) {
        return cycle_length <= girth;
    }

    virtual ~ace_evaluator() {}
};

bool greater_equal_componentwise(vector< int > const &bigger, vector< int > const &smaller) {
    if (bigger.size() != smaller.size()) {
        die("greater_equal_componentwise: vector sizes do not match: %d and %d", (int) bigger.size(), (int) smaller.size());
    }
    for (unsigned i = 0, i_max = (unsigned int)bigger.size(); i < i_max; ++i) {
        if (bigger[i] < smaller[i]) {
            return false;
        }
    }
    return true;
}

typedef void (*block_fun_t) (matrix< int > &, int, int, settings const &);

struct base_matrix_interruption_exception : std::exception {};

void block_bidiagonal(matrix< int > &mx, int offset, int height, settings const &s) {
    for (int i = 0; i < height; ++i) {
        mx(offset + i, offset + i) = 2;
        if (i + 1 < height) {
            mx(offset + i + 1, offset + i) = 2;
        }
    }

    // special column
    int span = height;
    if (s.can_select("special_column_span")) {
        s.select("special_column_span").cast_to(span);
    }
    int mid = ((offset + height - span) + (offset + height - 1)) / 2;
    mx(offset + height - span, offset + height - 1) = 2;
    mx(mid, offset + height - 1) = 1;
}

void block_unidiagonal(matrix< int > &mx, int offset, int height, settings const &) {
    for (int i = 0; i < height; ++i) {
        mx(offset + i, offset + i) = 2;
    }
}

map< string, block_fun_t > block_algos;

base_matrix random_base_matrix(
    int num_check_symbols,
    int code_max_length,
    vector< int > const &information_column_weights,
    int num_trials_base,
    int num_trials_column,
    settings const &checking_column_settings,
    settings const &row_constraints
) {
    bool debug = false;
    if (row_constraints.can_select("debug")) {
        row_constraints.select("debug").cast_to(debug);
    }

    // 0. If first time here, create all algorithms.
    if (block_algos.size() == 0) {
        block_algos.insert(make_pair(string("bidiagonal"),  &block_bidiagonal));
        block_algos.insert(make_pair(string("unidiagonal"), &block_unidiagonal));
    }

    // 1. Check and determine the block sizes.
    vector< int > block_sizes;
    vector< settings > block_settings;

    checking_column_settings.select("blocks").cast_to(block_settings);

    int sum = 0, n_zeros = 0;
    unsigned zero_idx = -1;
    for (unsigned i = 0, i_max = (unsigned int)block_settings.size(); i < i_max; ++i) {
        block_sizes.push_back(0);
        if (block_settings[i].can_select("height")) {
            block_settings[i].select("height").cast_to(block_sizes.back());
            sum += block_sizes.back();
        } else {
            ++n_zeros;
            zero_idx = i;
        }
    }
    if (n_zeros > 1) {
        die("Illegal base matrix configuration (checking column settings): more than one block without height");
    }
    if (n_zeros == 1 && sum >= num_check_symbols) {
        die("Illegal base matrix configuration (checking column settings): missing height is to be restored to %d", num_check_symbols - sum);
    }
    if (n_zeros == 0 && sum != num_check_symbols) {
        die("Illegal base matrix configuration (checking column settings): sum of heights is %d but should be %d", sum, num_check_symbols);
    }
    if (n_zeros == 1) {
        block_sizes[zero_idx] = num_check_symbols - sum;
    }

    vector< int > block_offsets(block_sizes.size());
    for (unsigned i = 0, i_max = (unsigned int)block_sizes.size(); i + 1 < i_max; ++i) {
        block_offsets[i + 1] = block_offsets[i] + block_sizes[i];
    }

    // 2. Initialize the diagonals.
    matrix< int > mx_init(num_check_symbols, code_max_length);
    for (unsigned i = 0, i_max = (unsigned int)block_settings.size(); i < i_max; ++i) {
        string block_type;
        block_settings[i].select("type").cast_to(block_type);
        map< string, block_fun_t >::const_iterator itr = block_algos.find(block_type);
        if (itr == block_algos.end()) {
            die("Unknown block type: %s", block_type.c_str());
        }
        itr->second(mx_init, block_offsets[i], block_sizes[i], block_settings[i]);
    }

    if (debug) {
        printf("--------------------------\n");
        for (int row = 0; row < num_check_symbols; ++row) {
            for (int col = 0; col < code_max_length; ++col) {
                printf("%d", mx_init(row, col));
            }
            printf("\n");
        }
        printf("--------------------------\n");
    }

    // 3. Reading and checking the column weight patches
    vector< int > checking_column_patches_weight(num_check_symbols);
    vector< int > checking_column_patches_maxind(num_check_symbols);
    vector< settings > patch_settings;
    if (checking_column_settings.can_select("patches")) {
        checking_column_settings.select("patches").cast_to(patch_settings);

        for (unsigned i = 0, i_max = (unsigned int)patch_settings.size(); i < i_max; ++i) {
            string col;
            int weight;
            vector< int > indices;
            patch_settings[i].select("col").cast_to(col);
            patch_settings[i].select("weight").cast_to(weight);
            read_intervals(col, std::back_inserter(indices));
            for (unsigned j = 0, j_max = (unsigned int)indices.size(); j < j_max; ++j) {
                int val = indices[j];
                if (val <= 0 || val > num_check_symbols) {
                    die("Illegal base matrix configuration (checking column settings): patch contains column %d for %d checking symbols", val, num_check_symbols);
                }
                if (checking_column_patches_weight[val - 1] == 0) {
                    checking_column_patches_weight[val - 1] = weight;
                } else if (checking_column_patches_weight[val - 1] != weight) {
                    die ("Illegal base matrix configuration (checking column settings): conflicting patches for column %d: weights %d and %d",
                        val, checking_column_patches_weight[val - 1], weight);
                }
            }
        }
    }
    int checking_column_mandatory_weight = 0;
    int checking_column_patches_sumweight = 0;
    for (int col = 0; col < num_check_symbols; ++col) {
        int count_nonzero = 0;
        for (int i = 0; i < num_check_symbols; ++i) {
            if (mx_init(i, col) != 0) {
                ++count_nonzero;
            }
        }
        checking_column_mandatory_weight += count_nonzero;
        for (unsigned i = 0, i_max = (unsigned int)block_sizes.size(); i < i_max; ++i) {
            if (checking_column_patches_maxind[col] + block_sizes[i] < col) {
                checking_column_patches_maxind[col] += block_sizes[i];
            } else {
                break;
            }
        }
        if (checking_column_patches_weight[col] != 0) {
            if (checking_column_patches_weight[col] < count_nonzero) {
                die("Illegal base matrix configuration (checking column settings): column %d patched to weight %d but already has weight %d",
                    col + 1, checking_column_patches_weight[col], count_nonzero);
            }
            checking_column_patches_weight[col] -= count_nonzero;
        }
        if (checking_column_patches_weight[col] > checking_column_patches_maxind[col]) {
            die("Illegal base matrix configuration (checking column settings): column %d patched to weight %d but has space only for %d",
                col + 1, checking_column_patches_weight[col], checking_column_patches_maxind[col]);
        }
        checking_column_patches_sumweight += checking_column_patches_weight[col];
    }

    // 4. Configuring row generator
    // 4.1. Creating the row weight limits
    vector< pair< int, int > > row_weight_limits;
    {
        // 4.1.1. Reading possible row weights. If not stated, use all weights
        vector< int > possible_row_weights;
        if (row_constraints.can_select("possible_weights")) {
            row_constraints.select("possible_weights").cast_to(possible_row_weights);
            for (int i = 0; i <= code_max_length; ++i) {
                row_weight_limits.push_back(make_pair(0, 0));
            }
            for (unsigned i = 0, i_max = (unsigned int)possible_row_weights.size(); i < i_max; ++i) {
                int value = possible_row_weights[i];
                if (value >= 0 && value <= code_max_length) {
                    row_weight_limits[value] = make_pair(0, num_check_symbols);
                }
            }
        } else {
            for (int i = 0; i <= code_max_length; ++i) {
                row_weight_limits.push_back(make_pair(0, i == 0 ? 0 : num_check_symbols));
            }
        }

        // 4.1.2. Reading and applying constraints.
        vector< vector< int > > row_weight_constraints;
        if (row_constraints.can_select("weight_constraints")) {
            row_constraints.select("weight_constraints").cast_to(row_weight_constraints);
            for (unsigned t = 0, t_max = (unsigned int)row_weight_constraints.size(); t < t_max; ++t) {
                vector< int > const &rwc = row_weight_constraints[t];
                if (rwc.size() != 3) {
                    die("Row constraints element has length %d, which is not 3", (int) (rwc.size()));
                }
                int index = rwc[0];
                if (0 <= index && index < code_max_length) {
                    row_weight_limits[index] = intersect_pairs(row_weight_limits[index], make_pair(rwc[1], rwc[2]));
                }
            }
        }
    }

    vector< int > real_information_column_weights(code_max_length - num_check_symbols);
    for (int w = 0, i = 0, w_max = (int) information_column_weights.size(); w < w_max; ++w) {
        for (int t = 0; t < information_column_weights[w]; ++t, ++i) {
            real_information_column_weights[i] = w;
        }
    }

    // 4.2, Evaluating sum of the weights
    int all_weights = checking_column_mandatory_weight + checking_column_patches_sumweight;
    for (unsigned i = 0, i_max = (unsigned int)real_information_column_weights.size(); i < i_max; ++i) {
        all_weights += real_information_column_weights[i];
    }

    // 4.2.1. Evaluating the minimum and maximum possible total row weights
    int min_weight_rowwise, max_weight_rowwise;
    bool satisfiable_weight_rowwise = true;
    {
        int weight_baseline = 0;
        int count_baseline = 0;
        for (unsigned w = 0, w_max = row_weight_limits.size(); w < w_max; ++w) {
            count_baseline += row_weight_limits[w].first;
            weight_baseline += row_weight_limits[w].first * w;
        }
        if (count_baseline > num_check_symbols) {
            satisfiable_weight_rowwise = false;
            min_weight_rowwise = 1;
            max_weight_rowwise = 0;
        } else {
            min_weight_rowwise = weight_baseline;
            int min_weight_count = count_baseline;
            for (unsigned w = 0, w_max = row_weight_limits.size(); w < w_max; ++w) {
                for (int t = row_weight_limits[w].first + 1; t <= row_weight_limits[w].second && min_weight_count < num_check_symbols; ++t) {
                    ++min_weight_count;
                    min_weight_rowwise += w;
                }
            }
            max_weight_rowwise = weight_baseline;
            int max_weight_count = count_baseline;
            for (unsigned w = row_weight_limits.size(); w > 0; --w) {
                for (int t = row_weight_limits[w - 1].first + 1; t <= row_weight_limits[w - 1].second && max_weight_count < num_check_symbols; ++t) {
                    ++max_weight_count;
                    max_weight_rowwise += w - 1;
                }
            }
        }
    }

    // 4.3. Running a knapsack solver on the row constraints.
    random_knapsack_solution knapsack(row_weight_limits, all_weights, num_check_symbols);
    if (!knapsack.has_solution()) {
        base_matrix bm;
        bm.solution_exists = false;

        if (!satisfiable_weight_rowwise) {
            printf("[base_matrix] Warning: row weight constraints are not satisfiable!\n");
        } else if (all_weights < min_weight_rowwise) {
            printf("[base_matrix] Warning: total column weight = %d, but minimum possible row weight = %d!\n", all_weights, min_weight_rowwise);
        } else if (all_weights > max_weight_rowwise) {
            printf("[base_matrix] Warning: total column weight = %d, but maximum possible row weight = %d!\n", all_weights, max_weight_rowwise);
        } else {
            printf("[base_matrix] Warning: knapsack found no solution!\n");
        }

        return bm;
    }

    // 4.4. Evaluating lower bounds on weights of particular rows
    vector< int > row_value_lower_bounds(num_check_symbols);
    for (int row = 0; row < num_check_symbols; ++row) {
        for (int col = 0; col < code_max_length; ++col) {
            if (mx_init(row, col)) {
                ++row_value_lower_bounds[row];
            }
        }
    }

    // 5. Making trials...
    base_matrix rv;
    rv.solution_exists = false;

    for (int trial_base = 0; trial_base < num_trials_base && (!rv.solution_exists || rv.score[1] > 1); ++trial_base) {
        try {
            console_exception_hook x_hook('x', "interrupts current base matrix generation", base_matrix_interruption_exception());
            if (debug) {
                std::ostringstream matrix_message;
                matrix_message << "\n=================\nCurrent matrix:";
                if (rv.solution_exists) {
                    matrix_message << "\n";
                    for (int row = 0; row < num_check_symbols; ++row) {
                        for (int col = 0; col < code_max_length; ++col) {
                            matrix_message << rv.data(row, col);
                        }
                        matrix_message << "\n";
                    }
                } else {
                    matrix_message << " none so far\n";
                }
                matrix_message << "=================\n";

                console_message_hook m_hook('p', "prints the current base matrix", matrix_message.str());
                console_check_hooks();
            }

            // 5.0. Matrix correlation statistics.
            int total_t2 = 0, total_t3 = 0;

            // 5.1. Synthesizing row weights to satisfy everyone.
            vector< int > row_values(num_check_symbols);
            int iterations = 0;
            // This could be better...
            do {
                knapsack.sample_solution(row_values);
                std::random_shuffle(row_values.begin(), row_values.end(), cp_random_instance);
            } while (!greater_equal_componentwise(row_values, row_value_lower_bounds) &&
                     ++iterations < HIGHLY_IMPROBABLE_INVERSE_PROBABILITY);

            if (iterations == HIGHLY_IMPROBABLE_INVERSE_PROBABILITY) {
                printf("Some non-working combination found, breaking current matrix generation\n");
                return rv;
            }

            for (int i = 0; i < num_check_symbols; ++i) {
                row_values[i] -= row_value_lower_bounds[i];
            }

            // 5.2. Filling the matrix's remainings, column by column
            matrix< int > mx = mx_init;
            // correlations(i, j) - how many ones are simultaneously in rows i and j
            matrix< int > correlations(num_check_symbols, num_check_symbols);
            vector< int > information_row_weights(num_check_symbols);

            for (int col = 0; col < code_max_length; ++col) {

                // 5.2.1. Collecting how much elements we need to choose, and from which region.
                int max_row_index = col < num_check_symbols ? checking_column_patches_maxind[col] : num_check_symbols;
                int how_much_bits = col < num_check_symbols ? checking_column_patches_weight[col] : real_information_column_weights[col - num_check_symbols];
                int trials_column = col < num_check_symbols ? 1 : num_trials_column;

                // 5.2.2. Determining valid row indices for our elements.
                vector< int > valid_indices;

                if (col == num_check_symbols && !can_assign_information_columns(real_information_column_weights, -1, row_values, valid_indices).first) {
                    if (debug) {
                        printf("[DEBUG] [base_matrix] Everything is bad from the very beginning\n");
                    }
                    break;
                }

                for (int i = 0; i < max_row_index; ++i) {
                    if (row_values[i] > 0) {
                        valid_indices.push_back(i);
                    }
                }
                if ((int) valid_indices.size() < how_much_bits) {
                    die("Something is broken in base matrix generation: no satisfiable assignment was found in the middle of the search");
                }

                // 5.2.3. Generating and testing random selections.
                vector< int > best_indices;
                int best_t2 = std::numeric_limits< int >::max();
                int best_t3 = std::numeric_limits< int >::max();
                for (int trial_column = 0; trial_column < trials_column && (best_t2 | best_t3) != 0; ++trial_column) {
                    console_check_hooks();

                    // 5.2.3.1. Random selection.
                    vector< int > selected_indices;
                    random_binom_realization((int)valid_indices.size(), how_much_bits, selected_indices);
                    for (int i = 0; i < how_much_bits; ++i) {
                        selected_indices[i] = valid_indices[selected_indices[i]];
                    }

                    if (!can_assign_information_columns(real_information_column_weights, max(0, col - num_check_symbols), row_values, selected_indices).first) {
                        if (debug) {
                            printf("[DEBUG] [base_matrix] col = %d, assignment failed, continuing\n", col);
                        }
                        if (trial_column + 1 < trials_column) {
                            continue;
                        }
                        if (debug) {
                            printf("[DEBUG] last attempt: using deterministic assignment\n");
                        }
                        vector< int > empty;
                        selected_indices = can_assign_information_columns(real_information_column_weights, max(0, col - num_check_symbols), row_values, empty).second;
                        if ((int) selected_indices.size() < how_much_bits) {
                            die("deterministic assignment fails! expected %d, found %d\n", how_much_bits, (int) selected_indices.size());
                        }
                        while ((int) selected_indices.size() > how_much_bits) {
                            selected_indices.pop_back();
                        }
                    }

                    // 5.2.3.2. Correlation heuristics.
                    int curr_t2 = 0;
                    int curr_t3 = 0;
                    for (int i = 0; i < how_much_bits; ++i) {
                        int sum = 0;
                        for (int j = i + 1; j < how_much_bits; ++j) {
                            sum += correlations(selected_indices[i], selected_indices[j]);
                        }
                        curr_t2 += sum;
                        curr_t3 += sum * (sum - 1) / 2;
                        int info_w = information_row_weights[selected_indices[i]];
                        curr_t3 += info_w * (info_w - 1) / 2 * (info_w - 2) / 3;
                    }

                    // 5.2.3.3. Updating the best.
                    if (curr_t3 < best_t3 || (curr_t3 == best_t3 && curr_t2 < best_t2)) {
                        best_t3 = curr_t3;
                        best_t2 = curr_t2;
                        best_indices = selected_indices;
                    }
                }

                if ((int)best_indices.size() < how_much_bits) {
                    if (debug) {
                        printf("[DEBUG] [base_matrix] col = %d, no satisfiable column found, breaking\n", col);
                    }
                    break;
                }

                // 5.2.4. Writing the new column and updating statistics
                for (int i = 0; i < how_much_bits; ++i) {
                    mx(best_indices[i], col) = 1;
                    --row_values[best_indices[i]];
                    if (col >= num_check_symbols) {
                        ++information_row_weights[best_indices[i]];
                    }
                }
                // 5.2.5. Updating more statistics, and decreasing row_values here
                for (int i = 0; i < num_check_symbols; ++i) {
                    if (mx(i, col) == 0) {
                        continue;
                    }
                    for (int j = i + 1; j < num_check_symbols; ++j) {
                        if (mx(j, col) > 0) {
                            ++correlations(i, j);
                        }
                    }
                }

                total_t2 += best_t2;
                total_t3 += best_t3;

                // 5.3. Computing ACE graph and its statistics
                if (col >= num_check_symbols) {
                    graph< int > ace_graph;
                    for (int j = 0; j <= col; ++j) {
                        int column_weight = 0;
                        for (int i = 0; i < num_check_symbols; ++i) {
                            if (mx(i, j) != 0) {
                                ++column_weight;
                            }
                        }
                        for (int i = 0; i < num_check_symbols; ++i) {
                            if (mx(i, j) != 0) {
                                ace_graph.add_bidi_edge(j, col + 1 + i, column_weight, column_weight);
                            }
                        }
                    }
                    ace_evaluator ace_ev(ace_graph);
                    int score0 = ace_ev.num_cycles;
                    int score1 = ace_ev.min_ace - 2 * ace_ev.girth;

                    if (col == num_check_symbols && ace_ev.girth < 4) {
                        if (debug) {
                            printf("[DEBUG] [base_matrix] girth at col = %d is %d\n", col, ace_ev.girth);
                        }
                        break; // continues with the next matrix
                    }

                    if (col + 1 == code_max_length && (!rv.solution_exists || score0 < rv.score[1])) {
                        printf("t=%d, ace=%d, num=%d, T2=%d, T3=%d \n", trial_base, score1, score0, total_t2, total_t3);
                        rv.score[0] = score1;
                        rv.score[1] = score0;
                        rv.score[2] = total_t2;
                        rv.score[3] = total_t3;
                        rv.solution_exists = true;
                        rv.data = mx;
                    }
                }
            }
        } catch (base_matrix_interruption_exception &ex) {
            break;
        }
    }

    if (rv.solution_exists) {
        // 6. Post-processing. As checking part is already done by generators,
        //                     we just put 2 to first rows in every information column.
        for (int i = num_check_symbols; i < code_max_length; ++i) {
            int j = 0;
            while (rv.data(j, i) == 0) {
                ++j;
            }
            rv.data(j, i) = 2;
        }
    }

    return rv;
}
