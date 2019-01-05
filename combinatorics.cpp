#include "commons_portable.h"
#include "combinatorics.h"
#include <algorithm>
#include <time.h>

#undef min
#undef max

using std::sort;
using std::vector;
using std::pair;
using std::make_pair;
using std::min;
using std::max;

/*************** divisor_table ***************/

global_divisor_table::global_divisor_table() {
    divisors.push_back(vector< int >(1, 1));
}

vector< int > const &global_divisor_table::operator () (int number) {
    if (number <= 0) {
        die("global_divisor_table is called with the argument %d", number);
    }
    while ((int) divisors.size() <= number) {
        int current = (int) divisors.size() + 1;
        divisors.push_back(vector< int >());
        vector< int > &divs = divisors.back();
        for (int div = 1; div * div <= current; ++div) {
            if (current % div == 0) {
                divs.push_back(div);
                int other = current / div;
                if (other != div) {
                    divs.push_back(other);
                }
            }
        }
        sort(divs.begin(), divs.end());
    }
    return divisors[number - 1];
}

global_divisor_table divisors;

/*************** binomial_table ***************/

global_binomial_table::global_binomial_table() {
    binoms.push_back(vector< long long >(1, 1));
}

long long evaluate_cached(vector< vector< long long > > &binoms, int n, int k) {
    if (binoms[n][k] == 0) {
        if (k == 0 || k == n) {
            binoms[n][k] = 1;
        } else {
            binoms[n][k] = evaluate_cached(binoms, n - 1, k) + evaluate_cached(binoms, n - 1, k - 1);
            if (binoms[n][k] < 0) {
                // This check is enough due to we are adding not multiplying
                die("binom(%d, %d) caused integer overflows", n, k);
            }
        }
    }
    return binoms[n][k];
}

long long global_binomial_table::operator () (int n, int k) {
    if (n < 0 || k < 0 || k > n) {
        die("binomial_table is called with the arguments n = %d, k = %d", n, k);
    }
    while ((int) binoms.size() <= n) {
        binoms.push_back(vector< long long >(binoms.size() + 1));
    }
    return evaluate_cached(binoms, n, k);
}

global_binomial_table binom;

/*************** gcd ***************/

int gcd(int a, int b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        int t = a % b;
        a = b;
        b = t;
    }
    return a > 0 ? a : -a;
}

/*************** addend_splitter ***************/

addend_splitter::addend_splitter(int n, int k): n(n), k(std::min(k, n)) {
    if (n <= 0 || k <= 0) {
        die("Illegal argument to addend_splitter(n = %d, k = %d)", n, k);
    }
}

long long addend_splitter::n_splits() const {
    return binom(n + k - 1, k - 1);
}

void addend_splitter::split(long long t, vector< int > &target) const {
    long long my_n_splits = n_splits();
    if (t < 0 || t >= my_n_splits) {
        die("Illegal number of split %d: must be >= 0 and < %d for n = %d, k = %d", t, my_n_splits, n, k);
    }
    target.clear();
    int nn = n;
    for (int kk = k; kk > 1; --kk) {
        int curr = 0;
        long long bc;
        while (t >= (bc = binom(nn - curr + kk - 2, kk - 2))) {
            t -= bc;
            ++curr;
        }
        target.push_back(curr);
        nn -= curr;
    }
    target.push_back(nn);
}

/*************** random_binom_realization ***************/

void random_binom_realization(int n, int w, vector< int > &target) {
    if (w > n) {
        die("random_binom_realization is called with (n = %d, w = %d)", n, w);
    }
    target.clear();
    for (int i = 0; i < n; ++i) {
        if (n - i - 1 < w || next_random_int(0, n - i) < w) {
            target.push_back(i);
            --w;
        }
    }
    if (w != 0) {
        die("random_binom_realization works wrong, w = %d != 0", w);
    }
}

/*************** random_knapsack_solution ***************/

random_knapsack_solution::random_knapsack_solution(
        vector< pair< int, int > > const &weight_constraints,
        int total_weight,
        int num_elements
) : construction_time((double) (clock()) / CLOCKS_PER_SEC),
    weight_constraints(weight_constraints),
    total_weight(total_weight),
    num_elements(num_elements),
    dp(
        weight_constraints.size() + 1,
        vector< vector< double > >(
            num_elements + 1,
            vector< double >(total_weight + 1)
        )
    )
{
    dp[0][0][0] = 1.0;
    for (unsigned w = 0, w_max = (unsigned int)weight_constraints.size(); w < w_max; ++w) {
        vector< vector< double > > const &prev = dp[w];
        vector< vector< double > >       &next = dp[w + 1];
        for (int delta = weight_constraints[w].first; delta <= weight_constraints[w].second; ++delta) {
            int dw = (int) (w) * delta;
            for (int elem = 0; elem + delta <= num_elements; ++elem) {
                vector< double > const &prev_inner = prev[elem];
                vector< double >       &next_inner = next[elem + delta];
                for (int wgt = 0; wgt + dw <= total_weight; ++wgt) {
                    next_inner[wgt + dw] += prev_inner[wgt];
                }
            }
        }
    }
    construction_time = (double) (clock()) / CLOCKS_PER_SEC - construction_time;
}

bool random_knapsack_solution::has_solution() const {
    return count_solutions() > 0;
}

double random_knapsack_solution::count_solutions() const {
    return dp[weight_constraints.size()][num_elements][total_weight];
}

bool random_knapsack_solution::sample_solution(vector< int > &weights) const {
    if (has_solution()) {
        if ((int) weights.size() != num_elements) {
            die("Requested vector has size %d, but needs to be %d", (int) weights.size(), num_elements);
        }
        unsigned weight_count = 0;
        int elements = num_elements;
        int curr_weight = total_weight;
        for (int w = (int) (weight_constraints.size()) - 1; w != -1; --w) {
            vector< double > counts;
            double count_sum = 0;
            for (int delta = weight_constraints[w].first; delta <= weight_constraints[w].second; ++delta) {
                if (elements - delta >= 0 && curr_weight - w * delta >= 0) {
                    counts.push_back(dp[w][elements - delta][curr_weight - w * delta]);
                    count_sum += counts.back();
                } else {
                    counts.push_back(0);
                }
            }
            if (count_sum != dp[w + 1][elements][curr_weight]) {
                die("DP is broken: backref sums = %lf, original DP value = %lf", count_sum, dp[w + 1][elements][curr_weight]);
            }
            double random_sampled = count_sum * next_random_01();
            int index = 0;
            while (index < (int) (counts.size()) && random_sampled - counts[index] > 0) {
                random_sampled -= counts[index];
                ++index;
            }
            index = min(index, (int) (counts.size()) - 1);
            for (int t = 0; t < index + weight_constraints[w].first; ++t) {
                weights[weights.size() - ++weight_count] = w;
            }
            elements -= index + weight_constraints[w].first;
            curr_weight -= w * (index + weight_constraints[w].first);
        }
        if (weight_count != weights.size() || elements != 0 || curr_weight != 0) {
            die("DP is broken: sample solution restoration resulted in weight_count = %u%u, elements = %d, curr_weight = %d",
                weight_count, weights.size(), elements, curr_weight);
        }
        return true;
    } else {
        return false;
    }
}
