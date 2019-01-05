#ifndef __COMBINATORICS_H
#define __COMBINATORICS_H

#include <vector>
#include <utility>

/*************** gcd ***************/

int gcd(int a, int b);

/*************** global_divisor_table ***************/

// A simple caching implementation of a divisor table.
// Don't create your own instances (unless multithreaded),
// use the single instance "divisors" instead.

class global_divisor_table {
private:
    std::vector< std::vector< int > > divisors;
public:
    global_divisor_table();

    // Returns a vector of divisors of the given "number",
    // sorted in the increasing order.
    std::vector< int > const &operator () (int number);
};

extern global_divisor_table divisors;

/*************** binomial_table ***************/

// A simple caching implementation of a binomial table.
// Don't create your own instances (unless multithreaded),
// use the single instance "binom" instead.

class global_binomial_table {
private:
    std::vector< std::vector< long long > > binoms;
public:
    global_binomial_table();

    // Returns the number of ways to choose "k" items
    // from "n" items, irrespective of their order.
    long long operator() (int n, int k);
};

extern global_binomial_table binom;

/*************** addend_splitter ***************/

class addend_splitter {
private:
    int n, k;
public:
    addend_splitter(int n, int k);

    long long n_splits() const;

    // Fills in the given vector with the "t"th lexicographical
    // split of "n" into "k" addends.
    void split(long long t, std::vector< int > &target) const;
};

/*************** random_binom_realization ***************/

void random_binom_realization(int n, int w, std::vector< int > &target);

/*************** random_knapsack_solution ***************/

class random_knapsack_solution {
    double construction_time;
    std::vector< std::pair< int, int > > weight_constraints;
    int total_weight;
    int num_elements;
    std::vector< std::vector< std::vector< double > > > dp;

public:
    random_knapsack_solution(
        std::vector< std::pair< int, int > > const &weight_constraints,
        int total_weight,
        int num_elements
    );

    bool   has_solution() const;
    double count_solutions() const;
    bool   sample_solution(std::vector< int > &weights) const;
};

#endif
