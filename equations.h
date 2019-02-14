#ifndef __EQUATIONS_H
#define __EQUATIONS_H

#include <limits>
#include <map>
#include <set>
#include <vector>

#include "data_structures.h"
#include "graph_spider.h"

struct equation_description {
    int first, second;
    int cycle_length;
};

struct result_tree_node {
    int  parent;
    int  edge_id;
    bool negative;
};

struct voltage_check_result {
    int              first_failed_equation;
    std::set< int >  failed_modulos;
    int              min_ok_modulo;
    int              max_nonok_modulo;

    bool check_modulo(int modulo) const;
};

struct equation_builder : graph_spider< vector_bag > {
private:
    int max_girth;
    std::map< vector_bag, equation_description > set2eq;

    std::vector< result_tree_node >     tree;
    std::vector< equation_description > equations;

    void remove_long_equations();
    void build_result();

protected:
    virtual void new_cycle_hook(int smaller_end, int larger_end, int cycle_length);
    virtual bool makes_sense_processing(int cycle_length);
public:
    equation_builder(graph< vector_bag > const &g, int max_girth = std::numeric_limits< int >::max());
    int n_equations() const;
    int n_vertices() const;
    int tree_girth() const;
    voltage_check_result check_voltages(std::vector< int > const &voltages,
                                        int max_modulo = std::numeric_limits< int >::max()) const;
	voltage_check_result check_voltages_and_coefs(std::vector< int > const &voltages,
		int max_modulo, std::vector< int > const &coef, int qmod ) const;
    virtual ~equation_builder();
};

#endif
