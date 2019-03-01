#include <algorithm>
#include <deque>
#include <map>
#include <set>
#include <utility>

#include "commons_portable.h"
#include "equations.h"
#include "data_structures.h"
#include "combinatorics.h"

using std::pair;
using std::set;
using std::vector;
using std::map;
using std::make_pair;

typedef map< vector_bag, equation_description >::iterator equations_iterator;

bool voltage_check_result::check_modulo(int modulo) const {
    return first_failed_equation == -1 &&
           modulo >= min_ok_modulo &&
           (modulo > max_nonok_modulo || failed_modulos.count(modulo) == 0);
}

equation_builder::equation_builder(graph< vector_bag > const &g, int max_girth)
    : graph_spider< vector_bag >(g)
    , max_girth(max_girth)
{
    run_spider();
    remove_long_equations();
    build_result();
}

bool equation_builder::makes_sense_processing(int cycle_length) {
    return max_girth > cycle_length;
}

void equation_builder::new_cycle_hook(int smaller_index, int greater_index, int cycle_length) {
    if (max_girth > cycle_length) {
        // Computing path sum
        vector_bag sum = nodes[greater_index].path_info - nodes[smaller_index].path_info;
        if (sum == edge_info_zero) {
            // If path sum is identity zero, update max girth
            max_girth = cycle_length;
        } else {
            // If path sum is not zero, create a new equation
            sum.normalize();
            equation_description new_equation;
            new_equation.first = smaller_index;
            new_equation.second = greater_index;
            new_equation.cycle_length = cycle_length;
            // Trying to insert it
            pair< vector_bag, equation_description > to_insert = make_pair(sum, new_equation);
            pair< equations_iterator, bool > if_inserted = set2eq.insert(to_insert);
            // If not inserted, something already was there, checking...
            if (!if_inserted.second) {
                // If our equation is shorter, we replace the existing one
                if (if_inserted.first->second.cycle_length > cycle_length) {
                    // This will be in the same place, so using hints to accelerate things
                    equations_iterator insertion_hint = if_inserted.first;
                    --insertion_hint;
                    set2eq.erase(if_inserted.first);
                    set2eq.insert(insertion_hint, to_insert);
                }
            }
        }
    }
}

void equation_builder::remove_long_equations() {
    for (equations_iterator itr = set2eq.begin(), itr_end = set2eq.end(); itr != itr_end; ) {
        equations_iterator copy = itr++;
        if (copy->second.cycle_length >= max_girth) {
            set2eq.erase(copy);
        }
    }
}

void equation_builder::build_result() {
    // 1. Mark active tree vertices
    vector< bit > active(nodes.size());
    // 1.1. Virtual root is active
    active[0] = true;
    // 1.2. For every equation go downwards from endpoints
    for (equations_iterator itr = set2eq.begin(), itr_end = set2eq.end(); itr != itr_end; ++itr) {
        equation_description const &eq = itr->second;
        for (int endpoint = 0; endpoint <= 1; ++endpoint) {
            int current = endpoint ? eq.first : eq.second;
            // If the vertex is already active, we have processed all the path towards root.
            while (!active[current]) {
                active[current] = true;
                current = nodes[current].parent;
            }
        }
    }
    // 2. Create remapper.
    vector< int > remapper(nodes.size(), -1);
    int n_saved_vertices = 0;
    for (int i = 0, i_max = (int) nodes.size(); i < i_max; ++i) {
        if (active[i]) {
            remapper[i] = n_saved_vertices++;
        }
    }
    // 3. Moving all equations to the "equations" variable
    for (equations_iterator itr = set2eq.begin(), itr_end = set2eq.end(); itr != itr_end; ++itr) {
        equations.push_back(itr->second);
        equation_description &eq = equations.back();
        eq.first = remapper[eq.first];
        eq.second = remapper[eq.second];
    }
    // 4. Moving the tree nodes to the "tree" variable
    for (int i = 0, i_max = (int) nodes.size(); i < i_max; ++i) {
        if (active[i]) {
            result_tree_node rtn;
            rtn.parent = remapper[nodes[i].parent];
            if (rtn.parent != 0) {
                vector< pair< int, int > > const &bag_contents = nodes[i].incoming_edge->info.data();
                if (bag_contents.size() != 1) {
                    die("I don't understand why edge info has size %d", (int) (bag_contents.size()));
                }
                if (bag_contents[0].second != -1 && bag_contents[0].second != 1) {
                    die("I don't understand why edge info has multiplicity of %d", bag_contents[0].second);
                }
                rtn.edge_id = bag_contents[0].first;
                rtn.negative = bag_contents[0].second == -1;
            } else {
                rtn.edge_id = 0;
                rtn.negative = false;
            }
            tree.push_back(rtn);
        }
    }
}

int equation_builder::n_equations() const {
    return (int) equations.size();
}

int equation_builder::n_vertices() const {
    return (int) tree.size() - 1;
}

int equation_builder::tree_girth() const {
    return max_girth;
}

equation_builder::~equation_builder() {}

voltage_check_result equation_builder::check_voltages
(
    vector< int > const &voltages,
    int max_modulo
) const {
    voltage_check_result rv;
    // 1. Filling up the nodes.
    vector< int > values(tree.size());
    for (int i = 1, i_max = (int) tree.size(); i < i_max; ++i) {
        result_tree_node const &curr = tree[i];
        if (curr.parent == 0) {
            // a root in one of the trees
            values[i] = 0;
        } else if (curr.negative) {
            values[i] = values[curr.parent] - voltages[curr.edge_id];
        } else {
            values[i] = values[curr.parent] + voltages[curr.edge_id];
        }
    }
    // 2. Testing the equations.
    int num_failed = 0;
    for (int i = 0, i_max = (int) equations.size(); i < i_max; ++i) {
        int sum = values[equations[i].first] - values[equations[i].second];
        if (sum == 0) {
            // An identity zero found.
            rv.first_failed_equation = num_failed;
            return rv;
        }
    }
    rv.first_failed_equation = -1;
    // 3. Collecting failing modulos.
    for (int i = 0, i_max = (int) equations.size(); i < i_max; ++i) {
        int sum = std::abs(values[equations[i].first] - values[equations[i].second]);
        vector< int > const &divs = divisors(sum);
        for (int j = 0, j_max = (int) divs.size(); j < j_max && divs[j] <= max_modulo; ++j) {
            rv.failed_modulos.insert(divs[j]);
        }
    }
    rv.max_nonok_modulo = rv.failed_modulos.size() > 0 ? *(rv.failed_modulos.rbegin()) : 1;
    rv.min_ok_modulo = 2;
    while (rv.failed_modulos.count(rv.min_ok_modulo) > 0) {
        ++rv.min_ok_modulo;
    }
    return rv;
}


voltage_check_result equation_builder::check_voltages_and_coefs
(
 vector< int > const &voltages,
 int max_modulo,
 vector< int > const &coef,
 int m_mod
 ) const 
{
	 voltage_check_result rv;
	 // 1. Filling up the nodes.
	 vector< int > values(tree.size());
	 vector< int > coef_values(tree.size());

	 for (int i = 1, i_max = (int) tree.size(); i < i_max; ++i) 
	 {
		 result_tree_node const &curr = tree[i];
		 if (curr.parent == 0) 
		 {
			 // a root in one of the trees
			 values[i] = 0;
			 coef_values[i] = 0;
		 } 
		 else 
			 if (curr.negative) 
			 {
				values[i] = values[curr.parent] - voltages[curr.edge_id];
				coef_values[i] = coef_values[curr.parent] - coef[curr.edge_id];
			 } 
			 else 
			 {
			    values[i] = values[curr.parent] + voltages[curr.edge_id];
				coef_values[i] = coef_values[curr.parent] + coef[curr.edge_id];
			 }
	 }

	 // 2. Testing the equations.
	 int num_failed = 0;
	 for (int i = 0, i_max = (int) equations.size(); i < i_max; ++i) 
	 {
		 int sum = values[equations[i].first] - values[equations[i].second];
		 if (sum == 0) 
		 {
			 // An identity zero found.
			 int coef_sum = coef_values[equations[i].first] - coef_values[equations[i].second];

			 if( coef_sum == 0 )
				rv.first_failed_equation = num_failed;
			 return rv;
		 }
	 }
	 rv.first_failed_equation = -1;
	 // 3. Collecting failing modulos.
	 for (int i = 0, i_max = (int) equations.size(); i < i_max; ++i) 
	 {
		 int sum = std::abs(values[equations[i].first] - values[equations[i].second]);
		 vector< int > const &divs = divisors(sum);
		 for (int j = 0, j_max = (int) divs.size(); j < j_max && divs[j] <= max_modulo; ++j) 
		 {
			 rv.failed_modulos.insert(divs[j]);
		 }
	 }
	 rv.max_nonok_modulo = rv.failed_modulos.size() > 0 ? *(rv.failed_modulos.rbegin()) : 1;
	 rv.min_ok_modulo = 2;
	 while (rv.failed_modulos.count(rv.min_ok_modulo) > 0) 
	 {
		 ++rv.min_ok_modulo;
	 }
	 return rv;
}
