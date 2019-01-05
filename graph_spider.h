#ifndef __GRAPH_SPIDER_H
#define __GRAPH_SPIDER_H

#include <vector>

#include "commons_portable.h"
#include "data_structures.h"

template<typename edge_info>
struct graph_spider_node {
    int parent;
    edge< edge_info > const *incoming_edge;
    edge_info path_info;

    graph_spider_node(int parent, edge< edge_info > const *incoming_edge, edge_info const &path_info)
        : parent(parent), incoming_edge(incoming_edge), path_info(path_info) {}
};

template<typename edge_info>
struct graph_spider {
    typedef graph_spider_node< edge_info > node_type;

protected:
    typedef typename graph< edge_info >::edge_type edge_type;

    graph< edge_info > g;
    std::vector< node_type > nodes;
    edge_info edge_info_zero;
    std::vector< std::vector< int > > vertex2prev, vertex2curr;

    // Notifies that a new cycle is detected.
    virtual void new_cycle_hook(int smaller_end, int larger_end, int cycle_length) = 0;
    // Asks if it is needed to process cycles of the given length.
    virtual bool makes_sense_processing(int cycle_length) = 0;

private:
    edge_type non_existing_edge;
    std::vector< edge_type > virtual_incoming_edges;

    void push_new_node(int parent, edge_type const *incoming_edge, edge_info const &path_info, int current_layer) {
        nodes.push_back(node_type(parent, incoming_edge, path_info));
        node_type const &new_node = nodes.back();
        int new_node_index = (int) nodes.size() - 1;
        if (new_node.parent != 0) {
            // Not a root node, so may check if it forms an equation.
            edge_type const &my_incoming = *new_node.incoming_edge;

            for (int selector = 0; selector <= 1; ++selector) {
                int cycle_length = 2 * current_layer - 1 + selector;
                // Checking all nodes with the same endpoint, which are in the (selector == 0 ? previous : current) layer.
                std::vector< int > const &other_same_end = (selector == 0 ? vertex2prev : vertex2curr)[my_incoming.target];
                for (int i = 0, i_max = (int) other_same_end.size(); i < i_max && makes_sense_processing(cycle_length); ++i) {
                    // If same edge, continue
                    if (new_node.incoming_edge->pair_id == nodes[other_same_end[i]].incoming_edge->pair_id) {
                        continue;
                    }
                    // Otherwise, running a hook
                    new_cycle_hook(other_same_end[i], new_node_index, cycle_length);
                }
            }

            // Inserting the new node into the current layer.
            if (makes_sense_processing(2 * current_layer)) {
                vertex2curr[my_incoming.target].push_back(new_node_index);
            }
        } else {
            // A non-root node is even useless to be insered in same_end.
        }
    }

    void expand_from_root(int root) {
        // Creating a root and adding it to nodes.
        int current_node = (int) nodes.size();
        // Cleaning vertex2prev and vertex2curr.
        for (int v = 0, v_max = g.n_vertices(); v < v_max; ++v) {
            vertex2prev[v].clear();
            vertex2curr[v].clear();
        }
        if (!makes_sense_processing(2)) {
            return;
        }

        push_new_node(0, &virtual_incoming_edges[root], edge_info_zero, 0);

        // For every layer (a set of node with the same distance to root)...
        for (int layer_number = 0; current_node < (int) nodes.size(); ++layer_number) {
            if (!makes_sense_processing(layer_number * 2 + 1)) {
                return;
            }
            // Moving vertex2curr into vertex2prev
            vertex2prev.swap(vertex2curr);
            for (int v = 0; v < g.n_vertices(); ++v) {
                vertex2curr[v].clear();
            }
            // Processing every node of the current layer
            int layer_end = (int) nodes.size();
            while (current_node < layer_end) {
                edge_type const &current_incoming = *nodes[current_node].incoming_edge;
                std::vector< edge_type > const &outgoing = g.outgoing_edges(current_incoming.target);
                // Checking all graph's outgoing edges for the corresponding vertex
                for (int ei = 0, ei_max = (int) outgoing.size(); ei < ei_max; ++ei) {
                    edge_type const &current_outgoing = outgoing[ei];
                    // If we are about to go backwards, don't do it!
                    // If we are trying to go below root, don't do it!
                    if (current_outgoing.pair_id == current_incoming.pair_id ||
                        current_outgoing.target < root) {
                        continue;
                    }
                    edge_info sum_infos = nodes[current_node].path_info + current_outgoing.info;
                    push_new_node(current_node, &current_outgoing, sum_infos, layer_number + 1);
                }
                ++current_node;
            }
        }
    }

public:
    graph_spider(graph< edge_info > const &g)
        : g(g)
        , edge_info_zero(edge_info())
        , vertex2prev(g.n_vertices())
        , vertex2curr(g.n_vertices())
        , non_existing_edge(g.n_vertices(), g.n_vertices(), edge_info(), -1)
    {
        for (int v = 0; v < g.n_vertices(); ++v) {
            virtual_incoming_edges.push_back(edge_type(g.n_vertices(), v, edge_info_zero, -1));
        }
    }

    void run_spider() {
        // Node 0 is the virtual node which is a parent for every root.
        nodes.push_back(node_type(0, &non_existing_edge, edge_info_zero));
        for (int root = 0; root < g.n_vertices(); ++root) {
            expand_from_root(root);
        }
    }

    virtual ~graph_spider() {}
};

#endif
