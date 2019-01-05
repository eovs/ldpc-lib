#ifndef __DATA_STRUCTURES_H
#define __DATA_STRUCTURES_H

#include <vector>

#include "commons_portable.h"

/*************** matrix< T > ***************/

template<typename T>
struct matrix {
private:
    std::vector< T > contents;
    int rows, cols;
public:
    matrix() : contents(0, T()), rows(0), cols(0) {}

    matrix(int rows, int cols, T const &init = T())
        : contents(rows * cols, init), rows(rows), cols(cols) {}

    int n_rows() const {
        return rows;
    }
    int n_cols() const {
        return cols;
    }
    T &operator() (int row, int col) {
#ifdef MATRIX_STRICT_CHECKS
        if (row >= rows || col >= cols || row < 0 || col < 0) {
            die("matrix out of bounds (mutable): requested (%d, %d), rows = %d, cols = %d", row, col, rows, cols);
        }
#endif
        return contents[row * cols + col];
    }
    T const &operator() (int row, int col) const {
#ifdef MATRIX_STRICT_CHECKS
        if (row >= rows || col >= cols || row < 0 || col < 0) {
            die("matrix out of bounds (immutable): requested (%d, %d), rows = %d, cols = %d", row, col, rows, cols);
        }
#endif
        return contents[row * cols + col];
    }

    bool operator == (matrix<T> const &that) const {
        if (n_rows() != that.n_rows() || n_cols() != that.n_cols()) {
            return false;
        }
        matrix<T> const &self = *this;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (self(i, j) != that(i, j)) {
                    return false;
                }
            }
        }
        return true;
    }
    bool operator != (matrix<T> const &that) const {
        return !(*this == that);
    }
};

/*************** nothing ***************/

// A "nothing" struct for an empty edge_info
struct nothing {
    inline nothing operator + (nothing const &) const { return *this; }
    inline nothing operator - (nothing const &) const { return *this; }
    inline nothing operator - () const { return *this; }
    inline nothing &operator += (nothing const &) { return *this; }
    inline nothing &operator -= (nothing const &) { return *this; }
    inline bool operator == (nothing const &) { return true; }
    inline bool operator != (nothing const &) { return false; }
};

/*************** edge ***************/

template<typename edge_info>
struct edge {
    int source, target;
    edge_info info;
    int pair_id;

    edge(int source, int target, edge_info const &info, int pair_id)
        : source(source), target(target), info(info), pair_id(pair_id) {}
};

/*************** graph ***************/

template<typename edge_info>
struct graph {
    typedef edge< edge_info > edge_type;
private:
    std::vector< std::vector< edge_type > > edges;
    int num_edges;
public:
    graph() : num_edges(0) {}

    void ensure_has_vertex(int vertex) {
        while ((int) edges.size() <= vertex) {
            edges.push_back(std::vector< edge_type >());
        }
    }

    void add_bidi_edge(int source, int target, edge_info const &fw, edge_info const &bw) {
        ensure_has_vertex(source);
        ensure_has_vertex(target);
        if (source == target) {
            die("graph.add_bidi_edge is called on equal source and target: %d", source);
        }
        int id = num_edges >> 1;
        edges[source].push_back(edge_type(source, target, fw, id));
        edges[target].push_back(edge_type(target, source, bw, id));
        num_edges += 2;
    }

    std::vector< edge_type > const &outgoing_edges(int vertex) const {
        return edges.at(vertex);
    }

    int vertex_outdegree(int vertex) const {
        return vertex < (int) edges.size() ? edges[vertex] : 0;
    }

    int n_vertices() const {
        return (int) edges.size();
    }

    int n_edges() const {
        return num_edges;
    }
};

/*************** vector_bag ***************/

class vector_bag {
    std::vector< std::pair< int, int > > contents;
    int compare_to(vector_bag const &that) const;
public:
    vector_bag();
    vector_bag(int single_value);
    vector_bag(vector_bag const &that);
    vector_bag &operator = (vector_bag const &that);
    vector_bag &operator += (vector_bag const &that);
    vector_bag &operator -= (vector_bag const &that);
    vector_bag operator + (vector_bag const &that) const;
    vector_bag operator - (vector_bag const &that) const;
    vector_bag operator - () const;

    void normalize();

    bool operator == (vector_bag const &that) const;
    bool operator != (vector_bag const &that) const;
    bool operator <  (vector_bag const &that) const;
    bool operator <= (vector_bag const &that) const;
    bool operator >  (vector_bag const &that) const;
    bool operator >= (vector_bag const &that) const;

    std::vector< std::pair< int, int > > const &data() const;
};

#endif
