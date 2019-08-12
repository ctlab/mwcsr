#ifndef SRC_SUBGRAPH_H
#define SRC_SUBGRAPH_H

#include <vector>

#include "DynamicGraph.h"
#include "index.h"
#include "definitions.h"
#include "graph.h"
#include "module.h"

namespace annealing {

    class Subgraph {
        dgraph::DynamicGraph dynamic_graph;
        const mwcsr::Graph graph;
        Index module_edges;
        Index boundary;
        Index module_vertices;
        vector<size_t> vdegree;

        double weight = 0;
        size_t n_vertices = 0;
        std::vector<dgraph::EdgeToken> tokens;
    public:
        explicit Subgraph(const mwcsr::Graph& supergraph);
        void add_vertex(size_t v);
        void add_edge(size_t e);
        bool remove_edge(size_t e);
        void remove_vertex(size_t v);

        size_t random_adjacent_edge(RandomEngine& re) const;
        size_t random_inner_edge(RandomEngine& re) const;
        size_t any_vertex() const;
        size_t size() const;
        size_t number_of_edges() const;
        double score() const;
        size_t boundary_size() const;
        bool contains_vertex(size_t v) const;
        const mwcsr::Edge& edge(size_t e) const;
        size_t degree(size_t v) const;
        Module get_snapshot() const;
    };
}

#endif //SRC_SUBGRAPH_H
