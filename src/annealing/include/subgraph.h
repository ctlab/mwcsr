#ifndef SRC_SUBGRAPH_H
#define SRC_SUBGRAPH_H

#include <vector>

#include "DynamicGraph.h"
#include "index.h"
#include "definitions.h"
#include "../../include/graph.h"
#include "module.h"

namespace annealing {

    class Subgraph {
        dgraph::DynamicGraph dynamic_graph;
        mwcsr::Graph graph;
        Index module_edges;
        Index boundary;
        Index module_vertices;
        vector<size_t> vdegree;
        vector<size_t> signal_utilization;

        double weight = 0;
        size_t n_vertices = 0;
        std::vector<dgraph::EdgeToken> tokens;
    public:
        explicit Subgraph(const mwcsr::Graph& supergraph);
        void add_vertex(size_t v);
        void add_edge(size_t e);
        bool remove_edge(size_t e);
        void remove_vertex(size_t v);

        const mwcsr::Edge& edge(size_t e) const;
        size_t degree(size_t v) const;
        size_t size() const;
        size_t number_of_edges() const;
        double score() const;
        size_t boundary_size() const;

        size_t any_vertex() const;
        size_t random_adjacent_edge(RandomEngine& re) const;
        size_t random_inner_edge(RandomEngine& re) const;

        double add_vertex_diff(size_t v) const;
        double remove_vertex_diff(size_t v) const;
        double add_edge_diff(size_t e) const;
        double remove_edge_diff(size_t e) const;

        bool contains_vertex(size_t v) const;
        Module get_snapshot() const;
    private:
        double diff(std::vector<size_t> signal, bool add) const;
        void signals_add(std::vector<size_t> signals);
        void signals_remove(std::vector<size_t> signals);
    };
}

#endif //SRC_SUBGRAPH_H
