#ifndef ANNEALING_GRAPH_H
#define ANNEALING_GRAPH_H

#include <cstddef>
#include <vector>
#include <unordered_map>
#include <memory>

namespace mwcsr {
    class Edge;

    class EdgeRep {
        bool to_delete = false;
        std::vector<size_t> signals;
    public:
        EdgeRep(const std::vector<size_t>& signals) :signals(signals){}
        friend class Edge;
    };

    class Edge {
        std::shared_ptr<EdgeRep> ptr;
        size_t v;
        size_t u;
        size_t id;
    public:
        Edge(size_t from, size_t to, const std::vector<size_t>& signals, size_t num);
        Edge(const Edge&) = default;
        Edge();
        bool operator==(const Edge& e) const;
        size_t opposite(size_t v) const;
        std::vector<size_t> edge_signals() const;
        size_t num() const;
        size_t from() const;
        size_t to() const;
        void remove();
        bool removed() const;
        void set_signals(const std::vector<size_t>& signals);
    };

    class Graph {
        std::vector<double> signal_weights;
        std::vector<std::vector<size_t>> v_signals;
        std::vector<std::vector<Edge>> adj;
        std::vector<Edge> edges;
        size_t m;

    public:
        Graph();
        Graph(const Graph& other) = default;
        Graph(size_t n, const std::vector<double>& signal_weights);
        void add_edge(size_t v, size_t u, std::vector<size_t> signal);
        void set_signals(size_t v, std::vector<size_t> signal);
        void remove_edge(size_t e);
        void remove_vertex(size_t v);

        std::vector<Edge> neighbours(size_t v);
        std::vector<size_t> vertex_signals(size_t v) const;
        size_t size() const;
        size_t edgeset_size() const;
        Edge& edge(size_t e);
        const Edge& const_edge(size_t e) const;
        size_t num_signals() const;
        double signal_weight(size_t signal) const;

        void absorb_vertex_signals();
    };

}

#endif //ANNEALING_GRAPH_H
