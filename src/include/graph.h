#ifndef ANNEALING_GRAPH_H
#define ANNEALING_GRAPH_H

#include <cstddef>
#include <vector>

namespace mwcsr {

    class Edge {
        size_t v;
        size_t u;
        size_t id;
        size_t signal;
    public:
        Edge(size_t from, size_t to, size_t signal, size_t num);

        bool operator==(const Edge& e) const;
        size_t opposite(size_t v) const;
        size_t edge_signal() const;
        size_t num() const;
        size_t from() const;
        size_t to() const;
    };

    class Graph {
        std::vector<double> signal_weights;
        std::vector<size_t> vertex_signals;
        std::vector<std::vector<Edge>> adj;
        std::vector<Edge> edges;
        size_t m;

    public:
        Graph();
        Graph(const Graph& other) = default;
        Graph(size_t n, const std::vector<double>& signal_weights);
        void add_edge(size_t v, size_t u, size_t signal);
        void set_signal(size_t v, size_t signal);

        const std::vector<Edge>& neighbours(size_t v) const;
        size_t vertex_signal(size_t v) const;
        size_t size() const;
        size_t edgeset_size() const;
        const Edge& edge(size_t e) const;
        size_t num_signals() const;
        double signal_weight(size_t signal) const;
    };

}

#endif //ANNEALING_GRAPH_H
