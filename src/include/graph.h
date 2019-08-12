#ifndef ANNEALING_GRAPH_H
#define ANNEALING_GRAPH_H

#include <cstddef>
#include <vector>

namespace mwcsr {

    class Edge {
        size_t v;
        size_t u;
        size_t id;
        double w;
    public:
        Edge(size_t from, size_t to, double weight, size_t num);

        bool operator==(const Edge& e) const;
        size_t opposite(size_t v) const;
        double weight() const;
        size_t num() const;
        size_t from() const;
        size_t to() const;
    };

    class Graph {
        std::vector<double> vertex_weights;
        std::vector<std::vector<Edge>> adj;
        std::vector<Edge> edges;
        size_t m;

    public:
        Graph();
        Graph(const Graph& other) = default;
        explicit Graph(size_t n);
        void add_edge(size_t v, size_t u, double weight);
        void set_weight(size_t v, double w);

        const std::vector<Edge>& neighbours(size_t v) const;
        size_t size() const;
        size_t edgeset_size() const;
        const Edge& edge(size_t e) const;
        double weight(size_t v) const;
    };

}

#endif //ANNEALING_GRAPH_H
