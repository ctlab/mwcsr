#ifndef ANNEALING_GRAPH_H
#define ANNEALING_GRAPH_H

#include <cstddef>
#include <vector>

namespace annealing {
    class Edge {
        size_t v;
        size_t u;
        size_t id;
        double w;
    public:
        Edge(size_t from, size_t to, double weight, size_t num);
        bool operator==(const Edge& e);
        size_t opposite(size_t v);
        double weight();
    };

    class Graph {
        std::vector<double> vertex_weights;
        std::vector<std::vector<Edge>> adj;
        size_t m;
    public:
        explicit Graph(size_t n);
        void add_edge(size_t v, size_t u, double weight);
        const std::vector<Edge>& neighbours(size_t v);
    };

}

#endif //ANNEALING_GRAPH_H
