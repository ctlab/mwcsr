#include <stdexcept>

#include "graph.h"

using std::vector;

namespace annealing {

    Edge::Edge(size_t from, size_t to, double weight, size_t num) :v(from), u(to), w(weight), id(num) {}

    double Edge::weight() {
        return w;
    }

    size_t Edge::opposite(size_t w) {
        if (w == v) {
            return u;
        } else if (w == u) {
            return v;
        } else {
            throw std::invalid_argument("both ends of the edge do not match to the argument");
        }
    }

    bool Edge::operator==(const Edge &e) {
        return e.id == id;
    }

    Graph::Graph(size_t n) {
        vertex_weights.resize(n);
        adj.resize(n, vector<Edge>());
    }

    void Graph::add_edge(size_t v, size_t u, double weight) {
        adj[v].emplace_back(v, u, weight, m++);
    }

    const vector<Edge> &Graph::neighbours(size_t v) {
        return adj[v];
    }
}
