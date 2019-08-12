#include <stdexcept>
#include <cmath>

#include "include/graph.h"

namespace mwcsr {

    using std::vector;

    Edge::Edge(size_t from, size_t to, double weight, size_t num) : v(from), u(to), w(weight), id(num) {}

    double Edge::weight() const {
        return w;
    }

    size_t Edge::opposite(size_t w) const {
        if (w == v) {
            return u;
        } else if (w == u) {
            return v;
        } else {
            throw std::invalid_argument("both ends of the edge do not match to the argument");
        }
    }

    bool Edge::operator==(const Edge& e) const {
        return e.id == id;
    }

    size_t Edge::num() const {
        return id;
    }

    size_t Edge::from() const {
        return v;
    }

    size_t Edge::to() const {
        return u;
    }

    Graph::Graph(size_t n) :m(0) {
        vertex_weights.resize(n);
        adj.resize(n, vector<Edge>());
    }

    void Graph::add_edge(size_t v, size_t u, double weight) {
        Edge e(v, u, weight, m++);
        adj[v].push_back(e);
        adj[u].push_back(e);
        edges.push_back(e);
    }

    const vector<Edge>& Graph::neighbours(size_t v) const {
        return adj[v];
    }

    size_t Graph::size() const {
        return adj.size();
    }

    size_t Graph::edgeset_size() const {
        return m;
    }

    const Edge& Graph::edge(size_t e) const {
        return edges[e];
    }

    double Graph::weight(size_t v) const {
        return vertex_weights[v];
    }

    Graph::Graph() {}

    void Graph::set_weight(size_t v, double w) {
        vertex_weights.at(v) = w;
    }

}
