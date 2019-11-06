#include <stdexcept>

#include "include/graph.h"

namespace mwcsr {

    using std::vector;

    Edge::Edge(size_t from, size_t to, size_t signal, size_t num) : v(from), u(to), id(num), signal(signal) {}

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

    size_t Edge::edge_signal() const {
        return signal;
    }

    Graph::Graph(size_t n, const std::vector<double>& signal_weights) :m(0), signal_weights(signal_weights) {
        vertex_signals.resize(n);
        adj.resize(n, vector<Edge>());
    }

    void Graph::add_edge(size_t v, size_t u, size_t signal) {
        Edge e(v, u, signal, m++);
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

    Graph::Graph() :m(0) {}

    void Graph::set_signal(size_t v, size_t signal) {
        vertex_signals.at(v) = signal;
    }

    size_t Graph::vertex_signal(size_t v) const {
        return vertex_signals.at(v);

    }

    size_t Graph::num_signals() const {
        return signal_weights.size();
    }

    double Graph::signal_weight(size_t signal) const {
        return signal_weights.at(signal);
    }

}
