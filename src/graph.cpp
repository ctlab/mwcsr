#include <stdexcept>
#include <algorithm>

#include "include/graph.h"

namespace mwcsr {

using std::vector;

Edge::Edge(size_t from, size_t to, const std::vector<size_t>& signals, size_t num) :ptr(std::make_shared<EdgeRep>(signals)),
                                                                            v(from), u(to), id(num) {}

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

std::vector<size_t> Edge::edge_signals() const {
    return ptr->signals;
}

void Edge::remove() {
    ptr->to_delete = true;
}

bool Edge::removed() const {
    return ptr->to_delete;
}

Edge::Edge() {}

void Edge::set_signals(const std::vector<size_t>& signals) {
    ptr->signals = signals;
}

Graph::Graph(size_t n, const std::vector<double>& signal_weights) :signal_weights(signal_weights), m(0) {
    v_signals.resize(n);
    adj.resize(n, vector<Edge>());
}

void Graph::add_edge(size_t v, size_t u, std::vector<size_t> signal) {
    size_t id = m++;
    Edge e(v, u, {signal}, id);
    adj[v].push_back(e);
    adj[u].push_back(e);
    edges.push_back(e);
}

vector<Edge> Graph::neighbours(size_t v) {
    std::vector<Edge> ret;
    std::copy_if(adj.at(v).begin(), adj.at(v).end(), std::back_inserter(ret), [](const Edge& e) {
        return !e.removed();
    });
    if (ret.size() != adj.at(v).size()) {
        adj.at(v) = ret;
    }
    return ret;
}

size_t Graph::size() const {
    return adj.size();
}

size_t Graph::edgeset_size() const {
    return m;
}

Edge& Graph::edge(size_t e) {
    return edges[e];
}

const Edge& Graph::const_edge(size_t e) const {
    return edges[e];
}

Graph::Graph() :m(0) {}

void Graph::set_signals(size_t v, std::vector<size_t> signal) {
    v_signals.at(v) = signal;
}

std::vector<size_t> Graph::vertex_signals(size_t v) const {
    return v_signals.at(v);

}

size_t Graph::num_signals() const {
    return signal_weights.size();
}

double Graph::signal_weight(size_t signal) const {
    return signal_weights.at(signal);
}

void Graph::remove_edge(size_t edge_num) {
    auto e = edge(edge_num);
    e.remove();
}

void Graph::remove_vertex(size_t v) {
    for (Edge e: adj[v]) {
        e.remove();
    }
    adj[v].clear();
}

void Graph::absorb_vertex_signals() {
    for (Edge e: edges) {
        auto ss = e.edge_signals();
        auto vs = vertex_signals(e.from());
        auto us = vertex_signals(e.to());
        ss.insert(ss.end(), vs.begin(), vs.end());
        ss.insert(ss.end(), us.begin(), us.end());
        std::sort(ss.begin(), ss.end());
        ss.erase(std::unique(ss.begin(), ss.end()), ss.end());
        e.set_signals(ss);
    }
    for (size_t v = 0; v < size(); v++) {
        v_signals[v].clear();
    }
}

}
