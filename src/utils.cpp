//
// Created by Alexander Loboda on 10/1/21.
//
#include "include/graph.h"
#include <Rcpp.h>

namespace {
    using Rcpp::NumericMatrix;
    using Rcpp::NumericVector;
    using Rcpp::List;
    using Rcpp::IntegerVector;
    using Rcpp::as;
}

namespace mwcsr {

mwcsr::Graph read_graph(List& instance) {
    auto edges = as<NumericMatrix>(instance["edgelist"]);
    auto signal_weights_r = as<NumericVector>(instance["signal_weights"]);
    auto vertex_signals = as<IntegerVector>(instance["vertex_signals"]);
    auto edge_signals = as<IntegerVector>(instance["edge_signals"]);
    std::vector<double> signal_weights(signal_weights_r.begin(), signal_weights_r.end());
    size_t n = as<IntegerVector>(instance["size"])[0];
    mwcsr::Graph g(n, signal_weights);

    for (size_t i = 0; i < (size_t)edges.nrow(); i++) {
        size_t v = edges(i, 0) - 1;
        size_t u = edges(i, 1) - 1;
        size_t signal = edge_signals[i];
        g.add_edge(v, u, {(size_t)signal});
    }

    for (size_t i = 0; i < (size_t)vertex_signals.size(); i++) {
        g.set_signals(i, {(size_t) vertex_signals[i]});
    }
    return g;
}
}
