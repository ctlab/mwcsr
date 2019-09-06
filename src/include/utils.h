#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include <Rcpp.h>
#include <iostream>

#include "graph.h"

using namespace Rcpp;

namespace mwcsr {

    Graph read_graph(List& instance) {
        auto edges = as<NumericMatrix>(instance["edgelist"]);
        auto signal_weights_r = as<NumericVector>(instance["signal_weights"]);
        auto vertex_signals = as<IntegerVector>(instance["vertex_signals"]);
        auto edge_signals = as<IntegerVector>(instance["edge_signals"]);
        std::vector<double> signal_weights(signal_weights_r.begin(), signal_weights_r.end());
        size_t n = as<IntegerVector>(instance["size"])[0];
        Graph g(n, signal_weights);

        for (size_t i = 0; i < edges.nrow(); i++) {
            size_t v = edges(i, 0) - 1;
            size_t u = edges(i, 1) - 1;
            size_t signal = edge_signals[i];
            g.add_edge(v, u, signal);
        }

        for (size_t i = 0; i < vertex_signals.size(); i++) {
            g.set_signal(i, vertex_signals[i]);
        }
        return g;
    }
}

#endif //SRC_UTILS_H
