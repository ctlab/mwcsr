#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include <Rcpp.h>

#include "graph.h"

using namespace Rcpp;

namespace mwcsr {

    Graph read_graph(List& instance) {
        auto edges = as<NumericMatrix>(instance["edgelist"]);
        auto edge_weights = as<NumericVector>(instance["edge_weights"]);
        auto vertex_weights = as<NumericVector>(instance["vertex_weights"]);
        size_t n = as<IntegerVector>(instance["size"])[0];
        Graph g(n);

        for (size_t i = 0; i < edges.nrow(); i++) {
            size_t v = edges(i, 0) - 1;
            size_t u = edges(i, 1) - 1;
            double weight = edge_weights[i];
            g.add_edge(v, u, weight);
        }

        for (size_t i = 0; i < vertex_weights.size(); i++) {
            g.set_weight(i, vertex_weights[i]);
        }
        return g;
    }
}

#endif //SRC_UTILS_H
