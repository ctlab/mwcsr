//
// Created by Alexander Loboda on 8/24/21.
//

#include <Rcpp.h>
#include <set>

#include "include/graph.h"
#include "include/utils.h"
#include "include/monitor.h"
#include "sgmwcs/include/relax_and_cut.h"

namespace {
    relax::Parameters from_r(Rcpp::List& params) {
        relax::Parameters ret;
        ret.beta_iterations = params["beta_iterations"];
        ret.heur_period = params["heur_iterations"];
        ret.iterations = params["max_iterations"];
        ret.sep_period = params["sep_iterations"];
        ret.verbose = params["verbose"];
        return ret;
    }

    using Rcpp::List;
    using mwcsr::Graph;
    using mwcsr::read_graph;
}

// [[Rcpp::export]]
List sgmwcs_solve(List& instance, List& solver_params) {
    constexpr int CHECK_USER_INTERRUPT_MILLIS = 100;

    Graph g = read_graph(instance);
    mwcsr::monitor int_monitor(Rcpp::checkUserInterrupt, CHECK_USER_INTERRUPT_MILLIS);

    List ret;
    relax::Parameters params = from_r(solver_params);
    relax::Solver solver(g, params, Rcpp::Rcout);
    solver.set_monitor(int_monitor);

    size_t best_vertex = 0;
    double best_weight = -std::numeric_limits<double>::infinity();

    for (size_t v = 0; v < g.size(); v++) {
        auto signals = g.vertex_signals(v);
        std::set<size_t> uniq(signals.begin(), signals.end());
        double weight = std::accumulate(uniq.begin(), uniq.end(), 0, [&g](double w, size_t i) -> double {
            return w + g.signal_weight(i);
        });

        if (weight > best_weight) {
            best_weight = weight;
            best_vertex = v;
        }
    }

    solver.solve();

    relax::Solution sol = solver.solution();

    if (sol.objective() > best_weight) {
        auto edges = sol.solution();
        ret["edges"] = IntegerVector(edges.begin(), edges.end()) + 1;
    } else {
        IntegerVector res{(int)best_vertex};
        ret["vertices"] = res + 1;
    }

    ret["lb"] = std::max(sol.objective(), best_weight);
    ret["ub"] = solver.upper_bound();

    return ret;
}
