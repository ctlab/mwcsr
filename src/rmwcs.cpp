#include <Rcpp.h>
#include <memory>
#include "Instance.h"
#include "SolverClassic.h"
#include "SolverCardinality.h"
#include "SolverBudget.h"

// [[Rcpp::export]]
Rcpp::List rmwcs_solve(Rcpp::List& network, Rcpp::List& params) {
    Instance instance(network);
    Parameters parameters(params);

    std::unique_ptr<SolverLag> solver;
    if (network.containsElementNamed("budget")) {
        solver.reset(new SolverBudget(instance, parameters));
    } else if (network.containsElementNamed("cardinality")) {
        solver.reset(new SolverCardinality(instance, parameters));
    } else {
        solver.reset(new SolverClassic(instance, parameters));
    }

    if (instance.nNodes > 0) {
        solver->solve();
    }

    int n = instance.nTrueNodes;
    auto solution = instance.incumbent;
    std::vector<unsigned> vertices;

    for (unsigned i = 0; i < n; i++) {
        if (solution[i]) {
            vertices.push_back(i + 1);
        }
    }

    Rcpp::List ret;
    ret["graph"] = Rcpp::IntegerVector(vertices.begin(), vertices.end());
    ret["ub"] = Rcpp::NumericVector::create(instance.bestBoundLag);
    return ret;
}
