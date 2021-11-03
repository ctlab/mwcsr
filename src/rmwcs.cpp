#include <Rcpp.h>
#include <memory>
#include "rmwcs/include/Instance.h"
#include "rmwcs/include/SolverClassic.h"
#include "rmwcs/include/SolverCardinality.h"
#include "rmwcs/include/SolverBudget.h"

// [[Rcpp::export]]
Rcpp::List rmwcs_solve(Rcpp::List& network, Rcpp::List& params) {
    constexpr int CHECK_USER_INTERRUPT_MILLIS = 100;

    Instance instance(network);
    Parameters parameters(params);

    mwcsr::monitor interruption_monitor(Rcpp::checkUserInterrupt, CHECK_USER_INTERRUPT_MILLIS);

    std::unique_ptr<SolverLag> solver;
    if (network.containsElementNamed("budget")) {
        solver.reset(new SolverBudget(instance, parameters, interruption_monitor));
    } else if (network.containsElementNamed("cardinality")) {
        solver.reset(new SolverCardinality(instance, parameters, interruption_monitor));
    } else {
        solver.reset(new SolverClassic(instance, parameters, interruption_monitor));
    }

    if (instance.nNodes > 0) {
        solver->solve();
    }

    int n = instance.nTrueNodes;
    auto solution = instance.incumbent;
    std::vector<unsigned> vertices;

    if (instance.incumbentFound) {
        for (int i = 0; i < n; i++) {
            if (solution[i]) {
                vertices.push_back(i + 1);
            }
        }
    }

    Rcpp::List ret;
    ret["graph"] = Rcpp::IntegerVector(vertices.begin(), vertices.end());
    ret["lb"] = Rcpp::NumericVector::create(instance.incumbentObjLag);
    ret["ub"] = Rcpp::NumericVector::create(instance.bestBoundLag);
    return ret;
}
