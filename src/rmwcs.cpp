#include <Rcpp.h>
#include <memory>
#include <instance/Instance.h>
#include <solverLag/SolverClassic.h>
#include <solverLag/SolverCardinality.h>
#include <solverLag/SolverBudget.h>

// [[Rcpp::export]]
Rcpp::IntegerVector rmwcs_solve(Rcpp::List& data, Rcpp::List& params) {
    Instance instance(params, data);

    std::unique_ptr<SolverLag> solver;
    if (data.containsElementNamed("budget")){
        solver.reset(new SolverBudget(instance, instance.params.maxIter));
    } else if (data.containsElementNamed("card")) {
        solver.reset(new SolverCardinality(instance, instance.params.maxIter));
    } else {
        solver.reset(new SolverClassic(instance, instance.params.maxIter));
    }

    if (instance.nNodes > 0) {
        solver->solve();
    }

    int n = instance.nTrueNodes;
    auto solution = instance.incumbent;
    std::vector<unsigned> vertices;

    for(unsigned i = 0; i < n; i++){
        if(solution[i]){
            vertices.push_back(i + 1);
        }
    }
    return Rcpp::IntegerVector(vertices.begin(), vertices.end());
}
