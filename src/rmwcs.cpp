#include <Rcpp.h>
#include <instance/Instance.h>
#include <solverLag/SolverClassic.h>

// [[Rcpp::export]]
Rcpp::IntegerVector rmwcs_solve(Rcpp::List& graph, Rcpp::List& params){
    Instance instance(params, graph);
    SolverClassic solver(instance, instance.params.maxIter);
    solver.solve();

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
