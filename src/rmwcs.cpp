#include <Rcpp.h>
#include <instance/Instance.h>
#include <solverLag/SolverClassic.h>

//' [[Rcpp::export]]
double rmwcs_solve(Rcpp::List graph, Rcpp::List params){
    Instance instance(params, graph);
    SolverClassic solver(instance, instance.params.maxIter);
    solver.solve();
    return instance.incumbentObjLag;
}