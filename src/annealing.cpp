#include <Rcpp.h>
#include <random>

#include "utils.h"
#include "simulated_annealing.h"
#include "fast_schedule.h"
#include "boltzmann_schedule.h"

using namespace Rcpp;
using namespace mwcsr;
using namespace annealing;

// [[Rcpp::export]]
List sa_solve(List& instance, List& solver) {
    Graph g = read_graph(instance);
    double initial = solver["initial_temperature"];
    double final = solver["final_temperature"];

    std::mt19937 re;
    SimulatedAnnealing sa(g, re);

    if (solver["schedule"] == "1") {
        FastSchedule cs(initial, final);
        sa.run(cs);
    } else {
        BoltzmannSchedule cs(initial, final);
        sa.run(cs);
    }
    auto vs = sa.vertices();
    auto es = sa.edges();
    vector<size_t> edges;
    std::transform(es.begin(), es.end(), edges.begin(), [](Edge e){return e.num();});
    List ret;
    ret["vertices"] = IntegerVector(vs.begin(), vs.end());
    ret["edges"] = IntegerVector(edges.begin(), edges.end());
    return ret;
}