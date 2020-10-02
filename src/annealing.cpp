#include <Rcpp.h>
#include <random>

#include "include/utils.h"
#include "include/simulated_annealing.h"
#include "include/fast_schedule.h"
#include "include/boltzmann_schedule.h"
#include "include/monitor.h"

using namespace Rcpp;
using namespace mwcsr;
using namespace annealing;

// [[Rcpp::export]]
List sa_solve(List& instance, List& solver) {
    constexpr int CHECK_USER_INTERRUPT_MILLIS = 100;

    Graph g = read_graph(instance);
    double initial = solver["initial_temperature"];
    double final = solver["final_temperature"];
    std::string schedule = solver["schedule"];    

    std::mt19937 re;
    SimulatedAnnealing sa(g, re);

    mwcsr::monitor int_monitor(Rcpp::checkUserInterrupt, CHECK_USER_INTERRUPT_MILLIS);
    if (schedule == "fast") {
        FastSchedule cs(initial, final);
        sa.run(cs, int_monitor);
    } else {
        BoltzmannSchedule cs(initial, final);
        sa.run(cs, int_monitor);
    }
    auto vs = sa.vertices();
    auto es = sa.edges();
    vector<size_t> edges(es.size(), 0);
    std::transform(es.begin(), es.end(), edges.begin(), [](Edge e){return e.num() + 1;});
    List ret;
    ret["vertices"] = IntegerVector(vs.begin(), vs.end());
    ret["edges"] = IntegerVector(edges.begin(), edges.end());
    return ret;
}
