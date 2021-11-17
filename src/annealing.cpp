#include <Rcpp.h>
#include <random>

#include "annealing/include/simulated_annealing.h"
#include "annealing/include/fast_schedule.h"
#include "annealing/include/boltzmann_schedule.h"
#include "include/monitor.h"
#include "include/utils.h"

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
    if (instance.containsElementNamed("warm_start_vertices")) {
        auto vertices = as<IntegerVector>(instance["warm_start_vertices"]);
        for (int v: vertices) {
            sa.add_vertex(v - 1);
        }
    }

    if (instance.containsElementNamed("warm_start_edges")) {
        auto edge_ids = as<IntegerVector>(instance["warm_start_edges"]);
        for (int eid: edge_ids) {
            sa.add_edge(eid - 1);
        }
    }

    if (instance.containsElementNamed("warm_start_weight")) {
        if (std::abs(sa.score() - as<NumericVector>(instance["warm_start_weight"])[0]) > 1e-6) {
            Rcpp::stop("The warm start solution has different weight. Must not happen.");
        }
    }

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
    // TODO check if that's ok
    std::transform(es.begin(), es.end(), edges.begin(), [](Edge e){return e.num() + 1;});
    List ret;
    ret["vertices"] = IntegerVector(vs.begin(), vs.end()) + 1;
    ret["edges"] = IntegerVector(edges.begin(), edges.end());
    return ret;
}
