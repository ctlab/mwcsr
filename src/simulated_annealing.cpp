#include "include/simulated_annealing.h"

namespace {
    double probability(double ds, double temp) {
        return std::min(1.0, exp(ds / temp));
    }
}

namespace annealing {

    SimulatedAnnealing::SimulatedAnnealing(const Graph& graph, RandomEngine& random_engine)
                                          :random_engine(random_engine),
                                           graph(graph),
                                           unif{random_engine},
                                           sub(graph) {}

    void SimulatedAnnealing::run(CoolingSchedule& schedule) {
        while (schedule.is_hot()) {
            temperature = schedule.next();
            strike();
        }
    }

    void SimulatedAnnealing::strike() {
        if (sub.size() == 0) {
            empty_module_step();
        } else {
            edge_step();
        }
        if (sub.score() > best_score) {
            best_score = sub.score();
            best = sub.get_snapshot();
        }
    }

    void SimulatedAnnealing::empty_module_step() {
        size_t v = uniform(graph.size());
        if (accepts(graph.weight(v))) {
            sub.add_vertex(v);
        }
    }

    bool SimulatedAnnealing::accepts(double diff) {
        double prob = probability(diff, temperature);
        return unif() < prob;
    }

    void SimulatedAnnealing::edge_step() {
        size_t bdr_sz = sub.boundary_size();
        size_t mdl_sz = sub.number_of_edges();

        if (sub.size() == 1) {
            ++mdl_sz;
        }

        size_t r = uniform(bdr_sz + mdl_sz);
        if (r < bdr_sz) {
            add_from_bdr();
        } else {
            remove_from_module();
        }
    }

    size_t SimulatedAnnealing::uniform(size_t n) {
        std::uniform_int_distribution<size_t> sampler(0, n - 1);
        return sampler(random_engine);
    }

    void SimulatedAnnealing::add_from_bdr() {
        size_t e = sub.random_adjacent_edge(random_engine);
        const Edge& edge = sub.edge(e);
        size_t v = edge.from();
        size_t u = edge.to();
        double diff = 0.0;
        if (!sub.contains_vertex(v)) {
            diff += graph.weight(v);
        }
        if (!sub.contains_vertex(u)) {
            diff += graph.weight(u);
        }
        diff += edge.weight();
        if (accepts(diff)) {
            sub.add_edge(e);
        }
    }

    void SimulatedAnnealing::remove_from_module() {
        if (sub.size() == 1) {
            sub.remove_vertex(sub.any_vertex());
            return;
        }

        size_t e = sub.random_inner_edge(random_engine);
        const Edge& edge = sub.edge(e);
        double diff = -edge.weight();
        size_t v = edge.from();
        size_t u = edge.to();
        if (sub.degree(v) == 1 && sub.degree(u) == 1) {
            if (unif() > 0.5) {
                std::swap(v, u);
            }
            diff -= graph.weight(u);
            if (accepts(diff)) {
                sub.remove_edge(e);
            }
            return;
        }

        if (sub.degree(v) == 1) {
            std::swap(v, u);
        }

        if (sub.degree(u) == 1) {
            diff -= graph.weight(u);
            if (accepts(diff)) {
                sub.remove_edge(e);
            }
            return;
        }
        sub.remove_edge(e);
    }

    vector<size_t> SimulatedAnnealing::vertices() {
        return best.vertices();
    }

    vector<Edge> SimulatedAnnealing::edges() {
        return best.edges();
    }

    StandardUniformDistribution::StandardUniformDistribution(RandomEngine &re) :re(re) {}

    double StandardUniformDistribution::operator()(){
        return unif(re);
    }
}
