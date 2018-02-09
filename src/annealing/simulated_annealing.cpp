#include "simulated_annealing.h"

namespace {
    double probability(double e, double es, double temp) {
        return std::min(1.0, exp(es - e / temp));
    }
}

namespace annealing {

    SimulatedAnnealing::SimulatedAnnealing(const Graph &graph, std::function<uint_fast32_t> &random_engine,
                                           CoolingSchedule& cooling_schedule)
                                          :graph(graph),
                                           temp_func(cooling_schedule),
                                           tokens{graph.edgeset_size()},
                                           random_engine(random_engine),
                                           unif{random_engine},
                                           dynamic_graph{graph.size()},
                                           module_vertices{graph.size(), false},
                                           module_edges{graph.edgeset_size(), random_engine},
                                           boundary{graph.edgeset_size(), random_engine} {}

    void SimulatedAnnealing::run() {

    }

    void SimulatedAnnealing::step() {
        temperature = temp_func.temperature();
        if (size == 0) {
            empty_module_step();
        } else {
            edge_step();
        }
    }

    void SimulatedAnnealing::add_vertex(size_t v) {
        size++;
        module_vertices[v] = true;
        for (Edge e: graph.neighbours(v)) {
            size_t u = e.opposite(v);
            size_t e_id = e.num();
            if (!(module_edges.contains(e_id) || boundary.contains(e_id))) {
                boundary.add(e_id);
            }
        }
        score += graph.weight(v);
    }

    void SimulatedAnnealing::add_edge(size_t e) {
        boundary.remove(e);
        module_edges.add(e);
        Edge edge = graph.edge(e);
        size_t v = edge.from();
        size_t u = edge.to();
        if (!module_vertices[v]) {
            add_vertex(v);
        }
        if (!module_vertices[u]) {
            add_vertex(u);
        }
        tokens[e] = std::move(dynamic_graph.add(v, u));
        score += edge.weight();
    }

    bool SimulatedAnnealing::remove_edge(size_t e) {
        Edge edge = graph.edge(e);
        size_t v = edge.from();
        size_t u = edge.to();
        dynamic_graph.remove(std::move(tokens[e]));
        size_t comp_size = dynamic_graph.component_size(v);
        if (comp_size < size - 1 && comp_size != 1) {
            tokens[e] = std::move(dynamic_graph.add(v, u));
            return false;
        }
        if (comp_size == size) {
            boundary.add(e);
            return true;
        }
        if (comp_size == size - 1) {
            remove_vertex(v);
        } else {
            remove_vertex(u);
        }
        score -= edge.weight();
        return true;
    }

    void SimulatedAnnealing::remove_vertex(size_t v) {
        size--;
        for (Edge e: graph.neighbours(v)) {
            if (boundary.contains(e.num())) {
                boundary.remove(e.num());
            }
        }
        module_vertices[v] = false;
        score -= graph.weight(v);
    }

    void SimulatedAnnealing::empty_module_step() {
        size_t v = uniform(graph.size());
        if (accepts(graph.weight(v))) {
            add_vertex(v);
        }
    }

    bool SimulatedAnnealing::accepts(double diff) {
        double prob = probability(score, score + diff, temperature);
        return unif() < prob;
    }

    void SimulatedAnnealing::edge_step() {
        size_t bdr_sz = boundary.size();
        size_t mdl_sz = module_edges.size();
        size_t r = uniform(bdr_sz + mdl_sz);
        if (r < bdr_sz) {
            add_from_bdr();
        } else {
            remove_from_module();
        }
    }

    size_t SimulatedAnnealing::uniform(size_t n) {
        std::uniform_int_distribution uniform_int_distribution(0, n - 1);
        return (size_t)uniform_int_distribution(random_engine);
    }

    void SimulatedAnnealing::add_from_bdr() {
        size_t e = boundary.random();
        const Edge& edge = graph.edge(e);
        size_t v = edge.from();
        size_t u = edge.to();
        double diff = 0.0;
        if (!module_vertices[v]) {
            diff += graph.weight(v);
        }
        if (!module_vertices[u]) {
            diff += graph.weight(u);
        }
        diff += edge.weight();
        if (accepts(diff)) {
            add_edge(e);
        }
    }

    void SimulatedAnnealing::remove_from_module() {
        size_t e = module_edges.random();
        const Edge& edge = graph.edge(e);
        // TODO TODO TODO TODO TODO TODO TODOOOOOO TODODODODO
    }

    StandardUniformDistribution::StandardUniformDistribution(RandomEngine &re) :re(re) {}

    double StandardUniformDistribution::operator()(){
        return unif(re);
    }
}