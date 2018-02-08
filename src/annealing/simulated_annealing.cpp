#include "simulated_annealing.h"

namespace annealing {

    SimulatedAnnealing::SimulatedAnnealing(const Graph &graph, std::function<uint_fast32_t> &random_engine)
                                          :graph(graph),
                                           score(0.0),
                                           size(0),
                                           tokens{graph.edgeset_size()},
                                           random_engine(random_engine),
                                           dynamic_graph{graph.size()},
                                           module_vertices{graph.size(), random_engine},
                                           module_edges{graph.edgeset_size(), random_engine},
                                           boundary{graph.edgeset_size(), random_engine} {}

    void SimulatedAnnealing::run() {

    }

    void SimulatedAnnealing::step() {

    }

    void SimulatedAnnealing::add_vertex(size_t v) {
        size++;
        module_vertices.add(v);
        for (Edge e: graph.neighbours(v)) {
            size_t u = e.opposite(v);
            size_t e_id = e.num();
            if (!(module_edges.contains(e_id) || boundary.contains(e_id))) {
                boundary.add(e_id);
            }
        }
    }

    void SimulatedAnnealing::add_edge(size_t e) {
        boundary.remove(e);
        module_edges.add(e);
        Edge edge = graph.edge(e);
        size_t v = edge.from();
        size_t u = edge.to();
        if (!module_vertices.contains(v)) {
            add_vertex(v);
        }
        if (!module_vertices.contains(u)) {
            add_vertex(u);
        }
        tokens[e] = std::move(dynamic_graph.add(v, u));
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
        return true;
    }

    void SimulatedAnnealing::remove_vertex(size_t v) {
        size--;
        for (Edge e: graph.neighbours(v)) {
            boundary.remove(e.num());
        }
        module_vertices.remove(v);
    }
}