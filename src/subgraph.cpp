#include "include/subgraph.h"
#include "include/graph.h"

namespace {
    using mwcsr::Graph;
    using mwcsr::Edge;
}

namespace annealing {

    Subgraph::Subgraph(const Graph &supergraph) :dynamic_graph{(unsigned)supergraph.size()},
                                                 graph(supergraph),
                                                 module_edges{supergraph.edgeset_size()},
                                                 boundary{supergraph.edgeset_size()},
                                                 module_vertices{supergraph.size()},
                                                 vdegree(std::vector<size_t>(supergraph.size(), 0)),
                                                 tokens{supergraph.edgeset_size()}{

    }

    size_t Subgraph::size() const {
        return n_vertices;
    }

    double Subgraph::score() const {
        return weight;
    }

    void Subgraph::add_vertex(size_t v) {
        n_vertices++;
        module_vertices.add(v);
        for (Edge e: graph.neighbours(v)) {
            size_t e_id = e.num();
            if (!(module_edges.contains(e_id) || boundary.contains(e_id))) {
                boundary.add(e_id);
            }
        }
        weight += graph.weight(v);
    }


    void Subgraph::add_edge(size_t e) {
        boundary.remove(e);
        module_edges.add(e);
        Edge edge = graph.edge(e);
        size_t v = edge.from();
        size_t u = edge.to();
        ++vdegree[v];
        ++vdegree[u];
        if (!module_vertices.contains(v)) {
            add_vertex(v);
        }
        if (!module_vertices.contains(u)) {
            add_vertex(u);
        }
        tokens[e] = std::move(dynamic_graph.add(v, u));
        weight += edge.weight();
    }

    bool Subgraph::remove_edge(size_t e) {
        Edge edge = graph.edge(e);
        auto v = edge.from();
        auto u = edge.to();
        dynamic_graph.remove(std::move(tokens[e]));
        size_t comp_size = dynamic_graph.component_size(v);
        if (comp_size < n_vertices - 1 && comp_size != 1) {
            tokens[e] = dynamic_graph.add(v, u);
            return false;
        }
        module_edges.remove(e);
        weight -= edge.weight();

        --vdegree[v];
        --vdegree[u];

        if (comp_size == n_vertices) {
            boundary.add(e);
            return true;
        }

        if (comp_size == n_vertices - 1) {
            remove_vertex(u);
        } else {
            remove_vertex(v);
        }

        return true;
    }

    void Subgraph::remove_vertex(size_t v) {
        n_vertices--;
        for (Edge e: graph.neighbours(v)) {
            if (boundary.contains(e.num())) {
                boundary.remove(e.num());
            }
        }
        module_vertices.remove(v);
        weight -= graph.weight(v);
    }

    size_t Subgraph::boundary_size() const {
        return boundary.size();
    }

    size_t Subgraph::any_vertex() const {
        return module_vertices.content()[0];
    }

    size_t Subgraph::random_inner_edge(RandomEngine &re) const {
        return module_edges.random(re);
    }

    size_t Subgraph::random_adjacent_edge(RandomEngine &re) const {
        return boundary.random(re);
    }

    const mwcsr::Edge& Subgraph::edge(size_t e) const {
        return graph.edge(e);
    }

    bool Subgraph::contains_vertex(size_t v) const {
        return module_vertices.contains(v);
    }

    size_t Subgraph::degree(size_t v) const {
        return vdegree.at(v);
    }

    Module Subgraph::get_snapshot() const {
        return Module(graph, module_vertices.content(), module_edges.content());
    }

    size_t Subgraph::number_of_edges() const {
        return module_edges.size();
    }
}
