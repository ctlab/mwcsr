#ifndef ANNEALING_SIMULATED_ANNEALING_H
#define ANNEALING_SIMULATED_ANNEALING_H

#include <cstddef>
#include <random>
#include <functional>
#include <vector>

#include "graph.h"
#include "index.h"
#include "definitions.h"

#include "../dgraph/DynamicGraph.h"

namespace annealing {
    using namespace dgraph;

    class SimulatedAnnealing {
        RandomEngine& random_engine;
        DynamicGraph dynamic_graph;
        const Graph& graph;
        Index module_vertices;
        Index module_edges;
        Index boundary;
        double score;
        size_t size;
        std::vector<EdgeToken> tokens;

    public:
        explicit SimulatedAnnealing(const Graph& graph, RandomEngine& random_engine);
        void run();

    private:
        void step();
        void add_vertex(size_t v);
        void add_edge(size_t e);
        bool remove_edge(size_t e);
        void remove_vertex(size_t v);
    };

}

#endif //ANNEALING_SIMULATED_ANNEALING_H
