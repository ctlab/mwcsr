#ifndef ANNEALING_SIMULATED_ANNEALING_H
#define ANNEALING_SIMULATED_ANNEALING_H

#include <cstddef>
#include <random>
#include <functional>
#include <vector>

#include "graph.h"
#include "index.h"

#include "../dgraph/DynamicGraph.h"

namespace annealing {
    class SimulatedAnnealing {
        std::function<uint_fast32_t>& random_engine;
    public:
        explicit SimulatedAnnealing(Graph graph, std::function<uint_fast32_t>& random_engine);
    };

}

#endif //ANNEALING_SIMULATED_ANNEALING_H
