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
#include "cooling_schedule.h"

namespace annealing {
    using namespace dgraph;

    class StandardUniformDistribution {
        RandomEngine& re;
        std::uniform_real_distribution unif;
    public:
        explicit StandardUniformDistribution(RandomEngine& re);
        double operator()();

    };

    class SimulatedAnnealing {
        RandomEngine& random_engine;
        StandardUniformDistribution unif;

        DynamicGraph dynamic_graph;
        const Graph& graph;
        Index module_edges;
        Index boundary;
        vector<bool> module_vertices;

        double score = 0;
        size_t size = 0;
        std::vector<EdgeToken> tokens;
        CoolingSchedule temp_func;
        double temperature;

        vector<size_t> best_module;
        double best_score = 0;

    public:
        explicit SimulatedAnnealing(const Graph& graph, RandomEngine& random_engine,
                                    CoolingSchedule& cooling_schedule);
        void run();

    private:
        void step();
        void empty_module_step();
        void edge_step();

        void add_vertex(size_t v);
        void add_edge(size_t e);
        bool remove_edge(size_t e);
        void remove_vertex(size_t v);

        bool accepts(double diff);

        size_t uniform(size_t n);

        void add_from_bdr();
        void remove_from_module();
    };

}

#endif //ANNEALING_SIMULATED_ANNEALING_H
