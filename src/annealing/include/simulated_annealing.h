#ifndef ANNEALING_SIMULATED_ANNEALING_H
#define ANNEALING_SIMULATED_ANNEALING_H

#include <cstddef>
#include <random>
#include <functional>
#include <vector>

#include "../../include/graph.h"
#include "index.h"
#include "definitions.h"
#include "module.h"
#include "subgraph.h"
#include "../../include/monitor.h"

#include "DynamicGraph.h"
#include "cooling_schedule.h"

namespace annealing {
    using dgraph::DynamicGraph;
    using dgraph::EdgeToken;

    class StandardUniformDistribution {
        RandomEngine& re;
        std::uniform_real_distribution<> unif;
    public:
        explicit StandardUniformDistribution(RandomEngine& re);

        double operator()();
    };

    class SimulatedAnnealing {
        RandomEngine random_engine;
        const Graph graph;
        StandardUniformDistribution unif;

        Subgraph sub;
        double temperature = 0;

        double best_score = 0;
        Module best;

    public:
        explicit SimulatedAnnealing(const Graph& graph, RandomEngine& random_engine);

        void run(CoolingSchedule& schedule, monitor monitor);
        vector<size_t> vertices();
        vector<Edge> edges();

        void add_vertex(int v);
        void add_edge(int eid);

        double score() const;
    private:
        void strike();

        void empty_module_step();
        void edge_step();
        bool accepts(double diff);
        size_t uniform(size_t n);
        void add_from_bdr();
        void remove_from_module();
    };

}

#endif //ANNEALING_SIMULATED_ANNEALING_H
