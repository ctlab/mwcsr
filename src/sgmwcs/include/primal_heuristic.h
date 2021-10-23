//
// Created by Alexander Loboda on 9/15/21.
//

#ifndef MWCSR_PRIMAL_HEURISTIC_H
#define MWCSR_PRIMAL_HEURISTIC_H

#include "../../include/graph.h"
#include "solution.h"

#include <functional>
#include <unordered_set>

namespace relax {

class PrimalHeuristic {
    mwcsr::Graph g;
    std::function<double(size_t)> weight;
    std::vector<size_t> active;
    std::vector<bool> current;
    std::unordered_map<size_t, size_t> visit;
    std::unordered_map<size_t, double> shortest_distance;
    std::unordered_map<size_t, int> previous;
    size_t iteration = 0;

    Solution best;
public:
    PrimalHeuristic(mwcsr::Graph g, std::function<double(size_t)> wf, std::vector<size_t> active,
                    std::vector<bool> current);
    Solution run_heuristic();
private:
    void run_from_point(size_t e);
    double wei(size_t e);
};

}

#endif //MWCSR_PRIMAL_HEURISTIC_H
