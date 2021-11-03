//
// Created by Alexander Loboda on 8/24/21.
//

#ifndef MWCSR_RELAX_AND_CUT_H
#define MWCSR_RELAX_AND_CUT_H

#include <memory>

#include "cut.h"
#include "../../include/graph.h"
#include "active_pool.h"
#include "../../include/monitor.h"
#include "solution.h"

namespace relax {

struct Parameters {
    double alpha = 2.0;
    unsigned beta_iterations = 5;
    unsigned iterations = 1000;
    unsigned heur_period = 10;
    unsigned sep_period = 10;
    unsigned max_age = 5;
    unsigned report_period = 10;
    bool verbose = false;
    bool debug = true;
};

class Solver {
    double lb;
    double ub;
    double current_bound;
    double subgradient_norm;
    unsigned no_improvement;

    mwcsr::Graph g;
    Cuts cuts;
    VariableFactory variable_factory;
    Solution best_solution;

    ActivePool edges;
    ActivePool active_variables;
    Parameters parameters;

    std::vector<Variable> signal_variables;
    std::vector<Variable> edge_variables;
    std::vector<Variable> variables;

    mwcsr::monitor monitor;
    std::ostream& os;
public:
    explicit Solver(const mwcsr::Graph& g, const Parameters&, std::ostream& out);

    void solve();
    Solution solution() const;
    double upper_bound() const;
    void set_monitor(const mwcsr::monitor& m);
private:
    void add_initial_cuts();
    void calculate_current_solution();

    void reset_variable_weights();
    double objective();

    std::vector<size_t> solution_subgraph();
    void separate_cuts(std::vector<size_t>);

    void check_previous_cuts();
    Solution primal_heuristic();
    void probing(double bound);
    void update_multipliers(double alpha);

    void print_stats(unsigned, double) const;
};

}

#endif //MWCSR_RELAX_AND_CUT_H
