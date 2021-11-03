//
// Created by Alexander Loboda on 8/24/21.
//

#include "include/relax_and_cut.h"
#include "include/components.h"
#include "include/primal_heuristic.h"
#include <ostream>

namespace {

template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    int size_s = std::snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
    if(size_s <= 0){ throw std::runtime_error( "Error during formatting." ); }
    auto size = static_cast<size_t>(size_s);
    std::unique_ptr<char[]> buf(new char[size]);
    std::snprintf(buf.get(), size, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

}

namespace relax {

Solver::Solver(const mwcsr::Graph& graph, const Parameters& params, std::ostream& os)
                                                                    : lb(0.0), ub(std::numeric_limits<double>::infinity()),
                                                                      subgradient_norm(0.0),
                                                                      no_improvement(0),
                                                                      g(graph),
                                                                      edges(g.edgeset_size()),
                                                                      active_variables(g.edgeset_size() + g.num_signals()),
                                                                      parameters(params),
                                                                      os(os) {
    g.absorb_vertex_signals();

    for (unsigned i = 0; i < g.num_signals(); i++) {
        signal_variables.push_back(variable_factory.take(g.signal_weight(i), "s" + std::to_string(i)));
    }

    for (unsigned i = 0; i < g.edgeset_size(); i++) {
        edge_variables.push_back(variable_factory.take(0.0, "e" + std::to_string(i)));
    }

    auto append = [this](const std::vector<Variable> vars) {
        variables.insert(variables.end(), vars.begin(), vars.end());
    };

    append(signal_variables);
    append(edge_variables);

    add_initial_cuts();
}

void Solver::add_initial_cuts() {
    std::vector<size_t> signal_cuts;
    std::vector<Cut> initial_cuts;

    for (size_t s = 0; s < g.num_signals(); s++) {
        signal_cuts.push_back(initial_cuts.size());
        initial_cuts.push_back(Cut({signal_variables[s]}, {}));
    }

    for (size_t e = 0; e < g.edgeset_size(); e++) {
        for (size_t s: g.edge(e).edge_signals()) {
            initial_cuts.push_back(Cut({edge_variables[e]}, {signal_variables[s]}));
            initial_cuts[signal_cuts[s]].rhs() += {edge_variables[e]};
        }
    }
    for (Cut& c: initial_cuts) {
        cuts.add(c);
    }
}

void Solver::solve() {
    double alpha = parameters.alpha;
    if (parameters.verbose) {
        os << "Solving SGMWCS problem.\n";
        os << "CV = currently violated\nCN = currently non violated\nCO = nonviolated cuts beyond maximum age\n";
        os << string_format("%10s%6s%10s%12s%7s%7s%7s%9s%9s\n",
                            "Iteration", "Alpha", "Obj", "Best Bound", "CV", "CN", "CO", "FixedTo0", "FixedTo1") << std::endl;
    }
    for (unsigned it = 0; it < parameters.iterations; it++) {
        monitor.check();

        calculate_current_solution();
        current_bound = objective();
        if (current_bound < ub) {
            ub = current_bound;
            no_improvement = 0;
        } else {
            ++no_improvement;
        }
        if (lb + EPS >= ub) {
            if (parameters.verbose) {
                print_stats(it, alpha);
            }
            return;
        }

        check_previous_cuts();
        if (it % parameters.sep_period == 0) {
            auto edges = solution_subgraph();
            separate_cuts(edges);
        }
        bool solution_found = false;
        if (it % parameters.heur_period == 0) {
            auto sol = primal_heuristic();
            if (sol.objective() > lb) {
                solution_found = true;
                lb = sol.objective();
                best_solution = sol;
            }
        }
        probing(current_bound);
        if (no_improvement >= parameters.beta_iterations) {
            no_improvement = 0;
            alpha /= 2;
        }
        if (parameters.verbose && (solution_found || it % parameters.report_period == 0)) {
            print_stats(it, alpha);
        }
        update_multipliers(alpha);
    }
}

void Solver::set_monitor(const mwcsr::monitor& m) {
    monitor = m;
}

void Solver::calculate_current_solution() {
    reset_variable_weights();
    cuts.calculate_variable_weights();
    for (auto v: active_variables.all_active()) {
        variables[v].setInstantValue();
    }
}

void Solver::reset_variable_weights() {
    for (auto v: active_variables.all_active()) {
        variables[v].reset_weight();
    }
}

double Solver::objective() {
    double obj = 0.0;
    for (auto v: active_variables.all_active()) {
        obj += variables[v].objective_part();
    }
    obj += cuts.objective_part();
    return obj;
}

std::vector<size_t> Solver::solution_subgraph() {
    std::vector<size_t> active_edges;
    auto active = edges.all_active();
    std::copy_if(active.begin(), active.end(), std::back_inserter(active_edges), [this](int edge) {
        return edge_variables.at(edge).instant_value() == 1;
    });
    return active_edges;
}

void Solver::separate_cuts(std::vector<size_t> edges) {
    std::vector<Component> components = Component::get_components(g, edges);
    for (int i = 0; i < (int)components.size() - 1; i++) {
        for (int j = i + 1; j < (int)components.size() - 1; j++) {
            auto& bigger = components[i];
            auto& smaller = components[i + 1];
            auto find_best = [this](const Component& c) -> size_t {
                auto ids = c.component_edges();
                return *std::max_element(ids.begin(), ids.end(), [this](size_t i, size_t j) {
                    return edge_variables[i].weight() < edge_variables[j].weight();
                });
            };
            size_t from_bigger = find_best(bigger);
            size_t from_smaller = find_best(smaller);
            Cut cut({edge_variables[from_bigger], edge_variables[from_smaller]}, {});
            cut.rhs() += 1;
            for (size_t e: components.at(i).component_edges()) {
                cut.rhs() += {edge_variables[e]};
            }
            if (cuts.add(cut)) {
                break;
            }
        }
    }
}

void Solver::check_previous_cuts() {
    subgradient_norm = cuts.check_previous(parameters.max_age);
}

void Solver::probing(double bound) {
    // check if variable can be fixed to zero or one, remove zero valued variable from everywhere
    for (auto i: active_variables.all_active()) {
        auto v = variables[i];
        if (v.instant_value() == 1) {
            if (bound - v.weight() < lb) {
                v.fix_value(1);
            }
        } else {
            if (bound + v.weight() < lb) {
                v.fix_value(0);
            }
        }
    }

    std::vector<Component> components = Component::get_components(g, edges.all_active());
    for (auto& c: components) {
        if (c.get_revenue() + EPS < lb) {
            for (size_t edge: c.component_edges()) {
                edge_variables.at(edge).fix_value(0);
            }
        }
    }

    cuts.try_fix();
    cuts.normalize();

    for (size_t e: edges.all_active()) {
        if (edge_variables[e].fixed()) {
            if (edge_variables[e].instant_value() == 0) {
                edges.remove(e);
                g.remove_edge(e);
            }
        }
    }

    for (size_t i: active_variables.all_active()) {
        auto v = variables[i];
        // TODO: exclude but consider
        if (v.fixed() && v.instant_value() == 0) {
            active_variables.remove(i);
        }
    }
}

void Solver::update_multipliers(double alpha) {
    double theta = alpha * (current_bound - lb) / subgradient_norm;
    cuts.step(theta);
}

Solution Solver::primal_heuristic() {
    auto f = [this](size_t e) -> double{
        return edge_variables.at(e).weight();
    };
    auto active = edges.all_active();
    std::vector<bool> current;
    std::transform(active.begin(), active.end(), std::back_inserter(current), [this](size_t e) {
        return edge_variables.at(e).instant_value() == 1;
    });
    PrimalHeuristic ph(g, f, active, current);
    return ph.run_heuristic();
}

Solution Solver::solution() const {
    return best_solution;
}

void Solver::print_stats(unsigned it, double alpha) const {
    // "Iteration", "Obj", "Best Bound", "CV", "CN", "CO", "FixedTo1", "FixedTo0");
    unsigned cv = 0, co = 0, cn = 0;
    for (size_t i = 0; i < cuts.size(); i++) {
        const auto& c = cuts.get_const(i);
        if (c.violated()) {
            cv++;
        } else if (c.non_violated_series() >= parameters.max_age && c.mutliplier() == 0.0 && c.subderivative() == 0.0) {
            co++;
        } else {
            cn++;
        }
    }
    unsigned f0 = 0, f1 = 0;
    for (Variable v: variables) {
        if (v.fixed()) {
            if (v.instant_value() == 0) {
                f0++;
            } else {
                f1++;
            }
        }
    }
    os << string_format("%10d%6.2f%10.3f%12.3f%7d%7d%7d%9d%9d", it, alpha, lb, ub, cv, cn, co, f0, f1) << std::endl;
}

double Solver::upper_bound() const {
    return ub;
}

}
