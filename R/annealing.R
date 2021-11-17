sa_class <- "simulated_annealing_solver"
schedules <- c("fast", "boltzmann")

check_sa_solver <- function(solver) {
    if (!inherits(solver, sa_class)) {
        stop("Function called with an invalid SA solver instance")
    }
}

#' @export
parameters.simulated_annealing_solver <- function(solver) {
    params(parameter("schedule", type = "mc", mc = schedules),
        parameter("initial_temperature", type = "float", positive = TRUE),
        parameter("final_temperature", type = "float", positive = TRUE),
        parameter("verbose", type = "logical"))
}

#' Construct an annealing solver
#'
#' Simulated annealing is a heuristic method of solving optimization problems.
#' Typically, it allows to find some good solution in a short time. This
#' implementation doesn't compute any upper bound on solution, so there is no
#' guarantee of optimality of solution provided.
#'
#' Algorithm maintains connected subgraph staring with empty subgraph.
#' Each iteration one random action is considered. It may be a removal of a
#' vertex or an edge which does not separate graph or addition of an extra vertex or
#' an edge neighboring existing graph. If the subgraph is empty all vertices
#' are considered as candidates to form a subgaph. After candidate is chosen two
#' subgraph scores are considered: for a new subgraph and for an old one. Simulated
#' annealing operates with a notion of temperature. The candidate action is
#' accepted with probability: p(S'|S) = exp(-E / T), where E is weight difference
#' between subgraphs and T is current temperature.
#'
#' Temperature is calculated in each iteration. in `mwcsr` there are two
#' temperature schedules supported. So called Boltmann annealing uses the formula:
#' T(k) = T0 / (ln(1 + k)), while in case of fast annealing this one is used:
#' T(k) = T0 / k, where k is iteration number.
#'
#' To tune the algorithm it is useful to realize how typical changes in the goal
#' function for single actions are distributed. Calculating the acceptance probabilities
#' at initial temperature and final temperatures may help to choose schedule and
#' temperatures.
#'
#' @seealso [rnc_solver] will probably be a better choice with minimal tuning necessary
#' @param schedule boltzmann annealing or fast annealing
#' @param initial_temperature initial value for the temperature
#' @param final_temperature final value for the temperature
#' @param verbose whether be verbose or not
#' @return An object of class `mwcs_solver`
#' @export
annealing_solver <- function(schedule = c("fast", "boltzmann"),
                             initial_temperature = 1.0,
                             final_temperature = 1e-6,
                             verbose = FALSE) {
    solver_ctor(c(sa_class, mwcs_solver_class))
}

#' @rdname solve_mwcsp
#' @order 5
#' @param warm_start warm start solution, an object of the class mwcsp_solution.
#' @export
solve_mwcsp.simulated_annealing_solver <- function(solver, instance, warm_start, ...) {
    if (!inherits(solver, sa_class)) {
        stop("Not a simulated annealing solver")
    }

    if (!(get_instance_type(instance)$type %in% c("MWCS", "GMWCS", "SGMWCS"))) {
        stop("Annealing solver supports only MWCS, GMWCS, SGMWCS problem instances")
    }

    signal_instance <- to_signal_instance(instance)

    inst_rep <- instance_from_graph(signal_instance)
    inst_rep[["vertex_signals"]] <- match(igraph::V(signal_instance)$signal, names(signal_instance$signals)) - 1
    inst_rep[["edge_signals"]] <- match(igraph::E(signal_instance)$signal, names(signal_instance$signals)) - 1
    inst_rep[["signal_weights"]] <- signal_instance$signals

    if (!missing(warm_start)) {
        ws_sol <- warm_start$warm_start_solution
        inst_rep[["warm_start_vertices"]] <- ws_sol$vertices
        inst_rep[["warm_start_edges"]] <- ws_sol$edges
        inst_rep[["warm_start_weight"]] <- ws_sol$weight
    }

    res <- sa_solve(inst_rep, solver)
    if (length(res$edges) == 0) {
        g <- igraph::induced_subgraph(instance, vids = res$vertices)
    } else {
        g <- igraph::subgraph.edges(instance, eids = res$edges)
    }
    weight <- get_weight(g)
    res$weight <- weight
    solution(g, weight, solved_to_optimality = FALSE, warm_start_solution = res)
}
