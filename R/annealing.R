sa_class <- "simulated_annealing_solver"
schedules <- c("fast", "boltzmann")

check_sa_solver <- function(solver) {
    if (!inherits(solver, sa_class)) {
        stop("Function called with an invalid SA solver instance")
    }
}

#' @export
parameters.simulated_annealing_solver <- function(solver) {
    params(parameter("normalization", type = "logical"),
        parameter("schedule", type = "mc", mc = schedules),
        parameter("initial_temperature", type = "float", positive = TRUE),
        parameter("final_temperature", type = "float", positive = TRUE),
        parameter("verbose", type = "logical"))
}

#' ctor for annealing solver
#'
#' Simulated annealing is a heuristic method of solving optimization problems.
#' Typically, it allows to find some good solution in a short time. This
#' implementation doesn't compute any upper bound on solution, so there is no
#' guarantee of optimality of solution provided.
#' @param normalization adjust all weights so that all possible changes have
#' mean square 1.0
#' @param schedule boltzmann annealing or fast annealing
#' @param initial_temperature initial value for the temperature
#' @param final_temperature final value for the temperature
#' @param verbose whether be verbose or not
#' @export
annealing_solver <- function(normalization = TRUE,
                             schedule = c("fast", "boltzmann"),
                             initial_temperature = 1.0,
                             final_temperature = 1e-6,
                             verbose = FALSE) {
    solver_ctor(c(sa_class, mwcs_solver_class))
}

normalize_weights <- function(instance) {
    weights <- c(V(instance)$weight, E(instance)$weight)
    support <- sum(weights[weights > 0]) / 2
    if (any(weights < -EPS)) {
        support <- min(support, -min(weights))
    }
    instance
}

#' Solve generalized maximum weight subgraph problem using simulated annealing
#'
#' @inheritParams solve_mwcsp.sgmwcs_solver
#' @param solver An annealing solver object.
#' @param instance An igraph instance to work with.
#' @return An object of class mwcsp_solution.
#' @export
solve_mwcsp.simulated_annealing_solver <- function(solver, instance,
                                                   signals = NULL, ...) {
    if (!inherits(solver, sa_class)) {
        stop("Not a simulated annealing solver")
    }
    instance <- normalize_signals(instance, signals)

    inst_rep <- instance_from_graph(instance$graph)
    inst_rep[["vertex_signals"]] <- as.integer(igraph::V(instance$graph)$signal) - 1
    inst_rep[["edge_signals"]] <- as.integer(igraph::E(instance$graph)$signal) - 1
    inst_rep[["signal_weights"]] <- instance$signals

    res <- sa_solve(inst_rep, solver)
    if (length(res$edges) == 0) {
        g <- igraph::induced_subgraph(instance$graph, vids = res$vertices)
    } else {
        g <- igraph::subgraph.edges(instance$graph, eids = res$edges)
    }
    weight <- sum(V(g)$weight) + sum(E(g)$weight)
    solution(g, weight, solved_to_optimality = NA)
}
