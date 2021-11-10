rnc_class <- "rnc_solver"

check_rnc_solver <- function(solver) {
    if (!inherits(solver, rnc_class)) {
        stop("Invalid relax&cut sgmwcs solver provided.")
    }
}

#' @export
parameters.rnc_solver <- function(solver) {
    params(parameter("max_iterations", type = "integer", positive = TRUE),
           parameter("beta_iterations", type = "integer", positive = TRUE),
           parameter("sep_iterations", type = "integer", positive = TRUE),
           parameter("heur_iterations", type = "integer", positive = TRUE),
           parameter("max_age", type = "integer", positive = TRUE),
           parameter("verbose", type = "logical"))
}

#' Construct relax-and-cut SGMWCS solver
#'
#' The solver is based on the same approach as rmwcs solver. Modifications to
#' the original scheme are introduced to tackle problems arising with introduction
#' of edge weights and signals. It is recommended to use rmwcs solver to solve MWCS
#' problems, while due to differences in primal heuristic it may be a good pratice
#' to run both solvers on the same problem.
#' @inheritParams rmwcs_solver
#' @return An object of class `mwcs_solver`
#' @seealso [rmwcs_solver]
#' @export
rnc_solver <- function(max_iterations = 1000L,
                       beta_iterations = 50L,
                       heur_iterations = 10L,
                       sep_iterations = 10L,
                       verbose = FALSE) {
    solver_ctor(c(rnc_class, mwcs_solver_class))
}

#' @rdname solve_mwcsp
#' @order 4
#' @export
solve_mwcsp.rnc_solver <- function(solver, instance, ...) {
    check_rnc_solver(solver)
    inst_type <- get_instance_type(instance)$type
    if (!(inst_type %in% c("MWCS", "GMWCS", "SGMWCS"))) {
        stop("Rnc solver supports only MWCS, GMWCS, SGMWCS problem instances")
    }

    signal_instance <- to_signal_instance(instance)

    inst_rep <- instance_from_graph(signal_instance)
    inst_rep[["vertex_signals"]] <- match(igraph::V(signal_instance)$signal, names(signal_instance$signals)) - 1
    inst_rep[["edge_signals"]] <- match(igraph::E(signal_instance)$signal, names(signal_instance$signals)) - 1
    inst_rep[["signal_weights"]] <- signal_instance$signals

    res <- sgmwcs_solve(inst_rep, solver)

    if (length(res$edges) == 0) {
        g <- igraph::induced_subgraph(instance, vids = res$vertices)
    } else {
        g <- igraph::subgraph.edges(instance, eids = res$edges)
    }

    weight <- get_weight(g)

    stopifnot(abs(weight - res$lb) < EPS)

    solution(g, weight, solved_to_optimality = abs(res$lb - res$ub) < EPS,
             upper_bound = res$ub)
}
