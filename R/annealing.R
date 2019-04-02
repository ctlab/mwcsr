sa_class <- "simulated_annealing_solver"
schedules <- c("fast", "boltzmann")

check_sa_solver <- function(solver) {
    if (!inherits(solver, sa_class)) {
        stop("Function called with an invalid SA solver instance")
    }
}

features.simulated_annealing_solver <- function(solver) {
    c(mwcs_class, gmwcs_class)
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
    weights <- c(instance$vertex_weights, instance$edge_weights)
    edges <- igraph::as_data_frame(instance$graph)
    for (i in 1:nrow(edges)) {
        with_endpoints <- instance$vertex_weights[as.integer(edges[i, ])]
        with_endpoints <- with_endpoints + instance$edge_weights[i]
        weights <- c(weights, with_endpoints)
    }
    weights <- c(weights, -weights)
    d <- sd(weights)
    instance$vertex_weights <- instance$vertex_weights / d
    instance$edge_weights <- instance$edge_weights / d
    instance
}

#' @export
solve_mwcsp.simulated_annealing_solver <- function(solver, instance) {
    if (solver$normalization) {
        instance <- normalize_weights(instance)
    }
    inst_rep <- to_list(instance)
    if (is.null(inst_rep$edge_weights)) {
        inst_rep$edge_weights <- rep(0, length(inst_rep$edgelist))
    }
    res <- sa_solve(inst_rep, solver)
    if (length(res$edges) == 0) {
        instance$solution <- igraph::induced_subgraph(instance$graph,
                                                      vids = res$vertices)
    } else {
        instance$solution <- igraph::subgraph.edges(instance$graph,
                                                    eids = res$edges)
    }
    instance
}
