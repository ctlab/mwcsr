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
    l <- list(parameter("normalization", type = "logical"),
        parameter("schedule", type = "mc", mc = schedules),
        parameter("initial_temperature", type = "float", positive = TRUE),
        parameter("final_temperature", type = "float", positive = TRUE),
        parameter("verbose", type = "logical"))
    names(l) <- lapply(l, function(x) x$name)
    l
}

#' ctor for annealing solver
#' @param normalize adjust all weights so that all possible changes have
#' mean square 1.0
#' @param boltzmann annealing or fast annealing
#' @param initial_temperature initial temperature
#' @param final_temperature final temperature
#' @export
annealing_solver <- function(normalization = TRUE,
                             schedule = c("fast", "boltzmann"),
                             initial_temperature = 1.0,
                             final_temperature = 1e-6,
                             verbose = FALSE) {
    x <- structure(list(), class = c(sa_class, mwcs_solver_class))
    params <- mget(names(formals()))
    do.call(set_parameters, c(list(solver = x), params))
}

#' @export
solve_mwcsp.simulated_annealing_solver <- function(solver, instance) {
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
