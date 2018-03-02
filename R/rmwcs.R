rmwcs_class <- "rmwcs_solver"

sep_methods <- c("strong", "fast")
subgradients <- c("classic", "average", "cft")

check_rmwcs_solver <- function(solver) {
    if (!inherits(solver, rmwcs_class)) {
        stop("Function called with an invalid rmwcs solver instance")
    }
}

#' @export
parameters.rmwcs_solver <- function(solver) {
    list(parameter("timelimit", type = "integer", positive =  TRUE,
                   is_null_possible = TRUE),
         parameter("max_iterations", type = "integer", positive = TRUE,
                   is_null_possible = TRUE),
         parameter("beta_iterations", type = "integer", positive = TRUE),
         parameter("separation", type = "mc", mc = sep_methods),
         parameter("sep_iterations", type = "integer", positive = TRUE),
         parameter("start_constraints", type = "logical"),
         parameter("pegging", type = "logical"),
         parameter("sep_iter_freeze", type = "integer", positive = TRUE),
         parameter("heur_iterations", type = "integer", positive = TRUE),
         parameter("subgradient", type = "mc", mc = subgradients),
         parameter("beta", type = "float"),
         parameter("verbose", type = "logical"))
}

#' Generates a rmwcs solver with corresponding parameters
#' @param timelimit Timelimit in seconds
#' @param max_iterations Maximum number of subgradient iterations
#' @param beta_iterations Number of nonimproving iterations until beta is halfed
#' @param separation Separation: "strong" of "fast"
#' @param sep_iterations Extending the life of non-violated inequalities
#' @param startcons Whether to add flow-conservation/degree cons at start
#' @param pegging Pegging (variable fixing)
#' @param sep_iter_freeze After how many iterations we are checking added ineqs
#' @param heur_iterations After how many iterations we are doing heuristics
#' @param subgradient Subgradient: "classic", "average", "cft"
#' @param beta Beta for subgradient
#' @param verbose Whether to print solving progress
#' @export
#' @import igraph
rmwcs <- function(timelimit = 1800L,
                  max_iterations = 1000L,
                  beta_iterations = 50L,
                  separation = "strong",
                  start_constraints = TRUE,
                  pegging = TRUE,
                  sep_iterations= 1L,
                  sep_iter_freeze = 1L,
                  heur_iterations = 1L,
                  subgradient = "classic",
                  beta = 2.0,
                  verbose = FALSE) {
    x <- structure(list(), class = c(rmwcs_class, mwcs_class))
    params <- mget(names(formals()))
    do.call(set_parameters(c(list(solver = x), params)))
}

#' @export
solve.rmwcs_solver <- function(solver, instance) {
    check_sa_solver(solver)
    check_mwcs(instance)
    instance_rep <- to_list(instance)
    features <- c(mwcs_class, budget_class, cardinality_class)
    check_features(instance, features)

    vs <- rmwcs_solve(instance_rep, solver)
    if (length(vs) > 1 && edge_problem) {
        vs <- vs[vs > length(V(g))]
        vs <- vs - length(V(g))
        subgraph.edges(g, eids = vs)
    } else {
        induced.subgraph(g, vids = vs)
    }
}
