gmwcs_class <- "gmwcs_solver"

check_gmwcs_solver <- function(solver) {
    if (!inherits(solver, gmwcs_class)) {
        stop("Function called with an invalid SA solver instance")
    }
}

#' @export
parameters.gmwcs_sovler <- function(solver) {
    l <- list(parameter("cplex_bin", type = "file"),
         parameter("cplex_jar", type = "file"),
         parameter("gmwcs_jar", type = "file"),
         parameter("threads", type = "integer", positive = TRUE),
         parameter("timelimit", type = "integer", positive = TRUE,
                   is_null_possible = TRUE),
         parameter("verbose", type = "logical"))
    names(l) <- lapply(l, function(x) x$name)
    l
}

#' Returns a solver object for GMWCS solver
#' @param cplex_bin a path to cplex library file
#' @param cplex_jar a path to cplex JNI jar file
#' @param gmwcs_jar a path to GMWCS solver jar file
#' @param threads number of threads for simultaneous computation
#' @param timelimits maximum number of seconds to solve the problem
#' @param verbose whether or not be verbose
#' @export
gmwcs_solver <- function (cplex_bin,
                          cplex_jar,
                          gmwcs_jar,
                          threads = parallel::detectCores(),
                          timelimit = NULL,
                          verbose = TRUE){
    x <- structure(list(), class = c(gmwcs_class, mwcs_class))
    params <- mget(names(formals()))
    do.call(set_parameters, c(list(solver = x), params))
}
