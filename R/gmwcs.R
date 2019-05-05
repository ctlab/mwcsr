gmwcs_class <- "gmwcs_solver"

check_gmwcs_solver <- function(solver) {
    if (!inherits(solver, gmwcs_class)) {
        stop("Function called with an invalid SA solver instance")
    }
}

#' @export
parameters.gmwcs_solver <- function(solver) {
    params(parameter("cplex_bin", type = "file"),
         parameter("cplex_jar", type = "file"),
         parameter("threads", type = "integer", positive = TRUE),
         parameter("timelimit", type = "integer", positive = TRUE,
                   is_null_possible = TRUE),
         parameter("memory", type = "char"),
         parameter("verbose", type = "logical"))
}

java_init <- function(solver) {
    solver$jar <- system.file("java", "gmwcs-solver.java", package = "mwcsr")
    params <- c(
        paste0("-Xmx", solver$memory),
        "-Xss16M"
    )
    .jinit(classpath = solver$jar, parameters = params, force.init = TRUE)
    solver
}

#' Returns a solver object for GMWCS solver
#'
#' This solver uses reformulation of MWCS problem in terms of mixed integer
#' programming. The later problem can be efficiently solved with
#' commercial optimization software. This solver uses CPLEX and requires
#' it to be installed.
#' @param cplex_bin a path to cplex library file
#' @param cplex_jar a path to cplex JNI jar file
#' @param threads number of threads for simultaneous computation
#' @param timelimit maximum number of seconds to solve the problem
#' @param memory maximum amount of memory(-Xmx flag)
#' @param verbose whether or not be verbose
#' @export
gmwcs_solver <- function (cplex_bin,
                          cplex_jar,
                          threads = parallel::detectCores(),
                          timelimit = NULL,
                          memory = "2G",
                          verbose = TRUE){
    if (!requireNamespace("rJava", quietly = TRUE)) {
        stop("Package \"rJava\" needed for this function to work. Please install it.",
            call. = FALSE)
    }
    java_init(solver_ctor(c(gmwcs_class, mwcs_solver_class)))
}

#' Solve generalized maximum weight subgraph problem using CPLEX solver
#' @inheritParams solve_mwcsp
#' @return An object of class mwcsp_solution.
#' @export
solve_mwcsp.gmwcs_solver <- function(solver, instance, ...) {
}
