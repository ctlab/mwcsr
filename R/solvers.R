mwcs_solver_class <- "mwcs_solver"
parameter_class <- "solver_parameter"

check_mwcs_solver <- function(x) {
    if (!inherits(x, mwcs_solver_class)) {
        stop("Not a MWCS solver")
    }
}

check_features <- function(instance, features) {
    extra_classes <- setdiff(class(instance), features)
    if (length(extra_classes) > 0) {
        stop(paste("The solver doesn't support these features:", extra_classes))
    }
}

#' Solves a MWCS instance
#' @param solver a solver object
#' @param instance an MWCS instance
#' @export
solve_mwcsp <- function(solver, instance) {
    check_mwcs(instance)
    check_mwcs_solver(solver)
    check_features(instance, features(solver))

    UseMethod("solve_mwcsp")
}

solve_mwcsp.default <- function(solver, instance) {
    stop("An abstract solver can't solve an instance")
}

features <- function(solver) UseMethod("features")

features.default <- function (solver) {
    stop("An abstract solver doesn't have features")
}

#' Sets values of specific parameters
#' @param solver a solver
#' @param ... listed parameter names and values assigned to them
#' @export
set_parameters <- function(solver, ...) {
    params <- parameters(solver)
    actual <- list(...)
    for (param_name in names(actual)) {
        param <- params[[param_name]]
        if (!param_name %in% names(params)) {
            warning(paste("Unknown parameter", param_name, "for the solver"))
        } else {
            solver[[param_name]] <- check_parameter(param, actual[[param_name]])
        }
    }
    solver
}

#' @export
parameters.default <- function(...) {
    stop("Abstact solver doesn't have parameters")
}

#' The method returns all parameters supported by specific solver
#' @param solver a solver object
#' @export
parameters <- function(solver) UseMethod("parameters")

#' Sets time limitation for a solver
#' @param x a variable name.
#' @param value a value to be assigned to x.
#' @export
`timelimit<-` <- function(x, value) {
    set_parameters(x, timelimit = value)
}
