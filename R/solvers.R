mwcs_class <- "mwcs_solver"
parameter_class <- "solver_parameter"

check_mwcs_solver <- function(x) {
    if (!inherits(x, mwcs_class)) {
        stop("Not a MWCS solver")
    }
}

check_features <- function(instance, features) {
    extra_classes <- setdiff(class(instance), features)
    if (length(extra_classes) > 0) {
        stop(paste("The solver doesn't support these features:", extra_classes))
    }
}

#' @export
solve <- function(solver, instance) UseMethod("solve")

#' @export
set_parameters <- function(solver, ...) {
    params <- parameters(solver)
    actual <- list(...)
    for (param_name in names(actual)) {
        param <- params[[param_name]]
        if (!param_name %in% names(params)) {
            stop(paste("Unknown parameter", param_name, "for the solver"))
        }
        check_parameter(param, actual[[param_name]])
        solver[[param_name]] <- actual[[param_name]]
    }
    solver
}

#' @export
parameters.default <- function(...) {
    stop("Abstact solver doesn't have parameters")
}

#' @export
parameters <- function(solver) UseMethod("parameters")

#' Sets time limitation for a solver
#' @param x a variable name.
#' @param value a value to be assigned to x.
#' @export
`timelimit<-` <- function(x, value) UseMethod("timelimit<-")
