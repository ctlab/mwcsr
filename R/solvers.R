mwcs_solver_class <- "mwcs_solver"
parameter_class <- "solver_parameter"

to_list <- function(instance) {
    l <- instance
    l$graph <- NULL
    class(l) <- NULL
    l$edgelist <- as_edgelist(instance$graph, names = FALSE)
    l$size <- length(V(instance$graph))
    l
}

good_values <- function(values) {
    values <- setNames(as.numeric(values), names(values))
    if (any(is.na(values))) {
        stop("NA weight presented or came from coercion to numeric type")
    }
    values
}

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

solver_ctor <- function(classes) {
    params <- as.list(parent.frame(1L))
    x <- structure(list(), class = classes)
    do.call(set_parameters, c(list(solver = x), params))
}

#' Solves a MWCS instance
#' @param solver a solver object
#' @param instance an MWCS instance
#' @export
solve_mwcsp <- function(solver, instance) {
    check_mwcs_solver(solver)

    if (!igraph:is_igraph(graph)) {
        stop("Not a graph object")
    }

    if(igraph::is_directed(graph)){
        stop("Not an undirected graph")
    }

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
`[<-.mwcs_solver` <- function(x, i = NULL, j = NULL, value) {
    params <- list(solver = x)
    params[[i]] <- value
    do.call(set_parameters, params)
}

#' @export
`[.mwcs_solver` <- function(x, i = NULL, j = NULL, value) {
    x[[i]]
}

#' @export
print.mwcs_solver <- function(x) {
    check_mwcs_solver(x)
    cat("MWCS Solver: ", class(x)[1], "\n\n")

    cat("Parameters:\n")
    get_param_name <- function(x) setNames(x$name, NULL)
    param_names <- sapply(parameters(x), get_param_name)
    params <- data.frame(name = param_names,
                         value = sapply(param_names, function(y) x[y]),
                         row.names = NULL)
    print(params, right = FALSE, row.names = FALSE)
    cat("\n")

    cat("Supported MWCS instance types: \n")
    cat(paste(sapply(features(x), function(x) paste0("  ", x)),
              collapse = "\n"))
    cat("\n")
    invisible(x)
}

#' @export
parameters.default <- function(...) {
    stop("Object is not an MWCS solver or is an abstract solver")
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
