mwcs_solver_class <- "mwcs_solver"
parameter_class <- "solver_parameter"

EPS <- 1e-8

instance_from_graph <- function(graph) {
    edgelist <- as_edgelist(graph, names = FALSE)
    list(edgelist = matrix(as.integer(edgelist), ncol = 2),
         size = length(V(graph)))
}

attr_sets <- list(V = list(list  = igraph::list.vertex.attributes,
                      query = igraph::vertex_attr,
                      name  = "nodes"),
                  E = list(list  = igraph::list.edge.attributes,
                      query = igraph::edge_attr,
                      name  = "edges"))

validate_attr_values <- function(graph, attr, type = "V", nonnegative = FALSE) {
    attr_set <- attr_sets[[type]]

    if (!(attr %in% attr_set$list(graph))) {
        return(paste0("Couldn't find ", attr_set$name, " attribute '", attr, "'"))
    }
    values <- as.numeric(attr_set$query(graph, attr))
    if (!all(is.finite(values))) {
        return(paste0("All values of ", attr_set$name, " attribute '", attr,
                        "' must be finite."))
    }
    if (nonnegative) {
        if (any(values < 0)) {
            return(paste0("All values of ", attr_set$name, " attribute '", attr,
                        "' must be nonnegative."))
        }
    }
    NULL
}

check_mwcs_solver <- function(x) {
    if (!inherits(x, mwcs_solver_class)) {
        stop("Not a MWCS solver")
    }
}

check_graph <- function(g) {
    if (!igraph::is.igraph(g)) {
        stop("Not a graph object")
    }
    if (igraph::is_directed(g)) {
        stop("Not an undirected graph")
    }
    if (igraph::any_multiple(g)) {
        stop("Multiple edges (edges connecting the same vertices) are not supported")
    }
}

solver_ctor <- function(classes) {
    params <- as.list(parent.frame(1L))
    x <- structure(list(), class = classes)
    do.call(set_parameters, c(list(solver = x), params))
}

#' Solves a MWCS instance.
#'
#' Generic function for solving MWCS instances using solvers collected in the package.
#'
#' MWCS instance here is represented as an undirected graph, an `igraph` object.
#' The package supports four types of instances: Simple MWCS, Generalized MWCS,
#' Budget MWCS, signal MWCS problems. All the necessary weights and costs are
#' passed by setting vertex and edge attributes. See [get_instance_type] to check
#' if the `igraph` object is a correct MWCS instance. For Simple MWCS problem
#' numeric vertex attribute `weight` must be set. For generalized version `weight`s
#' can be provided for edges. For budget version of the problem in addition to
#' vertex weights it is required that `igraph` object would have `budget_cost` vertex
#' attribute with positive numeric values.
#'
#' Signal MWCS instance is quite different. There is no `weight` attribute for
#' neither vertices nor edges. Instead, vertex and edge attribute `signal` should
#' be provided with signal names. A numeric vector containing weights for the signals
#' should be assigned to graph attribute `signals`.
#'
#' See vignette for description of the supported problems. See `igraph` package
#' documentation for more details about getting/setting necesasry attributes.
#'
#' @order 1
#' @examples
#'
#' library(igraph)
#'
#' # for a MWCS instance
#'
#' data(mwcs_example)
#' head(V(mwcs_example)$weight)
#'
#' # for a GMWCS instance
#' data(gmwcs_example)
#' head(E(gmwcs_example)$weight)
#'
#' # for a SGMWCS instance
#' data(sgmwcs_example)
#' head(V(sgmwcs_example)$signal)
#' head(E(sgmwcs_example)$signal)
#'
#' head(sgmwcs_example$signals)
#'
#' @param solver a solver object returned by rmwcs_solver, annealing_solver, rnc_solver or virgo_solver.
#' @param instance an MWCS instance, an igraph object with problem-related vertex, edge and graph attributes. See details.
#' @param ... other arguments to be passed.
#' @return An object of class `mwcsp_solution` consisting of resulting subgraph,
#' its weight and other information about solution provided.
#' @export
solve_mwcsp <- function(solver, instance, ...) {
    check_mwcs_solver(solver)
    check_graph(instance)
    inst_type <- get_instance_type(instance)
    if (!inst_type$valid) {
        stop(paste("Not a valid instance, call get_instance_type() to see errors"))
    }

    UseMethod("solve_mwcsp")
}

#' Sets values of specific parameters
#' @param solver a solver
#' @param ... listed parameter names and values assigned to them
#' @return The solver with parameters changed.
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
print.mwcs_solver <- function(x, ...) {
    check_mwcs_solver(x)
    cat("MWCS Solver: ", class(x)[1], "\n\n")

    cat("Parameters:\n")
    get_param_name <- function(x) setNames(x$name, NULL)
    param_names <- sapply(parameters(x), get_param_name)
    params <- data.frame(name = param_names,
                         value = sapply(param_names,
                            function(y) if(is.null(x[y])) "NULL" else x[y]),
                         row.names = NULL)
    print(params, right = FALSE, row.names = FALSE)
    cat("\n")

    invisible(x)
}

#' @export
parameters.default <- function(...) {
    stop("Object is not an MWCS solver or it is an abstract solver")
}

#' The method returns all parameters supported by specific solver
#' @param solver a solver object
#' @return A table containing parameter names and possible values for each parameter.
#' @export
parameters <- function(solver) UseMethod("parameters")

#' Sets time limitation for a solver
#' @param x a variable name.
#' @param value a value to be assigned to x.
#' @return The solver with new timelimit set.
#' @export
`timelimit<-` <- function(x, value) {
    set_parameters(x, timelimit = value)
}

check_signals <- function(instance) {
    if (!is.numeric(instance$signals) || !all(is.finite(instance$signals))) {
        return("`signals` attribute is not a vector of finite numbers")
    }

    if (any(duplicated(names(instance$signals)))) {
        return("Graph `signals` attribute has duplicated names")
    }

    if (!"signal" %in% list.vertex.attributes(instance)) {
        return("No `signal` attribute for nodes")
    }

    if (!all(V(instance)$signal %in% names(instance$signals))) {
        return("All node signals should be present in `signals` graph attribute")
    }

    if (!"signal" %in% list.edge.attributes(instance)) {
        return("No `signal` attribute for edges")
    }

    if (!all(E(instance)$signal %in% names(instance$signals))) {
        return("All edge signals should be present in `signals` graph attribute")
    }
    NULL
}

#' Check the type and the validity of an MWCS instance
#' @param instance `igraph` object, containing an instance to be checked
#' @return A list with members `type` containing the type of the instance,
#' `valid` -- boolean flag indicating whether the instance is valid or not,
#' `errors` -- a character vector containing the error messages
#' @examples
#' data(mwcs_example)
#' get_instance_type(mwcs_example)
#' @return A list with two fields: the type of the instance with which it will
#' be treated by solve_mwcsp function and boolean showing validness of the instance.
#' @export
get_instance_type <- function(instance) {
    check_graph(instance)
    res <- list(type="unknown", valid=FALSE, errors=NULL)
    if ("signals" %in% names(graph.attributes(instance))) {
        res$type <- "SGMWCS"
        res$errors <- check_signals(instance)
    } else if ("weight" %in% list.edge.attributes(instance) &&
               "weight" %in% list.vertex.attributes(instance)) {

        res$type <- "GMWCS"
        res$errors <- c(validate_attr_values(instance, "weight", "V"),
                        validate_attr_values(instance, "weight", "E"))
    } else if ("weight" %in% list.vertex.attributes(instance)) {
        if ("budget_cost" %in% list.vertex.attributes(instance)) {
            res$type <- "Budget MWCS"
            res$errors <- validate_attr_values(instance, "budget_cost", "V",
                                               nonnegative = TRUE)
        } else {
            res$type <- "MWCS"
        }
        res$errors <- c(res$errors, validate_attr_values(instance,
                                                         "weight", "V"))
    } else {
        res$errors <- "Can't determine type of the instance"
    }

    res$valid <- is.null(res$errors)

    res
}
