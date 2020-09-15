mwcs_solver_class <- "mwcs_solver"
parameter_class <- "solver_parameter"

EPS <- 1e-8

instance_from_graph <- function(graph) {
    edgelist <- as_edgelist(graph, names = FALSE)
    list(edgelist = matrix(as.integer(edgelist), ncol = 2),
         size = length(V(graph)))
}

vertex_attr_set <- list(list  = igraph::list.vertex.attributes,
                        query = igraph::vertex_attr,
                        name  = "vertex")

edge_attr_set <- list(list  = igraph::list.edge.attributes,
                      query = igraph::edge_attr,
                      name  = "edge")

attr_values <- function(graph, attr, type = "V", nonnegative = FALSE) {
    attr_set <- NULL
    if (type == "V") {
        attr_set <- vertex_attr_set
    } else {
        attr_set <- edge_attr_set
    }
    if (!(attr %in% attr_set$list(graph))) {
        stop(paste0("Can't find ", attr_set$name, " attribute '", attr, "'"))
    }
    values <- as.numeric(attr_set$query(graph, attr))
    if (any(is.na(values))) {
        stop(paste0("NA value found for the ", attr_set$name, " attribute '",
                    attr, "' or it came from coercion to numeric type"))
    }
    if (nonnegative) {
        if (any(values < 0)) {
            stop(paste0("All values of ", attr_set$name, " attribute '", attr,
                        "' must be nonnegative"))
        }
    }
    values
}

check_mwcs_solver <- function(x) {
    if (!inherits(x, mwcs_solver_class)) {
        stop("Not a MWCS solver")
    }
}

solver_ctor <- function(classes) {
    params <- as.list(parent.frame(1L))
    x <- structure(list(), class = classes)
    do.call(set_parameters, c(list(solver = x), params))
}

#' Solves a MWCS instance
#' @param solver a solver object.
#' @param instance an MWCS instance.
#' @param ... other arguments passed to other methods.
#' @export
solve_mwcsp <- function(solver, instance, ...) {
    check_mwcs_solver(solver)

    if (!igraph::is_igraph(instance)) {
        stop("Not a graph object")
    }

    if(igraph::is_directed(instance)){
        stop("Not an undirected graph")
    }

    UseMethod("solve_mwcsp")
}

solve_mwcsp.default <- function(solver, instance) {
    stop("An abstract solver can't solve an instance")
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
#' @export
parameters <- function(solver) UseMethod("parameters")

#' Sets time limitation for a solver
#' @param x a variable name.
#' @param value a value to be assigned to x.
#' @export
`timelimit<-` <- function(x, value) {
    set_parameters(x, timelimit = value)
}


#' Check the type and the validity of an MWCS instance
#' @param instance `igraph` object, containing an instance to be checked
#' @return a list with members `type` containing the type of the instance,
#' `valid` -- boolean flag indicating whether the instance is valid or not,
#' `errors` -- a character vector containing the error messages
#' @examples
#' data(mwcs_example)
#' get_instance_type(mwcs_example)
#' @export
get_instance_type <- function(instance) {
    res <- list(type="unknown", valid=FALSE, errors=character())
    if ("signals" %in% names(graph.attributes(instance))) {
        res$type <- "SGMWCS"

        if (!is.numeric(instance$signals) || !all(is.finite(instance$signals))) {
            res$errors <- "`signals` attribute is not a vector of finite numbers"
            return(res)
        }

        if (!"signal" %in% list.vertex.attributes(instance)) {
            res$errors <- "No `signal` attribute for nodes"
            return(res)
        }

        if (!all(V(instance)$signal %in% names(instance$signals))) {
            res$errors <- "All node signals should be present in `signals` graph attribute"
            return(res)
        }

        if (!"signal" %in% list.edge.attributes(instance)) {
            res$errors <- "No `signal` attribute for edges"
            return(res)
        }

        if (!all(E(instance)$signal %in% names(instance$signals))) {
            res$errors <- "All edge signals should be present in `signals` graph attribute"
            return(res)
        }

        res$valid <- TRUE
        return(res)
    }

    if ("weight" %in% list.edge.attributes(instance) &&
        "weight" %in% list.vertex.attributes(instance)) {
        res$type <- "GMWCS"

        if (!all(is.finite(V(instance)$weight))) {
            res$errors <- "Nodes `weight` attribute is not a vector of finite numbers"
            return(res)
        }

        if (!all(is.finite(E(instance)$weight))) {
            res$errors <- "Edges `weight` attribute is not a vector of finite numbers"
            return(res)
        }


        res$valid <- TRUE
        return(res)
    }

    if ("weight" %in% list.vertex.attributes(instance)) {
        res$type <- "MWCS"

        if (!all(is.finite(V(instance)$weight))) {
            res$errors <- "Nodes `weight` attribute is not a vector of finite numbers"
            return(res)
        }

        res$valid <- TRUE
        return(res)
    }

    res$errors <- "Can't determine type of the instance"
    return(res)

}
