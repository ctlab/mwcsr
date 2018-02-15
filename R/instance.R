#' ctor for mwcs_instance
#' @param graph a graph
#' @param parse_vertex_weights wether or not parse vertex attribute "weight"
#' @param parse_edge_weights wether or not parse edge attribute "weight"
#' @param parse_budgets wether or not parse vertex attribute "budget"
#' @export
mwcs_instance <- function(graph,
                          parse_vertex_weights = TRUE,
                          parse_edge_weights = FALSE,
                          parse_budgets = FALSE) {
    obj <- structure(list(graph = graph), class = "mwcs_instance")
    if (!igraph::is_igraph(graph)) {
        stop("Not a graph object")
    }

    if(igraph::is_directed(graph)){
        stop("Not an undirected graph")
    }

    if (parse_vertex_weights) {
        if (!("weight" %in% list.vertex.attributes(graph))) {
            stop("Couldn't parse vertex weights: no such vertex attribute")
        }
        vertex_weights(obj) <- V(graph)$weight
    } else {
        vertex_weights(obj) <- rep(0, length(V(graph)))
    }

    if (parse_edge_weights) {
        if (!("weight") %in% list.edge.attributes(graph)) {
            stop("Couldn't parse edge weights: no such edge attribute")
        }
        edge_weights(obj) <- E(graph)$weight
    }

    if (parse_budgets) {
        if (!("budget") %in% list.vertex.attributes(graph)) {
            stop("Couldn't parse budgets: no such vertex attribute")
        }
        budgets(obj) <- V(graph)$budget
    }
    obj$solved_to_optimality <- FALSE
    obj
}

check_mwcs <- function(x) {
    if (!inherits(x, "mwcs_instance")) {
        stop("Not a MWCS instance")
    }
}

size <- function(instance) {
    if (!inherits(x, "mwcs_instance")) {
        stop("Not a MWCS instance")
    }
    len(V(x$graph))
}

good_values <- function(values) {
    values <- setNames(as.numeric(values), names(values))
    if (any(is.na(values))) {
        stop("NA weight presented or came from coercion to numeric type")
    }
    values
}

set_parameter <- function(x, value, parameter, sub_op, initial = FALSE) {
    check_mwcs(x)
    value <- good_values(value)

    elements <- sub_op(x$graph)

    if (initial) {
        x[[parameter]] <- rep(0.0, length(elements))
    }

    if (length(names(value)) > 0) {
        diff <- setdiff(names(value), names(elements))
        if (length(diff) != 0) {
            stop(paste("These elemenents are not presented in the graph: ",
                       names(diff)))
        }
        positions <- match(names(value), names(elements))
        x[[parameter]][positions] <- value
    } else {
        if (length(value) > length(elements)) {
            stop("Too many values to assign")
        }
        x[[parameter]][1:length(elements)] <- value
    }
    x
}

#' assignment operator for vertex weights
#' @param x a variable name.
#' @param value a value to be assigned to x.
#' @export
`vertex_weights<-` <- function(x, value) {
    set_parameter(x, value, parameter = "vertex_weights", sub_op = igraph::V)
}

#' assignment operator for edge weights
#' @param x a variable name.
#' @param value a value to be assigned to x.
#' @export
`edge_weights<-` <- function (x, value) {
    x <- set_parameter(x, value, parameter = "edge_weights", sub_op = igraph::E,
                  initial = inherits(x, "gmwcs_instance"))
    if (!inherits(x, "gmwcs_instance")) {
        class(x) <- c("gmwcs_instance", class(x))
    }
    x
}

#' assignment operator for budgets
#' @param x a variable name.
#' @param value a value to be assigned to x.
#' @export
`budgets<-` <- function(x, value) {
    x <- set_parameter(x, value, parameter = "budgets", sub_op = igraph::V,
                  initial = inherits(x, "budget_mwcs_instance"))
    if (!inherits(x, "budget_mwcs_instance")) {
        class(x) <- c("budget_mwcs_instance", class(x))
    }
    x
}

#' assignment operator for the root
#' @param x a variable name.
#' @param value a value to be assigned to x.
#' @export
`root<-` <- function(x, value) {
    check_mwcs(x)
    if (is.null(value)) {
        class(x) <- class(x)[class(x) != "rooted_mwcs_instance"]
        return(x)
    }
    if (is.integer(value)) {
        if (value < 1 | value > length(V(x$graph))) {
            stop("No such vertex in the graph")
        }
        x$root <- value
    } else if (is.character(value)) {
        if (!(value %in% names(V(x$graph)))) {
            stop("No such vertex in the graph")
        }
        x$root <- which(names(V(x$graph)) == value)
    } else {
        stop("Argument should be integer of vertex name")
    }
    if (!inherits(x, "rooted_mwcs_instance")) {
        class(x) <- c("rooted_mwcs_instance", class(x))
    }
    x
}

#' equality operator for mwcs instances
#' @param x,y mwcs instances
#' @export
`==.mwcs_instance` <- function(x, y) {
    setequal(class(x), class(y)) &
    setequal(names(x), names(y)) &
    all.equal(x[names(y)], y[names(y)])
}
