mwcs_class <- "mwcs_instance"
budget_class <- "budget_mwcs_instance"
rooted_class <- "rooted_mwcs_instance"
gmwcs_class <- "gmwcs_instance"
cardinality_class <- "cardinality_mwcs_instance"

remove_class <- function(x, class) {
    class(x) <- class(x)[class(x) != class]
    x
}

to_list <- function(instance) {
    l <- instance
    l$graph <- NULL
    class(l) <- NULL
    l$edgelist <- as_edgelist(instance$graph, names = FALSE)
    l$size <- length(V(instance$graph))
    l
}

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
    obj <- structure(list(graph = graph, upper_bound = NA, solution = NULL),
                     class = mwcs_class)
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
        if (!("budget_cost") %in% list.vertex.attributes(graph)) {
            stop("Couldn't parse budgets: no such vertex attribute")
        }
        budget_costs(obj) <- V(graph)$budget_cost
    }
    obj$solved_to_optimality <- FALSE
    obj
}

#' Returns the best known estimation on the upper bound of the score
#' @param instance an instance of the MWCS problem
#' @export
upper_bound <- function(instance) {
    instance$ub
}

#' Returns the most weighted connected subgraph
#' @param instance an instance of MWCS problem
#' @export
solution <- function(instance) {
    instance$solution
}

check_mwcs <- function(x) {
    if (!inherits(x, mwcs_class)) {
        stop("Not a MWCS instance")
    }
}

size <- function(instance) {
    check_mwcs(instance)
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
            warning(paste("No such vertices in graph:",
                       do.call(paste, c(list(sep = ", "), lapply(diff, list)))))
        }
        positions <- match(setdiff(names(value), diff), names(elements))
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
    check_mwcs(x)
    if (is.null(value)) {
        x$edge_weights <- NULL
        return(remove_class(x, gmwcs_class))
    }
    x <- set_parameter(x, value, parameter = "edge_weights", sub_op = igraph::E,
                  initial = inherits(x, gmwcs_class))
    if (!inherits(x, gmwcs_class)) {
        class(x) <- c("gmwcs_instance", class(x))
    }
    x
}

remove_budget_class <- function(x) {
    x <- remove_class(x, budget_class)
    x$budget <- NULL
    x$budget_costs <- NULL
    x
}

#' assignment operator for budget costs
#' @param x a variable name.
#' @param value a value to be assigned to x.
#' @export
`budget_costs<-` <- function(x, value) {
    check_mwcs(x)
    if (inherits(x, cardinality_class)) {
        stop("Budget and cardinality features are mutually exclusive")
    }
    if (is.null(value)) {
        return(remove_budget_class(x))
    }
    x <- set_parameter(x, value, parameter = "budget_costs", sub_op = igraph::V,
                  initial = inherits(x, budget_class))
    if (!inherits(x, budget_class)) {
        class(x) <- c(budget_class, class(x))
        x$budget <- Inf
    }
    x
}

#' assignment operator for the budget
#' @param x a variable name.
#' @param value a value to be assigned to x.
#' @export
`budget<-` <- function(x, value) {
    check_mwcs(x)
    if (inherits(x, cardinality_class)) {
        stop("Budget and cardinality features are mutually exclusive")
    }
    if (is.null(value)) {
        return(remove_budget_class(x))
    }
    value <- as.numeric(value)
    if (is.na(value)) {
        stop("Budget mustn't be a NA")
    }
    if (!inherits(x, budget_class)) {
        budget_costs(x) <- 0
    }
    x$budget <- value
    x
}

#' assignment operator for the maximum cardinality of solution
#' @export
#' @param x a variable name.
#' @param value a value to be assigned to x.
`max_cardinality<-` <- function(x, value) {
    check_mwcs(x)
    if (inherits(x, budget_class)) {
        stop("Budget and cardinality features are mutually exclusive")
    }
    if (is.null(value)) {
        x$cardinality <- NULL
        return(remove_class(x, cardinality_class))
    }
    value <- as.integer(value)
    if (is.na(value) | value < 0) {
        stop("Invalid argument")
    }
    if (!inherits(x, cardinality_class)) {
        class(x) <- c(cardinality_class, class(x))
    }
    x$cardinality <- value
    x
}

#' assignment operator for the root
#' @param x a variable name.
#' @param value a value to be assigned to x.
#' @export
`root<-` <- function(x, value) {
    check_mwcs(x)
    if (is.null(value)) {
        x$root <- NULL
        return(remove_class(x, rooted_class))
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
    if (!inherits(x, rooted_class)) {
        class(x) <- c(rooted_class, class(x))
    }
    x
}

#' equality operator for mwcs instances
#' @param x,y mwcs instances
#' @export
`==.mwcs_instance` <- function(x, y) {
    check_mwcs(x)
    check_mwcs(y)
    setequal(class(x), class(y)) &
    setequal(names(x), names(y)) &
    all.equal(x[names(y)], y[names(y)])
}
