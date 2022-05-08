solution <- function(graph, weight, solved_to_optimality = FALSE, ...) {
    stopifnot((igraph::vcount(graph) == 0) || igraph::is_connected(graph))
    obj <- c(list(graph = graph, weight = weight,
                  solved_to_optimality = solved_to_optimality), list(...))
    structure(obj, class = "mwcsp_solution")
}

#' Calculate weight of the solution. MWCS, GMWCS and SGMWCS instances are supported
#' @param solution Either `mwcsp_solution` or `igraph`` object representing the solution
#' @return Weight of the subgraph
#' @export
#' @importFrom methods is
get_weight <- function(solution) {
    if (is(solution, "mwcsp_solution")) {
        return(get_weight(solution$graph))
    }
    type <- get_instance_type(solution)$type
    if (type == 'SGMWCS') {
        active_signals <- unique(c(V(solution)$signal, E(solution)$signal))
        scores <- c(solution$signals)
        return (sum(scores[active_signals]))
    }
    else if (type == 'GMWCS') {
        return (sum(c(V(solution)$weight, E(solution)$weight)))
    } else if (type == 'MWCS') {
        return (sum(c(V(solution)$weight)))
    } else {
        stop('unexpected graph type. Expected one of: MWCS, GMWCS, SGMWCS')
    }
}
