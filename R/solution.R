solution <- function(graph, weight, solved_to_optimality = FALSE, ...) {
    obj <- c(list(graph = graph, weight = weight,
                  solved_to_optimality = solved_to_optimality), list(...))
    structure(obj, class = "mwcsp_solution")
}
