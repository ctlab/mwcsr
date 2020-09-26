#' @export
parameters.scip_stp <- function(solver) {
    params(parameter("cplex_bin", type = "file", is_null_possible=TRUE),
         parameter("cplex_jar", type = "file", is_null_possible=TRUE),
         parameter("threads", type = "integer", positive = TRUE),
         parameter("timelimit", type = "integer", positive = TRUE,
                   is_null_possible = TRUE),
         parameter("penalty", type = "float", nonnegative = TRUE, is_null_possible=FALSE),
         parameter("memory", type = "char"),
         parameter("mst", type = "logical"),
         parameter("log", type = "integer", is_null_possible=TRUE))
}


write_stp <- function(g, stp_file) {
    edges <- igraph::as_edgelist(g, names = FALSE)
    nodes <- cbind("T ", seq_along(V(g)), V(g)$weight)
    n <- length(V(g))
    append <- function(s) {
        write(s, file=stp_file, append=TRUE)

    }
    write("33D32945 STP File, STP Format Version 1.0\n", file=stp_file, append=FALSE)
    append("SECTION Graph")
    append(paste0("Nodes ", n))
    append(paste0("Edges ", length(edges)/2))
    utils::write_table(cbind("E ", edges, file=stp_file, append=TRUE))
    append("END\nSECTION Terminals")
    append(paste0("Terminals ", n))
    utils::write_table(nodes, file=stp_file, append=TRUE)
    append("END\nEOF")
}

run_solver <- function(solver, instance) {
    graph_file <- tempfile("graph", fileext = "stp")
    write_stp(instance, file.path(graph_file))
    solver$run_main()
    
}