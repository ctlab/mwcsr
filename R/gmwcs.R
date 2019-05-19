gmwcs_class <- "gmwcs_solver"

check_gmwcs_solver <- function(solver) {
    if (!inherits(solver, gmwcs_class)) {
        stop("Function called with an invalid SA solver instance")
    }
}

#' @export
parameters.gmwcs_solver <- function(solver) {
    params(parameter("cplex_bin", type = "file"),
         parameter("cplex_jar", type = "file"),
         parameter("threads", type = "integer", positive = TRUE),
         parameter("timelimit", type = "integer", positive = TRUE,
                   is_null_possible = TRUE),
         parameter("memory", type = "char"),
         parameter("verbose", type = "logical"))
}

java_init <- function(solver) {
    params <- c(
        paste0("-Xmx", solver$memory),
        "-Xss64M"
    )
    rJava::.jpackage("mwcsr")
    rJava::.jaddClassPath(solver$cplex_jar)
    rJava::J("java.lang.System")$load(solver$cplex_bin)
    withCallingHandlers(rJava::new(rJava::J("ilog.cplex.IloCplex")),
                        error = function(e) {
                            message("CPLEX solver cannot be initilized.")
                            message("Please check all the paths to CPLEX jar and binary file")
                            stop(e)
                        })
    solver
}

#' Returns a solver object for GMWCS solver
#'
#' This solver uses reformulation of MWCS problem in terms of mixed integer
#' programming. The later problem can be efficiently solved with
#' commercial optimization software. This solver uses CPLEX and requires
#' it to be installed.
#' @param cplex_bin a path to cplex library file
#' @param cplex_jar a path to cplex JNI jar file
#' @param threads number of threads for simultaneous computation
#' @param timelimit maximum number of seconds to solve the problem
#' @param memory maximum amount of memory(-Xmx flag)
#' @param verbose whether or not be verbose
#' @export
gmwcs_solver <- function (cplex_bin,
                          cplex_jar,
                          threads = parallel::detectCores(),
                          timelimit = NULL,
                          memory = "2G",
                          verbose = TRUE){
    if (!requireNamespace("rJava", quietly = TRUE)) {
        stop("Package \"rJava\" is required for this function to work. Please install it.",
            call. = FALSE)
    }
    java_init(solver_ctor(c(gmwcs_class, mwcs_solver_class)))
}

write_files <- function(g, nodes_file, edges_file) {
    node_weights <- attr_values(g, "weight", "V")
    edge_weights <- attr_values(g, "weight", "E")

    edges <- igraph::as_edgelist(g, names = FALSE)
    edges <- cbind(edges, edge_weights)
    vertices <- cbind(1:length(V(g)), node_weights)
    utils::write.table(vertices, file = nodes_file, quote = FALSE, sep = "\t",
                       row.names = FALSE, col.names = FALSE)
    utils::write.table(edges, file = edges_file, quote = FALSE, sep = "\t",
                       row.names = FALSE, col.names = FALSE)
}

#' Solve generalized maximum weight subgraph problem using CPLEX solver
#' @inheritParams solve_mwcsp
#' @return An object of class mwcsp_solution.
#' @export
solve_mwcsp.gmwcs_solver <- function(solver, instance, ...) {
    if (!inherits(solver, gmwcs_class)) {
        stop("Not a gmwcs solver")
    }
    java_init(solver)
    if (!"weight" %in% list.vertex.attributes(instance)) {
        warning("No `weight` vertex attribute. Setting weights to zero")
        V(instance)$weight <- 0
    }
    if (!"weight" %in% list.edge.attributes(instance)) {
        warning("No `weight` edge attribute. Setting weights to zero")
        E(instance)$weight <- 0
    }

    nodes_file <- tempfile()
    edges_file <- tempfile()
    write_files(instance, nodes_file, edges_file)
    on.exit(file.remove(nodes_file))
    on.exit(file.remove(edges_file))
    rJava::J("ru.ifmo.ctddev.gmwcs.Main")$main(rJava::.jarray("-h"))
}
