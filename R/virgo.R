virgo_java_class = "ru.itmo.ctlab.virgo.Main"
virgo_class = "virgo_solver"

#' @export
parameters.virgo_solver <- function(solver) {
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

init_solver <- function(solver) {
    solver_name <- "virgo-solver.jar"
    solver_jar <- system.file("java", solver_name, package="mwcsr")
    command <- paste0("-XX:ActiveProcessorCount=", solver$threads)
    if (is.null(solver$cplex_bin) || is.null(solver$cplex_jar)) {
        command <- c(command, "-cp", solver_jar, virgo_java_class)
    } else {
        command <- c(command, sprintf("-Djava.library.path=%s", solver$cplex_bin),
                    "-cp", paste(solver_jar, solver$cplex_jar, sep=':'),
                    virgo_java_class
                    )

    }
    solver$run_main <- function(cli_args, loglevel = 2) {
        exit_code <- NULL
        if (loglevel == 0) {
            exit_code <- system2("java", c(command, cli_args), stdout = FALSE)
        } else {
            exit_code <- system2("java", c(command, cli_args))
        }
        if (exit_code != 0) {
            stop("Failed to run java-based solver.")
        }
    }
    solver
}

find_cplex_jar <- function(cplex_dir) {
    files <- list.files(cplex_dir, pattern = "cplex.jar", recursive=T, full.names = TRUE)
    if (length(files) > 0) {
        return(files[1])
    } else {
        return(NULL)
    }
}

find_cplex_bin <- function(cplex_dir) {
    files <- list.files(cplex_dir, pattern = "libcplex\\d+.(so|jnilib|dll)", recursive=T, full.names = TRUE)
    if (length(files) > 0) {
        return(dirname(files[1]))
    } else {
        return(NULL)
    }
}


#' Construct a virgo solver
#'
#' This solver uses reformulation of MWCS problem in terms of mixed integer
#' programming. The later problem can be efficiently solved with
#' commercial optimization software. Exact version of solver uses CPLEX and requires
#' it to be installed. CPLEX 12.7.1 or higher is required.
#'
#' The solver currently does not support repeated negative signals, i.e. every
#' negative signal should be present only once among all edges and vertices.
#'
#' You can access solver directly using `run_main` function. See example.
#' @param cplex_dir a path to dir containing cplex_bin and cplex_jar,
#'        setting this to NULL sets `mst`` param to `TRUE`
#' @param cplex_bin a path to cplex binary dir
#' @param cplex_jar a path to cplex jar file
#' @param threads number of threads for simultaneous computation
#' @param timelimit maximum number of seconds to solve the problem
#' @param memory maximum amount of memory(-Xmx flag)
#' @param log verbosity level
#' @param penalty additional edge penalty for graph edges
#' @param mst whether to use approximate MST solver, no CPLEX files required with this parameter
#'            is set to `TRUE`
#' @return An object of class `mwcs_solver`.
#' @references Loboda A., Artyomov M., and Sergushichev A. (2016)
#' "Solving generalized maximum-weight connected subgraph problem for network enrichment analysis"
#'  doi:10.1007/978-3-319-43681-4_17
#' @export
#' @examples
#' data("sgmwcs_small_instance")
#' approx_vs <- virgo_solver(mst=TRUE, threads = 1)
#' approx_vs$run_main("-h")
#' sol <- solve_mwcsp(approx_vs, sgmwcs_small_instance)
#' \dontrun{
#' vs <- virgo_solver(cplex_dir='/path/to/cplex')
#' sol <- solve_mwcsp(approx_vs, sgmwcs_example)
#' }
virgo_solver <- function (cplex_dir,
                          threads = parallel::detectCores(),
                          timelimit = NULL,
                          penalty = 0.0,
                          memory = "2G",
                          log = 0,
                          cplex_bin=NULL,
                          cplex_jar=NULL,
                          mst=FALSE) {

    if (missing(cplex_dir)) {
        if (!mst && (is.null(cplex_bin) || is.null(cplex_jar))) {
            stop("Either provide `cplex_dir` paramter or both `cplex_bin` and `cplex_jar`")
        }
        cplex_dir<-NULL
    } else {
        if (is.null(cplex_dir)) {
            mst <- TRUE
        }
    }


    if (is.null(cplex_jar) && !mst) {
        cplex_jar <- find_cplex_jar(cplex_dir)
        if (is.null(cplex_jar)) {
            stop(paste0("Could not find `cplex.jar` in ", cplex_dir))
        }
    }

    if (is.null(cplex_bin) && !mst) {
        cplex_bin <- find_cplex_bin(cplex_dir)
        if (is.null(cplex_bin)) {
            stop(paste0("Could not find libcplex files in ", cplex_dir))
        }
    } else if (!mst && !dir.exists(cplex_bin)) {
        # cplex_bin should be directory, not the lib-file itself
        cplex_bin <- dirname(cplex_bin)
    }

    rm(cplex_dir) # don't need this parameter anymore

    init_solver(solver_ctor(c(virgo_class,mwcs_solver_class)))
}

write_files <- function(g, nodes_file, edges_file, signals_file, signals) {
    edges <- igraph::as_edgelist(g, names = FALSE)
    write_tbl <- function(x, file, rn) {
        utils::write.table(x, file = file, quote = FALSE, sep = "\t",
                       row.names = rn, col.names = FALSE)
    }
    if (is.null(signals)) {
        node_weights <- V(g)$weight
        if ("weight" %in% list.edge.attributes(g)) {
            edge_weights <- E(g)$weight
        } else {
            edge_weights <- rep(0, ecount(g))
        }

        edges <- cbind(edges, edge_weights)
        vertices <- cbind(seq_along(V(g)), node_weights)
    } else {
        node_signals <- igraph::vertex_attr(g, "signal")
        edge_signals <- igraph::edge_attr(g, "signal")
        edges <- cbind(edges, edge_signals)
        vertices <- cbind(seq_along(V(g)), node_signals)
        write_tbl(signals, signals_file, TRUE)
    }
    write_tbl(vertices, nodes_file, FALSE)
    write_tbl(edges, edges_file, FALSE)
}

cli_args <- function(solver, nodes.file, edges.file, signals.file = NULL, stats.file = NULL) {
    loglevel <- max(1, solver$log)
    params <- c("-n", nodes.file, "-e", edges.file, "-m", solver$threads, "-l", loglevel)
    if (!is.null(solver$timelimit)) {
        params <- c(params, "-t", solver$timelimit)
    }
    if (!is.null(signals.file)) {
        params <- c(params, "-s", signals.file)
        params <- c(params, "-type", "sgmwcs")
    } else {
        params <- c(params, "-type", "gmwcs")
    }
    if (!is.null(solver$penalty)) {
        params <- c(params, "-p", solver$penalty)
    }
    if (!is.null(stats.file)) {
        params <- c(params, "-f", stats.file)
    }
    if (solver$mst) {
        params <- c(params, "-mst")
    }
    params
}

run_solver <- function(solver, instance, sgmwcs, signals = NULL) {
    graph_dir <- tempfile("graph")
    dir.create(graph_dir, showWarnings=FALSE)
    edges.file <- file.path(graph_dir, "edges.txt")
    nodes.file <- file.path(graph_dir, "nodes.txt")
    signals.file <- if (sgmwcs) file.path(graph_dir, "signals.txt") else NULL
    stats.file <- file.path(graph_dir, "stats.tsv")

    write_files(instance, nodes.file, edges.file, signals.file, signals)
    args <- cli_args(solver, nodes.file, edges.file, signals.file, stats.file)
    solver$run_main(args, solver$log)

    nodes <- utils::read.table(paste0(nodes.file, ".out"), comment.char = "#")
    edges <- utils::read.table(paste0(edges.file, ".out"), comment.char = "#")

    if (!sgmwcs) {
        nodes <- nodes[nodes[, 2] != "n/a", ]
        edges <- edges[edges[, 3] != "n/a", ]
    }

    stats <- utils::read.table(stats.file, header = T)
    mwcs <- (if (nrow(nodes) <= 1) {
        induced_subgraph(instance, as.integer(nodes[, 1]))
    } else {
        eids <- get.edge.ids(instance, t(edges[,1:2]))
        subgraph.edges(instance, eids)
    })
    list(mwcs=mwcs, stats=stats)
}

#' @rdname solve_mwcsp
#' @order 2
#' @export
#'
solve_mwcsp.virgo_solver <- function(solver, instance, ...) {
    inst_type <- get_instance_type(instance)
    if (inst_type$type == "SGMWCS" && inst_type$valid) {
        return(solve_sgmwcs(solver, instance, ...))
    } else if (inst_type$type %in% c("GMWCS", "MWCS") && inst_type$valid) {
        return(solve_gmwcs(solver, instance, ...))
    } else {
        msg <- "Not a valid MWCS, GMWCS or SGMWCS instance"
        if (inst_type$type != "unknown") {
            msg <- paste0(msg, sprintf("\nThe instance looks like %s", inst_type$type))
            if (!inst_type$valid) {
                msg <- paste0(msg, ", but there was an error:")
                msg <- paste0(c(msg, inst_type$errors), collapse = "\n")
            }
        }
        stop(msg)
    }
}

solve_sgmwcs <- function(solver, instance, ...) {
    signals <- instance$signals

    neg_signals <- names(which(signals < 0))
    if (any(table(c(V(instance)$signal, E(instance)$signal))[neg_signals] > 1, na.rm=TRUE)) {
        stop("Instances with repeated negative signals are not supported")
    }

    V(instance)$index <- seq_len(vcount(instance))
    E(instance)$index <- seq_len(ecount(instance))

    sol <- run_solver(solver, instance, sgmwcs = TRUE,
                      data.frame(score=instance$signals, row.names=names(instance$signals)))
    mwcs <- sol$mwcs
    sigs <- union(V(mwcs)$signal, E(mwcs)$signal)
    weight <- sum(signals[sigs])
    if (length(E(mwcs)) > 0) {
        ret <- igraph::subgraph.edges(instance, E(mwcs)$index)
    } else {
        ret <- igraph::induced_subgraph(instance, V(mwcs)$index)
    }
    solution(ret, weight, sol$stats$isOpt == 1, stats=sol$stats)
}

solve_gmwcs <- function(solver, instance, ...) {
    inst_type <- get_instance_type(instance)
    if (!inst_type$type %in% c("GMWCS", "MWCS") || !inst_type$valid) {
        stop("Not a valid GMWCS instance")
    }

    sol <- run_solver(solver, instance, sgmwcs = FALSE)
    mwcs <- sol$mwcs

    weight <- sum(V(mwcs)$weight, E(mwcs)$weight)
    solution(mwcs, weight, sol$stats$isOpt == 1, stats=sol$stats)
}
