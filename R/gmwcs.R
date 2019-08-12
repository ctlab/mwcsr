gmwcs_class <- "gmwcs_solver"
sgmwcs_class <- "sgmwcs_solver"

gmwcs_java_class <- "ru.ifmo.ctddev.gmwcs.Main"
sgmwcs_java_class <- "ru.itmo.ctlab.sgmwcs.Main"

check_gmwcs_solver <- function(solver) {
    if (!inherits(solver, gmwcs_class)) {
        stop("Function called with an invalid gmwcs solver instance")
    }
}

check_sgmwcs_solver <- function(solver) {
    if (!inherits(solver, sgmwcs_class)) {
        stop("Function called with an invalid sgmwcs solver instance")
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

#' @export
parameters.sgmwcs_solver <- function(solver) {
    parameters.gmwcs_solver(solver)
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

check_rjava <- function() {
    if (!requireNamespace("rJava", quietly = TRUE)) {
        stop("Package \"rJava\" is required for this function to work. Please install it.",
            call. = FALSE)
    }
}

#' Construct a GMWCS solver
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
    check_rjava()
    java_init(solver_ctor(c(gmwcs_class, mwcs_solver_class)))
}

#' @rdname gmwcs_solver
#' @export
sgmwcs_solver <- function (cplex_bin,
                           cplex_jar,
                           threads = parallel::detectCores(),
                           timelimit = NULL,
                           memory = "2G",
                           verbose = TRUE){
    check_rjava()
    java_init(solver_ctor(c(sgmwcs_class, mwcs_solver_class)))
}

write_files <- function(g, nodes_file, edges_file, signals_file, signals) {
    node_weights <- attr_values(g, "weight", "V")
    edge_weights <- attr_values(g, "weight", "E")

    edges <- igraph::as_edgelist(g, names = FALSE)
    write_tbl <- function(x, file) {
        utils::write.table(x, file = file, quote = FALSE, sep = "\t",
                       row.names = FALSE, col.names = FALSE)
    }
    if (is.null(signals)) {
        edges <- cbind(edges, edge_weights)
        vertices <- cbind(1:length(V(g)), node_weights)
    } else {
        node_signals <- igraph::vertex_attr(g, "signal")
        edge_signals <- igraph::edge_attr(g, "signal")
        edges <- cbind(edges, edge_signals)
        vertices <- cbind(1:length(V(g)), node_signals)
        signals <- data.frame(singal = names(signals), weights = signals)
        write_tbl(signals, signals_file)
    }
    write_tbl(vertices, nodes_file)
    write_tbl(edges, edges_file)
}

check_attr <- function(instance, attr, default, hint = NULL) {
    if (! attr %in% igraph::list.vertex.attributes(instance)) {
        warning(paste0("No `", attr, "` vertex attribute. Setting to ", default),
                call. = FALSE)
        instance <- igraph::set.vertex.attribute(instance, attr, default)
    }
    vertex_values <- igraph::get.vertex.attribute(instance, attr)
    if (any(is.na(vertex_values))) {
        warning(paste0("NA values provided for `", attr, "` vertex attribute. Setting to ",
                       ifelse(is.null(hint), default, hint)), call. = FALSE)
        vertex_values[is.na(vertex_values)] <- default
        instance <- igraph::set.vertex.attribute(instance, attr, value = vertex_values)
    }
    if (! attr %in% igraph::list.edge.attributes(instance)) {
        warning(paste0("No `", attr, "` edge attribute. Setting to ", default),
                call. = FALSE)
        instance <- igraph::set.edge.attribute(instance, attr, default)
    }
    edge_values <- get.vertex.attribute(instance, attr)
    if (any(is.na(edge_values))) {
        warning(paste0("NA values provided for `", attr, "` edge attribute. Setting to ",
                       ifelse(is.null(hint), default, hint)), call. = FALSE)
        edge_values[is.na(edge_values)] <- default
        instance <- igraph::set.vertex.attribute(instance, attr, value = edge_values)
    }
    instance
}

cli_args <- function(solver, nodes_file, edges_file, signals_file = NULL) {
    params <- c("-n", nodes_file, "-e", edges_file, "-m", solver$threads)
    if (!is.null(solver$timelimit)) {
        params <- c(params, "-t", solver$timelimit)
    }
    if (!is.null(signals_file)) {
        params <- c(params, "-s", signals_file)
    }
    params
}

run_solver <- function(solver, instance, sgmwcs, signals = NULL) {
    java_init(solver)

    nodes_file <- tempfile()
    edges_file <- tempfile()
    signals_file <- if (sgmwcs) tempfile() else NULL
    on.exit(file.remove(nodes_file, edges_file))
    if (sgmwcs) {
        on.exit(file.remove(signals_file))
    }
    write_files(instance, nodes_file, edges_file, signals_file, signals)
    args <- cli_args(solver, nodes_file, edges_file, signals_file)
    args <- rJava::.jarray(args)
    rJava::J(if (sgmwcs) sgmwcs_java_class else gmwcs_java_class)$main(args)

    nodes <- utils::read.table(paste0(nodes_file, ".out"), comment.char = "#")
    edges <- utils::read.table(paste0(edges_file, ".out"), comment.char = "#")

    if (!sgmwcs) {
        nodes <- nodes[nodes[, 2] != "n/a", ]
        edges <- edges[edges[, 3] != "n/a", ]
    }

    if (nrow(edges) == 0) {
        igraph::induced_subgraph(instance, as.integer(nodes[, 1]))
    } else {
        eids <- apply(edges, 1, function(x) get.edge.ids(instance, x[1:2]))
        igraph::subgraph.edges(instance, eids)
    }
}

#' @rdname solve_mwcsp
#' @param signals named vector of weights for signals
#' @export
solve_mwcsp.sgmwcs_solver <- function(solver, instance, signals = NULL, ...) {
    if (!inherits(solver, sgmwcs_class)) {
        stop("Not a sgmwcs solver")
    }

    original <- instance

    if (is.null(signals)) {
        instance <- check_attr(instance, "weight", 0)
        V(instance)$signal <- paste0("V", 1:igraph::vcount(instance))
        E(instance)$signal <- paste0("E", 1:igraph::ecount(instance))
        signals <- setNames(c(V(instance)$weight, E(instance)$weight),
                            c(V(instance)$signal, E(instance)$signal))
    }

    instance <- check_attr(instance, "signal", NA, hint = "a default signal with zero weight.")
    observed_signals <- setdiff(union(V(instance)$signal, E(instance)$signal), NA)
    non_assigned_signals <- setdiff(observed_signals, names(signals))
    if (length(non_assigned_signals) != 0) {
        warning(paste0("No weight assigned for the following signals: ",
                       paste0(non_assigned_signals, collapse = ", ")))
        V(instance)$signal[V(instance)$signal %in% non_assigned_signals] <- NA
        E(instance)$signal[E(instance)$signal %in% non_assigned_signals] <- NA
    }

    #renaming signals to S1..Sn
    map <- stats::setNames(paste0("S", 1:(length(signals) + 1)),
                           c(names(signals), "NA"))
    mapSignal <- function(x) {
        x <- sapply(x, function(x) ifelse(is.na(x), map[length(map)], map[x]))
    }
    V(instance)$signal <- mapSignal(V(instance)$signal)
    E(instance)$signal <- mapSignal(E(instance)$signal)
    names(signals) <- mapSignal(names(signals))
    if (any(is.na(V(instance)$signal)) || any(is.na(E(instance)$signal))) {
        signals <- c(signals, setNames(0, map[length(map)]))
    }

    V(instance)$id <- 1:length(V(instance))
    E(instance)$id <- 1:length(E(instance))
    mwcs <- run_solver(solver, instance, sgmwcs = TRUE, signals)
    weight <- sum(signals[unique(union(V(mwcs)$signal, E(mwcs)$signal))])
    if (length(E(mwcs)) > 0) {
        ret <- igraph::subgraph.edges(original, E(mwcs)$id)
    } else {
        ret <- igraph::induced_subgraph(original, V(mwcs)$id)
    }
    solution(ret, weight, FALSE)
}

#' @rdname solve_mwcsp
#' @export
solve_mwcsp.gmwcs_solver <- function(solver, instance, ...) {
    if (!inherits(solver, gmwcs_class)) {
        stop("Not a gmwcs solver")
    }
    instance <- check_attr(instance, "weight", 0)

    mwcs <- run_solver(solver, instance, sgmwcs = FALSE)

    weight <- sum(V(mwcs)$weight, E(mwcs)$weight)
    solution(mwcs, weight, FALSE)
}
