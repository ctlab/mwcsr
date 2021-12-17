scipjack_class = "scipjack_solver"

#' @export
parameters.scipjack_solver <- function(solver) {
    params(
        parameter("scipstp_bin", type="file", is_null_possible=FALSE),
        parameter("config_file", type="file", is_null_possible=TRUE))
}

write_stp <- function(g, stp_file) {
    edges <- igraph::as_edgelist(g, names = FALSE)
    nodes <- cbind("T", seq_along(V(g)), format(V(g)$weight, scientific=F, trim=T))
    n <- length(V(g))
    append <- function(s) {
        write(s, file=stp_file, append=TRUE)
    }
    write("33D32945 STP File, STP Format Version 1.0\n", file=stp_file,
         append=FALSE)

    append("SECTION Comments ")
    append("Problem \"Maximum Node Weight Connected Subgraph\"")
    append("END\n\n")

    append("SECTION Graph")
    append(paste0("Nodes ", n))
    append(paste0("Edges ", length(edges)/2))
    utils::write.table(cbind("E", edges), file=stp_file, quote=FALSE, sep = " ",
            col.names=FALSE, row.names=FALSE, append=TRUE)
    append("END\n\nSECTION Terminals")
    append(paste0("Terminals ", n))
    utils::write.table(nodes, file=stp_file, quote = FALSE, sep = " ",
            col.names=FALSE, row.names=FALSE, append=TRUE)
    append("END\n\nEOF")
}

append_output_file <- function(config_path, output_file) {
    write(paste0("stp/logfile=\"", output_file, "\""), file=config_path, append=TRUE)
}

find_output <- function(config_path) {
    lines <- readLines(config_path)
    idx <- grep('stp/logfile=*', lines)
    if (length(idx) == 0)
        path <- NULL
    else {
        path <- strsplit(lines[idx], '=')[[1]][2]
        path <- gsub("^\\s+|\\s+$", "", path)
    }
    path
}

run_scip <- function(solver, instance) {
    graph_dir <- tempfile("graph")
    dir.create(graph_dir, showWarnings=FALSE)
    default_output_file <- file.path(graph_dir, "out.stp")

    graph_file <- file.path(graph_dir, "in.stp")
    write_stp(instance, graph_file)

    config_path <- solver$config_file
    config_copy <- file.path(graph_dir, 'scip_config.s')
    if (!is.null(config_path)) {
        file.copy(config_path, config_copy)
    } else {
        file.create(config_copy)
    }

    output_file <- find_output(config_copy)
    if (is.null(output_file)) {
        append_output_file(config_copy, default_output_file)
        output_file <- default_output_file
    }

    system2(solver$scipstp_bin, c("-f", graph_file, "-s", config_copy), )

    lines <- readLines(output_file)
    nodes <- lines[grepl("^V \\d+", lines)]
    nodes <- as.integer(unlist(sub("V (\\d+)", '\\1', nodes)))

    induced_subgraph(instance, nodes)
}

run_scip_solver <- function(solver, instance) {
    V(instance)$index <- seq_len(vcount(instance))
    E(instance)$index <- seq_len(ecount(instance))
    mwcs <- run_scip(solver, instance)
    if (length(E(mwcs)) > 0) {
        ret <- igraph::subgraph.edges(instance, E(mwcs)$index)
    } else {
        ret <- igraph::induced_subgraph(instance, V(mwcs)$index)
    }
    solution(ret, get_weight(ret), TRUE, stats=NULL)
}

#' Construct a SCIP-jack solver
#'
#' This solver requires STP extension of \href{https://scipopt.org/#scipoptsuite}{SCIP-jack} solver.
#' To use this class you first need to download and build `SCIP-jack` and
#' `SCIPSTP` application.
#'
#' You can access solver directly using `run_scip` function. See example.
#' @param scipstp_bin path to `scipstp binary`.
#' @param config_file scipstp-formatted file. Parameters list is accessible
#' at  \href{https://www.scipopt.org/doc-6.0.2/html/PARAMETERS.php}{Official SCIP website}.
#'
#' @export
#' @references Rehfeldt D., Koch T. (2019)
#' "Combining NP-Hard Reduction Techniques and Strong Heuristics in an Exact Algorithm for the Maximum-Weight Connected Subgraph Problem."
#' \doi{10.1137/17M1145963}
#' @examples
#' \dontrun{
#' data("bionet_example")
#' scip <- scipjack_solver(scipstp_bin='/path/to/scipoptsuite/build/bin/applications/scipstp')
#' sol <- solve_mwcsp(scip, bionet_example)
#' }

scipjack_solver <- function(scipstp_bin,
                            config_file=NULL) {
    solver_ctor((c(scipjack_class, mwcs_solver_class)))
}

#' @rdname solve_mwcsp
#' @order 5
#' @export
solve_mwcsp.scipjack_solver <- function(solver, instance, ...) {
    inst_type <- get_instance_type(instance)
    if (inst_type$type == "MWCS" && inst_type$valid) {
        return(run_scip_solver(solver, instance))
    } else {
        stop("Not a valid MWCS instance")
    }
}
