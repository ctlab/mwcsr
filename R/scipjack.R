scipjack_class = "scipjack_solver"

#' @export
parameters.scipjack_solver <- function(solver) {
    params(
        parameter("scipstp_bin", type="file", is_null_possible=FALSE),
        parameter("output_dir", type="file", is_null_possible=FALSE))
}

write_stp <- function(g, stp_file) {
    edges <- igraph::as_edgelist(g, names = FALSE)
    nodes <- cbind("T", seq_along(V(g)), V(g)$weight)
    n <- length(V(g))
    append <- function(s) {
        write(s, file=stp_file, append=TRUE)
    }
    write("33D32945 STP File, STP Format Version 1.0\n", file=stp_file,
         append=FALSE)
    
    append("SECTION Comments ")
    append("Problem \"Maximum Node Weight Connected Subgraph\"")
    append("END\n")

    append("SECTION Graph")
    append(paste0("Nodes ", n))
    append(paste0("Edges ", length(edges)/2))
    utils::write.table(cbind("E", edges), file=stp_file, quote=FALSE, sep = " ",
            col.names=FALSE, row.names=FALSE, append=TRUE)
    append("END\nSECTION Terminals")
    append(paste0("Terminals ", n))
    utils::write.table(nodes, file=stp_file, quote = FALSE,
            col.names=FALSE, row.names=FALSE, sep = " ", append=TRUE)
    append("END\nEOF")
}

gen_config_file <- function(config_path, output_file) {
    write(paste0("stp/logfile=\"", output_file, "\""), file=config_path, append=FALSE)
}

run_scip <- function(solver, instance) {
    graph_dir <- tempfile("graph")
    dir.create(graph_dir, showWarnings=FALSE)

    graph_file <- file.path(graph_dir, "in.stp")
    write_stp(instance, graph_file)

    config_path <- file.path(graph_dir, "scipconfig.s")
    output_file <- file.path(graph_dir, "out.stp")
    gen_config_file(config_path, output_file)

    system2(solver$scipstp_bin, c("-f", graph_file, "-s", config_path))

    lines <- readLines(output_file)
    nodes <- lines[grepl("^V \\d+", lines)]
    # edges <- lines[grepl("^E \\d+ \\d+", lines)]
    # edges <- as.integer(unlist(strsplit(sub("E (\\d+) (\\d+)", '\\1 \\2', edges), ' ')))
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

#' @export
scipjack_solver <- function(scipstp_bin,
                            output_dir = NULL) {
    output_dir <- tempfile("graph")
    dir.create(output_dir, showWarnings=FALSE)
    solver_ctor((c(scipjack_class, mwcs_solver_class)))
}

#' @rdname solve_mwcsp
#' @export
#'
solve_mwcsp.scipjack_solver <- function(solver, instance, ...) {
    inst_type <- get_instance_type(instance)
    if (inst_type$type == "MWCS" && inst_type$valid) {
        return(run_scip_solver(solver, instance))
    } else {
        stop("Not a valid MWCS instance")
    }
}