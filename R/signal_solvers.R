to_signal_instance <- function(instance) {
    signal_instance <- NULL
    inst_type <- get_instance_type(instance)

    if (inst_type$type == "SGMWCS" && inst_type$valid) {
        signal_instance <- instance
    } else if (inst_type$type == "GMWCS" && inst_type$valid) {
        signal_instance <- normalize_sgmwcs_instance(instance,
                                                     nodes.group.by = NULL,
                                                     edges.group.by = NULL)
    } else if (inst_type$type == "MWCS" && inst_type$valid) {
        inst2 <- instance
        E(inst2)$weight <- 0
        signal_instance <- normalize_sgmwcs_instance(inst2,
                                                     nodes.group.by = NULL,
                                                     edges.group.by = NULL)
    } else {
        stop("Not a valid MWCS, GMWCS or SGMWCS instance")
    }
    signal_instance
}


#' Helper function to convert an `igraph` object into a proper SGMWCS instance
#'
#' This function generates new `igraph` object with additional `signals` field added.
#' The way the instance is constructed is defined by the function parameters.
#' Nodes and edges are grouped separately, grouping columns are defined
#' by `nodes.group.by` and `edges.group.by` arguments. `group.only.positive` param specifies
#' whether to group only positive-weighted (specified by `nodes/edges.weight.column`) nodes and edges.
#'
#' @param g Graph to convert
#' @param nodes.weight.column Nodes column name (e.g. weight, score, w) for scoring
#' @param edges.weight.column Edges column name for scoring
#' @param nodes.group.by Nodes grouping column (e.g. signal, group, class)
#' @param edges.group.by Edges grouping column
#' @param group.only.positive Whether to group only positive-scored nodes/edges#'
#' @return An `igraph` object with proper attributes set.
#' @export
#' @examples
#' data("gatom_example")
#' normalize_sgmwcs_instance(gatom_example)
#'
#' @importFrom stats setNames
normalize_sgmwcs_instance <- function(g,
                                      nodes.weight.column="weight",
                                      edges.weight.column="weight",
                                      nodes.group.by="signal",
                                      edges.group.by="signal",
                                      group.only.positive=TRUE) {
    check_graph(g)
    if ("signals" %in% names(graph.attributes(g))) {
        warning("Input graph alredy looks like SGMWCS instance, doing nothing")
        return(g)
    }

    nt <- as_data_frame(g, what="vertices")
    if (!nodes.weight.column %in% colnames(nt)) {
        stop(sprintf("No %s node attribute present", nodes.weight.column))
    }

    if (!all(is.finite(nt[[nodes.weight.column]]))) {
        stop(sprintf("Not all node weights are finite numbers"))
    }

    if (!is.null(nodes.group.by)) {
        if (all(nodes.group.by %in% colnames(nt))) {
            signalParts <- nt[, nodes.group.by, drop=F]
            signalParts[is.na(signalParts)] <- paste0("__n", seq_len(sum(is.na(signalParts))))
            nt$signal <- do.call(paste, c(signalParts, sep="\r"))
            if (group.only.positive) {
                nt$signal <- paste(nt$signal, ifelse(nt[[nodes.weight.column]] > 0,
                                                     "",
                                                     seq_len(nrow(nt))), sep="\r")
            }
        } else {
            stop(sprintf("Can't collapse nodes, not all fields are present: %s",
                      paste0(setdiff(nodes.group.by, colnames(nt)), collapse=", ")))
        }

    } else {
        nt$signal <- paste0("sn_", seq_len(nrow(nt)))
    }
    nt <- nt[, c("signal", nodes.weight.column)]
    colnames(nt) <- c("signal", "weight")


    et <- as_data_frame(g, what="edges")
    if (!edges.weight.column %in% colnames(et)) {
        stop(sprintf("No %s edge attribute present", edges.weight.column))
    }

    if (!all(is.finite(et[[edges.weight.column]]))) {
        stop(sprintf("Not all edge weights are finite numbers"))
    }

    if (!is.null(edges.group.by)) {
        if (all(edges.group.by %in% colnames(et))) {
            signalParts <- et[, edges.group.by, drop=F]
            signalParts[is.na(signalParts)] <- paste0("__e", seq_len(sum(is.na(signalParts))))
            et$signal <- do.call(paste, c(signalParts, sep="\r"))

            if (group.only.positive) {
                et$signal <- paste(et$signal, ifelse(et[[edges.weight.column]] > 0,
                                                     "",
                                                     seq_len(nrow(et))), sep="\r")
            }
        } else {
            stop(sprintf("Can't collapse edges, not all fields present: %s",
                      paste0(setdiff(edges.group.by, colnames(et)), collapse=", ")))
        }

    } else {
        et$signal <- paste0("se_", seq_len(nrow(et)))
    }

    et <- et[, c("signal", edges.weight.column)]
    colnames(et) <- c("signal", "weight")

    st <- rbind(nt, et)

    st <- unique(st)
    if (any(duplicated(st$signal))) {
        stop(sprintf("Multiple weights are present for signal %s",
                     st$signal[which(duplicated(st$signal))[1]]))
    }

    # renaming signals
    old_signals <- st$signal
    st$signal <- paste0("s", seq_len(nrow(st)))
    old2new <- setNames(st$signal, old_signals)
    nt$signal <- old2new[nt$signal]
    et$signal <- old2new[et$signal]

    ret <- g
    ret$signals <- setNames(st$weight, st$signal)
    V(ret)$signal <- nt$signal
    E(ret)$signal <- et$signal
    ret
}
