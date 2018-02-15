#' Generates a rmwcs solver with corresponding parameters
#' @param timelimit Timelimit in seconds
#' @param max_iterations Maximum number of subgradient iterations
#' @param beta_iterations Number of nonimproving iterations until beta is halfed
#' @param separation Separation: "strong" of "fast"
#' @param max_age Extending the life of non-violated inequalities
#' @param startcons Whether to add flow-conservation/degree cons at start
#' @param pegging Pegging (variable fixing)
#' @param sep_iterations After how many iterations we are separating
#' @param sep_iter_freeze After how many iterations we are checking added ineqs
#' @param heur_iterations After how many iterations we are doing heuristics
#' @param subgradient Subgradient: "classic", "average", "cft"
#' @param beta Beta for subgradient
#' @param verbose Whether to print solving progress
#' @export
#' @import igraph
rmwcs <- function(timelimit = 1800L,
                  max_iterations = 1000L,
                  beta_iterations = 50L,
                  separation = "strong",
                  max_age = 1000L,
                  startcons = TRUE,
                  pegging = TRUE,
                  sep_iterations= 1L,
                  sep_iter_freeze = 1L,
                  heur_iterations = 1L,
                  subgradient = "classic",
                  beta = 2.0,
                  verbose = FALSE) {

    args <- mget(names(formals()), sys.frame(sys.nframe()))
    fs <- formals()
    for(arg in names(args)){
        if(is.null(args[[arg]]) || is.na(args[[arg]])){
            stop(paste0("Argument '", arg, "' mustn't be null or NA"))
        }
        args[[arg]] <- methods::as(args[[arg]], class(fs[[arg]]))
    }

    sep_methods <- c("strong", "fast")
    subgradients <- c("classic", "average", "cft")

    match.arg(arg = separation, choices = sep_methods)
    match.arg(arg = subgradient, choices = subgradients)
    args$separation <- pmatch(args$separation, sep_methods)
    args$subgradient <- pmatch(args$subgradient, subgradients)

    function(g, max_cardinality, budget) {
        if(!is_igraph(g)){
            stop("Not a graph object")
        }
        if(is_directed(g)){
            stop("Not an undirected graph")
        }
        if(!("score" %in% list.vertex.attributes(g))){
            stop("Please provide vertex attribute 'score'")
        }
        if(!missing(max_cardinality) && !missing(budget)){
            stop("You cannot set both cardinallity and budget restrictions")
        }
        if(!missing(budget)){
            if(!("cost" %in% list.vertex.attributes(g))){
                warning("No costs provided. Setting to ones")
                g <- set.vertex.attribute(g, name = "cost", value = 1)
            }
        }
        scores <- as.numeric(V(g)$score)
        costs <- as.numeric(V(g)$cost)
        if(any(is.na(scores)) || (!missing(budget) && any(is.na(costs)))){
            stop("Invalid score or cost")
        }
        edge_list = as_edgelist(g, names = FALSE)
        edge_problem = FALSE
        if ("score" %in% list.edge.attributes(g)) {
            edge_scores <- as.numeric(E(g)$score)
            stopifnot(all(edge_scores <= 0))
            edge_list <- cbind(edge_list, edge_scores)
            edge_problem = TRUE
        }

        instance <- list(edgelist = edge_list, scores = scores)

        if(!missing(max_cardinality))
            instance$card = as.integer(max_cardinality)

        if(!missing(budget)){
            instance$budget = budget
            instance$costs = costs
        }

        vs <- rmwcs_solve(instance, args)
        if (length(vs) > 1 && edge_problem) {
            vs <- vs[vs > length(V(g))]
            vs <- vs - length(V(g))
            subgraph.edges(g, eids = vs)
        } else {
            induced.subgraph(g, vids = vs)
        }
    }
}

#' ctor for mwcs solver
#' @export
mwcs.solver <- function() {

}

#' ctor for gmwcs solver
#' @export
gmwcs.solver <- function() {

}

#' ctor for java solver
#' @export
java_solver <- function(cplex_jar,
                        cplex_bin,
                        solver_jar,
                        timelimit,
                        root,
                        threads,
                        unrooted_time,
                        rooted_time) {
    rJava::.jinit(classpath = c(cplex_jar, solver_jar), parameters = "-Xmx4m")

    function(g) {
        args <- rJava::.jarray(c("-h"))
        rJava::J("ru.ifmo.ctddev.gmwcs.Main")$main(rJava::.jarray("-h"))
    }
}
