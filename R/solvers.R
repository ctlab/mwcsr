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
                  beta_iterations = 5L,
                  separation = "strong",
                  max_age = 10L,
                  startcons = TRUE,
                  pegging = TRUE,
                  sep_iterations= 10L,
                  sep_iter_freeze = 50L,
                  heur_iterations = 10L,
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

    function(g, max_cardinallity, budget) {
        if(!is_igraph(g)){
            stop("Not a graph object")
        }
        if(!("score" %in% list.vertex.attributes(g))){
            stop("Please provide vertex attribute 'score'")
        }
        scores <- V(g)$score
        scores <- as.numeric(scores)
        if(any(is.na(scores))){
            stop("Invalid scores")
        }

        instance <- list(edgelist = as_edgelist(g), scores = scores)
        vs <- rmwcs_solve(instance, args)
        induced.subgraph(g, vids = vs)
    }
}
