#' @export
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
                  beta = 2.0) {

    args <- mget(names(formals()), sys.frame(sys.nframe()))
    print(args)
    fs <- formals()
    for(arg in names(args)){
        if(is.null(args[[arg]]) || is.na(args[[arg]])){
            stop(paste0("Argument '", arg, "' mustn't be null or NA"))
        }
        args[[arg]] <- as(args[[arg]], class(fs[[arg]]))
    }

    sep_methods <- c("strong", "fast")
    subgradients <- c("classic", "average", "cft")

    match.arg(arg = separation, choices = sep_methods)
    match.arg(arg = subgradient, choices = subgradients)
    args$sep_methods <- pmatch(args$separation, sep_methods)
    args$subgradinent <- pmatch(args$subgradient, subgradients)

    function(g, max_cardinallity, budget) {
        scores <- V(g)$scores
        stopifnot(length(scores) == length(V(g)))
        scores <- as.numeric(scores)
        instance <- list(edgelist = as_edgelist(g), scores = scores)
        rmwcs_solve(instance, args)
    }
}
