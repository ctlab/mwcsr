rmwcs_class <- "rmwcs_solver"

sep_methods <- c("strong", "fast")
subgradients <- c("classic", "average", "cft")

check_rmwcs_solver <- function(solver) {
    if (!inherits(solver, rmwcs_class)) {
        stop("Function called with an invalid rmwcs solver instance")
    }
}

#' @export
parameters.rmwcs_solver <- function(solver) {
    params(parameter("timelimit", type = "integer", positive =  TRUE,
                   is_null_possible = TRUE),
         parameter("max_iterations", type = "integer", positive = TRUE,
                   is_null_possible = TRUE),
         parameter("max_age", type = "integer", positive = TRUE),
         parameter("beta_iterations", type = "integer", positive = TRUE),
         parameter("separation", type = "mc", mc = sep_methods),
         parameter("sep_iterations", type = "integer", positive = TRUE),
         parameter("start_constraints", type = "logical"),
         parameter("pegging", type = "logical"),
         parameter("sep_iter_freeze", type = "integer", positive = TRUE),
         parameter("heur_iterations", type = "integer", positive = TRUE),
         parameter("subgradient", type = "mc", mc = subgradients),
         parameter("beta", type = "float"),
         parameter("verbose", type = "logical"))
}

#' Generate a rmwcs solver
#'
#' The method is based on relax-and-cut approach and allows to solve
#' Maximum Weight Subgraph Probleam and its budget and cardinality variants.
#' By constructing lagrangian
#' relaxation of MWCS problem necessary graph connectivity constraints are introduced
#' in the objective function giving upper bound on the weight of the optimal
#' solution. On the other side, primal heuristic uses individul contribution
#' of the variables to lagrangian relaxation to find possible solution of the initial
#' problem. The relaxation is then optimized by using iterative subgradient method.
#'
#' One iteration of algorithm includes solving lagrangian relaxation and updating
#' lagrange multipliers. It may also contain cuts (or connectivity constraints) separation process, run of
#' heuristic method, variable fixing routine. The initial step size for
#' subgradient method can be passed as `beta` argument. If there is no improvement in
#' upper bound in consequtive `beta_iterations` iterations the step size is
#' halved. There are three possible strategies for updating multipliers. See the references
#' section for the article where differences are discussed.
#'
#' Some initial cuts are added at the start of the algorithm if `start_constraints`
#' is set to \code{TRUE}. Other constraints are separated on the fly and are
#' unaffected in the next `sep_iter_freeze` iterations of the subgradient mehod.
#' Then the corresponding lagrange mutipliers are updated from iteration to iteration.
#' Aging procedure for cuts is incorporated in the algorithm meaning constraint multipliers
#' are updated for non-violated cuts for up to `max_age` iterations from
#' the point where a cut was violated last time. There are two separation methods
#' implemented: fast and strong, where tha latter is supposed to minimize number of
#' variables used in generated constraint while in the former case there is no need to explore
#' whole graph to construct a constraint.
#'
#' A variant of MST approximation of PCSTP is used as Primal Heuristic.
#' See references for more details.
#'
#' The frequences
#' of running separation process and heuristic are specified in
#' `sep_iterations` and `heur_iterations`.
#'
#' @param timelimit Timelimit in seconds
#' @param max_iterations Maximum number of iterations
#' @param start_constraints Whether to add flow-conservation/degree constraints at start
#' @param separation Separation: "strong" or "fast"
#' @param sep_iterations Frequency of separating cuts (in iterations)
#' @param beta_iterations Number of nonimproving iterations until beta is halved
#' @param pegging variable fixing
#' @param max_age number of iterations in aging procedure for non-violated  cuts
#' @param sep_iter_freeze Number of iterations when a newly separated cut is anaffected by subgradient algorithm.
#' @param heur_iterations Frequency of calling heuristic method (in iterations)
#' @param subgradient Subgradient: "classic", "average", "cft"
#' @param beta Initial step size of subgradient algorithm
#' @param verbose Should the solving progress and stats be printed?
#' @return An object of class `mwcs_solver`.
#' @references Álvarez-Miranda E., Sinnl M. (2017)
#' "A Relax-and-Cut framework for large-scale maximum weight connected subgraph problems"
#' \doi{10.1016/j.cor.2017.05.015}
#' @export
#' @import igraph
rmwcs_solver <- function(timelimit = 1800L,
                  max_iterations = 1000L,
                  beta_iterations = 5L,
                  separation = "strong",
                  start_constraints = TRUE,
                  pegging = TRUE,
                  max_age = 10,
                  sep_iterations= 10L,
                  sep_iter_freeze = 50L,
                  heur_iterations = 10L,
                  subgradient = "classic",
                  beta = 2.0,
                  verbose = FALSE) {
    solver_ctor(c(rmwcs_class, mwcs_solver_class))
}

#' @rdname solve_mwcsp
#' @param max_cardinality integer maximum number of vertices in solution.
#' @param budget numeric maximum budget of solution.
#' @order 3
#' @export
solve_mwcsp.rmwcs_solver <- function(solver, instance, max_cardinality = NULL,
                                     budget = NULL, ...) {
    check_rmwcs_solver(solver)
    inst_type <- get_instance_type(instance)$type
    if (!(inst_type %in% c("MWCS", "Budget MWCS"))) {
        stop("Rmwcs solver supports only MWCS and Budget MWCS problem instances")
    }

    if (!is.null(max_cardinality) && !is.null(budget)) {
        stop("One of the arguments 'max_cardinality' and 'budget' must be NULL")
    }

    inst_type <- get_instance_type(instance)
    if (!inst_type$type %in% c("MWCS", "Budget MWCS") && !inst_type$valid) {
        stop("Instance is not a valid MWCS nor Budget MWCS instance")
    }

    instance <- igraph::simplify(instance)

    solver$separation <- pmatch(solver$separation, sep_methods) - 1
    solver$subgradient <- pmatch(solver$subgradient, subgradients) - 1

    instance_rep <- instance_from_graph(instance)
    instance_rep$vertex_weights <- V(instance)$weight

    if (!is.null(budget)) {
        instance_rep$budget_cost <- V(instance)$budget_cost
        budget <- as.numeric(budget)
        instance_rep$budget <- budget
    }

    if (!is.null(max_cardinality)) {
        max_cardinality <- as.integer(max_cardinality)
        stopifnot(max_cardinality > 0)
        instance_rep$cardinality <- max_cardinality
    }

    vs <- rmwcs_solve(instance_rep, solver)

    subgraph <- igraph::induced_subgraph(instance, vs$graph)
    weight <- sum(instance_rep$vertex_weights[vs$graph])
    opt <- isTRUE(all.equal(weight, vs$ub))
    solution(subgraph, weight, opt, upper_bound = vs$ub, incumbent = vs$lb,
             solved_to_optimality = abs(vs$lb - vs$ub) < EPS)
}
