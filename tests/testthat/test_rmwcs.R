context("rmwcs solver")

data("GAM")
set.seed(42L)

test_that("rmwcs solver can be constructed with specific parameters", {
    solver <- rmwcs()
    expect_equal(solver$separation, "strong")
    timelimit(solver) <- 20
    expect_equal(solver$timelimit, 20)
})

test_that("rmwcs solver works on specific test", {
    solver <- rmwcs()
    g <- make_ring(5) %>% set.vertex.attribute("weight", value = 1:-3)
    instance <- mwcs_instance(g, parse_vertex_weights = TRUE)
    instance <- solve_mwcsp(solver, instance)
    expect(instance$solved_to_optimality)
})

test_that("rmwcs solver doesn't crash on simple graphs", {
    solver <- rmwcs()
    size <- 10
    cardinality <- 3

    test_graph <- function(make) {
        scores <- runif(size) - 0.5;
        g <- make(size) %>% set.vertex.attribute("weight", value = scores)
        instance <- mwcs_instance(g)
        solve_mwcsp(solver, instance)

        max_cardinality(instance) <- cardinality
        instance <- solve_mwcsp(solver, instance)
        expect_lte(length(V(solution(instance))), 3)

        g <- set.vertex.attribute(g, name = "budget_cost", value = runif(size))
        instance <- mwcs_instance(g, parse_vertex_weights = FALSE,
                                  parse_budgets = TRUE)
        vertex_weights(instance) <- 1
        budget(instance) <- cardinality
        s <- solve_mwcsp(solver, instance)

        expect_lt(length(V(solution(s))), cardinality + 1)
    }

    test_graph(make_ring)
    test_graph(function(x) make_star(x, mode = "undirected"))
    test_graph(make_full_graph)
})

test_that("rmwcs solver builds connected solutions on GAM instances", {
    solver <- rmwcs()
    for (graph in GAM) {
        instance <- mwcs_instance(graph)
        sol <- solve_mwcsp(solver, instance)
        expect(is.connected(solution(sol)))
    }
})
