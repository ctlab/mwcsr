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
    s <- solve_mwcsp(solver, g)
    expect(s$solved_to_optimality)
})

test_that("rmwcs solver doesn't crash on simple graphs", {
    solver <- rmwcs()
    size <- 10
    card <- 3

    test_graph <- function(make) {
        scores <- runif(size) - 0.5;
        g <- make(size) %>% set.vertex.attribute("weight", value = scores)
        solve_mwcsp(solver, g)

        solution <- solve_mwcsp(solver, g, max_cardinality = card)
        expect_lte(length(V(solution$graph)), 3)

        g <- set.vertex.attribute(g, name = "budget_cost", value = runif(size))
        V(g)$weight <- 1
        solution <- solve_mwcsp(solver, g, budget = card)

        expect_lte(sum(V(solution$graph)$budget_cost), card)
    }

    test_graph(make_ring)
    test_graph(function(x) make_star(x, mode = "undirected"))
    test_graph(make_full_graph)
})

test_that("rmwcs solver builds connected solutions on GAM instances", {
    solver <- rmwcs()
    for (graph in GAM) {
        sol <- solve_mwcsp(solver, graph)
        expect(is.connected(sol$graph))
    }
})
