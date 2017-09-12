context("rmwcs solver")

data("GAM")
set.seed(42L)

test_that("rmwcs solver works on specific test", {
    solve <- rmwcs()
    g <- make_ring(5) %>% set.vertex.attribute("score", value = 1:-3)
    gp <- solve(g)
    expect_equal(sum(V(gp)$score), 1)
})

test_that("rmwcs solver doesn't crash on simple graphs", {
    solve <- rmwcs()
    size <- 10
    cardinality <- 3

    test_graph <- function(make) {
        scores <- runif(size) - 0.5;
        g <- make(size) %>% set.vertex.attribute("score", value = scores)
        solve(g)
        solve(g, max_cardinality = cardinality)
        g <- set.vertex.attribute(g, name = "cost", value = runif(size))
        solve(g, budget = cardinality)
    }

    test_graph(make_ring)
    test_graph(function(x) make_star(x, mode = "undirected"))
    test_graph(make_full_graph)
})

test_that("rmwcs solver builds connected solutions on GAM instances", {
    solve <- rmwcs()
    for (instance in GAM) {
        sol <- solve(instance)
        expect_true(is.connected(sol))
    }
})
