context("rmwcs solver")

data("GAM")
set.seed(20170905L)

test_that("rmwcs solver works on specific test", {
    solve <- rmwcs(verbose = TRUE)
    g <- make_ring(5) %>% set.vertex.attribute("score", value = 1:-3)
    gp <- solve(g)
    expect_equal(sum(V(gp)$score), 1)
})

test_that("rmwcs solver doesn't crash on simple graphs", {
    solve <- rmwcs()
    size <- 10

    test_graph <- function(make) {
        scores <- runif(size) - 0.5;
        solve(make(size) %>% set.vertex.attribute("score", value = scores))
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
