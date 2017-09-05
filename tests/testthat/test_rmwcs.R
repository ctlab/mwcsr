context("rmwcs solver")

data("GAM")
set.seed(20170905L)

test_that("rmwcs solver doesn't crash on simple graphs", {
    solve <- rmwcs()
    size <- 10

    test_graph <- function(make) {
        scores <- runif(size) - 0.5;
        solve(make(size) %>% set.vertex.attribute("score", value = scores))
    }

    test_graph(make_ring)
    test_graph(make_star)
    test_graph(make_full_graph)
})

test_that("rmwcs solver works fine on GAM instances", {
    solve <- rmwcs()
    for (instance in GAM) {
        sol <- solve(instance)
        expect_true(is.connected(sol))
    }
})
