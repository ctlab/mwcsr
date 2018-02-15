context("mwcs instances")

test_that("mwcs_instance dumb constructor works", {
    g <- igraph::make_ring(2)
    instance <- mwcs_instance(g, parse_vertex_weights = FALSE)
    expect_s3_class(instance, "mwcs_instance")
})

test_that("weights on mwcs assign correctly", {
    g <- igraph::make_star(3, mode = "undirected")
    V(g)$name <- letters[3:1]

    instance <- mwcs_instance(g, parse_vertex_weights = FALSE)
    vertex_weights(instance) <- setNames(5.0, "a")
    expect_equal(instance$vertex_weights, c(0.0, 0.0, 5.0))

    vertex_weights(instance) <- 1
    expect_equal(instance$vertex_weights, c(1.0, 1.0, 1.0))

    vertex_weights(instance) <- setNames(c(1.0, 2.0, 3.0), letters[1:3])
    expect_equal(instance$vertex_weights, c(3.0, 2.0, 1.0))
})

test_that("parsing of edge weights and budgets works", {
    g <- igraph::make_ring(4)
    edge_weights <- seq(1, 2, length.out = 4)
    vertex_weights <- seq(1, 0, length.out = 4)
    budgets <- seq(0, 1, length.out = 4)
    E(g)$weight <- edge_weights
    V(g)$weight <- vertex_weights
    V(g)$budget <- budgets
    instance <- mwcs_instance(g, parse_edge_weights = TRUE, parse_budgets = TRUE)

    expected <- mwcs_instance(g, parse_vertex_weights = FALSE)
    class(expected) <- c("budget_mwcs_instance", "gmwcs_instance", "mwcs_instance")
    expected$vertex_weights <- vertex_weights
    expected$edge_weights <- edge_weights
    expected$budgets <- budgets

    expect_setequal(class(instance), class(expected))
    expect_setequal(names(instance), names(expected))
    expect_equal(instance[names(expected)], expected[names(expected)])
})

test_that("error is throwing when trying to assign bad values", {
    g <- igraph::make_ring(4)

    V(g)$weight <- c(1, 2, 3, "")
    expect_error(mwcs_instance(g))

    V(g)$weight <- c(1, 2, 3, 4)
    instance <- mwcs_instance(g)

    expect_error(budgets(instance) <- setNames(5, "a"))
})

test_that("root setting works", {
    g <- make_ring(3)
    V(g)$name <- letters[1:3]

    instance <- mwcs_instance(g, parse_vertex_weights = FALSE)
    root(instance) <- 1L
    expect_is(instance, "rooted_mwcs_instance")

    expect_error(root(instance) <- 0L)
    expect_error(root(instance) <- NA)

    root(instance) <- "a"
    expect_equal(instance$root, 1)
})
