library(igraph)

test_that("one vertex best solution is returned", {
    solver <- rnc_solver(max_iterations = 100)
    g <- igraph::make_graph(c("A", "B"), directed = FALSE)
    V(g)$weight <- c(5, 5)
    E(g)$weight <- -7
    solution <- solve_mwcsp(solver, g)
    expect_equal(solution$weight, 5)
})

test_that("sgmwcs rnc solver works on specific test", {
    solver <- rnc_solver(max_iterations = 100)
    g <- igraph::make_ring(5)
    V(g)$weight <- 1:-3
    E(g)$weight <- 1
    solution <- solve_mwcsp(solver, g)
    expect_equal(solution$weight, 2)
    expect_equal(solution$solved_to_optimality, TRUE)
    expect_true(igraph::is_connected(solution$graph))
})

test_that("sgmwcs rnc solver does not crush on GAM instances", {
    solver <- rnc_solver(max_iterations = 100)
    solution <- solve_mwcsp(solver, gam_example)
    expect_gte(length(V(solution$graph)), 0)
    expect_true(igraph::is_connected(solution$graph))
})

test_that("sgmwcs rnc solver gives good solution for a GAM instance", {
    rnc <- rnc_solver(max_iterations = 100)
    solution <- solve_mwcsp(rnc, gmwcs_example)
    expect_gte(get_weight(solution$graph), 200)
    expect_true(igraph::is_connected(solution$graph))
})
