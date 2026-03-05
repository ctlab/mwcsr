test_that("SA solver works on specific test", {
    solver <- annealing_solver(initial_temperature = 10)
    g <- igraph::make_ring(5)
    V(g)$weight <- 1:-3
    E(g)$weight <- 1
    solution <- solve_mwcsp(solver, g)
    expect_equal(solution$weight, 2)
})

test_that("The SA solver does not crush on GAM instances", {
    solver <- annealing_solver()
    solution <- solve_mwcsp(solver, gam_example)
    expect_gte(length(V(solution$graph)), 0)
})

test_that("The SA solver gives a good solution for a GAM instance", {
    solver <- annealing_solver(schedule = "boltzmann",
                               initial_temperature = 2.0, final_temperature = 0.125)
    solution <- solve_mwcsp(solver, gam_example)
    expect_gt(solution$weight, 200)
})

test_that("SA solver handles graphs with self-loops", {
    # Regression test: self-loop edges caused a crash because DynamicGraph
    # ignores self-loops (null token), making vertices appear isolated in the
    # dynamic connectivity structure even though a self-loop edge remained in
    # module_edges. This led to orphan edges and eventual assertion failure.
    solver <- annealing_solver()
    g <- igraph::make_ring(3)
    g <- igraph::add_edges(g, c(1L, 1L))  # self-loop on vertex 1
    V(g)$weight <- c(2, -1, -1)
    E(g)$weight <- c(1, 1, 1, 0.5)
    expect_no_error(solve_mwcsp(solver, g))
})
