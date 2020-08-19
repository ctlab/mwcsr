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
    for (graph in GAM) {
        solution <- solve_mwcsp(solver, graph)
        expect_gte(length(V(solution$graph)), 0)
    }
})

test_that("The SA solver gives good solution for a GAM instance", {
    solver <- annealing_solver(normalization = FALSE, schedule = "boltzmann",
                               initial_temperature = 2.0, final_temperature = 0.125)
    g <- GAM[[1]]
    solution <- solve_mwcsp(solver, g)
    expect_gt(solution$weight, 500)
    # expect_gt(solution$weight, 0)
})
