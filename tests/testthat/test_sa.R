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
