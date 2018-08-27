context("simulated annealing solver")

test_that("SA solver works on specific test", {
    solver <- annealing_solver(initial_temperature = 10)
    g <- make_ring(5) %>% set.vertex.attribute("weight", value = 1:-3)
    instance <- mwcs_instance(g, parse_vertex_weights = TRUE)
    edge_weights(instance) <- 1
    instance <- solve_mwcsp(solver, instance)
    expect(length(V(instance$solution)) != 0, "pu")
})

test_that("The SA solver does not crush on GAM instances", {
    solver <- annealing_solver(initial_temperature = 10)
    for (graph in GAM) {
        instance <- mwcs_instance(graph)
        instance <- solve_mwcsp(solver, instance)
        expect_gte(length(V(solution(instance))), 0)
    }
})

test_that("The SA solver gives good solution for a GAM instance", {
    solver <- annealing_solver(normalization = FALSE, schedule = "boltzmann",
                               initial_temperature = 2.0, final_temperature = 0.125)
    g <- GAM[[1]]
    V(g)$name <- as.character(1:length(V(g)))
    instance <- mwcs_instance(g)
    instance <- solve_mwcsp(solver, instance)
    weight <- sum(as.numeric(V(induced_subgraph(g, V(solution(instance))$name))$weight))
    expect_gte(weight, 500)
})
