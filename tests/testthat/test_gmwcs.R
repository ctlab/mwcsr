data("GAM")
set.seed(42L)

cplex_lib <- "./cplex.so"
cplex_jar <- "./cplex.jar"

gmwcs_instance_edges <- data.frame(from = c(1, 1, 2, 3, 4, 1),
                                   to = c(2, 3, 3, 4, 5, 5),
                                   weight = c(4.0, 7.0, 5.0, 1.0, -2.0, -1.5))
gmwcs_instance <- igraph::graph_from_data_frame(gmwcs_instance_edges,
                                                directed = FALSE)
igraph::V(gmwcs_instance)$weight <- c(-3.0, -5.0, 0.0, 2.0, 1.0)

sgmwcs_edges <- data.frame(from = c(1, 2, 2, 3, 4, 5, 1),
                           to = c(2, 3, 4, 4, 6, 6, 5),
                           signal = c("S6", "S2", "S7", "S2", "S8", "S9", "S10"))
sgmwcs_instance <- igraph::graph_from_data_frame(sgmwcs_edges,
                                                 directed = FALSE)
igraph::V(sgmwcs_instance)$signal <- c("S1", "S3", "S4", "S5", "S1", "S1")
signals <- stats::setNames(c(7.0, -20.0, 40.0, 15.0, 8.0, 3.0, -7.0, -10.0,
                             -2.0, -15.3, 1.0),
                           paste0("S", 1:11))

test_that("gmwcs_solver works", {
    if (!file.exists(cplex_jar) || !file.exists(cplex_lib)) {
        skip("No CPLEX available")
    }
    solver <- gmwcs_solver(cplex_lib, cplex_jar)
    solution <- solve_mwcsp(solver, gmwcs_instance)
    expect_equal(solution$weight, 11)
})

test_that("sgmwcs_solver works", {
    if (!file.exists(cplex_jar) || !file.exists(cplex_lib)) {
        skip("No CPLEX available")
    }
    solver <- sgmwcs_solver(cplex_lib, cplex_jar)
    solution <- solve_mwcsp(solver, gmwcs_instance)
    expect_equal(solution$weight, 11)
})
