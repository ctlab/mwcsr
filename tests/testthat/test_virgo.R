cplex_dir <- Sys.getenv("CPLEX_HOME")


gmwcs_instance_edges <- data.frame(from = c(1, 1, 2, 3, 4, 1),
                                   to = c(2, 3, 3, 4, 5, 5),
                                   weight = c(4.0, 7.0, 5.0, 1.0, -2.0, -1.5))
gmwcs_instance <- igraph::graph_from_data_frame(gmwcs_instance_edges,
                                                directed = FALSE)
igraph::V(gmwcs_instance)$weight <- c(-3.0, -5.0, 0.0, 2.0, 1.0)

mwcs_instance <- igraph::remove.edge.attribute(gmwcs_instance, "weight")

sgmwcs_edges <- data.frame(from = c(1, 2, 2, 3, 4, 5, 1),
                           to = c(2, 3, 4, 4, 6, 6, 5),
                           signal = c("S6", "S2", "S7", "S12", "S8", "S9", "S10"))
sgmwcs_instance <- igraph::graph_from_data_frame(sgmwcs_edges,
                                                 directed = FALSE)
igraph::V(sgmwcs_instance)$signal <- c("S1", "S3", "S4", "S5", "S1", "S1")
sgmwcs_instance$signals <- stats::setNames(c(7.0, -20.0, 40.0, 15.0, 8.0, 3.0, -7.0, -10.0,
                             -2.0, -15.3, 1.0, -20),
                           paste0("S", 1:12))

test_that("virgo_solver works on MWCS", {
    if (!file.exists(cplex_dir)) {
        skip("No CPLEX available")
    }
    solver <- virgo_solver(cplex_dir=cplex_dir)
    solution <- solve_mwcsp(solver, mwcs_instance)
    expect_equal(solution$weight, 3)
})


test_that("virgo_solver works on GMWCS", {
    if (!file.exists(cplex_dir)) {
        skip("No CPLEX available")
    }
    solver <- virgo_solver(cplex_dir=cplex_dir)
    solution <- solve_mwcsp(solver, gmwcs_instance)
    expect_equal(solution$weight, 11)
})

test_that("virgo_solver works on SGMWCS", {
    if (!file.exists(cplex_dir)) {
        skip("No CPLEX available")
    }
    solver <- virgo_solver(cplex_dir=cplex_dir)
    solution <- solve_mwcsp(solver, sgmwcs_instance)
    expect_equal(solution$weight, 51)
})


test_that("heuristic virgo_solver works on MWCS", {
    solver <- virgo_solver(cplex_dir=NULL)
    solution <- solve_mwcsp(solver, mwcs_instance)
    expect_equal(solution$weight, 3)
})


test_that("heuristic virgo_solver works on GMWCS", {
    solver <- virgo_solver(cplex_dir=NULL)
    solution <- solve_mwcsp(solver, gmwcs_instance)
    expect_gt(solution$weight, 0)
    expect_lte(solution$weight, 11)
})

test_that("heuristic virgo_solver works on SGMWCS", {
    solver <- virgo_solver(cplex_dir=NULL)
    solution <- solve_mwcsp(solver, sgmwcs_instance)
    expect_equal(solution$weight, 51)
})

test_that("virgo solver does not supported repeated negative signals", {
    sgmwcs_edges <- data.frame(from = c(1, 2, 2, 3, 4, 5, 1),
                               to = c(2, 3, 4, 4, 6, 6, 5),
                               signal = c("S6", "S2", "S7", "S2", "S8", "S9", "S10"))
    sgmwcs_instance <- igraph::graph_from_data_frame(sgmwcs_edges,
                                                     directed = FALSE)
    igraph::V(sgmwcs_instance)$signal <- c("S1", "S3", "S4", "S5", "S1", "S1")
    sgmwcs_instance$signals <- stats::setNames(c(7.0, -20.0, 40.0, 15.0, 8.0, 3.0, -7.0, -10.0,
                                                 -2.0, -15.3, 1.0),
                                               paste0("S", 1:11))

    solver <- virgo_solver(cplex_dir=NULL)
    expect_error(solution <- solve_mwcsp(solver, sgmwcs_instance))

})

test_that("heuristic virgo_solver works on SGMWCS", {
    solver <- virgo_solver(cplex_dir=NULL)
    si <- sgmwcs_instance
    si$signals <- c(si$signals, "neg"=-1)
    solution <- solve_mwcsp(solver, si)
    expect_true(!is.null(solution))
})

