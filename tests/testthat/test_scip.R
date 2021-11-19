scip_bin <- Sys.getenv("SCIPSTP_BIN")

test_that("SCIP solver does not crash on GAM instances", {
    if (!file.exists(scip_bin)) {
        skip("SCIP is not available")
    }
    solver <- scipjack_solver(scip_bin)
    for (graph in GAM) {
        solution <- solve_mwcsp(solver, graph)
        expect_gte(length(V(solution$graph)), 0)
    }
})

test_that("SCIP solver supports config file", {
    if (!file.exists(scip_bin)) {
        skip("SCIP is not available")
    }
    scip_path <- system.file('ext', 'scip_config.s')
    solver <- scipjack_solver(scip_bin, scip_path)
    for (graph in GAM) {
        solution <- solve_mwcsp(solver, graph)
        expect_gte(length(V(solution$graph)), 0)
    }
})
