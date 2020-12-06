scip_bin <- Sys.getenv("SCIPSTP_BIN")

test_that("SCIP solver does not crash on GAM instances", {
    solver <- scip_solver(scip_bin)
    for (graph in GAM) {
        solution <- solve_mwcsp(solver, graph)
        expect_gte(length(V(solution$graph)), 0)
    }
})

test_that("SCIP solver supports config file", {
    scip_path <- system.file('ext', 'scip_config.s')
    solver <- scip_solver(scip_bin, scip_path)
    for (graph in GAM) {
        solution <- solve_mwcsp(solver, graph)
        expect_gte(length(V(solution$graph)), 0)
    }
})
