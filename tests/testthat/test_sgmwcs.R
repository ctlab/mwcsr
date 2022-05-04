library(igraph)

test_that("normalize_sgmwcs_instance works", {
    et <- data.frame(from = c(1, 2, 2),
                     to = c(2, 3, 4),
                     signal = c("S5", "S6", "S7"),
                     weight = c(-1, -1, -1))
    nt <- data.frame(node=c(1, 2, 3, 4),
                     signal=c("S1", "S2", "S3", "S4"),
                     weight = c(-1, -1, -1, -1))
    g <- graph_from_data_frame(et, directed = FALSE, vertices = nt)
    instance <- normalize_sgmwcs_instance(g)
    expect_true("signals" %in% names(graph.attributes(instance)))
    expect_equal(length(instance$signals), 7)

    expect_warning(normalize_sgmwcs_instance(instance))

    # absent attributes
    expect_error(normalize_sgmwcs_instance(g, nodes.weight.column = "score"))
    expect_error(normalize_sgmwcs_instance(g, edges.weight.column = "score"))
    expect_error(normalize_sgmwcs_instance(g, nodes.group.by = "mz"))
    expect_error(normalize_sgmwcs_instance(g, edges.group.by = "mz"))

    et1 <- et
    et1$weight[1] <- Inf
    g1 <- graph_from_data_frame(et1, directed = FALSE, vertices = nt)
    expect_error(normalize_sgmwcs_instance(g1))

    nt1 <- nt
    nt1$weight[1] <- Inf
    g1 <- graph_from_data_frame(et, directed = FALSE, vertices = nt1)
    expect_error(normalize_sgmwcs_instance(g1))

    nt1 <- nt
    nt1$weight <- c(1, 2, 3, 4)
    nt1$signal <- c("S1", "S1", "S3", "S4")
    g1 <- graph_from_data_frame(et, directed = FALSE, vertices = nt1)
    expect_error(normalize_sgmwcs_instance(g1))

})

# METNET-T-57
test_that("normalize_sgmwcs_instance handles NA signals as unique", {
    et <- data.frame(from = c(1, 2, 2),
                     to = c(2, 3, 4),
                     signal = c(NA, NA, NA),
                     weight = c(-1, -2, 3))
    nt <- data.frame(node=c(1, 2, 3, 4),
                     signal=c(NA, NA, NA, "S4"),
                     weight = c(-10, 1, 4, -1))
    g <- graph_from_data_frame(et, directed = FALSE, vertices = nt)
    instance <- normalize_sgmwcs_instance(g)
    expect_true(get_instance_type(instance)$valid)
})

test_that("get_instance_type fails when signal names are repeated", {
    et <- data.frame(from = c(1, 2, 2),
                     to = c(2, 3, 4),
                     signal = c("S5", "S6", "S7"),
                     weight = c(-1, -1, -1))
    nt <- data.frame(node=c(1, 2, 3, 4),
                     signal=c("S1", "S2", "S3", "S4"),
                     weight = c(-1, -1, -1, -1))
    g <- graph_from_data_frame(et, directed = FALSE, vertices = nt)
    instance <- normalize_sgmwcs_instance(g)
    instance$signals <- c(instance$signals, instance$signals[1])
    expect_true(!get_instance_type(instance)$valid)
})
