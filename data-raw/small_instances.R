library(igraph)
library(devtools)


gmwcs_instance_edges <- data.frame(from = c(1, 1, 2, 3, 4, 1),
                                   to = c(2, 3, 3, 4, 5, 5),
                                   weight = c(4.0, 7.0, 5.0, 1.0, -2.0, -1.5))
gmwcs_small_instance <- igraph::graph_from_data_frame(gmwcs_instance_edges,
                                                directed = FALSE)
igraph::V(gmwcs_small_instance)$weight <- c(-3.0, -5.0, 0.0, 2.0, 1.0)

mwcs_small_instance <- igraph::remove.edge.attribute(gmwcs_small_instance, "weight")

sgmwcs_edges <- data.frame(from = c(1, 2, 2, 3, 4, 5, 1),
                           to = c(2, 3, 4, 4, 6, 6, 5),
                           signal = c("S6", "S2", "S7", "S12", "S8", "S9", "S10"))
sgmwcs_small_instance <- igraph::graph_from_data_frame(sgmwcs_edges,
                                                 directed = FALSE)
igraph::V(sgmwcs_small_instance)$signal <- c("S1", "S3", "S4", "S5", "S1", "S1")
sgmwcs_small_instance$signals <- stats::setNames(c(7.0, -20.0, 40.0, 15.0, 8.0, 3.0, -7.0, -10.0,
                             -2.0, -15.3, 1.0, -20),
                           paste0("S", 1:12))

use_data(mwcs_small_instance, overwrite = TRUE)
use_data(gmwcs_small_instance, overwrite = TRUE)
use_data(sgmwcs_small_instance, overwrite = TRUE)
