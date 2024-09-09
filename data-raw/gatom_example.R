library(igraph)
library(devtools)

load("data/gatom_example.rda")

example <- gatom_example

for (attr in setdiff(list.vertex.attributes(example), c("name", "label", "weight", "signal"))) {
    example <- delete_vertex_attr(example, attr)
}

for (attr in setdiff(list.edge.attributes(example), c("name", "label", "weight", "signal"))) {
    example <- delete_edge_attr(example, attr)
}

sgmwcs_example <- mwcsr::normalize_sgmwcs_instance(example)
sgmwcs_example <- delete_vertex_attr(sgmwcs_example, "weight")
sgmwcs_example <- delete_edge_attr(sgmwcs_example, "weight")

example <- delete_vertex_attr(example, "signal")
example <- delete_edge_attr(example, "signal")

mwcs_example <- delete_edge_attr(example, "weight")
gmwcs_example <- example

use_data(mwcs_example, overwrite = TRUE)
use_data(gmwcs_example, overwrite = TRUE)
use_data(sgmwcs_example, overwrite = TRUE)

