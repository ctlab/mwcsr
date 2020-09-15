library(igraph)
library(devtools)

load(system.file("inst/gatom_example.rda", package="mwcsr"))

example <- gatom_example

for (attr in setdiff(list.vertex.attributes(example), c("name", "label", "weight", "signal"))) {
    example <- remove.vertex.attribute(example, attr)
}

for (attr in setdiff(list.edge.attributes(example), c("name", "label", "weight", "signal"))) {
    example <- remove.edge.attribute(example, attr)
}

sgmwcs_example <- mwcsr::normalize_sgmwcs_instance(example)
sgmwcs_example <- remove.vertex.attribute(sgmwcs_example, "weight")
sgmwcs_example <- remove.edge.attribute(sgmwcs_example, "weight")

example <- remove.vertex.attribute(example, "signal")
example <- remove.edge.attribute(example, "signal")

mwcs_example <- remove.edge.attribute(example, "weight")
gmwcs_example <- example

use_data(mwcs_example, overwrite = TRUE)
use_data(gmwcs_example, overwrite = TRUE)
use_data(sgmwcs_example, overwrite = TRUE)

