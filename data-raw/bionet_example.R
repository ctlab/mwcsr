library(BioNet)
library(DLBCL)
data(dataLym)
data(interactome)
library(igraph)
library(devtools)
pvals <- cbind(t = dataLym$t.pval, s = dataLym$s.pval)
rownames(pvals) <- dataLym$label
pval <- aggrPvals(pvals, order = 2, plot = FALSE)
subnet <- subNetwork(dataLym$label, interactome)
subnet <- rmSelfLoops(subnet)
subnet
fb <- fitBumModel(pval, plot = FALSE)
scores <- scoreNodes(subnet, fb, fdr = 0.001)
bionet_example <- igraph.from.graphNEL(subnet, weight = FALSE)
V(bionet_example)$weight <- scores

use_data(bionet_example, overwrite=TRUE)

