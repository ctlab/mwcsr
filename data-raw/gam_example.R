library(igraph)
library(devtools)
library(stringr)

regex_filter <- function(lines, pattern, cols){
    ms <- as.matrix(str_match(lines, pattern)[, cols])
    ms[!is.na(ms[, 1]), , drop = F]
}

readSTP <- function(filename){
    lines <- readLines(filename)
    n <- regex_filter(lines, "Nodes\\s*([\\d]+)", 2)
    edge_list <- regex_filter(lines, "E\\s*([\\d]+)\\s*([\\d]+)", 2:3)
    scores <- regex_filter(lines, "T\\s*([\\d]+)\\s*(.+)", 2:3)

    make_empty_graph(n = n, directed = F) %>%
        add_edges(t(edge_list)) %>%
        set.vertex.attribute("weight", scores[, 1], as.numeric(scores[, 2]))
}

dir <- "inst/extdata/"

download.file("http://dimacs11.zib.de/instances/MWCS-GAM.zip", "gam.zip")
unzip("gam.zip", junkpaths = T, exdir = dir)

instance <- list.files(dir)[1]
path <- paste0(dir, instance)

gam_example <- readSTP(path)

use_data(gam_example, overwrite = T)
