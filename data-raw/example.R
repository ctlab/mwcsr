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
        set.vertex.attribute("score", scores[, 1], scores[, 2])
}

dir <- "inst/extdata/"

download.file("http://dimacs11.zib.de/instances/MWCS-GAM.zip", "gam.zip")
unzip("gam.zip", junkpaths = T, exdir = dir)

names <- list.files(dir)

# Take every 5th instance to reduce size
names <- names[seq(1, length(names), 5)]

suppressWarnings(
    GAM <- sapply(names, function (x) readSTP(paste0(dir, x)), simplify = F)
)

use_data(GAM, overwrite = T)
