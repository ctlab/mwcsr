[![codecov](https://codecov.io/github/ctlab/mwcsr/branch/master/graphs/badge.svg)](https://app.codecov.io/gh/ctlab/mwcsr/branch/master)

# mwcsr

A package for solving maximum weight connected subgraph (MWCS) problem
and its variants.

Supported MWCS variants are:

-   classic (simple) MWCS, where only vertices are weighted;
-   budget MWCS, where vertices are parametrized by costs and overall
    budget is limited;
-   generalized MWCS (GMWCS), where both vertices and edges are
    weighted;
-   signal generalized MWCS (SGMWCS), where both vertices and edges are
    marked with weighted “signals”, and a weight of a subgraph is
    calculated as a sum of weights of its unique signals.

Currently, four solvers are supported:

-   heuristic relax-and-cut solver `rmwcs_solver` for MWCS and Budget
    MWCS;
-   heuristic relax-and-cut solver `rnc_solver` for MWCS/GMWCS/SGMWCS;
-   heuristic simulated annealing solver `annealing_solver` for
    MWCS/GMWCS/SGMWCS;
-   exact (if CPLEX library is available) or heuristic (without CPLEX)
    solver `virgo_solver` for MWCS/GMWCS/SGMWCS.

## Installation

The package can be installed from GitHub using `devtools`:

``` r
library(devtools)
install_github("ctlab/mwcsr")
```

## Quick start

Load `mwcsr`, as well as `igraph` package, which contains functions for
graph manipulations.

``` r
library(mwcsr)
library(igraph)
```

Let’s load an example instance of MWCS problem. The instance is a simple
`igraph` object with `weight` vertex attribute.

``` r
data("mwcs_example")
print(mwcs_example)
```

    ## IGRAPH f86c56f UN-- 194 209 -- 
    ## + attr: name (v/c), label (v/c), weight (v/n), label (e/c)
    ## + edges from f86c56f (vertex names):
    ##  [1] C00022_2--C00024_0  C00022_0--C00024_1  C00025_0--C00026_0 
    ##  [4] C00025_1--C00026_1  C00025_2--C00026_2  C00025_4--C00026_4 
    ##  [7] C00025_7--C00026_7  C00024_1--C00033_0  C00024_0--C00033_1 
    ## [10] C00022_0--C00041_0  C00022_1--C00041_1  C00022_2--C00041_2 
    ## [13] C00036_0--C00049_0  C00036_1--C00049_1  C00036_2--C00049_2 
    ## [16] C00036_4--C00049_4  C00037_1--C00065_0  C00037_0--C00065_1 
    ## [19] C00022_0--C00074_5  C00022_1--C00074_6  C00022_2--C00074_7 
    ## [22] C00024_0--C00083_0  C00024_1--C00083_1  C00026_1--C00091_0 
    ## + ... omitted several edges

``` r
summary(V(mwcs_example)$weight)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## -0.7379 -0.7379  1.9294  5.9667  7.2931 38.1546

Now let us initialize a heuristic relax-and-cut MWCS solver
(Alvarez-Miranda and Sinnl, 2017):

``` r
rcsolver <- rmwcs_solver()
```

Now we can use this solver to solve the example instance:

``` r
m <- solve_mwcsp(rcsolver, mwcs_example)
print(m$graph)
```

    ## IGRAPH 4fe8482 UN-- 162 164 -- 
    ## + attr: name (v/c), label (v/c), weight (v/n)
    ## + edges from 4fe8482 (vertex names):
    ##  [1] C00022_0--C00024_1  C00022_0--C00041_0  C00022_0--C00074_5 
    ##  [4] C00022_0--C00149_0  C00022_1--C00041_1  C00022_1--C00074_6 
    ##  [7] C00022_1--C00149_2  C00022_2--C00024_0  C00022_2--C00041_2 
    ## [10] C00022_2--C00074_7  C00022_2--C00149_1  C00024_0--C00033_1 
    ## [13] C00024_0--C00222_0  C00024_1--C00033_0  C00024_1--C00222_2 
    ## [16] C00025_0--C00026_0  C00025_1--C00026_1  C00025_2--C00026_2 
    ## [19] C00025_4--C00026_4  C00025_7--C00026_7  C00026_0--C00091_1 
    ## [22] C00026_0--C00311_1  C00026_1--C00091_0  C00026_1--C00311_0 
    ## + ... omitted several edges

``` r
print(m$weight)
```

    ## [1] 1178.432

## Using exact CPLEX-based Virgo solver

The `mwcsr` package also provide and interface to exact CPLEX-based
Virgo solver  
(<https://github.com/ctlab/virgo-solver>) which can be used to solve
MWCS, GMWCS and SGMWCS instances to provable optimality.

To setup this solver CPLEX libraries has to be available. CPLEX can be
downloaded from the official web-site:
<https://www.ibm.com/products/ilog-cplex-optimization-studio>. Free
licence can be obtained for academic purposes.

First, initialize `cplex_dir` variable to contain path to CPLEX
libraries (for example, /opt/ibm/ILOG/CPLEX_Studio129/).

``` r
cplex_dir <- '<replace with path to cplex folder>'
```

Then initialize the solver:

``` r
vsolver <- virgo_solver(cplex_dir=cplex_dir)
```

And run it on the same MWCS instance:

``` r
m <- solve_mwcsp(vsolver, mwcs_example)
print(m$graph)
```

    ## IGRAPH 9e2522d UN-- 164 166 -- 
    ## + attr: name (v/c), label (v/c), weight (v/n), label (e/c)
    ## + edges from 9e2522d (vertex names):
    ##  [1] C00022_2--C00024_0  C00022_0--C00024_1  C00025_0--C00026_0 
    ##  [4] C00025_1--C00026_1  C00025_2--C00026_2  C00025_4--C00026_4 
    ##  [7] C00025_7--C00026_7  C00024_1--C00033_0  C00024_0--C00033_1 
    ## [10] C00022_0--C00041_0  C00022_1--C00041_1  C00022_2--C00041_2 
    ## [13] C00036_0--C00049_0  C00036_1--C00049_1  C00036_2--C00049_2 
    ## [16] C00036_4--C00049_4  C00037_1--C00065_0  C00037_0--C00065_1 
    ## [19] C00022_0--C00074_5  C00022_1--C00074_6  C00022_2--C00074_7 
    ## [22] C00026_1--C00091_0  C00042_2--C00091_0  C00026_0--C00091_1 
    ## + ... omitted several edges

While the solution is a bit different its weight is the same as found
before. The solutions differs only in zero-weight vertices.

``` r
get_weight(m$graph)
```

    ## [1] 1178.432

However, Virgo guarantees that the result is optimal, unless the solver
was interrupted on time limit.

``` r
m$solved_to_optimality
```

    ## [1] TRUE

Next, consider a GMWCS instance which additionally has edge weights:

``` r
data("gmwcs_example")
gmwcs_example
```

    ## IGRAPH f86c56f UNW- 194 209 -- 
    ## + attr: name (v/c), label (v/c), weight (v/n), label (e/c), weight
    ## | (e/n)
    ## + edges from f86c56f (vertex names):
    ##  [1] C00022_2--C00024_0 C00022_0--C00024_1 C00025_0--C00026_0 C00025_1--C00026_1
    ##  [5] C00025_2--C00026_2 C00025_4--C00026_4 C00025_7--C00026_7 C00024_1--C00033_0
    ##  [9] C00024_0--C00033_1 C00022_0--C00041_0 C00022_1--C00041_1 C00022_2--C00041_2
    ## [13] C00036_0--C00049_0 C00036_1--C00049_1 C00036_2--C00049_2 C00036_4--C00049_4
    ## [17] C00037_1--C00065_0 C00037_0--C00065_1 C00022_0--C00074_5 C00022_1--C00074_6
    ## [21] C00022_2--C00074_7 C00024_0--C00083_0 C00024_1--C00083_1 C00026_1--C00091_0
    ## [25] C00042_2--C00091_0 C00026_0--C00091_1 C00042_1--C00091_1
    ## + ... omitted several edges

``` r
summary(E(gmwcs_example)$weight)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## -1.2715 -0.5762 -0.1060  0.5452  1.0458  8.5829

The same solver can be used to solve this instance:

``` r
m <- solve_mwcsp(vsolver, gmwcs_example)
print(m$graph)
```

    ## IGRAPH c12fb8d UNW- 182 181 -- 
    ## + attr: name (v/c), label (v/c), weight (v/n), label (e/c), weight
    ## | (e/n)
    ## + edges from c12fb8d (vertex names):
    ##  [1] C00022_2--C00024_0  C00022_0--C00024_1  C00025_0--C00026_0 
    ##  [4] C00025_1--C00026_1  C00025_2--C00026_2  C00025_4--C00026_4 
    ##  [7] C00025_7--C00026_7  C00024_1--C00033_0  C00024_0--C00033_1 
    ## [10] C00022_0--C00041_0  C00022_1--C00041_1  C00022_2--C00041_2 
    ## [13] C00036_0--C00049_0  C00036_1--C00049_1  C00036_2--C00049_2 
    ## [16] C00036_4--C00049_4  C00037_1--C00065_0  C00037_0--C00065_1 
    ## [19] C00022_0--C00074_5  C00022_1--C00074_6  C00022_2--C00074_7 
    ## + ... omitted several edges

``` r
get_weight(m$graph)
```

    ## [1] 1296.41

Finally, let consider an SGMWCS instance. The weights of nodes and edges
are defined not directly, but through the `signals` attribute:

``` r
data("sgmwcs_example")
sgmwcs_example
```

    ## IGRAPH f86c56f UN-- 194 209 -- 
    ## + attr: signals (g/n), name (v/c), label (v/c), signal (v/c), label
    ## | (e/c), signal (e/c)
    ## + edges from f86c56f (vertex names):
    ##  [1] C00022_2--C00024_0 C00022_0--C00024_1 C00025_0--C00026_0 C00025_1--C00026_1
    ##  [5] C00025_2--C00026_2 C00025_4--C00026_4 C00025_7--C00026_7 C00024_1--C00033_0
    ##  [9] C00024_0--C00033_1 C00022_0--C00041_0 C00022_1--C00041_1 C00022_2--C00041_2
    ## [13] C00036_0--C00049_0 C00036_1--C00049_1 C00036_2--C00049_2 C00036_4--C00049_4
    ## [17] C00037_1--C00065_0 C00037_0--C00065_1 C00022_0--C00074_5 C00022_1--C00074_6
    ## [21] C00022_2--C00074_7 C00024_0--C00083_0 C00024_1--C00083_1 C00026_1--C00091_0
    ## [25] C00042_2--C00091_0 C00026_0--C00091_1 C00042_1--C00091_1
    ## + ... omitted several edges

``` r
str(V(sgmwcs_example)$signal)
```

    ##  chr [1:194] "s1" "s1" "s1" "s2" "s3" "s4" "s4" "s4" "s4" "s4" "s5" "s5" ...

``` r
str(E(sgmwcs_example)$signal)
```

    ##  chr [1:209] "s103" "s104" "s105" "s106" "s107" "s108" "s109" "s110" "s110" ...

``` r
head(sgmwcs_example$signals)
```

    ##        s1        s2        s3        s4        s5        s6 
    ##  5.008879 -0.737898 -0.737898 20.112627 19.890279  2.069292

Again, we can solve this instance with Virgo solver:

``` r
m <- solve_mwcsp(vsolver, sgmwcs_example)
print(m$graph)
```

    ## IGRAPH f81ae9a UN-- 40 39 -- 
    ## + attr: signals (g/n), name (v/c), label (v/c), signal (v/c), index
    ## | (v/n), label (e/c), signal (e/c), index (e/n)
    ## + edges from f81ae9a (vertex names):
    ##  [1] C00022_0--C00024_1  C00025_1--C00026_1  C00024_1--C00033_0 
    ##  [4] C00022_0--C00041_0  C00036_1--C00049_1  C00037_1--C00065_0 
    ##  [7] C00022_0--C00074_5  C00026_1--C00091_0  C00042_2--C00091_0 
    ## [10] C00042_2--C00091_51 C00111_6--C00118_2  C00117_6--C00119_4 
    ## [13] C00042_2--C00122_1  C00022_0--C00149_0  C00122_1--C00149_0 
    ## [16] C00036_1--C00149_1  C00122_1--C00149_1  C00111_6--C00184_0 
    ## [19] C00117_6--C00199_1  C00118_2--C00231_1  C00199_1--C00231_1 
    ## + ... omitted several edges

``` r
get_weight(m$graph)
```

    ## [1] 270.7676

## Running Virgo heuristics without CPLEX

In case CPLEX is not available, Virgo solver can be run in the heuristic
mode. Just set `cplex_dir` parameter to `NULL`:

``` r
vhsolver <- virgo_solver(cplex_dir=NULL)
```

While the results are not optimal, sometimes they can be enough for
practical applications:

``` r
m <- solve_mwcsp(vhsolver, mwcs_example)
get_weight(m$graph)
```

    ## [1] 1174.737

``` r
m$solved_to_optimality
```

    ## [1] FALSE

``` r
m <- solve_mwcsp(vhsolver, gmwcs_example)
get_weight(m$graph)
```

    ## [1] 1276.772

``` r
m <- solve_mwcsp(vhsolver, sgmwcs_example)
get_weight(m$graph)
```

    ## [1] 227.3432
