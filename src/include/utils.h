#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include <Rcpp.h>

#include "graph.h"

using namespace Rcpp;

namespace mwcsr {

Graph read_graph(List& instance);

}

#endif //SRC_UTILS_H
