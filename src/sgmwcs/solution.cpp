//
// Created by Alexander Loboda on 9/15/21.
//

#include "include/solution.h"

relax::Solution::Solution() :obj(0.0){}

relax::Solution::Solution(std::vector<size_t> elements, double obj) :edges(elements), obj(obj){}

void relax::Solution::add_edge(size_t e) {
    edges.push_back(e);
}

void relax::Solution::set_obj(double value) {
    obj = value;
}

double relax::Solution::objective() const {
    return obj;
}

std::vector<size_t> relax::Solution::solution() const {
    return edges;
}

void relax::Solution::add_obj(double value) {
    obj += value;
}
