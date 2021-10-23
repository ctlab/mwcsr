//
// Created by Alexander Loboda on 9/3/21.
//

#include "include/active_pool.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace relax {

ActivePool::ActivePool(size_t n) :active(n), positions(n) {
    std::iota(active.begin(), active.end(), 0);
    std::iota(positions.begin(), positions.end(), 0);
}

void ActivePool::remove(size_t k) {
    if(active.empty()) {
        throw std::logic_error("Removing from empty list");
    }
    size_t pos = positions.at(k);
    size_t latest = active[active.size() - 1];
    positions[latest] = pos;
    active[pos] = latest;
    active.pop_back();
}

std::vector<size_t> ActivePool::all_active() {
    return active;
}
}