//
// Created by Alexander Loboda on 9/3/21.
//

#ifndef MWCSR_ACTIVE_POOL_H
#define MWCSR_ACTIVE_POOL_H

#include <vector>
#include <cstdlib>

namespace relax {

class ActivePool {
    std::vector<size_t> positions;
    std::vector<size_t> active;
public:
    explicit ActivePool(size_t n);
    void remove(size_t k);
    std::vector<size_t> all_active();
};

}

#endif //MWCSR_ACTIVE_POOL_H
