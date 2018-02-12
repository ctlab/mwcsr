#ifndef ANNEALING_INDEX_H
#define ANNEALING_INDEX_H

#include <cstddef>
#include <vector>
#include <random>

#include "definitions.h"

namespace annealing {

    class Index {
        std::vector<size_t> map;
        std::vector<size_t> data;
        std::vector<bool> exists;
        RandomEngine& re;
        size_t n;
    public:
        Index(size_t n, RandomEngine& re);
        void add(size_t v);
        void remove(size_t v);
        bool contains(size_t v);
        size_t random();
        size_t operator()(size_t);
        size_t size();
        std::vector<size_t> content();
    };

}

#endif //ANNEALING_INDEX_H
