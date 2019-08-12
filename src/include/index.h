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
        size_t n;
    public:
        explicit Index(size_t n);
        void add(size_t v);
        void remove(size_t v);
        bool contains(size_t v) const;
        size_t random(RandomEngine& re) const;
        size_t operator()(size_t) const;
        size_t size() const;
        std::vector<size_t> content() const;
    };

}

#endif //ANNEALING_INDEX_H
