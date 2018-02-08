#ifndef ANNEALING_INDEX_H
#define ANNEALING_INDEX_H

#include <cstddef>
#include <vector>
#include <random>
#include <functional>

namespace annealing {

    class Index {
        std::vector<size_t> map;
        std::vector<size_t> data;
        std::vector<bool> exists;
        std::function<uint_fast32_t> &re;
        size_t n;
    public:
        Index(size_t n, std::function<uint_fast32_t> &re);
        void add(size_t v);
        void remove(size_t v);
        bool contains(size_t v);
        size_t random();
        size_t operator()(size_t);
    };

}

#endif //ANNEALING_INDEX_H
