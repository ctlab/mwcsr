#include <utility>
#include <stdexcept>

#include "include/index.h"

namespace annealing {

    Index::Index(size_t n) :n(n) {
        map.resize(n, 0);
        exists.resize(n, false);
    }

    void Index::add(size_t v) {
        map[v] = data.size();
        exists[v] = true;
        data.push_back(v);
    }

    void Index::remove(size_t v) {
        if (!exists[v]) {
            throw std::invalid_argument("removing non-existing element of an index");
        }
        size_t pos = map[v];
        size_t last = data.size() - 1;
        std::swap(data[pos], data[last]);
        exists[v] = false;
        map[data[pos]] = pos;
        data.pop_back();
    }

    bool Index::contains(size_t v) const {
        return v < n && exists[v];
    }

    size_t Index::random(RandomEngine& re) const {
        if (data.empty()) {
            throw std::logic_error("random element of an index when it is empty");
        }
        std::uniform_int_distribution<size_t> unif(0, data.size() - 1);
        return data[unif(re)];
    }

    size_t Index::operator()(size_t v) const {
        return map.at(v);
    }

    size_t Index::size() const {
        return data.size();
    }

    std::vector<size_t> Index::content() const {
        return data;
    }
}
