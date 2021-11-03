//
// Created by Alexander Loboda on 9/1/21.
//

#ifndef MWCSR_CUT_H
#define MWCSR_CUT_H

#include <vector>
#include <functional>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <assert.h>

#include "variable.h"
#include "active_pool.h"

namespace relax {

class Cut;

}

template<>
struct std::hash<relax::Cut>;

namespace relax {

class Cut {
    unsigned age;
    double lambda;
    double subgradient;
    VariableSum left;
    VariableSum right;

public:
    Cut(std::initializer_list<Variable> lhs, std::initializer_list<Variable> rhs);
    Cut(const Cut& other);
    Cut& operator=(const Cut& other);

    Cut(Cut&& other) noexcept;
    VariableSum& lhs();
    VariableSum& rhs();
    friend typename ::std::hash<Cut>;
    void swap(Cut& other);
    friend bool operator==(const Cut& lhs, const Cut& rhs);

    void calculate_variable_weights() const;
    double objective_part() const;
    void calculate_subgradient();
    double subderivative() const;
    unsigned non_violated_series() const;
    unsigned update_age();

    std::vector<Variable> variables() const;

    bool violatable();
    bool violated() const;

    void free();
    void step(double d);

    bool try_fix() const;
    double mutliplier() const;

    friend void swap(relax::Cut& lhs, relax::Cut& rhs) noexcept;
    friend std::ostream& operator<<(std::ostream& os, const Cut&);
};

}

template<>
struct std::hash<relax::Cut> {
    std::size_t operator()(relax::Cut const& s) const noexcept {
        size_t h1 = std::hash<relax::VariableSum>{}(s.left);
        size_t h2 = std::hash<relax::VariableSum>{}(s.right);
        return (h1 * 0x1f1f1f1f) ^ h2;
    }

};

namespace relax {

template<typename T>
class Index {
    std::vector<T> list;
    std::unordered_map<T, size_t> positions;
public:
    Index() = default;

    std::vector<T> elements() const {
        return list;
    }

    decltype(list.cbegin()) begin() const {
        return list.cbegin();
    }

    decltype(list.cend()) end() const {
        return list.cend();
    }

    void remove(const T& obj) {
        auto iter = positions.find(obj);
        if (iter == positions.end()) {
            throw std::invalid_argument("Removing non-existing element from index.");
        }
        size_t pos = iter->second;
        positions.erase(iter);
        swap(list[pos], list[list.size() - 1]);
        if (pos != list.size() - 1) {
            positions[list[pos]] = pos;
        }
        list.pop_back();
    }

    void add(const T& obj) {
        list.push_back(obj);
        positions[obj] = list.size() - 1;
        assert(exists(obj));
    }

    T& operator[](size_t i) {
        return list.at(i);
    }

    const T& operator[](size_t i) const {
        return list.at(i);
    }

    size_t size() const {
        return list.size();
    }

    bool exists(const relax::Cut& cut) const {
        return positions.find(cut) != positions.end();
    }

    void replace(size_t i, const T& cut) {
        if (exists(cut)) {
            remove(cut);
        } else {
            T& c = list.at(i);
            assert(exists(c));
            positions.erase(c);
            list.at(i) = cut;
            positions[cut] = i;
        }
    }
};

}

namespace relax {

class Cuts {
    Index<Cut> cuts;
public:
    bool add(const Cut& cut);
    void remove(size_t i);
    const Cut& operator[](size_t i) const;
    bool exists(const Cut&);
    size_t size() const;

    void calculate_variable_weights() const;
    double objective_part();

    void normalize();
    const Cut& get_const(size_t i) const;

    friend std::ostream& operator<<(std::ostream& os, const Cuts&);

    double check_previous(unsigned int i);

    void step(double d);

    void try_fix();
};

}

#endif //MWCSR_CUT_H
