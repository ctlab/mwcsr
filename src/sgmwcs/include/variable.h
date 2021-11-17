//
// Created by Alexander Loboda on 9/1/21.
//

#ifndef MWCSR_VARIABLE_H
#define MWCSR_VARIABLE_H

#include <vector>
#include <memory>
#include <algorithm>
#include <string>
#include <ostream>
#include <functional>

namespace relax {

constexpr double EPS = 1e-6;

class Variable {
    struct VariableRep {
        int lb;
        int ub;
        int instant_value;
        double weight;
        double obj_weight;
        std::string name;

        VariableRep(double weight, std::string name);
    };

    std::shared_ptr<VariableRep> ptr;
    size_t id;

    Variable(int id, double obj_weight, std::string name);

    friend class VariableFactory;

    friend ::std::hash<Variable>;
public:
    Variable(const Variable&) = default;

    void fix_value(int val);

    friend bool operator==(const Variable& a, const Variable& b);

    void reset_weight();

    int instant_value() const;

    void setInstantValue();

    void append_prize(double d) const;

    double objective_part();

    double weight();

    double upper_bound() const;

    double lower_bound() const;

    bool fixed() const;

    friend std::ostream& operator<<(std::ostream& os, const Variable& v) {
        os << v.ptr->name;
        return os;
    }
};

class VariableFactory {
    int latest_id;
public:
    VariableFactory();

    VariableFactory(const VariableFactory&) = default;

    Variable take(double obj_weight = 0.0, std::string name = "");
};

class VariableSum {
    std::vector<Variable> summands;
    double load = 0.0;
public:
    explicit VariableSum(const Variable& variable);

    VariableSum(const VariableSum&) = default;

    VariableSum(VariableSum&& other) noexcept = default;

    VariableSum(std::initializer_list<Variable> lst);

    VariableSum& operator=(VariableSum&& other) noexcept = default;

    void swap(VariableSum& other);

    double get_const_part() const;

    VariableSum& operator+=(VariableSum);

    VariableSum& operator+=(double);

    void calculate_variable_weights(double sign) const;

    double subgradient_part();

    std::vector<Variable> variables() const;

    friend ::std::hash<VariableSum>;

    friend bool operator==(const VariableSum& a, const VariableSum& b);

    friend std::ostream& operator<<(std::ostream& os, const VariableSum&);

    double upper_bound() const;

    double lower_bound() const;

    double instant_value() const;
};

}

namespace std {
template<>
struct hash<relax::Variable> {
    std::size_t operator()(relax::Variable const& s) {
        return std::hash<size_t>{}(s.id);
    }
};

template<>
struct hash<relax::VariableSum> {
    std::size_t operator()(relax::VariableSum const& s) {
        size_t seed = 0;
        std::hash<relax::Variable> hasher;
        std::vector<size_t> hashes;
        std::transform(s.summands.cbegin(), s.summands.cend(), std::back_inserter(hashes),
                       [&hasher](relax::Variable const& v) {
                           return hasher(v);
                       });
        hashes.push_back(std::hash<double>{}(s.load));
        for (auto h: hashes) {
            seed ^= h + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

}

#endif //MWCSR_VARIABLE_H
