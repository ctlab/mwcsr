//
// Created by Alexander Loboda on 9/1/21.
//

#include "include/variable.h"

#include <stdexcept>
#include <iostream>

namespace relax {

Variable::VariableRep::VariableRep(double weight, std::string name) :lb(0), ub(1), instant_value(0.0),
                                                                     obj_weight(weight), name(name) {}

VariableSum::VariableSum(const Variable& s) :summands({s}), load(0.0) { }

void VariableSum::swap(VariableSum& other) {
    summands.swap(other.summands);
    std::swap(load, other.load);
}

Variable::Variable(int id, double obj_weight, std::string name) :ptr(std::make_shared<VariableRep>(obj_weight, name)),
                                                                 id(id) {}

void Variable::fix_value(int val) {
    if (val < ptr->lb || val > ptr->ub) {
        throw std::range_error("Trying to fix variable with value that is out of possible range.");
    }
    ptr->lb = val;
    ptr->ub = val;
    ptr->instant_value = val;
}

VariableFactory::VariableFactory() :latest_id(0){}

Variable VariableFactory::take(double obj_weight, std::string name) {
    return Variable(latest_id++, obj_weight, name);
}

bool operator==(const Variable& a, const Variable& b) {
    return a.id == b.id;
}

void Variable::reset_weight() {
    ptr->weight = ptr->obj_weight;
}

void Variable::append_prize(double d) const {
    ptr->weight += d;
}

void Variable::setInstantValue() {
    if (ptr->weight > 0) {
        ptr->instant_value = ptr->ub;
    } else {
        ptr->instant_value = ptr->lb;
    }
}

double Variable::objective_part() {
    return ptr->instant_value * ptr->weight;
}

int Variable::instant_value() const {
    return ptr->instant_value;
}

double Variable::weight() {
    return ptr->weight;
}

double Variable::upper_bound() const {
    return ptr->ub;
}

double Variable::lower_bound() const {
    return ptr->lb;
}

bool Variable::fixed() const {
    return lower_bound() == upper_bound();
}

bool operator==(const VariableSum& lhs, const VariableSum& rhs) {
    return lhs.summands == rhs.summands && lhs.load == rhs.load;
}

VariableSum::VariableSum(std::initializer_list<Variable> lst) :load(0.0) {
    for (Variable v: lst) {
        if (v.fixed()) {
            load += v.instant_value();
        } else {
            summands.push_back(v);
        }
    }
}

VariableSum& VariableSum::operator+=(VariableSum other) {
    summands.insert(summands.end(), other.summands.begin(), other.summands.end());
    return *this;
}

VariableSum& VariableSum::operator+=(double x) {
    load += x;
    return *this;
}

void VariableSum::calculate_variable_weights(double sign) const {
    for (auto& v: summands) {
        v.append_prize(sign);
    }
}

double VariableSum::get_const_part() const {
    return load;
}

double VariableSum::subgradient_part() {
    double ret = 0.0;
    for (auto s: summands) {
        ret += s.instant_value();
    }
    return ret;
}

double VariableSum::upper_bound() const {
    double ret = load;
    for (Variable v: summands) {
        ret += v.upper_bound();
    }
    return ret;
}

double VariableSum::lower_bound() const {
    double ret = load;
    for (Variable v: summands) {
        ret += v.lower_bound();
    }
    return ret;
}

double VariableSum::instant_value() const {
    double ret = load;
    for (Variable v: summands) {
        ret += v.instant_value();
    }
    return ret;
}

std::vector<Variable> VariableSum::variables() const {
    return summands;
}

std::ostream& operator<<(std::ostream& os, const VariableSum& vs) {
    if (vs.summands.size() == 0) {
        os << 0;
        return os;
    }
    for (size_t i = 0; i < vs.summands.size(); i++) {
        if (i != 0) {
            os << " + ";
        }
        os << vs.summands.at(i);
    }
    if (vs.load != 0) {
        os << " + " + std::to_string(vs.load);
    }
    return os;
}

}
