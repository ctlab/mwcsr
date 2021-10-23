//
// Created by Alexander Loboda on 9/1/21.
//

#include "include/cut.h"
#include <iostream>

namespace relax {

bool Cuts::add(const Cut& cut) {
    if (exists(cut)) {
        return false;
    }
    cuts.add(cut);
    return true;
}

void Cuts::remove(size_t i) {
    cuts.remove(cuts[i]);
}

const Cut& Cuts::operator[](size_t i) const {
    return cuts[i];
}

bool Cuts::exists(const Cut& c) {
    return cuts.exists(c);
}

size_t Cuts::size() const {
    return cuts.size();
}

void Cuts::calculate_variable_weights() const {
    for (auto& cut: cuts) {
        cut.calculate_variable_weights();
    }
}

double Cuts::objective_part() {
    double ret = 0.0;
    for (auto& c: cuts) {
        ret += c.objective_part();
    }
    return ret;
}

void Cuts::normalize() {
    for (size_t i = 0; i < cuts.size(); i++) {
        const auto& c = cuts[i];
        auto vars = c.variables();
        if (std::any_of(vars.begin(), vars.end(), [](const Variable& v) {return v.fixed();})) {
            Cut cut(c);
            cuts.replace(i, cut);
        }
    }
}

const Cut& Cuts::get_const(size_t i) const {
    return cuts[i];
}

void Cut::swap(Cut& other) {
    left.swap(other.left);
    right.swap(other.right);
    std::swap(age, other.age);
    std::swap(subgradient, other.subgradient);
    std::swap(lambda, other.lambda);
}

bool operator==(const Cut& lhs, const Cut& rhs) {
    return lhs.left == rhs.left && lhs.right == rhs.right;
}

Cut::Cut(std::initializer_list<Variable> lhs, std::initializer_list<Variable> rhs) : left(lhs), right(rhs),
                                                                                     age(0), subgradient(0.), lambda(0.) {}

Cut::Cut(Cut&& other) noexcept :left(std::move(other.left)), right(std::move(other.right)),
                                age(other.age), subgradient(other.subgradient), lambda(other.lambda) {}

void Cut::calculate_variable_weights() const {
    if (lambda < EPS) {
        return;
    }

    left.calculate_variable_weights(-lambda);
    right.calculate_variable_weights(lambda);
}

double Cut::objective_part() const {
    return lambda * (right.get_const_part() - left.get_const_part());
}

void Cut::calculate_subgradient() {
    subgradient = right.get_const_part() - left.get_const_part();
    subgradient -= left.subgradient_part();
    subgradient += right.subgradient_part();
}

bool Cut::violatable() {
    return left.upper_bound() > right.lower_bound();
}

bool Cut::violated() const {
    return left.instant_value() > right.instant_value();
}

void Cut::free() {
    subgradient = 0.0;
}

void Cut::step(double theta) {
    lambda = std::max(0.0, lambda - theta * subgradient);
}

double Cut::subderivative() const {
    return subgradient;
}

std::vector<Variable> Cut::variables() const {
    auto ret = left.variables();
    auto other = right.variables();
    ret.insert(ret.end(), other.begin(), other.end());
    return ret;
}

VariableSum& Cut::lhs() {
    return left;
}

VariableSum& Cut::rhs() {
    return right;
}

bool Cut::try_fix() const {
    // TODO: include it in the code
    if(left.lower_bound() == right.upper_bound()) {
        for (Variable v: left.variables()) {
            if (!v.fixed()) {
                v.fix_value(0);
            }
        }
        for (Variable v: right.variables()) {
            if (!v.fixed()) {
                v.fix_value(1);
            }
        }
        return true;
    }
    return false;
}

Cut::Cut(const Cut& other) :left({}), right({}) {
    subgradient = other.subgradient;
    lambda = other.lambda;
    age = other.age;
    auto copy = [](VariableSum& curr, const VariableSum& other) {
        curr += other.get_const_part();
        for (Variable v : other.variables()) {
            if (v.fixed()) {
                curr += v.instant_value();
            } else {
                curr += {v};
            }
        }
    };
    copy(left, other.left);
    copy(right, other.right);
}

Cut& Cut::operator=(const Cut& other) {
    Cut c(other);
    swap(c);
    return *this;
}

const double Cut::mutliplier() const {
    return lambda;
}

unsigned Cut::non_violated_series() const {
    return age;
}

unsigned Cut::update_age() {
    if (violated()) {
        age = 0;
    } else {
        age++;
    }
    return age;
}

void swap(Cut& lhs, Cut& rhs) noexcept {
    lhs.swap(rhs);
}

std::ostream& operator<<(std::ostream& os, const Cut& c) {
    os << c.left << " ≤ " << c.right << "\t(λ = " << c.lambda << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Cuts& cuts) {
    os << "Cuts: " << std::endl;
    for (size_t i = 0; i < cuts.size(); i++) {
        auto& c = cuts.get_const(i);
        os << "\t" << c << std::endl;
    }
    return os;
}

double Cuts::check_previous(unsigned max_age) {
    double subgradient_norm = 0.0;
    for (size_t i = 0; i < cuts.size(); i++) {
        auto& c = cuts[i];
        c.calculate_subgradient();
        if (!c.violatable()) {
            c.free();
            remove(i);
            continue;
        }
        auto age = c.update_age();
        if (c.mutliplier() == 0.0 && age > max_age && !c.violated()) {
            c.free();
        }
        subgradient_norm += c.subderivative() * c.subderivative();
    }
    return subgradient_norm;
}

void Cuts::step(double theta) {
    for (size_t i = 0; i < cuts.size(); i++) {
        cuts[i].step(theta);
    }

}

void Cuts::try_fix() {
    for (auto& cut: cuts) {
        cut.try_fix();
    }
}


}