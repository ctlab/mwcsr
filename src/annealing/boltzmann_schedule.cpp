#include <cmath>
#include <iostream>

#include "include/boltzmann_schedule.h"

namespace annealing {

    BoltzmannSchedule::BoltzmannSchedule(double t0, double min) : CoolingSchedule(t0, min) {}

    double BoltzmannSchedule::next() {
        current = t0 / (1.0 + (std::log(++k)));
        return current;
    }
}
