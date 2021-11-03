#include "include/cooling_schedule.h"

namespace annealing {
    
    CoolingSchedule::CoolingSchedule(double t0, double min) :min(min), t0(t0), k(0) {
        current = std::numeric_limits<double>::infinity();
    }

    bool CoolingSchedule::is_hot() {
        return current >= min;
    }

    double CoolingSchedule::temperature() {
        return current;
    }

}
