#include "include/cooling_schedule.h"

namespace annealing {
    
    CoolingSchedule::CoolingSchedule(double t0, double min) :k(0), t0(t0), min(min) {
        current = std::numeric_limits<double>::infinity();
    }

    bool CoolingSchedule::is_hot() {
        return current >= min;
    }

    double CoolingSchedule::temperature() {
        return current;
    }

}
