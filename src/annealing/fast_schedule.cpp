#include "include/fast_schedule.h"

namespace annealing {
    FastSchedule::FastSchedule(double t0, double min) : CoolingSchedule(t0, min) {}

    double FastSchedule::next() {
        current = t0 / ++k;
        return current;
    }
}
