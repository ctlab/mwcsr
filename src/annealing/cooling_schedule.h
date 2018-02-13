#ifndef ANNEALING_COOLINGSCHEDULE_H
#define ANNEALING_COOLINGSCHEDULE_H

#include <limits>

namespace annealing {

    class CoolingSchedule {
    protected:
        double min;
        double t0;
        double current;
        unsigned k;
    public:
        CoolingSchedule(double t0, double min);
        double temperature();
        virtual double next() = 0;
        bool is_hot();
    };

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

#endif //ANNEALING_COOLINGSCHEDULE_H
