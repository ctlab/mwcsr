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

}

#endif //ANNEALING_COOLINGSCHEDULE_H
