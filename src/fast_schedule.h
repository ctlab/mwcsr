#ifndef ANNEALING_FAST_SCHEDULE_H
#define ANNEALING_FAST_SCHEDULE_H

#include "cooling_schedule.h"

namespace annealing {

    class FastSchedule : public CoolingSchedule {
    public:
        FastSchedule(double t0, double min);
        double next() override;
    };

}

#endif //ANNEALING_FAST_SCHEDULE_H
