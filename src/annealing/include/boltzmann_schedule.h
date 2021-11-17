#ifndef ANNEALING_BOLTZMANN_ANNEALING_H
#define ANNEALING_BOLTZMANN_ANNEALING_H

#include "cooling_schedule.h"

namespace annealing {

    class BoltzmannSchedule : public CoolingSchedule {
    public:
        BoltzmannSchedule(double t0, double min);
        double next() override;
    };

}

#endif //ANNEALING_BOLTZMANN_ANNEALING_H
