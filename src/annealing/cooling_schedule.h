#ifndef ANNEALING_COOLINGSCHEDULE_H
#define ANNEALING_COOLINGSCHEDULE_H

namespace annealing {

    class CoolingSchedule {
        double t0;
        size_t k;
    public:
        virtual double temperature() = 0;
        virtual bool is_hot() = 0;
    };

}

#endif //ANNEALING_COOLINGSCHEDULE_H
