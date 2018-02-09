#ifndef ANNEALING_COOLINGSCHEDULE_H
#define ANNEALING_COOLINGSCHEDULE_H


class CoolingSchedule {
public:
    virtual double temperature() = 0;
};


#endif //ANNEALING_COOLINGSCHEDULE_H
