//
// Created by Alexander Loboda on 10/1/20.
//

#include <chrono>

#ifndef SRC_INTERRUPTIONMONITOR_H
#define SRC_INTERRUPTIONMONITOR_H

namespace mwcsr {

class monitor {
    typedef std::chrono::steady_clock clock_t;
    clock_t::time_point last_check;
    std::chrono::milliseconds interval;
    std::function<void()> callable;
public:
    monitor(std::function<void()> callable, int millis_interval);
    void check();
};

}

#endif //SRC_INTERRUPTIONMONITOR_H
