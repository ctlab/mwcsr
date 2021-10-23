//
// Created by Alexander Loboda on 10/1/20.
//

#include <Rcpp.h>

#include <utility>

#include "include/monitor.h"

namespace mwcsr {

monitor::monitor(std::function<void()> callable, int millis_interval)
                                    :interval(millis_interval), callable(std::move(callable)), stab(false) {
    last_check = clock_t::now();
}

void monitor::check() {
    if (!stab) {
        auto tp = clock_t::now();
        auto diff = tp - last_check;
        if (diff > interval) {
            last_check = tp;
            callable();
        }
    }
}

monitor::monitor() :interval(0), stab(true) {
    last_check = clock_t::now();
}

}