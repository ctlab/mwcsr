/*
 * ProgramOptions.h
 *
 * class for handling the program options
 *
 *  Created on: Mar 28, 2015
 *      Author: markus
 */

#ifndef PROGRAMOPTIONS_H_
#define PROGRAMOPTIONS_H_

#include <Rcpp.h>

struct Parameters {
    /// Input file (instance)
    int problem = 0;
    int budget = 0;
    int solver = 0;
    int maxIter = 1000;
    int betaIter = 5;
    bool integer = false;
    int separation = 0;
    int maxAge = 10;
    int outputlag = 1;
    bool startcons = true;
    bool pegging = true;
    int sepIterFreeze = 50;
    int sepIter = 10;
    int heurIter = 10;
    int subgradient = 0;
    double beta = 2.0;
    int timelimit = 1800;
    int setting = 0;

    explicit Parameters(Rcpp::List);
    ~Parameters() = default;
};

#endif /* PROGRAMOPTIONS_H_ */
