/*
 * SolverClassic.h
 *
 *  Created on: Aug 14, 2015
 *      Author: markus
 */

#ifndef SOLVERLAG_SOLVERCLASSIC_H_
#define SOLVERLAG_SOLVERCLASSIC_H_

#include "SolverLag.h"

class SolverClassic : public SolverLag {
    double calculateCurrentSolutionPegging();
    double calculateCurrentSolution(bool save);
    bool primalHeuristic();
    int lagrangianPegging();
    int addInitCuts();

public:
    SolverClassic(Instance& instance, Parameters& params, mwcsr::monitor monitor);
    // SolverClassic(const SolverClassic&);

    ~SolverClassic();
};

#endif /* SOLVERLAG_SOLVERCLASSIC_H_ */
