/*
 * SolverClassic.h
 *
 *  Created on: Aug 14, 2015
 *      Author: markus
 */

#ifndef SOLVERLAG_SOLVERCARDINALITY_H_
#define SOLVERLAG_SOLVERCARDINALITY_H_

#include "SolverLag.h"

class SolverCardinality : public SolverLag {
    double calculateCurrentSolution(bool save);
    bool primalHeuristic();
    int lagrangianPegging();
    int addInitCuts();

    double weightLast;
    double weightOutside;

public:
    SolverCardinality(Instance& instance, Parameters& params, mwcsr::monitor monitor);

    ~SolverCardinality();
};

#endif /* SOLVERLAG_SOLVERCLASSIC_H_ */
