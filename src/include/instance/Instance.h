/*
 * Instance.h
 *
 *  Created on: Mar 29, 2015
 *      Author: markus
 */

#ifndef INSTANCE_INSTANCE_H_
#define INSTANCE_INSTANCE_H_

#include <Rcpp.h>
#include <vector>
#include "utility/Parameters.h"

using namespace std;

class Instance {
    int calculateComponents();

    int preprocessing();

    int degreeOneTest();
    int uselessComponentsTest();
    void rebuildDatastructures();
    int degreeZeroTest();

    vector<bool> nodesToRemove;

  public:
    Instance(List, List);
    ~Instance() = default;

    struct cut {
        double rhsConst = 0;
        vector<int> lhs;
        vector<int> rhs;
    };

    double transformInternalValue(double value) const;

    Parameters params;

    vector<double> myPrizes;
    vector<double> myBudgetCost;

    vector<bool> realTerminals;
    vector<int> myTerminals;
    vector<bool> trueTerminals;
    vector<int> myTrueTerminals;
    vector<vector<int>> adjList;

    vector<int> map;

    vector<int> componentArray;
    vector<vector<int>> components;
    vector<double> maxRevenueInComponent;
    vector<int> componentFixed;

    unsigned nNodes;
    int nEdges;
    int nTerminals;
    int nTrueNodes;
    int nTrueEdges;

    int nComponents;

    double maxPrize;
    double minWeight;
    double sumPrizes;

    double budget;

    int nRealTerminals;

    int nFlowNodes;
    int nFlowArcs;
    std::pair<int, int> *flowArcs;

    vector<int> fixedToZero;
    vector<int> fixedToOne;
    vector<cut> myCuts;

    bool incumbentFound = false;
    vector<bool> incumbent;

    int nFixedZero = 0;
    int nFixedOne = 0;
    int iterationsLag = -1;
    double runtimeLag = -1;
    double bestBoundLag = -999999999999;
    double incumbentObjLag = 999999999999;
    double gapLag = -1;
    int solSize = 0;

    void readInstance(List vector);
};

#endif /* INSTANCE_INSTANCE_H_ */
