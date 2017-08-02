/*
 * Instance.h
 *
 *  Created on: Mar 29, 2015
 *      Author: markus
 */

#ifndef INSTANCE_INSTANCE_H_
#define INSTANCE_INSTANCE_H_

//#include <ogdf/basic/Graph.h>
#include <vector>

using namespace std;

class Instance {

    void readMWCS();
    void readPCSTP();
    void readDilkina();
    void readSTP();
    void readMWCSBudget();

    int calculateComponents();

    int preprocessing();

    int degreeOneTest();
    int uselessComponentsTest();
    void rebuildDatastructures();
    int degreeZeroTest();

    vector<bool> nodesToRemove;

  public:
    Instance();
    ~Instance();

    struct cut {
        double rhsConst = 0;
        vector<int> lhs;
        vector<int> rhs;
    };

    double transformInternalValue(double value) const;

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

    int nNodes;
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
    double runtimeCplex = -1;
    double bestBoundCplex = -999999999999;
    double incumbentObjCplex = 999999999999;
    double gapCplex = -1;
    double roottime = -1;
    double rootgap = -1;
    double rootbound = -1;
    int nBBNodes = -1;
    int cplexStatus = -1;
    int solSize = 0;
};

#endif /* INSTANCE_INSTANCE_H_ */
