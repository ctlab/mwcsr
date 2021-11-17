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
#include "Parameters.h"

class Instance {
    int calculateComponents();

    int preprocessing();

    int degreeOneTest();
    int uselessComponentsTest();
    void rebuildDatastructures();
    int degreeZeroTest();

    std::vector<bool> nodesToRemove;
    void readInstance(Rcpp::List& instance);
    void readEdges(Rcpp::IntegerMatrix&);
    void addEdge(unsigned, unsigned);
    void initStructures(unsigned);

public:
    explicit Instance(Rcpp::List&);
    ~Instance() = default;

    struct cut {
        double rhsConst = 0;
        std::vector<int> lhs;
        std::vector<int> rhs;
    };

    double transformInternalValue(double value) const;

    std::vector<double> myPrizes;
    std::vector<double> myBudgetCost;

    std::vector<bool> realTerminals;
    std::vector<int> myTerminals;
    std::vector<bool> trueTerminals;
    std::vector<int> myTrueTerminals;
    std::vector<std::vector<int>> adjList;

    std::vector<int> map;

    std::vector<int> componentArray;
    std::vector<std::vector<int>> components;
    std::vector<double> maxRevenueInComponent;
    std::vector<int> componentFixed;

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
    int cardCons;

    int nRealTerminals;

    std::vector<int> fixedToZero;
    std::vector<int> fixedToOne;
    std::vector<cut> myCuts;

    bool incumbentFound = false;
    std::vector<bool> incumbent;

    int nFixedZero = 0;
    int nFixedOne = 0;
    int iterationsLag = -1;
    double runtimeLag = -1;
    double bestBoundLag = std::numeric_limits<double>::infinity();
    double incumbentObjLag = 0;
    double gapLag = -1;
    int solSize = 0;

    void findSimpleSolution();
};

#endif /* INSTANCE_INSTANCE_H_ */
