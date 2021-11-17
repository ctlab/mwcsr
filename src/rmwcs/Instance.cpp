/*
 * Instance.cpp
 *
 *  Created on: Mar 29, 2015
 *      Author: markus
 */

#include "include/Instance.h"

#include <queue>

using Rcpp::NumericVector;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::IntegerMatrix;
using Rcpp::as;

using std::vector;
using std::queue;

Instance::Instance(List& network)
        : components{vector<vector<int>>()},
          maxRevenueInComponent{vector<double>()}, nComponents{0}, maxPrize{0},
          minWeight{std::numeric_limits<double>::max()}, sumPrizes{0},
          budget(std::numeric_limits<double>::infinity()), nRealTerminals{0} {

    readInstance(network);

    nTrueNodes = nNodes;
    nTrueEdges = nEdges;

    nComponents = calculateComponents();

    findSimpleSolution();
    preprocessing();
    rebuildDatastructures();
    nComponents = calculateComponents();

    fixedToOne = vector<int>(nNodes, 0);
    fixedToZero = vector<int>(nNodes, 0);
}

void Instance::rebuildDatastructures() {

    vector<double> myPrizes2;
    vector<double> myBudgetCost2;
    vector<bool> realTerminals2;
    vector<int> myTerminals2;
    vector<vector<int>> adjList2;

    vector<int> backMap = vector<int>(nNodes);

    int newNNodes = 0;

    for (unsigned i = 0; i < nNodes; ++i) {
        if (!nodesToRemove[i]) {
            backMap[i] = newNNodes;
            map.push_back(i);
            adjList2.emplace_back();
            myPrizes2.push_back(myPrizes[i]);
            myBudgetCost2.push_back(myBudgetCost[i]);
            realTerminals2.push_back(realTerminals[i]);
            if (realTerminals[i])
                myTerminals2.push_back(newNNodes);
            newNNodes++;
        }
    }

    int newNEdges = 0;
    for (unsigned i = 0; i < nNodes; ++i) {
        if (!nodesToRemove[i]) {
            for (unsigned j = 0; j < adjList[i].size(); ++j) {
                if (!nodesToRemove[adjList[i][j]]) {
                    adjList2[backMap[i]].push_back(backMap[adjList[i][j]]);
                    newNEdges++;
                }
            }
        }
    }

    newNEdges /= 2;

    myPrizes = myPrizes2;
    myBudgetCost = myBudgetCost2;

    realTerminals = realTerminals2;
    myTerminals = myTerminals2;
    adjList = adjList2;

    nNodes = newNNodes;
    nEdges = newNEdges;
    nTerminals = myTerminals.size();

    trueTerminals = vector<bool>(nNodes, false);
}

int Instance::preprocessing() {
    int numberRemoved = 0;
    numberRemoved += uselessComponentsTest();
    numberRemoved += degreeOneTest();
    numberRemoved += degreeZeroTest();

    return numberRemoved;
}

int Instance::uselessComponentsTest() {
    int numberRemoved = 0;

    for (int i = 0; i < nComponents; ++i) {
        if (maxRevenueInComponent[i] < maxPrize) {
            numberRemoved += components[i].size();
            for (unsigned j = 0; j < components[i].size(); ++j) {
                nodesToRemove[components[i][j]] = true;
            }
        }
    }

    return numberRemoved;
}

int Instance::degreeOneTest() {
    int numberRemoved = 0;

    vector<int> toRemove;
    do {
        toRemove.clear();

        for (unsigned n = 0; n < nNodes; ++n) {
            if ((adjList[n].size() == 1 && !realTerminals[n])) {
                toRemove.push_back(n);
            }
        }

        for (unsigned i = 0; i < toRemove.size(); ++i) {

            int node = toRemove[i];
            int adjacentNode = adjList[node][0];
            adjList[node].clear();

            if (adjList[adjacentNode].size() == 0)
                continue;

            unsigned j = 0;
            for (; j < adjList[adjacentNode].size(); ++j) {
                if (adjList[adjacentNode][j] == node) {
                    break;
                }
            }
            adjList[adjacentNode].erase(adjList[adjacentNode].begin() + j);
        }

        numberRemoved += toRemove.size();
    } while (toRemove.size() > 0);

    return numberRemoved;
}

int Instance::degreeZeroTest() {
    int numberRemoved = 0;

    for (unsigned n = 0; n < nNodes; ++n) {
        if (adjList[n].size() == 0 && !nodesToRemove[n]) {
            nodesToRemove[n] = true;
            numberRemoved++;
        }
    }
    return numberRemoved;
}

int Instance::calculateComponents() {
    int numberOfComponents = 0;

    componentArray = vector<int>(nNodes, 0);
    components.clear();
    maxRevenueInComponent.clear();

    for (unsigned n = 0; n < nNodes; n++) {
        if (componentArray[n] == 0) {
            vector<int> componentHelper;
            double revInComp = 0.0;
            numberOfComponents++;
            componentArray[n] = numberOfComponents;
            if (myPrizes[n] > 0)
                revInComp += myPrizes[n];
            componentHelper.push_back(n);
            queue<int> myQueue;
            myQueue.push(n);
            while (!myQueue.empty()) {
                int m = myQueue.front();
                myQueue.pop();
                for (int k : adjList[m]) {
                    if (componentArray[k] == 0) {
                        componentArray[k] = numberOfComponents;
                        componentHelper.push_back(k);
                        if (myPrizes[k] > 0)
                            revInComp += myPrizes[k];
                        myQueue.push(k);
                    }
                }
            }
            components.push_back(componentHelper);
            maxRevenueInComponent.push_back(revInComp);
            componentFixed.push_back(0);
        }
    }

    return numberOfComponents;
}

double Instance::transformInternalValue(double value) const {
    return value;
}

void Instance::initStructures(unsigned nNodes) {
    realTerminals = vector<bool>(nNodes, false);
    trueTerminals = vector<bool>(nNodes, false);
    myTerminals = vector<int>(nNodes, -1);
    adjList = vector<vector<int>>(nNodes);
    nodesToRemove = vector<bool>(nNodes, false);
    myPrizes = vector<double>(nNodes);
    myBudgetCost = vector<double>(nNodes, 1);

    for (unsigned i = 0; i < nNodes; i++){
        adjList.emplace_back();
        myTerminals[i] = i;
    }
}

void Instance::addEdge(unsigned v, unsigned u) {
    adjList[v].push_back(u);
    adjList[u].push_back(v);
}

void Instance::readEdges(IntegerMatrix& edges) {
    nEdges = edges.nrow();
    bool edge_problem = edges.ncol() == 3;

    if (edge_problem) {
        for (int i = 0; i < nEdges; i++) {
            int from = edges(i, 0) - 1;
            int to = edges(i, 1) - 1;
            double weight = edges(i, 2);
            myPrizes[i + nNodes] = weight;
            addEdge(from, i + nNodes);
            addEdge(i + nNodes, to);
        }
        nNodes += nEdges;
        nEdges *= 2;
    } else {
        for (int i = 0; i < nEdges; i++) {
            int from = edges(i, 0) - 1;
            int to = edges(i, 1) - 1;
            addEdge(from, to);
        }
    }
}

void Instance::readInstance(List& instance) {
    auto edges = as<IntegerMatrix>(instance["edgelist"]);
    auto scores = as<NumericVector>(instance["vertex_weights"]);

    nNodes = static_cast<unsigned>(scores.size());
    unsigned m = static_cast<unsigned>(edges.nrow());
    initStructures(edges.ncol() == 3 ? nNodes + m : nNodes);

    NumericVector costs(myBudgetCost.begin(), myBudgetCost.end());
    if (instance.containsElementNamed("costs")) {
        costs = as<NumericVector>(instance["costs"]);
    }
    if (instance.containsElementNamed("budget")) {
        budget = instance["budget"];
    }
    if (instance.containsElementNamed("cardinality")) {
        cardCons = instance["cardinality"];
    }

    for (unsigned i = 0; i < nNodes; i++) {
        myPrizes[i] = scores[i];
        myBudgetCost[i] = costs[i];
        if (myPrizes[i] > 0) {
            realTerminals[i] = true;
        }
        if (myPrizes[i] > maxPrize) {
            maxPrize = myPrizes[i];
        }
    }

    readEdges(edges);
    nEdges *= 2;
}

void Instance::findSimpleSolution() {
    for (size_t i = 0; i < nNodes; i++) {
        if (myPrizes[i] > 0 && myBudgetCost[i] < budget) {
            solSize = 1;
            incumbent = std::vector<bool>(nNodes, false);
            incumbent[i] = true;
            incumbentObjLag = myPrizes[i];
            incumbentFound = true;
        }
    }
}
