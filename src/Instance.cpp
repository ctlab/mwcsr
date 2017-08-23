/*
 * Instance.cpp
 *
 *  Created on: Mar 29, 2015
 *      Author: markus
 */

#include "instance/Instance.h"

#include <fstream>
#include <queue>

using Rcpp::NumericVector;
using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::as;

Instance::Instance(List parameters, List network)
    : components{vector<vector<int>>()},
      maxRevenueInComponent{vector<double>()}, nComponents{0}, maxPrize{0},
      minWeight{std::numeric_limits<double>::max()}, sumPrizes{0},
      nRealTerminals{0}, params(Parameters(parameters)) {

    readInstance(network);

    cout << "Instance reading finished" << endl;
    cout << "Nodes: \t \t \t" << nNodes << endl;
    cout << "Edges: \t \t \t" << nEdges << endl;
    cout << "TerminalsRead: \t \t" << nTerminals << endl;

    nTrueNodes = nNodes;
    nTrueEdges = nEdges;

    nComponents = calculateComponents();
    cout << "Components:  \t \t " << nComponents << endl;
    if (params.problem == 1)
        cout << "Budget:  \t \t " << budget << endl;
    if (params.problem == 2)
        cout << "Cardinality:  \t \t " << params.cardCons << endl;
    // exit(0);
    if (preprocessing() > 0) {
        rebuildDatastructures();
        cout << "Size after preprocessing" << endl;
        cout << "Nodes: \t \t \t" << nNodes << endl;
        cout << "Edges: \t \t \t" << nEdges << endl;
        cout << "Terminals:  \t \t" << nTerminals << endl;
        nComponents = calculateComponents();
        cout << "Components:  \t \t " << nComponents << endl;
    }

    fixedToOne = vector<int>(nNodes, 0);
    fixedToZero = vector<int>(nNodes, 0);

    nFlowNodes = 2 * nNodes;
    nFlowArcs = 2 * nEdges + nNodes;

    flowArcs = new std::pair<int, int>[nFlowArcs];
    int counter = 0;
    for (int i = 0; i < nNodes; ++i) {
        flowArcs[counter] = std::make_pair(i, i + nNodes);
        counter++;
    }

    for (int i = 0; i < nNodes; ++i) {
        for (int j : adjList[i]) {
            flowArcs[counter] = std::make_pair(i + nNodes, j);
            counter++;
        }
    }
}

/*void Instance::readMWCS() {

    std::ifstream infile(params.file);
    std::string line;

    std::string magicNumber;
    std::getline(infile, line);
    std::istringstream iss(line);
    iss.precision(16);

    iss >> magicNumber;

    if (magicNumber != "33D32945") {
        cout << "Not in STP Format" << endl;
        exit(0);
    }

    minWeight = 0;
    int srcID, dstID;
    int edgeCounter = 0;
    int termID;
    double prize;
    int termCounter = 0;
    int count = 0;
    int countRealTerminals = 0;
    while (std::getline(infile, line)) {
        if (line.empty())
            continue;

        iss.clear();
        iss.str(line);
        std::string token;

        iss >> token;
        // cerr<<count<<" "<<token<<" "<<line<<" "<<iss.str()<<endl;
        count++;
        if (token == "Nodes") {

            iss >> nNodes;
            adjList = vector<vector<int>>(nNodes);
            adjList = vector<vector<int>>(nNodes);
            realTerminals = vector<bool>(nNodes, false);
            trueTerminals = vector<bool>(nNodes, false);
            myPrizes = vector<double>(nNodes);
            myBudgetCost = vector<double>(nNodes, 1);
            nodesToRemove = vector<bool>(nNodes);
            for (int i = 0; i < nNodes; i++) {
                adjList[i] = vector<int>();
            }

        } else if (token == "Edges") {
            iss >> nEdges;
            // cerr<<nEdges<<endl;

        } else if (token == "E") {
            iss >> srcID >> dstID;
            // cerr<<edgeCounter<<endl;
            adjList[srcID - 1].push_back(dstID - 1);
            adjList[dstID - 1].push_back(srcID - 1);
            // myCost[edgeCounter]=0;
            edgeCounter++;
            edgeCounter++;
        } else if (token == "Terminals") {
            iss >> nTerminals;
            myTerminals = vector<int>(nTerminals, -1);
        } else if (token == "T") {
            iss >> termID >> prize;
            myTerminals[termCounter] = termID - 1;
            myPrizes[termID - 1] = prize;
            // cout.precision(16);
            // cout<<(termID-1)<<" "<<prize<<endl;

            if (myPrizes[termID - 1] > 0) {
                realTerminals[termID - 1] = true;
                countRealTerminals++;
            }

            if (myPrizes[termID - 1] > maxPrize) {
                maxPrize = myPrizes[termID - 1];
            }

            termCounter++;
        }
    }

    if (params.problem == 2 && params.cardMultiplier >= 1) {
        int nPositive = 0;
        for (unsigned i = 0; i < myPrizes.size(); ++i) {
            if (myPrizes[i] > 0) {
                nPositive++;
            }
        }
        params.cardCons = floor(nPositive * params.cardMultiplier / 100.0);
    }

    // cout<<"Real Terminals \t"<<countRealTerminals<<endl;
}

void Instance::readMWCSBudget() {

    std::ifstream infile(params.file);
    std::string line;

    std::string magicNumber;
    std::getline(infile, line);
    std::istringstream iss(line);
    iss.precision(16);

    iss >> magicNumber;

    if (magicNumber != "33D32945") {
        cout << "Not in STP Format" << endl;
        exit(0);
    }

    minWeight = 0;
    int srcID, dstID;
    int edgeCounter = 0;
    int termID;
    double prize;
    int termCounter = 0;
    int count = 0;
    int countRealTerminals = 0;
    int myMaxPrize = 0;
    while (std::getline(infile, line)) {
        if (line.empty())
            continue;

        iss.clear();
        iss.str(line);
        std::string token;

        iss >> token;
        // cerr<<count<<" "<<token<<" "<<line<<" "<<iss.str()<<endl;
        count++;
        if (token == "Nodes") {

            iss >> nNodes;
            adjList = vector<vector<int>>(nNodes);
            // cerr<<nNodes<<endl;
            adjList = vector<vector<int>>(nNodes);
            realTerminals = vector<bool>(nNodes, false);
            trueTerminals = vector<bool>(nNodes, false);
            myPrizes = vector<double>(nNodes);
            myBudgetCost = vector<double>(nNodes);
            nodesToRemove = vector<bool>(nNodes);
            for (int i = 0; i < nNodes; i++) {
                adjList[i] = vector<int>();
            }

        } else if (token == "Edges") {
            iss >> nEdges;
            // cerr<<nEdges<<endl;

        } else if (token == "E") {
            iss >> srcID >> dstID;
            // cerr<<edgeCounter<<endl;
            adjList[srcID - 1].push_back(dstID - 1);
            adjList[dstID - 1].push_back(srcID - 1);
            // myCost[edgeCounter]=0;
            edgeCounter++;
            edgeCounter++;
        } else if (token == "Terminals") {
            iss >> nTerminals;
            myTerminals = vector<int>(nTerminals, -1);
        } else if (token == "T") {
            iss >> termID >> prize >> myMaxPrize;
            myTerminals[termCounter] = termID - 1;
            myPrizes[termID - 1] = prize;
            myBudgetCost[termID - 1] = myMaxPrize;
            // cout.precision(16);
            // cout<<(termID-1)<<" "<<prize<<endl;
            if (myPrizes[termID - 1] > 0) {
                realTerminals[termID - 1] = true;
                countRealTerminals++;
            }

            if (myPrizes[termID - 1] > maxPrize) {
                maxPrize = myPrizes[termID - 1];
            }

            termCounter++;
        }
    }

    if (params.problem == 2 && params.cardMultiplier >= 1) {
        int nPositive = 0;
        for (unsigned i = 0; i < myPrizes.size(); ++i) {
            if (myPrizes[i] > 0) {
                nPositive++;
            }
        }
        params.cardCons = floor(nPositive * params.cardMultiplier / 100.0);
    }

    if (params.problem == 1 && params.budgetMultiplier >= 1) {
        double totalCost = 0;
        for (int i = 0; i < nNodes; ++i) {
            totalCost += myBudgetCost[i];
            // cout<<myBudgetCost[i]<<" "<<myPrizes[i]<<endl;
        }

        budget = floor(totalCost * params.budgetMultiplier / 100.0);
    }

    cout << "Real Terminals \t" << countRealTerminals << endl;
}*/

void Instance::rebuildDatastructures() {

    vector<double> myPrizes2;
    vector<double> myBudgetCost2;
    vector<bool> realTerminals2;
    vector<int> myTerminals2;
    vector<vector<int>> adjList2;

    vector<int> backMap = vector<int>(nNodes);

    int newNNodes = 0;

    for (int i = 0; i < nNodes; ++i) {
        if (!nodesToRemove[i]) {
            backMap[i] = newNNodes;
            map.push_back(i);
            adjList2.push_back(vector<int>());
            myPrizes2.push_back(myPrizes[i]);
            myBudgetCost2.push_back(myBudgetCost[i]);
            realTerminals2.push_back(realTerminals[i]);
            if (realTerminals[i])
                myTerminals2.push_back(newNNodes);
            newNNodes++;
        }
    }

    int newNEdges = 0;
    for (int i = 0; i < nNodes; ++i) {
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

    // componentArray.init(G);
}

int Instance::preprocessing() {
    // return 0;

    int numberRemoved = 0;
    numberRemoved += uselessComponentsTest();
    numberRemoved += degreeOneTest();
    numberRemoved += degreeZeroTest();

    return numberRemoved;
}

int Instance::uselessComponentsTest() {
    int numberRemoved = 0;

    for (int i = 0; i < nComponents; ++i) {
        // cerr<<maxRevenueInComponent[i]<<" "<<maxPrize<<endl;
        if (maxRevenueInComponent[i] <= maxPrize) {
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

        for (int n = 0; n < nNodes; ++n) {
            if ((adjList[n].size() == 1 && !realTerminals[n])) {
                toRemove.push_back(n);
            }
            /*if((adjList[n].size()==1 && realTerminals[n]) &&
            params.problem==0)
            {
                    toRemove.push_back(n);
                    int adjacentNode=adjList[n][0];
                    myPrizes[adjacentNode]+=myPrizes[n];
                    if(myPrizes[adjacentNode]>0)
                            realTerminals[adjacentNode]=true;
            }*/
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
            // cout<<j<<" "<<adjList[adjacentNode].size()<<endl;
            // cout<<node<<" "<<adjList[adjacentNode][j]<<endl;
            adjList[adjacentNode].erase(adjList[adjacentNode].begin() + j);
        }

        numberRemoved += toRemove.size();
    } while (toRemove.size() > 0);

    return numberRemoved;
}

int Instance::degreeZeroTest() {
    int numberRemoved = 0;

    for (int n = 0; n < nNodes; ++n) {
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

    for (int n = 0; n < nNodes; n++) {
        if (componentArray[n] == 0) {
            // cerr<<n<<endl;
            vector<int> componentHelper;
            double revInComp = 0.0;
            numberOfComponents++;
            componentArray[n] = numberOfComponents;
            // cout<<n<<endl;
            // cout<<myPrizes[n]<<endl;
            if (myPrizes[n] > 0)
                revInComp += myPrizes[n];
            componentHelper.push_back(n);
            queue<int> myQueue;
            myQueue.push(n);
            // cerr<<"component "<<n<<" "<<myCurrentLabel<<endl;
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

    /*cout<<"number of components "<<numberOfComponents<<endl;
    for(int i=0;i<numberOfComponents;++i)
    {
            cout<<"size: \t"<<components[i].size()<<"\t
    "<<maxRevenueInComponent[i]<<endl;
    }*/

    return numberOfComponents;
}

double Instance::transformInternalValue(double value) const {
    return value;
}

void Instance::readInstance(List instance) {
    auto edges = as<NumericMatrix>(instance["edgelist"]);
    auto scores = as<NumericVector>(instance["scores"]);
    nNodes = static_cast<unsigned>(scores.size());
    nEdges = static_cast<unsigned>(edges.nrow());
    realTerminals = vector<bool>(nNodes, false);
    trueTerminals = vector<bool>(nNodes, false);
    myTerminals = vector<int>(nNodes, -1);
    adjList = vector<vector<int>>(nNodes);
    nodesToRemove = vector<bool>(nNodes, false);

    for (unsigned i = 0; i < nNodes; i++){
        adjList.emplace_back();
        myPrizes[i] = scores[i];
        myTerminals[i] = i;
        if(scores[i] > 0){
            realTerminals[i] = true;
        }
        if(scores[i] > maxPrize){
            maxPrize = scores[i];
        }
    }

    for(unsigned i = 0; i < nEdges; i++){
        int from = edges(i, 0);
        int to = edges(i, 1);
        adjList[from].push_back(to);
        adjList[to].push_back(from);
    }
}
