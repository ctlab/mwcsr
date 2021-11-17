/*
 * SolverCardinality.cpp
 *
 *  Created on: Aug 14, 2015
 *      Author: markus
 */

#include "include/SolverCardinality.h"
#include <queue>
#include <stack>

using std::vector;
using std::priority_queue;

using Rcpp::Rcout;

SolverCardinality::SolverCardinality(Instance& instance, Parameters& params, mwcsr::monitor monitor)
        : SolverLag(instance, params, monitor), weightLast{0.0}, weightOutside{0.0} {}

SolverCardinality::~SolverCardinality() {}

int SolverCardinality::addInitCuts() {
    int nCuts = 0;
    // return 0;
    for (unsigned i = 0; i < instance.nNodes; ++i) {
        cut myCut;
        myCut.lambda = 0;
        myCut.rhs = vector<nodevaluepair>();
        nodevaluepair n;
        n.id = i;
        n.value = 2.0;

        if (instance.myPrizes[i] > 0)
            n.value = 1.0;

        myCut.rhsConst = 0;
        myCut.rhs.push_back(n);

        myCut.lhs = vector<nodevaluepair>();

        for (int j : instance.adjList[i]) {
            nodevaluepair n;
            n.id = j;
            n.value = 1.0;
            myCut.lhs.push_back(n);
        }
        sort(myCut.lhs.begin(), myCut.lhs.end());
        myCut.myHash = myCut.myHasher(myCut.lhs, myCut.rhs, instance.nNodes);
        myCutHash.insert(myCut.myHash);

        myCuts.push_back(myCut);
        nCuts++;
    }

    return nCuts;
}

bool SolverCardinality::primalHeuristic() {

    bool improved = false;

    unsigned iter = 1;
    // if (myComponents.size()<iter)
    //	iter=myComponents.size();

    // if(iter==0)
    //	iter=1;

    // if(myComponents.size()==1)
    //	Rcout<<"component size 1!!!"<<"\n";

    vector<int> myHeurTerminals;
    vector<bool> myHeurTerminalsBool = vector<bool>(instance.nNodes, false);

    vector<double> lpValue;

    for (unsigned i = 0; i < instance.nNodes; ++i) {
        // Rcout<<(double)sumSolution[i]/iterations<<"\n";
        // lpValue.push_back((double)sumSolution[i]/iterations);

        if (currentSolution[i]) {
            lpValue.push_back(1);
            if (instance.realTerminals[i]) {
                myHeurTerminalsBool[i] = true;
                myHeurTerminals.push_back(i);
            }
        } else {
            lpValue.push_back(0);
        }
    }

    for (unsigned myIter = 0; myIter < iter; ++myIter) {
        int startID = 0;
        double bestVal = 0.0;
        double obj = 0;

        if (myComponents.size() > 0) {
            for (unsigned i = 0; i < myComponents[myIter].components.size();
                 ++i) {
                int c = myComponents[myIter].components[i];
                // Rcout<<c<<" "<<realPrizes[c]<<" "<<instance.myPrizes[c]<<"\n";

                if (myHeurTerminalsBool[c] && realPrizes[c] > bestVal) {
                    startID = c;
                    bestVal = realPrizes[c];
                }
            }
        } else {
            if (myIter > 0)
                break;

            for (unsigned i = 0; i < instance.nNodes; ++i) {
                int c = i;
                if (myHeurTerminalsBool[c] && realPrizes[c] > bestVal) {
                    startID = c;
                    bestVal = realPrizes[c];
                }
            }
        }

        vector<int> inComponentBool = vector<int>(instance.nNodes, false);
        vector<int> inComponent;
        vector<int> pred = vector<int>(instance.nNodes, -1);
        vector<bool> extracted = vector<bool>(instance.nNodes, false);
        vector<double> distance =
                vector<double>(instance.nNodes, std::numeric_limits<double>::max());
        vector<int> hop =
                vector<int>(instance.nNodes, std::numeric_limits<int>::max());

        priority_queue<nodevaluepair, std::vector<nodevaluepair>,
                std::greater<nodevaluepair>>
                myPQueue;
        vector<int> myBestSol =
                vector<int>(instance.nNodes, instance.nNodes + 1);

        myBestSol[startID] = 0;
        inComponentBool[startID] = true;
        distance[startID] = 0;
        pred[startID] = startID;
        nodevaluepair n;
        hop[startID] = 0;
        n.id = startID;
        n.value = 0;
        obj += instance.myPrizes[startID];
        myPQueue.push(n);

        // Rcout<<startID<<" "<<obj<<" "<<myIter<<"\n";

        // Rcout<<obj<<"\n";

        double best = obj;
        int numTerms = 0;
        int numBest = 0;

        int currentSize = 1;
        int bestSize = 1;

        while (!myPQueue.empty() && currentSize < instance.cardCons) {
            nodevaluepair k2 = myPQueue.top();
            myPQueue.pop();
            int k = k2.id;

            if (extracted[k])
                continue;

            extracted[k] = true;

            // Rcout<<"k "<<k<<" "<<distance[k]<<"
            // "<<myHeurTerminalsBool[k]<<"\n";

            if (!myHeurTerminalsBool[k] || inComponentBool[k]) {
                double toAdd = -instance.myPrizes[k] * (1 - lpValue[k]);
                if (toAdd <= 0)
                    toAdd = epsOpt;
                // Rcout<<toAdd<<"\n";

                for (int l : instance.adjList[k]) {
                    if (distance[l] > distance[k] + toAdd + epsOpt ||
                        (distance[l] == distance[k] + toAdd &&
                         hop[l] > hop[k] + 1))
                        // if (hop[l]>hop[k]+(1-lpValue[k]))
                    {
                        // cerr<<distance[l]<<" "<<distance[k] +
                        // (realPrizes[l]>0?realPrizes[l]:0)<<"\n";

                        distance[l] = distance[k] + toAdd;
                        hop[l] = hop[k] + 1;
                        pred[l] = k;

                        extracted[l] = false;
                        nodevaluepair lNv;
                        lNv.id = l;
                        lNv.value = distance[l];
                        myPQueue.push(lNv);
                        // Rcout<<l<<" "<<k<<" "<<distance[l]<<"
                        // "<<distance[k]<<" "<<currentSolution[l]<<"
                        // "<<lNv.value<<"\n";
                    }

                    if (distance[l] < 0) {
                        Rcout << l << " " << distance[l] << " " << distance[k]
                              << " " << toAdd << " " << inComponentBool[k]
                              << "\n";
                        Rf_error("Primal heuristic bug.");
                    }
                }
            } else {
                int currentNode = k;

                bool fit = true;
                int add = 0;
                while (!inComponentBool[currentNode]) {
                    // Rcout<<currentNode<<" "<<pred[currentNode]<<"\n";

                    add++;
                    if (currentSize + add > instance.cardCons) {
                        fit = false;
                        break;
                    }
                    currentNode = pred[currentNode];
                }

                // Rcout<<add<<" "<<currentSize<<" "<<fit<<"\n";

                if (fit) {
                    currentNode = k;
                    while (!inComponentBool[currentNode]) {
                        // cerr<<currentNode<<"\n";
                        inComponentBool[currentNode] = true;
                        hop[currentNode] = 0;
                        distance[currentNode] = 0;
                        extracted[currentNode] = false;
                        myBestSol[currentNode] = numTerms;
                        nodevaluepair nv;
                        nv.id = currentNode;
                        nv.value = distance[currentNode];
                        myPQueue.push(nv);
                        obj += instance.myPrizes[currentNode];
                        currentNode = pred[currentNode];
                        currentSize++;
                        // Rcout<<"pushed"<<"\n";
                    }
                    numTerms++;
                    if (obj > best) {
                        best = obj;
                        numBest = numTerms;
                        bestSize = currentSize;
                        // Rcout<<"best"<<obj<<" "<<currentSize<<"\n";
                    }
                }
            }
        }

        // Rcout<<obj<<" "<<currentSize<<"\n";
        if (best > incumbentObj && false) {
            obj = 0.0;
            currentSize = 0;
            // Rcout<<"best "<<best<<" "<<numBest<<"\n";
            inComponentBool = vector<int>(instance.nNodes, false);

            myPQueue = priority_queue<nodevaluepair, std::vector<nodevaluepair>,
                    std::greater<nodevaluepair>>();

            pred = vector<int>(instance.nNodes, -1);
            extracted = vector<bool>(instance.nNodes, false);
            distance = vector<double>(instance.nNodes,
                                      std::numeric_limits<double>::max());
            hop = vector<int>(instance.nNodes, std::numeric_limits<int>::max());

            int c = startID;
            inComponentBool[c] = true;
            distance[c] = 0;
            pred[c] = c;
            nodevaluepair n;
            n.id = c;
            n.value = 0;
            obj += instance.myPrizes[c];
            myPQueue.push(n);
            hop[c] = 0;

            int countTerm = 0;
            currentSize = 1;
            while (!myPQueue.empty() && numBest > 0 &&
                   currentSize < instance.cardCons) {
                nodevaluepair k2 = myPQueue.top();
                myPQueue.pop();
                int k = k2.id;

                if (extracted[k])
                    continue;

                extracted[k] = true;

                // Rcout<<"k "<<k<<" "<<distance[k]<<"
                // "<<myHeurTerminalsBool[k]<<"\n";

                if (!myHeurTerminalsBool[k] || inComponentBool[k]) {
                    double toAdd = -instance.myPrizes[k] * (1 - lpValue[k]);
                    if (toAdd <= 0)
                        toAdd = epsOpt;
                    // Rcout<<toAdd<<"\n";

                    for (int l : instance.adjList[k]) {
                        if (distance[l] > distance[k] + toAdd + epsOpt ||
                            (distance[l] == distance[k] + toAdd &&
                             hop[l] > hop[k] + 1))
                            // if (hop[l]>hop[k]+(1-lpValue[k]))
                        {
                            // cerr<<distance[l]<<" "<<distance[k] +
                            // (realPrizes[l]>0?realPrizes[l]:0)<<"\n";

                            distance[l] = distance[k] + toAdd;
                            pred[l] = k;
                            hop[l] = hop[k] + 1;
                            extracted[l] = false;
                            nodevaluepair lNv;
                            lNv.id = l;
                            lNv.value = distance[l];
                            myPQueue.push(lNv);
                            // Rcout<<l<<" "<<k<<" "<<distance[l]<<"
                            // "<<distance[k]<<" "<<currentSolution[l]<<"
                            // "<<lNv.value<<"\n";
                        }

                        if (distance[l] < 0) {
                            Rcout << l << " " << distance[l] << " "
                                  << distance[k] << " " << toAdd << " "
                                  << inComponentBool[k] << "\n";
                            exit(1);
                        }
                    }
                } else {
                    int currentNode = k;

                    bool fit = true;
                    int add = 0;
                    while (!inComponentBool[currentNode]) {
                        // Rcout<<currentNode<<" "<<pred[currentNode]<<"\n";

                        add++;
                        if (currentSize + add > instance.cardCons) {
                            fit = false;
                            break;
                        }
                        currentNode = pred[currentNode];
                    }

                    // Rcout<<add<<" "<<currentSize<<" "<<fit<<"\n";

                    if (fit) {
                        currentNode = k;

                        // Rcout<<k<<" "<<myHeurTerminalsBool[k]<<"
                        // "<<instance.myPrizes[k]<<" "<<distance[k]<<"\n";
                        while (!inComponentBool[currentNode]) {
                            // cerr<<currentNode<<"\n";
                            inComponentBool[currentNode] = true;
                            distance[currentNode] = 0;
                            extracted[currentNode] = false;
                            nodevaluepair nv;
                            nv.id = currentNode;
                            nv.value = distance[currentNode];
                            myPQueue.push(nv);
                            hop[currentNode] = 0;
                            // Rcout<<instance.myPrizes[currentNode]<<"\n";
                            obj += instance.myPrizes[currentNode];
                            currentNode = pred[currentNode];
                            currentSize++;
                            // Rcout<<"pushed"<<"\n";
                        }
                        countTerm++;
                        // Rcout<<countTerm<<" "<<numBest<<" "<<obj<<"\n";
                        if (countTerm == numBest) {
                            // Rcout<<"inbetween"<<obj<<"\n";
                            break;
                        }
                    }
                }
            }
            // Rcout<<"best"<<obj<<"\n";
        }

        // post-processing
        for (unsigned i = 0; i < instance.nNodes; ++i) {
            inComponentBool[i] = 0;
            if (myBestSol[i] < numBest)
                inComponentBool[i] = 1;
            // Rcout<<myBestSol[i]<<" "<<numBest<<"\n";
        }
        currentSize = bestSize;

        obj = best;
        vector<nodevaluepair> toAdd;

        do {
            toAdd.clear();
            // Rcout<<"before postprocessing "<<obj<<" "<<currentSize<<"\n";

            for (unsigned n = 0; n < instance.nNodes; ++n) {
                if (currentSize >= instance.cardCons)
                    break;
                if (!inComponentBool[n])
                    continue;
                for (int j : instance.adjList[n]) {
                    if (!inComponentBool[j] && instance.myPrizes[j] > 0) {
                        nodevaluepair nv;
                        nv.id = j;
                        nv.value = instance.myPrizes[j];
                        toAdd.push_back(nv);
                    }
                }
            }

            sort(toAdd.begin(), toAdd.end(), std::greater<nodevaluepair>());

            // Rcout<<"before"<<obj<<" "<<currentSize<<"\n";
            // Rcout<<"before"<<obj<<"\n";

            for (unsigned i = 0; i < toAdd.size(); ++i) {
                int node = toAdd[i].id;
                if (!inComponentBool[node]) {
                    inComponentBool[node] = true;
                    obj += instance.myPrizes[node];
                }
                currentSize++;
                if (currentSize >= instance.cardCons)
                    break;
            }
        } while (toAdd.size() > 0);

        // Rcout<<obj<<" "<<currentSize<<"\n";

        // Rcout<<toAdd.size()<<" "<<obj<<"\n";

        if (obj > incumbentObj) {
            incumbentObj = obj;
            improved = true;
            for (unsigned i = 0; i < instance.nNodes; ++i) {
                incumbent[i] = inComponentBool[i];
            }

            // Rcout<<"improved "<<obj<<"\n";
        }
    }
    return improved;
}

int SolverCardinality::lagrangianPegging() {
    int nFixed = 0;
    double boundPegging = 0;
    // Rcout<<"pegging"<<"\n";

    vector<int> fixToZero;
    vector<int> fixToOne;

    for (unsigned i = 0; i < instance.nNodes; ++i) {
        if (fixedToZero[i] || fixedToOne[i] ||
            fabs(realPrizes[i] - epsOpt) < epsOpt)
            continue;

        //	Rcout<<i<<" "<<fixedToZero[i]<<" "<<fixedToOne[i]<<"\n";

        if (!currentSolution[i]) {
            boundPegging = currentBound + realPrizes[i] - weightLast;
            // Rcout<<"boundPegging"<<boundPegging<<" "<<incumbentObj<<"
            // "<<realPrizes[i]<<" "<<currentSolution[i]<<" "<<i<<"\n";

            if (boundPegging < incumbentObj) {
                fixToZero.push_back(i);
                nFixed++;
            }
        } else if (currentSolution[i]) {
            boundPegging = currentBound - realPrizes[i] + weightOutside;
            // Rcout<<"boundPegging"<<boundPegging<<"\n";
            if (boundPegging < incumbentObj) {
                fixToOne.push_back(i);
                nFixed++;
            }
        }
    }
    // Rcout<<fixToZero.size()<<" "<<fixToOne.size()<<"\n";
    for (int i : fixToZero) {
        fixedToZero[i] = true;
        instance.nFixedZero++;
        for (int j : instance.adjList[i]) {
            unsigned k = 0;
            for (; k < instance.adjList[j].size(); ++k) {
                if (instance.adjList[j][k] == i)
                    break;
            }
            instance.adjList[j].erase(instance.adjList[j].begin() + k);
        }
        instance.adjList[i].clear();

        // Rcout<<i<<"\n";
    }

    for (int i : fixToOne) {
        fixedToOne[i] = true;
        instance.nFixedOne++;
        // Rcout<<i<<"\n";
    }

    return nFixed;
}

double SolverCardinality::calculateCurrentSolution(bool save) {

    savedObj = calculateReducedCosts();

    double obj = savedObj;

    vector<nodevaluepair> myNV;
    int solutionSize = 0;

    int nFixed = 0;
    for (unsigned i = 0; i < instance.nNodes; ++i) {
        currentSolution[i] = 0;

        if (fixedToOne[i]) {
            // nFixed++;
            currentSolution[i] = 1.0;
            obj += realPrizes[i];
            solutionSize++;
            if (save && currentSolution[i])
                sumSolution[i]++;
            if (solutionSize > instance.cardCons) {
                Rcout << "strange" << "\n";
            }
            continue;
        } else if (fixedToZero[i]) {
            currentSolution[i] = 0;
            nFixed++;
            continue;
        }
        nodevaluepair nv;
        nv.id = i;
        nv.value = realPrizes[i];
        myNV.push_back(nv);
    }

    // Rcout<<nFixed<<"\n";
    sort(myNV.begin(), myNV.end(), std::greater<nodevaluepair>());

    // Rcout<<"obj before"<<obj<<"\n";

    unsigned i = 0;
    for (; i < myNV.size(); ++i) {
        nodevaluepair nv = myNV[i];
        if (realPrizes[nv.id] <= 0)
            break;
        currentSolution[nv.id] = 1.0;
        obj += realPrizes[nv.id];
        solutionSize++;
        // Rcout<<nv.id<<" "<<realPrizes[nv.id]<<" "<<solutionSize<<"
        // "<<instance.myPrizes[nv.id]<<"\n";

        if (save && currentSolution[nv.id])
            sumSolution[nv.id]++;
        if (solutionSize >= instance.cardCons)
            break;
    }
    if (i < myNV.size())
        weightLast = myNV[i].value;
    else
        weightLast = 0;
    if ((i + 1) < myNV.size())
        weightOutside = myNV[i + 1].value;
    else
        weightOutside = 0;

    if (weightLast < 0)
        weightLast = 0;
    if (weightOutside < 0)
        weightOutside = 0;

    // Rcout<<"obj "<<obj<<"\n";

    return obj;
}
