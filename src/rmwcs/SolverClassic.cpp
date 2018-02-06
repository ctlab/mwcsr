/*
 * SolverClassic.cpp
 *
 *  Created on: Aug 14, 2015
 *      Author: markus
 */

#include "SolverClassic.h"
#include <queue>
#include <stack>

using std::vector;
using std::priority_queue;

using Rcpp::Rcout;

SolverClassic::SolverClassic(Instance& instance, Parameters& params)
        : SolverLag(instance, params) {
}

SolverClassic::~SolverClassic() {}

int SolverClassic::addInitCuts() {
    int nCuts = 0;
    // return 0;
    for (int i = 0; i < instance.nNodes; ++i) {
        cut myCut;
        myCut.lambda = 0;
        myCut.rhs = vector<nodevaluepair>();
        nodevaluepair n;
        n.id = i;
        n.value = 2.0;

        if (instance.realTerminals[i])
            n.value = 1.0;

        myCut.rhsConst = 0;
        myCut.rhs.push_back(n);

        myCut.lhs = vector<nodevaluepair>();
        // double minRhs=999;
        bool allpositive = true;

        for (int j : instance.adjList[i]) {
            nodevaluepair n;
            n.id = j;
            n.value = 1.0;
            // if(instance.myPrizes[j]<minRhs)
            //	minRhs=instance.myPrizes[j];
            if (instance.myPrizes[j] < 0)
                allpositive = false;

            myCut.lhs.push_back(n);
        }
        // if(minRhs<0 && instance.realTerminals[i])
        //	myCut.lambda=-minRhs;

        // if(minRhs<0 && !instance.realTerminals[i])
        //	myCut.lambda=-2*minRhs;
        // myCut.myHash=myCut.myHasher(myCut.lhs);
        // Rcout<<myCut.myHash<<"\n";

        if (!allpositive) {
            sort(myCut.lhs.begin(), myCut.lhs.end());
            // myCut.lambda=1e-5;
            // myCut.frozen=false;
            myCuts.push_back(myCut);
            myCut.myHash =
                    myCut.myHasher(myCut.lhs, myCut.rhs, instance.nNodes);
            myCutHash.insert(myCut.myHash);
        }
        nCuts++;
    }

    return nCuts;
}

bool SolverClassic::primalHeuristic() {

    bool improved = false;

    unsigned iter = 1;

    // if(myComponents.size()==1)
    //	Rcout<<"component size 1!!!"<<"\n";

    vector<int> myHeurTerminals;
    vector<bool> myHeurTerminalsBool = vector<bool>(instance.nNodes, false);

    vector<double> lpValue;

    double revInHeurTerms = 0.0;

    for (int i = 0; i < instance.nNodes; ++i) {
        if (currentSolution[i]) {
            lpValue.push_back(1);
            if (instance.realTerminals[i]) {
                myHeurTerminalsBool[i] = true;
                myHeurTerminals.push_back(i);
                revInHeurTerms += instance.myPrizes[i];
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

                if (myHeurTerminalsBool[c] && realPrizes[c] > bestVal) {
                    startID = c;
                    bestVal = realPrizes[c];
                }
            }
        } else {
            if (myIter > 0)
                break;

            for (int i = 0; i < instance.nNodes; ++i) {
                int c = i;
                if (myHeurTerminalsBool[c] && realPrizes[c] > bestVal) {
                    startID = c;
                    bestVal = realPrizes[c];
                }
            }
        }

        priority_queue<nodevaluepair, std::vector<nodevaluepair>,
                std::greater<nodevaluepair>>
                myPQueue;
        vector<int> inComponentBool = vector<int>(instance.nNodes, 0);
        vector<int> myBestSol = vector<int>(instance.nNodes, instance.nNodes + 1);
        vector<int> pred = vector<int>(instance.nNodes, -1);
        vector<int> extracted = vector<int>(instance.nNodes, 0);
        vector<double> distance = vector<double>(instance.nNodes, std::numeric_limits<double>::max());

        inComponentBool[startID] = true;
        myBestSol[startID] = -1;
        distance[startID] = 0;
        pred[startID] = startID;
        nodevaluepair n;
        n.id = startID;
        n.value = 0;
        obj += instance.myPrizes[startID];
        myPQueue.push(n);

        // Rcout<<obj<<"\n";

        double best = obj;
        int numTerms = 0;
        int numBest = 0;
        int nHeurTerm = myHeurTerminals.size();

        while (!myPQueue.empty() && numTerms < nHeurTerm) {
            if (revInHeurTerms + obj < incumbentObj)
                break;

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
                    if (distance[l] > distance[k] + toAdd + epsOpt) {
                        // cerr<<distance[l]<<" "<<distance[k] +
                        // (realPrizes[l]>0?realPrizes[l]:0)<<"\n";

                        distance[l] = distance[k] + toAdd;
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
                // Rcout<<k<<" "<<myHeurTerminalsBool[k]<<"
                // "<<instance.myPrizes[k]<<" "<<distance[k]<<"\n";
                revInHeurTerms -= instance.myPrizes[k];

                while (!inComponentBool[currentNode]) {
                    // cerr<<currentNode<<"\n";
                    inComponentBool[currentNode] = true;
                    myBestSol[currentNode] = numTerms;
                    distance[currentNode] = 0;
                    extracted[currentNode] = false;
                    nodevaluepair nv;
                    nv.id = currentNode;
                    nv.value = distance[currentNode];
                    myPQueue.push(nv);
                    obj += instance.myPrizes[currentNode];
                    currentNode = pred[currentNode];
                    // Rcout<<"pushed"<<"\n";
                }
                numTerms++;
                if (obj > best) {
                    best = obj;
                    numBest = numTerms;
                    // myBestSol=inComponentBool;
                    // Rcout<<"inbetween"<<instance.transformInternalValue(obj)<<"\n";
                }
            }
        }

        // post-processing
        for (int i = 0; i < instance.nNodes; ++i) {
            inComponentBool[i] = 0;
            if (myBestSol[i] < numBest)
                inComponentBool[i] = 1;
            // Rcout<<myBestSol[i]<<" "<<numBest<<"\n";
        }
        // inComponentBool=myBestSol;
        obj = best;
        vector<int> toAdd;

        do {
            toAdd.clear();
            for (int n = 0; n < instance.nNodes; ++n) {
                // adjList[n]=vector<int>();
                if (!inComponentBool[n])
                    continue;
                for (int j : instance.adjList[n]) {
                    // Rcout<<j<<" ";
                    if (!inComponentBool[j] && instance.myPrizes[j] > 0) {
                        toAdd.push_back(j);
                        // Rcout<<instance.myPrizes[j]<<"\n";
                    }
                }
            }

            for (unsigned i = 0; i < toAdd.size(); ++i) {
                int node = toAdd[i];
                if (!inComponentBool[node]) {
                    inComponentBool[node] = true;
                    obj += instance.myPrizes[node];
                }
            }
        } while (toAdd.size() > 0);

        if (obj > incumbentObj) {
            incumbentObj = obj;
            improved = true;
            double test = 0.0;
            for (int i = 0; i < instance.nNodes; ++i) {
                incumbent[i] = inComponentBool[i];
                // if(incumbent[i])Rcout<<i<<"\n";
                test += instance.myPrizes[i] * incumbent[i];
            }
            // Rcout<<test<<"\n";
            if (fabs(obj - test) > 1e-6) {
                Rf_error("Assertion failed on test incumbent solution");
            }

            // Rcout<<"improved "<<obj<<"\n";
        }
    }
    return improved;
}

int SolverClassic::lagrangianPegging() {
    int nFixed = 0;
    double boundPegging = 0;
    // Rcout<<"pegging"<<"\n";

    vector<int> fixToZero;
    vector<int> fixToOne;

    for (int i = 0; i < instance.nNodes; ++i) {
        if (fixedToZero[i] || fixedToOne[i])
            continue;

        //	Rcout<<i<<" "<<fixedToZero[i]<<" "<<fixedToOne[i]<<"\n";

        if (!currentSolution[i]) {
            boundPegging = currentBound + realPrizes[i];
            // Rcout<<"boundPegging"<<boundPegging<<" "<<incumbentObj<<"
            // "<<realPrizes[i]<<" "<<currentSolution[i]<<" "<<i<<"\n";

            if (boundPegging + epsInt < incumbentObj) {
                // Rcout<<"ZERO"<<i<<" "<<boundPegging<<" "<<incumbentObj<<"
                // "<<realPrizes[i] <<" "<<currentBound<<"\n";
                fixToZero.push_back(i);
                nFixed++;
            }
        } else if (currentSolution[i]) {
            boundPegging = currentBound - realPrizes[i];
            if (boundPegging + epsInt < incumbentObj) {
                // Rcout<<"ONE"<<i<<" "<<boundPegging<<" "<<incumbentObj<<"
                // "<<realPrizes[i] <<" "<<currentBound<<"\n";

                fixToOne.push_back(i);
                nFixed++;
            }
        }
    }
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

double SolverClassic::calculateCurrentSolution(bool save) {

    savedObj = calculateReducedCosts();

    double obj = savedObj;

    for (int i = 0; i < instance.nNodes; ++i) {
        currentSolution[i] = 0;

        if (fixedToOne[i]) {
            currentSolution[i] = 1.0;
            obj += realPrizes[i];
        } else if (fixedToZero[i]) {
            currentSolution[i] = 0;
        } else if (realPrizes[i] > 0.0) {
            currentSolution[i] = 1.0;
            obj += realPrizes[i];
        } else if (realPrizes[i] <= 0.0) {
            currentSolution[i] = 0;
        }

        if (save && currentSolution[i])
            sumSolution[i]++;
    }

    return obj;
}
