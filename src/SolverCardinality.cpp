/*
 * SolverCardinality.cpp
 *
 *  Created on: Aug 14, 2015
 *      Author: markus
 */

#include "solverLag/SolverCardinality.h"
#include <queue>
#include <stack>

// using namespace std;

SolverCardinality::SolverCardinality(Instance &_instance, int _maxIterations)
    : SolverLag(_instance, _maxIterations), weightLast{0.0}, weightOutside{
                                                                 0.0} {
}

SolverCardinality::~SolverCardinality() {}

int SolverCardinality::addInitCuts() {
    int nCuts = 0;
    // return 0;
    for (int i = 0; i < instance.nNodes; ++i) {
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
    //	cout<<"component size 1!!!"<<endl;

    vector<int> myHeurTerminals;
    vector<bool> myHeurTerminalsBool = vector<bool>(instance.nNodes, false);

    vector<double> lpValue;

    for (int i = 0; i < instance.nNodes; ++i) {
        // cout<<(double)sumSolution[i]/iterations<<endl;
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
                // cout<<c<<" "<<realPrizes[c]<<" "<<instance.myPrizes[c]<<endl;

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

        // cout<<startID<<" "<<obj<<" "<<myIter<<endl;

        // cout<<obj<<endl;

        double best = obj;
        int numTerms = 0;
        int numBest = 0;

        int currentSize = 1;
        int bestSize = 1;

        while (!myPQueue.empty() && currentSize < instance.params.cardCons) {
            nodevaluepair k2 = myPQueue.top();
            myPQueue.pop();
            int k = k2.id;

            if (extracted[k])
                continue;

            extracted[k] = true;

            // cout<<"k "<<k<<" "<<distance[k]<<"
            // "<<myHeurTerminalsBool[k]<<endl;

            if (!myHeurTerminalsBool[k] || inComponentBool[k]) {
                double toAdd = -instance.myPrizes[k] * (1 - lpValue[k]);
                if (toAdd <= 0)
                    toAdd = epsOpt;
                // cout<<toAdd<<endl;

                for (int l : instance.adjList[k]) {
                    if (distance[l] > distance[k] + toAdd + epsOpt ||
                        (distance[l] == distance[k] + toAdd &&
                         hop[l] > hop[k] + 1))
                    // if (hop[l]>hop[k]+(1-lpValue[k]))
                    {
                        // cerr<<distance[l]<<" "<<distance[k] +
                        // (realPrizes[l]>0?realPrizes[l]:0)<<endl;

                        distance[l] = distance[k] + toAdd;
                        hop[l] = hop[k] + 1;
                        pred[l] = k;

                        extracted[l] = false;
                        nodevaluepair lNv;
                        lNv.id = l;
                        lNv.value = distance[l];
                        myPQueue.push(lNv);
                        // cout<<l<<" "<<k<<" "<<distance[l]<<"
                        // "<<distance[k]<<" "<<currentSolution[l]<<"
                        // "<<lNv.value<<endl;
                    }

                    if (distance[l] < 0) {
                        cout << l << " " << distance[l] << " " << distance[k]
                             << " " << toAdd << " " << inComponentBool[k]
                             << endl;
                        exit(1);
                    }
                }
            } else {
                int currentNode = k;

                bool fit = true;
                int add = 0;
                while (!inComponentBool[currentNode]) {
                    // cout<<currentNode<<" "<<pred[currentNode]<<endl;

                    add++;
                    if (currentSize + add > instance.params.cardCons) {
                        fit = false;
                        break;
                    }
                    currentNode = pred[currentNode];
                }

                // cout<<add<<" "<<currentSize<<" "<<fit<<endl;

                if (fit) {
                    currentNode = k;
                    while (!inComponentBool[currentNode]) {
                        // cerr<<currentNode<<endl;
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
                        // cout<<"pushed"<<endl;
                    }
                    numTerms++;
                    if (obj > best) {
                        best = obj;
                        numBest = numTerms;
                        bestSize = currentSize;
                        // cout<<"best"<<obj<<" "<<currentSize<<endl;
                    }
                }
            }
        }

        // cout<<obj<<" "<<currentSize<<endl;
        if (best > incumbentObj && false) {
            obj = 0.0;
            currentSize = 0;
            // cout<<"best "<<best<<" "<<numBest<<endl;
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
                   currentSize < instance.params.cardCons) {
                nodevaluepair k2 = myPQueue.top();
                myPQueue.pop();
                int k = k2.id;

                if (extracted[k])
                    continue;

                extracted[k] = true;

                // cout<<"k "<<k<<" "<<distance[k]<<"
                // "<<myHeurTerminalsBool[k]<<endl;

                if (!myHeurTerminalsBool[k] || inComponentBool[k]) {
                    double toAdd = -instance.myPrizes[k] * (1 - lpValue[k]);
                    if (toAdd <= 0)
                        toAdd = epsOpt;
                    // cout<<toAdd<<endl;

                    for (int l : instance.adjList[k]) {
                        if (distance[l] > distance[k] + toAdd + epsOpt ||
                            (distance[l] == distance[k] + toAdd &&
                             hop[l] > hop[k] + 1))
                        // if (hop[l]>hop[k]+(1-lpValue[k]))
                        {
                            // cerr<<distance[l]<<" "<<distance[k] +
                            // (realPrizes[l]>0?realPrizes[l]:0)<<endl;

                            distance[l] = distance[k] + toAdd;
                            pred[l] = k;
                            hop[l] = hop[k] + 1;
                            extracted[l] = false;
                            nodevaluepair lNv;
                            lNv.id = l;
                            lNv.value = distance[l];
                            myPQueue.push(lNv);
                            // cout<<l<<" "<<k<<" "<<distance[l]<<"
                            // "<<distance[k]<<" "<<currentSolution[l]<<"
                            // "<<lNv.value<<endl;
                        }

                        if (distance[l] < 0) {
                            cout << l << " " << distance[l] << " "
                                 << distance[k] << " " << toAdd << " "
                                 << inComponentBool[k] << endl;
                            exit(1);
                        }
                    }
                } else {
                    int currentNode = k;

                    bool fit = true;
                    int add = 0;
                    while (!inComponentBool[currentNode]) {
                        // cout<<currentNode<<" "<<pred[currentNode]<<endl;

                        add++;
                        if (currentSize + add > instance.params.cardCons) {
                            fit = false;
                            break;
                        }
                        currentNode = pred[currentNode];
                    }

                    // cout<<add<<" "<<currentSize<<" "<<fit<<endl;

                    if (fit) {
                        currentNode = k;

                        // cout<<k<<" "<<myHeurTerminalsBool[k]<<"
                        // "<<instance.myPrizes[k]<<" "<<distance[k]<<endl;
                        while (!inComponentBool[currentNode]) {
                            // cerr<<currentNode<<endl;
                            inComponentBool[currentNode] = true;
                            distance[currentNode] = 0;
                            extracted[currentNode] = false;
                            nodevaluepair nv;
                            nv.id = currentNode;
                            nv.value = distance[currentNode];
                            myPQueue.push(nv);
                            hop[currentNode] = 0;
                            // cout<<instance.myPrizes[currentNode]<<endl;
                            obj += instance.myPrizes[currentNode];
                            currentNode = pred[currentNode];
                            currentSize++;
                            // cout<<"pushed"<<endl;
                        }
                        countTerm++;
                        // cout<<countTerm<<" "<<numBest<<" "<<obj<<endl;
                        if (countTerm == numBest) {
                            // cout<<"inbetween"<<obj<<endl;
                            break;
                        }
                    }
                }
            }
            // cout<<"best"<<obj<<endl;
        }

        // post-processing
        for (int i = 0; i < instance.nNodes; ++i) {
            inComponentBool[i] = 0;
            if (myBestSol[i] < numBest)
                inComponentBool[i] = 1;
            // cout<<myBestSol[i]<<" "<<numBest<<endl;
        }
        currentSize = bestSize;

        obj = best;
        vector<nodevaluepair> toAdd;

        do {
            toAdd.clear();
            // cout<<"before postprocessing "<<obj<<" "<<currentSize<<endl;

            for (int n = 0; n < instance.nNodes; ++n) {
                if (currentSize >= instance.params.cardCons)
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

            // cout<<"before"<<obj<<" "<<currentSize<<endl;
            // cout<<"before"<<obj<<endl;

            for (unsigned i = 0; i < toAdd.size(); ++i) {
                int node = toAdd[i].id;
                if (!inComponentBool[node]) {
                    inComponentBool[node] = true;
                    obj += instance.myPrizes[node];
                }
                currentSize++;
                if (currentSize >= instance.params.cardCons)
                    break;
            }
        } while (toAdd.size() > 0);

        // cout<<obj<<" "<<currentSize<<endl;

        // cout<<toAdd.size()<<" "<<obj<<endl;

        if (obj > incumbentObj) {
            incumbentObj = obj;
            improved = true;
            for (int i = 0; i < instance.nNodes; ++i) {
                incumbent[i] = inComponentBool[i];
            }

            // cout<<"improved "<<obj<<endl;
        }
    }
    return improved;
}

int SolverCardinality::lagrangianPegging() {
    int nFixed = 0;
    double boundPegging = 0;
    // cout<<"pegging"<<endl;

    vector<int> fixToZero;
    vector<int> fixToOne;

    for (int i = 0; i < instance.nNodes; ++i) {
        if (fixedToZero[i] || fixedToOne[i] ||
            fabs(realPrizes[i] - epsOpt) < epsOpt)
            continue;

        //	cout<<i<<" "<<fixedToZero[i]<<" "<<fixedToOne[i]<<endl;

        if (!currentSolution[i]) {
            boundPegging = currentBound + realPrizes[i] - weightLast;
            // cout<<"boundPegging"<<boundPegging<<" "<<incumbentObj<<"
            // "<<realPrizes[i]<<" "<<currentSolution[i]<<" "<<i<<endl;

            if (boundPegging < incumbentObj) {
                fixToZero.push_back(i);
                nFixed++;
            }
        } else if (currentSolution[i]) {
            boundPegging = currentBound - realPrizes[i] + weightOutside;
            // cout<<"boundPegging"<<boundPegging<<endl;
            if (boundPegging < incumbentObj) {
                fixToOne.push_back(i);
                nFixed++;
            }
        }
    }
    // cout<<fixToZero.size()<<" "<<fixToOne.size()<<endl;
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

        // cout<<i<<endl;
    }

    for (int i : fixToOne) {
        fixedToOne[i] = true;
        instance.nFixedOne++;
        // cout<<i<<endl;
    }

    return nFixed;
}

double SolverCardinality::calculateCurrentSolution(bool save) {

    savedObj = calculateReducedCosts();

    double obj = savedObj;

    vector<nodevaluepair> myNV;
    int solutionSize = 0;

    int nFixed = 0;
    for (int i = 0; i < instance.nNodes; ++i) {
        currentSolution[i] = 0;

        if (fixedToOne[i]) {
            // nFixed++;
            currentSolution[i] = 1.0;
            obj += realPrizes[i];
            solutionSize++;
            if (save && currentSolution[i])
                sumSolution[i]++;
            if (solutionSize > instance.params.cardCons) {
                cout << "strange" << endl;
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

    // cout<<nFixed<<endl;
    sort(myNV.begin(), myNV.end(), std::greater<nodevaluepair>());

    // cout<<"obj before"<<obj<<endl;

    unsigned i = 0;
    for (; i < myNV.size(); ++i) {
        nodevaluepair nv = myNV[i];
        if (realPrizes[nv.id] <= 0)
            break;
        currentSolution[nv.id] = 1.0;
        obj += realPrizes[nv.id];
        solutionSize++;
        // cout<<nv.id<<" "<<realPrizes[nv.id]<<" "<<solutionSize<<"
        // "<<instance.myPrizes[nv.id]<<endl;

        if (save && currentSolution[nv.id])
            sumSolution[nv.id]++;
        if (solutionSize >= instance.params.cardCons)
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

    // cout<<"obj "<<obj<<endl;

    return obj;
}
