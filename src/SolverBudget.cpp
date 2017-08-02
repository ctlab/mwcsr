/*
 * SolverBudget.cpp
 *
 *  Created on: Mar 29, 2015
 *      Author: markus
 */

#include "utility/ProgramOptions.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <solverLag/SolverBudget.h>
#include <stack>
//#include "combo.h"

SolverBudget::SolverBudget(Instance &_instance, int _maxIterations)
    : SolverLag(_instance, _maxIterations), M{vector<vector<double>>(
                                                instance.nNodes)} {

    for (int j = 0; j < instance.nNodes; j++) {
        if (instance.trueTerminals[j]) {
            //	cout<<j<<endl;
            fixedToOne[j] = true;
        }

        if (instance.myBudgetCost[j] > instance.budget) {
            fixedToZero[j] = true;
        }
    }

    for (int j = 0; j < instance.nNodes; j++) {
        M[j] = vector<double>(instance.budget + 1, 0);
    }
}

SolverBudget::~SolverBudget() {}

int SolverBudget::addInitCuts() {
    int nCuts = 0;
    // return 0;
    for (int i = 0; i < instance.nNodes; ++i) {
        cut myCut;
        // myCut.rhsConst=0;
        myCut.lambda = 0;
        myCut.rhs = vector<nodevaluepair>();
        nodevaluepair n;
        n.id = i;
        n.value = 1.0;
        if (instance.myPrizes[i] <= 0)
            n.value = 2.0;

        myCut.rhs.push_back(n);

        myCut.lhs = vector<nodevaluepair>();
        for (int j : instance.adjList[i]) {
            nodevaluepair n;
            n.id = j;
            n.value = 1.0;
            myCut.lhs.push_back(n);
            // cout<<j<<" ";
        }
        // cout<<endl;
        sort(myCut.lhs.begin(), myCut.lhs.end());
        // myCut.lambda=1e-5;
        // myCut.frozen=false;
        myCut.myHash = myCut.myHasher(myCut.lhs, myCut.rhs, instance.nNodes);
        myCutHash.insert(myCut.myHash);
        // myCut.myHash=myCut.myHasher(myCut.lhs);
        myCuts.push_back(myCut);
        nCuts++;
    }

    return nCuts;
}

double SolverBudget::calculateCurrentSolution(bool saveSol) {
    myBound = calculateReducedCosts();
    // cout<<myBound<<endl;
    // incumbentObj=180;
    // cout<<"myBound"<<myBound<<endl;

    /*
    for(int n=0;n<instance.nNodes;++n)
    {
            cout<<n<<" "<<realPrizes[n]<<" "<<instance.myPrizes[n]<<endl;
    }*/

    long C = instance.budget;
    /*
    item* myItems=new item[instance.nNodes];
    for(int j = 0; j <instance.nNodes; j++)
    {
            myItems->p=realPrizes[j];
            myItems->w=instance.myBudgetCost[j];
    }
    long myCobj=combo(&myItems[0],&myItems[instance.nNodes], C, 0, 99999999, 1,
    0);

    cout<<myCobj<<endl;*/
    vector<int> jf = vector<int>(instance.nNodes, -1);

    int jc = 0;
    for (int j = 0; j < instance.nNodes; j++) {
        // fill(M[j].begin(),M[j].end(),0.0);
        // if(saveSol)
        currentSolution[j] = 0;
        if (fixedToOne[j]) {
            // if(saveSol)
            currentSolution[j] = 1;
            C -= instance.myBudgetCost[j];
            myBound += realPrizes[j];
            // cout<<j<<endl;
            continue;
        } else if (realPrizes[j] <= 0.0 || fixedToZero[j]) {
            continue;
        }

        else if (realPrizes[j] > 0.0 && instance.myBudgetCost[j] == 0) {
            currentSolution[j] = 1;
            myBound += realPrizes[j];
            continue;
        }
        jf[jc] = j;
        ++jc;
    }

    // cout<<"myBound"<<myBound<<" budget left "<<C<<" "<<jc<<endl;

    int trueNodes = jc;

    /*
    for(int j = 0; j <trueNodes; j++)
    {
            cout<<realPrizes[jf[j]]<<"y("<<jf[j]<<")+ ";
    }
    cout<<endl;

    for(int j = 0; j <trueNodes; j++)
    {
            cout<<instance.myBudgetCost[jf[j]]<<"y("<<jf[j]<<")+ ";
    }
    cout<<endl;*/
    // cout<<C<<endl;
    vector<double> myEntries = vector<double>((C + 1) * trueNodes, 0);
    for (int j = 0; j < trueNodes; j++) {
        for (int i = 1; i <= C; i++) {
            if (j > 0) {
                myEntries[trueNodes * i + j] = myEntries[trueNodes * i + j - 1];
                if (instance.myBudgetCost[jf[j]] <= i) {
                    double value =
                        myEntries[trueNodes *
                                      (i - instance.myBudgetCost[jf[j]]) +
                                  j - 1] +
                        realPrizes[jf[j]];
                    if (value > myEntries[trueNodes * i + j])
                        myEntries[trueNodes * i + j] = value;
                }
            } else {
                if (instance.myBudgetCost[jf[j]] <= i) {
                    myEntries[trueNodes * i + j] = realPrizes[jf[j]];
                }
            }
        }
    }

    int j2 = trueNodes - 1;
    int weight2 = C;
    double test = 0;

    while (j2 >= 0 && weight2 >= 0) {
        if ((j2 == 0 && myEntries[j2 + weight2 * trueNodes] > 0) ||
            (j2 > 0 &&
             myEntries[j2 + weight2 * trueNodes] !=
                 myEntries[j2 - 1 + weight2 * trueNodes])) {
            // cout<<j<<endl;
            currentSolution[jf[j2]] = 1;

            if (saveSol) {
                sumSolution[jf[j2]]++;
            }
            test += realPrizes[jf[j2]];

            // test+=realPrizes[jf[j2]];
            // cout<<j<<" "<<jf[j]<<" "<<realPrizes[jf[j]]<<"
            // "<<instance.myBudgetCost[jf[j]]<<" "<<weight<<endl;
            weight2 = weight2 - instance.myBudgetCost[jf[j2]];
        }
        --j2;
    }

    // cout<<myBound<<" "<<myEntries[myEntries.size()-1]<<" "<<test<<" "<<C<<"
    // "<<weight2<<endl;

    /*

            for(int i = 1; i <= C; i++){
                    //cout<<i<<endl;
                    for(int j = 0; j <trueNodes; j++){
                            if(j > 0){
                                    M[j][i] = M[j-1][i];
                                    if (instance.myBudgetCost[jf[j]] <= i)
                                            M[j][i] = max(M[j][i],
       M[j-1][i-instance.myBudgetCost[jf[j]]]+realPrizes[jf[j]]);
                            }
                            else{
                                    M[j][i] = 0;
                                    if(instance.myBudgetCost[jf[j]] <= i)
                                            M[j][i] = max(M[j][i],
       realPrizes[jf[j]]);
                            }
                    }
            }

            int j=trueNodes-1;
            int weight=C;

            test=0;
            while(j>=0 && weight>=0)
            {
                    if ((j==0 && M[j][weight]>0) || (j>0 && M[j][weight] !=
       M[j-1][weight]))
                    {
                            //cout<<j<<endl;
                            currentSolution[jf[j]]=1;

                            if(saveSol)
                            {
                                    sumSolution[jf[j]]++;
                            }
                            test+=realPrizes[jf[j]];
                            //cout<<j<<" "<<jf[j]<<" "<<realPrizes[jf[j]]<<"
       "<<instance.myBudgetCost[jf[j]]<<" "<<weight<<endl;
                            weight = weight-instance.myBudgetCost[jf[j]];
                    }
                    --j;
            }
            cout<<myEntries[myEntries.size()-1]<<" "<<M[trueNodes-1][C]<<endl;*/

    /*
    double largest=0;
    for(double e:myEntries)
    {
            if(e>largest)
                    largest=e;
    }
    cout<<largest<<endl;*/

    // cout<<myEntries[myEntries.size()-1]<<endl;
    myBound += myEntries[myEntries.size() - 1];
    // cout<<"test "<<M[trueNodes-1][C]<<" "<<test<<endl;
    return myBound;
}

int SolverBudget::lagrangianPegging() {
    int nFixed = 0;
    vector<nodevaluepair> myItems;

    double objHeur = 0.0;

    for (int i = 0; i < instance.nNodes; ++i) {
        // fill(M[j].begin(),M[j].end(),0.0);
        // if(saveSol)
        // cout<<lagrangePrizes.cost[n] <<endl;

        if (instance.myBudgetCost[i] == 0 && realPrizes[i] >= 0.0) {
            objHeur += realPrizes[i];
            continue;
        }

        if (realPrizes[i] <= 0.0 || fixedToZero[i]) {
            continue;
        }

        nodevaluepair nv;
        nv.id = i;
        nv.value = realPrizes[i] / (instance.myBudgetCost[i]);
        myItems.push_back(nv);
    }

    sort(myItems.begin(), myItems.end(), std::greater<nodevaluepair>());
    double budgetIn = 0.0;

    vector<double> inItems = vector<double>(instance.nNodes, 0.0);
    vector<int> order;

    for (nodevaluepair nv : myItems) {
        // cout<<budgetIn<<" "<<objHeur<<" "<<nv.value<<endl;
        int n = nv.id;
        if (budgetIn + instance.myBudgetCost[n] <= instance.budget) {
            budgetIn += instance.myBudgetCost[n];
            objHeur += realPrizes[n];
            inItems[n] = 1;
            //	cout<<instance.myBudgetCost[n]<<" "<<n<<endl;
            order.push_back(n);
        } else {
            double remainder = instance.budget - budgetIn;
            double fraction = remainder / instance.myBudgetCost[n];
            budgetIn += remainder;
            //	cout<<remainder<<" "<<fraction<<" "<<instance.myCost[n]<<endl;
            objHeur += realPrizes[n] * fraction;
            // cout<<instance.myBudgetCost[n]<<" "<<n<<" "<<fraction<<endl;
            inItems[n] = fraction;
            order.push_back(n);
            break;
        }
    }

    //	cout<<budgetIn<<endl;
    // cout<<objHeur<<" "<<objKnapsack<<" "<<objSpanning<<"
    // "<<currentBound<<endl;
    for (int n = 0; n < instance.nNodes; ++n) {
        if (inItems[n] > 0)
            continue;

        if (fixedToZero[n] > 0)
            continue;

        if (instance.myBudgetCost[n] == 0)
            continue;

        // cout<<" "<<objHeur<<" "<<realPrizes[n]<<endl;
        double addedObj = objHeur + realPrizes[n];
        // cout<<addedObj<<endl;
        double additionalCost = instance.myBudgetCost[n];
        int myIndex = 1;
        int last = order[order.size() - myIndex];

        while (additionalCost > 0) {
            // cout<<last<<" "<<myIndex<<" "<<order.size()<<"
            // "<<additionalCost<<" "<<instance.myBudgetCost[last]<<"
            // "<<inItems[last]<<" "<<addedObj<<" "<<realPrizes[last]<<endl;

            if (additionalCost - instance.myBudgetCost[last] * inItems[last] >=
                0) {
                addedObj -= realPrizes[last] * inItems[last];
                additionalCost -= instance.myBudgetCost[last] * inItems[last];
            } else {
                double fraction = additionalCost / instance.myBudgetCost[last];
                addedObj -= realPrizes[last] * fraction;
                additionalCost -= instance.myBudgetCost[last] * fraction;
                break;
            }
            myIndex++;
            last = order[order.size() - myIndex];
            //			cout<<instance.myBudgetCost[last]<<endl;
            //			cout<<inItems[last]<<endl;
        }
        if (addedObj > objHeur + 0.0001) {
            cout << "bug " << addedObj << " " << objHeur << " "
                 << additionalCost << endl;
            exit(1);
        }

        // cout<<addedObj<<" "<<incumbentObj<<endl;

        if (myBound + addedObj < incumbentObj) {
            cout << "fixing " << realPrizes[n] << " " << myBound + addedObj
                 << " " << incumbentObj << " " << myBound + objHeur << endl;
            fixedToZero[n] = true;
            instance.nFixedZero++;
            nFixed++;
        }
    }
    return nFixed;
}

bool SolverBudget::primalHeuristic() {

    bool improvement = false;

    unsigned iter = 3;
    // if (myComponents.size()<iter)
    //	iter=myComponents.size();

    // if(iter==0)
    iter = 1;

    vector<int> myHeurTerminals;
    vector<bool> myHeurTerminalsBool = vector<bool>(instance.nNodes, false);

    vector<double> lpSolution;

    for (int i = 0; i < instance.nNodes; ++i) {
        if (currentSolution[i]) {
            lpSolution.push_back(1);
            if (instance.realTerminals[i]) {
                myHeurTerminalsBool[i] = true;
                myHeurTerminals.push_back(i);
            }
        } else {
            lpSolution.push_back(0);
        }
    }

    for (unsigned myIter = 0; myIter < iter; ++myIter) {
        int root = 0;

        double incObj = 0;
        double budget = 0;

        double bestVal = 0;

        if (myComponents.size() > 0) {
            for (unsigned i = 0; i < myComponents[0].components.size(); ++i) {
                int c = myComponents[0].components[i];

                if (myHeurTerminalsBool[c] && realPrizes[c] > bestVal) {
                    root = c;
                    bestVal = realPrizes[c];
                    // inComponentBool[c]=true;
                }
            }
        } else {
            for (int i = 0; i < instance.nNodes; ++i) {
                int c = i;
                if (myHeurTerminalsBool[c] && realPrizes[c] > bestVal) {
                    root = c;
                    bestVal = realPrizes[c];
                }
            }
        }

        vector<int> inComponentBool = vector<int>(instance.nNodes, false);

        priority_queue<nodevaluepair, vector<nodevaluepair>,
                       greater<nodevaluepair>>
            myPQueue;
        vector<int> pred = vector<int>(instance.nNodes, -1);
        vector<bool> extracted = vector<bool>(instance.nNodes, false);
        vector<double> distance =
            vector<double>(instance.nNodes, std::numeric_limits<double>::max());
        vector<int> myBestSol =
            vector<int>(instance.nNodes, instance.nNodes + 1);

        incObj += instance.myPrizes[root];
        budget += instance.myBudgetCost[root];
        extracted[root] = false;
        myBestSol[root] = 0;
        nodevaluepair nv;
        nv.id = root;
        nv.value = 0;
        myPQueue.push(nv);
        pred[root] = root;
        inComponentBool[root] = true;
        distance[root] = 0;
        // int badAdd=0;

        int numTerms = 0;
        int numBest = 0;
        int bestBudget = 0;
        double best = 0;

        while (!myPQueue.empty()) {
            nodevaluepair nv = myPQueue.top();
            myPQueue.pop();

            int k = nv.id;
            // cout<<k<<" "<<distance[k]<<budget<<" "<<extracted[k]<<"
            // "<<pred[k]<<endl;

            if (extracted[k])
                continue;

            extracted[k] = true;
            if (!myHeurTerminalsBool[k] || inComponentBool[k]) {
                double toAdd = -instance.myPrizes[k] * (1 - lpSolution[k]);
                if (toAdd <= 0)
                    toAdd = epsOpt;
                // cout<<toAdd<<endl;

                for (int l : instance.adjList[k]) {
                    if (fixedToZero[l])
                        continue;

                    if (pred[l] == k) {
                        // cout<<l<<" "<<pred[l]<<" "<<k<<endl;
                        continue;
                    }

                    if (budget + instance.myBudgetCost[l] > instance.budget)
                        continue;
                    // cout<<distance[l]<<" "<< distance[k]<<"
                    // "<<instance.myBudgetCost[l]<<" "<<lpSolution[l]<<"
                    // "<<l<<" "<<k<<endl;
                    if (distance[l] >
                        distance[k] +
                            instance.myBudgetCost[l] * (1 - lpSolution[l]) +
                            epsOpt)
                    // if (distance[l] > distance[k] +
                    // toAdd*(1-lpSolution[l])+epsOpt)
                    {
                        // cerr<<distance[l]<<" "<<distance[k]+
                        // instance.myBudgetCost[l]*(1-lpSolution[l])+epsOpt<<"
                        // "<<l<<" "<<k<<endl;

                        distance[l] =
                            distance[k] +
                            instance.myBudgetCost[l] * (1 - lpSolution[l]) +
                            epsOpt;
                        // distance[l] = distance[k] +
                        // toAdd*(1-lpSolution[l])+epsOpt;
                        pred[l] = k;

                        extracted[l] = false;
                        nodevaluepair lNv;
                        lNv.id = l;
                        lNv.value = distance[l];
                        myPQueue.push(lNv);
                    }
                }
            } else {
                // cout<<"k "<<k<<endl;
                int currentNode = k;

                bool fit = true;
                double add = 0;
                while (!inComponentBool[currentNode]) {
                    // cout<<currentNode<<"
                    // "<<instance.myBudgetCost[currentNode]<<"
                    // "<<pred[currentNode]<<endl;

                    add += instance.myBudgetCost[currentNode];
                    if (budget + add > instance.budget) {
                        fit = false;
                        break;
                    }
                    currentNode = pred[currentNode];
                }

                currentNode = k;
                // cout<<fit<<" "<<k<<" "<<budget<<" "<<add<<endl;

                if (fit) {
                    while (!inComponentBool[currentNode]) {
                        // cout<<currentNode<<endl;
                        inComponentBool[currentNode] = true;
                        incObj += instance.myPrizes[currentNode];
                        budget += instance.myBudgetCost[currentNode];
                        //	cout<<"incObj "<<incObj<<endl;
                        myBestSol[currentNode] = numTerms;
                        distance[currentNode] = 0;
                        extracted[currentNode] = false;
                        nodevaluepair nv;
                        nv.id = currentNode;
                        // cout<<"dist2 "<<distance[currentNode]<<endl;
                        nv.value = distance[currentNode];
                        myPQueue.push(nv);
                        currentNode = pred[currentNode];
                    }

                    numTerms++;
                    if (incObj > best) {
                        best = incObj;
                        numBest = numTerms;
                        bestBudget = budget;
                        // cout<<"best"<<obj<<" "<<currentSize<<endl;
                    }

                    // cout<<incObj<<" "<<budget<<" "<<instance.budget<<endl;
                }
            }
        }

        // cout<<best<<" "<<incObj<<" "<<numBest<<" "<<bestBudget<<endl;

        for (int i = 0; i < instance.nNodes; ++i) {
            inComponentBool[i] = 0;
            if (myBestSol[i] < numBest)
                inComponentBool[i] = 1;
            // cout<<myBestSol[i]<<" "<<numBest<<endl;
        }
        budget = bestBudget;

        incObj = best;

        vector<nodevaluepair> toAdd;

        do {
            toAdd.clear();

            // cout<<"before postprocessing "<<obj<<" "<<currentSize<<endl;

            for (int n = 0; n < instance.nNodes; ++n) {
                if (budget >= instance.budget)
                    break;
                if (!inComponentBool[n])
                    continue;
                for (int j : instance.adjList[n]) {
                    if (!inComponentBool[j] && instance.myPrizes[j] > 0 &&
                        budget + instance.myBudgetCost[j] <= instance.budget) {
                        nodevaluepair nv;
                        nv.id = j;
                        // nv.value=instance.myPrizes[j];
                        nv.value =
                            instance.myPrizes[j] / instance.myBudgetCost[j];
                        toAdd.push_back(nv);
                    }
                }
            }

            sort(toAdd.begin(), toAdd.end(), std::greater<nodevaluepair>());

            // cout<<"before"<<obj<<" "<<currentSize<<endl;
            // cout<<"before"<<obj<<endl;

            for (unsigned i = 0; i < toAdd.size(); ++i) {
                // cout<<"here"<<endl;
                int node = toAdd[i].id;

                if (budget + instance.myBudgetCost[node] > instance.budget) {
                    continue;
                }
                if (!inComponentBool[node]) {
                    inComponentBool[node] = true;
                    incObj += instance.myPrizes[node];
                    budget += instance.myBudgetCost[node];
                }
                if (budget >= instance.budget)
                    break;
            }
        } while (toAdd.size() > 0);

        bool foundAll = true;
        for (int i : instance.myTrueTerminals) {
            if (!inComponentBool[i]) {
                foundAll = false;
                break;
            }
        }

        if (incObj > incumbentObj && foundAll) {
            incumbentObj = incObj;
            for (int i = 0; i < instance.nNodes; ++i) {
                incumbent[i] = inComponentBool[i];
                if (incumbent[i])
                    cout << i << endl;
                // test+=instance.myPrizes[i]*incumbent[i];
            }

            improvement = true;
        }
    }

    return improvement;
}

#ifdef TEST
bool SolverBudget::primalHeuristic() {
    bool improvement = false;
    int root = instance.myTrueTerminals[0];

    vector<double> lpSolution = vector<double>(instance.nNodes, 0);
    for (int i = 0; i < instance.nNodes; ++i) {
        lpSolution[i] = sumSolution[i] / double(iterations + 1);
    }

    queue<int> myQueue;
    // myQueue.push(root);
    vector<bool> inComponent = vector<bool>(instance.nNodes, false);
    vector<int> myComponent;
    // myComponent.push_back(root);
    // inComponent[root]=true;

    double incObj = instance.myPrizes[root];
    double budget = instance.myBudgetCost[root];

    int toFind = instance.myTrueTerminals.size() - 1;

    // prim-i to connect terminals
    priority_queue<nodevaluepair, vector<nodevaluepair>, greater<nodevaluepair>>
        myPQueue;
    vector<int> pred = vector<int>(instance.nNodes, -1);
    vector<bool> extracted = vector<bool>(instance.nNodes, false);
    vector<double> distance =
        vector<double>(instance.nNodes, std::numeric_limits<double>::max());
    extracted[root] = false;
    nodevaluepair nv;
    nv.id = root;
    nv.value = 0;
    myPQueue.push(nv);
    myQueue.push(root);
    pred[root] = root;
    inComponent[root] = true;
    distance[root] = 0;
    // int badAdd=0;

    if (toFind > 0) {
        while (!myPQueue.empty() && toFind) {
            nodevaluepair nv = myPQueue.top();
            myPQueue.pop();

            int k = nv.id;
            // cout<<k<<" "<<distance[k]<<budget<<" "<<extracted[k]<<"
            // "<<pred[k]<<endl;

            if (extracted[k])
                continue;

            extracted[k] = true;
            if (!(instance.trueTerminals[k]) || inComponent[k]) {
                for (int l : instance.adjList[k]) {
                    if (fixedToZero[l])
                        continue;

                    if (pred[l] == k) {
                        // cout<<l<<" "<<pred[l]<<" "<<k<<endl;
                        continue;
                    }

                    if (budget + instance.myBudgetCost[l] > instance.budget)
                        continue;
                    // cout<<distance[l]<<" "<< distance[k]<<"
                    // "<<instance.myBudgetCost[l]<<" "<<lpSolution[l]<<"
                    // "<<l<<" "<<k<<endl;
                    if (distance[l] >
                        distance[k] +
                            instance.myBudgetCost[l] * (1 - lpSolution[l]) +
                            epsOpt) {
                        // cerr<<distance[l]<<" "<<distance[k]+
                        // instance.myBudgetCost[l]*(1-lpSolution[l])+epsOpt<<"
                        // "<<l<<" "<<k<<endl;

                        distance[l] =
                            distance[k] +
                            instance.myBudgetCost[l] * (1 - lpSolution[l]) +
                            epsOpt;
                        pred[l] = k;

                        extracted[l] = false;
                        nodevaluepair lNv;
                        lNv.id = l;
                        lNv.value = distance[l];
                        myPQueue.push(lNv);
                    }
                }
            } else if (instance.trueTerminals[k]) {
                // cout<<"k "<<k<<endl;
                int currentNode = k;

                bool fit = true;
                double add = 0;
                while (!inComponent[currentNode]) {
                    // cout<<currentNode<<"
                    // "<<instance.myBudgetCost[currentNode]<<"
                    // "<<pred[currentNode]<<endl;

                    add += instance.myBudgetCost[currentNode];
                    if (budget + add > instance.budget) {
                        fit = false;
                        // break;
                    }
                    currentNode = pred[currentNode];
                }

                currentNode = k;
                // cout<<fit<<" "<<k<<" "<<budget<<" "<<add<<endl;

                if (fit) {
                    toFind--;
                    while (!inComponent[currentNode]) {
                        // cout<<currentNode<<endl;
                        inComponent[currentNode] = true;
                        incObj += instance.myPrizes[currentNode];
                        budget += instance.myBudgetCost[currentNode];
                        //	cout<<"incObj "<<incObj<<endl;

                        distance[currentNode] = 0;
                        extracted[currentNode] = false;
                        nodevaluepair nv;
                        nv.id = currentNode;
                        // cout<<"dist2 "<<distance[currentNode]<<endl;
                        nv.value = distance[currentNode];
                        myPQueue.push(nv);
                        myQueue.push(currentNode);
                        currentNode = pred[currentNode];
                        // cout<<currentNode<<endl;
                        // cout<<"pushed"<<endl;
                    }
                }
                /*else
                {
                        //return false;
                        //badAdd++;
                        //if(badAdd==100)
                        //	break;
                        add=0;
                        while(!inComponent[currentNode])
                        {

                                add+=instance.myBudgetCost[currentNode];
                                distance[currentNode]=std::numeric_limits<double>::max();
                                extracted[currentNode]=true;
                                //cout<<currentNode<<"
                "<<extracted[currentNode]<<endl;


                                extracted[currentNode]=false;
                                nodevaluepair nv;
                                nv.id=currentNode;
                                nv.value=distance[currentNode];
                                myPQueue.push(nv);
                                currentNode=pred[currentNode];

                                if(budget+add>instance.budget)
                                {
                                        break;
                                }
                        }
                }*/
            }
        }
    }

    // exit(1);

    // cout<<toFind<<" "<<budget<<" "<<incObj<<endl;
    if (toFind > 0) {
        return false;
    }

    priority_queue<nodevaluepair, vector<nodevaluepair>, less<nodevaluepair>>
        myPQueue2;
    while (!myQueue.empty()) {
        int k = myQueue.front();
        myQueue.pop();
        // forall_adj_edges(e,k)
        for (int li : instance.adjList[k]) {
            if (budget + instance.myBudgetCost[li] > instance.budget)
                continue;

            if (fixedToZero[li])
                continue;

            // cout<<k<<" "<<instance.myBudgetCost[li]<<endl;

            if (!inComponent[li]) {
                // cout<<li<<" "<<instance.myBudgetCost[li]<<endl;
                nodevaluepair nv;
                nv.id = li;
                nv.value = realPrizes[li] / instance.myBudgetCost[li];
                myPQueue2.push(nv);
            }
        }
    }

    while (!myPQueue2.empty()) {
        nodevaluepair nv = myPQueue2.top();
        myPQueue2.pop();

        int k = nv.id;
        // cout<<k<<" "<<distance[k]<<budget<<endl;

        if (inComponent[k])
            continue;

        if (budget + instance.myBudgetCost[k] > instance.budget)
            continue;

        inComponent[k] = true;
        incObj += instance.myPrizes[k];
        budget += instance.myBudgetCost[k];

        for (int li : instance.adjList[k]) {
            if (fixedToZero[li])
                continue;
            if (budget + instance.myBudgetCost[li] > instance.budget)
                continue;

            if (!inComponent[li]) {
                nodevaluepair nv;
                nv.id = li;
                nv.value = realPrizes[li] / instance.myBudgetCost[li];
                myPQueue2.push(nv);
            }
        }
    }

    // cout<<incObj<<" "<<toFind<<" "<<budget<<endl;

    if (toFind == 0 && incObj > incumbentObj) {
        incumbentObj = incObj;
        improvement = true;
    }

    double testCost = 0.0;
    double testBudget = 0.0;

    for (int i = 0; i < instance.nNodes; ++i) {
        if (inComponent[i]) {
            testCost += instance.myPrizes[i];
            testBudget += instance.myBudgetCost[i];
        }
    }

    // cout<<testCost<<" "<<testBudget<<endl;

    return improvement;
}
#endif
