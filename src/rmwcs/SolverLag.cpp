/*
 * SolverLag.cpp
 *
 *  Created on: Mar 29, 2015
 *      Author: markus
 */

#include <stack>
#include <chrono>
#include <queue>
#include <utility>

#include "include/Parameters.h"
#include "include/SolverLag.h"

using std::vector;
using std::list;
using std::queue;
using std::max;

using Rcpp::Rcout;

namespace chrono = std::chrono;

bool cutToRemove(SolverLag::cut x) { return x.toRemove; }

SolverLag::SolverLag(Instance& instance, Parameters& params, mwcsr::monitor int_monitor)
        : instance(instance), params(params), int_monitor(std::move(int_monitor)), myCuts{list<cut>()},
          myNewCuts{list<cut>()},
          realPrizes{vector<double>(instance.nNodes, 0)},
          currentSolution{vector<double>(instance.nNodes, 0)},
          previousSolution{vector<double>(instance.nNodes, 0)}, sumSolution{vector<int>(instance.nNodes, 0)},
          incumbent{vector<bool>(
                instance.nNodes, 0)}, dualIncumbent{vector<int>(instance.nNodes, 0)},
          labels{vector<int>(
                instance.nNodes, 0)}, fixedToZero{vector<int>(instance.nNodes, 0)},
          fixedToOne{vector<int>(
                instance.nNodes, 0)}, incumbentObj{instance.incumbentObjLag}, subgradientSquared{1},
          subgradientNorm{0}, directionPrevSquared{0}, alpha{params.beta},
          noImprov{0}, numberOfComponents{0},
          bestBound{std::numeric_limits<double>::max()},
          currentBound{std::numeric_limits<double>::max()},
          previousBound{std::numeric_limits<double>::max()},
          bestBoundCFT{std::numeric_limits<double>::max()}, worstBoundCFT{std::numeric_limits<double>::lowest()},
          counterCFT{0}, maxIterations{params.maxIter}, iterations{0},
          sepIter(params.sepIter), sepIterFreeze(params.sepIterFreeze),
          savedObj(0.0), runtime(0.0) {}

SolverLag::~SolverLag() {}

int SolverLag::solve() {
    solveSubgradient(maxIterations);
    if (params.solver == 1) {
        writeFixingToInstance();
        writeSolutionToInstance();
    }
    return 1;
}

int SolverLag::writeCutsToInstance() {
    int nCuts = 0;
    instance.myCuts.clear();
    for (cut c : myCuts) {
        Instance::cut myCut;
        if (c.lambda > 0) {
            for (nodevaluepair i : c.lhs) {
                myCut.lhs.push_back(i.id);
            }
            for (nodevaluepair i : c.rhs) {
                myCut.rhs.push_back(i.id);
            }
            myCut.rhsConst = c.rhsConst;
            nCuts++;
            instance.myCuts.push_back(myCut);
        }
    }

    return nCuts;
}

int SolverLag::writeFixingToInstance() {

    for (unsigned i = 0; i < instance.nNodes; ++i) {
        instance.fixedToOne[i] = fixedToOne[i];
        instance.fixedToZero[i] = fixedToZero[i];
    }

    return 1;
}

int SolverLag::writeSolutionToInstance() {
    instance.incumbent = vector<bool>(instance.nNodes, false);

    for (unsigned i = 0; i < instance.nNodes; ++i) {
        instance.incumbent[i] = incumbent[i];
    }
    instance.incumbentFound = true;

    return 0;
}

int SolverLag::solveSubgradient(int maxIterations) {
    chrono::time_point <std::chrono::system_clock> startTime = chrono::system_clock::now();
    iterations = 0;

    double eps = 1e-10;

    int nFixed = 0;

    if (params.integer) {
        eps = 0;
    }

    bool boundImprov = false;
    addInitCuts();

    while (iterations < maxIterations && sqrt(subgradientSquared) > epsOpt) {
        int_monitor.check();
        boundImprov = false;
        subgradientSquared = 0.0;
        currentBound = calculateCurrentSolution(true);

        if (currentBound < bestBound) {
            bestBound = currentBound;
            boundImprov = true;
            if (params.solver == 1) {
                writeCutsToInstance();
            }
        }

        if (params.subgradient != 2) {
            if (boundImprov) {
                noImprov = 0;
                bestBound = currentBound;

                for (unsigned i = 0; i < instance.nNodes; ++i) {
                    dualIncumbent[i] = currentSolution[i];
                }
            } else
                noImprov++;
        }

        if (params.subgradient == 2) {
            if (boundImprov) {
                bestBound = currentBound;
                for (cut& c : myCuts) {
                    c.bestLambda = c.lambda;
                }
            }

            double delta = bestBound - incumbentObj;
            if (bestBound + delta < currentBound) {
                noImprov++;
            } else {
                noImprov = 0;
            }
        }

        double bestBoundCheck = bestBound;
        if (params.integer) {
            bestBoundCheck = floor(bestBoundCheck);
        }

        int numberOfCuts = createCuts(iterations);

        if (iterations % params.heurIter == 0) {
            primalHeuristic();
        }

        nFixed = 0;
        if ((boundImprov) && params.pegging > 0) {
            for (unsigned c = 0; c < instance.components.size(); ++c) {
                if (instance.maxRevenueInComponent[c] < incumbentObj &&
                    !instance.componentFixed[c]) {
                    instance.componentFixed[c] = true;
                    for (int j : instance.components[c]) {
                        if (!fixedToZero[j]) {
                            fixedToZero[j] = 1;
                            instance.nFixedZero++;
                        }
                    }
                }
            }

            nFixed = lagrangianPegging();
        }

        if (bestBoundCheck <= incumbentObj + eps) {
            break;
        }

        if (params.outputlag) {
            Rcout << std::setprecision(9) << "iteration: \t" << iterations
                  << "\t lagrangian bound: \t"
                  << instance.transformInternalValue(bestBound)
                  << "\t current bound: \t "
                  << instance.transformInternalValue(currentBound)
                  << "\t incumbent: \t "
                  << instance.transformInternalValue(incumbentObj)
                  << "\t number of active cuts: \t" << numberOfCuts << "\n";
        }
        if (nFixed)
            myCuts.erase(
                    std::remove_if(myCuts.begin(), myCuts.end(), cutToRemove),
                    myCuts.end());
        upgradeMultipliers();
        for (unsigned i = 0; i < instance.nNodes; ++i) {
            previousSolution[i] = currentSolution[i];
        }

        iterations++;
    }

    if (params.outputlag) {
        Rcout << std::setprecision(9) << "iteration: \t" << iterations
              << "\t lagrangian bound: \t"
              << instance.transformInternalValue(bestBound)
              << "\t incumbent: \t "
              << instance.transformInternalValue(incumbentObj) << "\n";
    }

    chrono::time_point <std::chrono::system_clock> endTime =
            chrono::system_clock::now();
    chrono::duration<double> elapsedSeconds = endTime - startTime;

    runtime = elapsedSeconds.count();
    writeStatistics();

    return 1;
}

double SolverLag::calculateReducedCosts() {
    double obj = 0.0;
    for (unsigned n = 0; n < instance.nNodes; ++n) {
        realPrizes[n] = instance.myPrizes[n];
    }

    for (cut& c : myCuts) {
        if (c.frozen)
            continue;

        if (c.lambda == 0)
            continue;
        for (nodevaluepair n : c.lhs) {
            realPrizes[n.id] += (n.value * c.lambda);
        }
        for (nodevaluepair n : c.rhs) {
            realPrizes[n.id] -= (n.value * c.lambda);
        }
        obj += c.rhsConst * c.lambda;
    }

    return obj;
}

void SolverLag::writeStatistics() {
    instance.bestBoundLag = bestBound;
    instance.incumbentObjLag = incumbentObj;
    instance.iterationsLag = iterations;
    instance.runtimeLag = runtime;

    instance.incumbent = vector<bool>(instance.nTrueNodes, false);
    instance.solSize = 0;
    for (unsigned i = 0; i < instance.nNodes; ++i) {
        if (incumbent[i]) {
            instance.solSize++;
            instance.incumbent[instance.map[i]] = true;
        }
    }

    if (instance.gapLag < epsOpt)
        instance.gapLag = 0;
}

int SolverLag::createCuts(int iter) {
    int numberOfCuts = 0;
    int cnt = 0;

    if (iter % sepIter == 0 && iter > 0)
        cnt += separateCuts();
    if (iter % sepIterFreeze == 0) {
        for (cut& c : myCuts) {
            c.frozen = false;
        }
        numberOfCuts += cnt;
    }

    numberOfCuts += checkPreviousCuts(true);

    return numberOfCuts;
}

void SolverLag::updateMultipliersSherali() {
    if (noImprov > params.betaIter) {
        noImprov = 0;
        alpha /= 2;
        currentBound = bestBound;;
        for (unsigned n = 0; n < instance.nNodes; ++n) {
            currentSolution[n] = dualIncumbent[n];
        }
        subgradientSquared = 0;
        checkPreviousCuts(false);
        for (cut& c : myCuts) {
            if (c.frozen)
                continue;
            c.directionPrevious = 0;
        }
    }
    directionPrevSquared = 0.0;
    for (cut& c : myCuts) {
        if (c.frozen)
            continue;
        directionPrevSquared += (c.directionPrevious * c.directionPrevious);
    }

    double sigma = 0;

    if (directionPrevSquared > epsOpt) {
        sigma = sqrt(subgradientSquared) / sqrt(directionPrevSquared);
    }

    double directionSquared = 0.0;
    for (cut& c : myCuts) {
        if (c.frozen)
            continue;
        c.direction = c.subgradient + sigma * c.directionPrevious;
        c.directionPrevious = c.direction;
        directionSquared += (c.direction * c.direction);
    }

    if (directionSquared < epsOpt) {
        directionSquared = subgradientSquared;
        for (cut& c : myCuts) {
            if (c.frozen)
                continue;
            c.direction = c.subgradient;
        }
    }

    double theta = alpha * (currentBound - incumbentObj) / (directionSquared);

    for (cut& c : myCuts) {
        if (c.frozen)
            continue;
        c.lambda = std::max(0.0, c.lambda - theta * c.direction);
    }
}

void SolverLag::updateMultipliersLucena() {
    if (noImprov > params.betaIter) {
        noImprov = 0;
        alpha /= 2;
    }

    double theta = alpha * (currentBound - incumbentObj) / (subgradientSquared);

    for (cut& c : myCuts) {
        if (c.frozen)
            continue;
        c.lambda = std::max(0.0, c.lambda - theta * c.subgradient);
    }
}

void SolverLag::updateMultipliersCFT() {
    if (noImprov >= params.betaIter) {
        noImprov = 0;
        alpha /= 2;
        for (cut& c : myCuts) {
            c.lambda = c.bestLambda;
        }
    }

    double theta = alpha * (currentBound - incumbentObj) / (subgradientSquared);

    for (cut& c : myCuts) {
        if (c.frozen)
            continue;
        c.lambda = std::max(0.0, c.lambda - theta * c.subgradient);
    }
}

void SolverLag::upgradeMultipliers() {
    if (params.subgradient == 2) {
        updateMultipliersCFT();
        return;
    }

    if (params.subgradient == 0) {
        updateMultipliersLucena();
        return;
    }

    if (params.subgradient == 1) {
        updateMultipliersSherali();
        return;
    }
}

int SolverLag::checkPreviousCuts(bool addCuts) {
    int numberOfCuts = 0;

    for (cut& c : myCuts) {
        if (c.frozen)
            continue;

        c.violated = true;
        c.subgradient = calculateSubgradientCuts(c);

        double lhs = 0.0;
        unsigned nFixedToZero = 0;
        for (nodevaluepair n : c.lhs) {
            lhs += n.value * currentSolution[n.id];
            if (fixedToZero[n.id])
                nFixedToZero++;
        }

        if (nFixedToZero == c.lhs.size() && c.rhs.size() == 1 && addCuts) {
            for (nodevaluepair n : c.rhs) {
                if (!fixedToZero[n.id]) {
                    fixedToZero[n.id] = true;

                    for (int j : instance.adjList[n.id]) {
                        unsigned k = 0;
                        for (; k < instance.adjList[j].size(); ++k) {
                            if (instance.adjList[j][k] == n.id)
                                break;
                        }
                        instance.adjList[j].erase(instance.adjList[j].begin() +
                                                  k);
                    }

                    instance.adjList[n.id].clear();
                }
            }
            c.toRemove = true;
            c.subgradient = 0;
        }

        double rhs = -c.rhsConst;
        nFixedToZero = 0;

        for (nodevaluepair n : c.rhs) {
            rhs += n.value * currentSolution[n.id];
            if (fixedToZero[n.id])
                nFixedToZero++;
        }

        if (nFixedToZero == c.rhs.size() + c.rhsConst && addCuts) {
            c.subgradient = 0;
            c.toRemove = true;
        }

        if (lhs < rhs) {
            numberOfCuts++;
            c.age = 0;
        } else {
            c.violated = false;
            c.age++;

            if (c.lambda == 0 && c.subgradient > 0 && c.age > params.maxAge) {
                c.subgradient = 0;
            }
        }
        subgradientSquared += (c.subgradient * c.subgradient);
    }

    if (addCuts) {
        for (cut& c : myNewCuts) {
            if (!c.frozen) {
                c.subgradient = calculateSubgradientCuts(c);
                c.directionPrevious = calculateSubgradientCutsPrevious(c);
                subgradientSquared += (c.subgradient * c.subgradient);
            }
            myCuts.push_back(c);
        }
    }

    return numberOfCuts;
}

double SolverLag::calculateSubgradientCuts(const cut& myCut) {
    double subg = myCut.rhsConst;
    for (nodevaluepair n : myCut.lhs) {
        subg += (n.value * currentSolution[n.id]);
    }

    for (nodevaluepair n : myCut.rhs) {
        subg -= (n.value * currentSolution[n.id]);
    }

    return subg;
}

double SolverLag::calculateSubgradientCutsPrevious(const cut& myCut) {
    double subg = myCut.rhsConst;
    for (nodevaluepair n : myCut.lhs) {
        subg += (n.value * previousSolution[n.id]);
    }

    for (nodevaluepair n : myCut.rhs) {
        subg -= (n.value * previousSolution[n.id]);
    }

    return subg;
}

int SolverLag::setVariableFixing(const vector<int>& toZero,
                                 const vector<int>& toOne) {
    int numberFixed = toZero.size() + toOne.size();

    for (unsigned i = 0; i < toZero.size(); ++i) {
        fixedToZero[toZero[i]] = true;
    }

    for (unsigned i = 0; i < toOne.size(); ++i) {
        fixedToOne[toOne[i]] = true;
    }

    return numberFixed;
}

void SolverLag::initCuts(list<SolverLag::cut>& cuts) {
    myCuts = cuts;
    for (cut& c : myCuts) {
        c.subgradient = 0;
        c.direction = 0;
        c.directionPrevious = 0;
        c.age = 0;
    }
}

int SolverLag::separateCuts() {
    int numberOfCuts = 0;
    myNewCuts.clear();

    int myCurrentLabel = 1;
    std::fill(labels.begin(), labels.end(), 0);
    myComponents.clear();

    for (unsigned n = 0; n < instance.nNodes; ++n) {
        if (currentSolution[n] == 1 && instance.realTerminals[n] &&
            labels[n] == 0) {

            CompStruct myComponentHelper;
            vector<bool> componentBoundary =
                    vector<bool>(instance.nNodes, false);
            vector<int> componentBoundaryIndexed;
            vector<bool> inComponent = vector<bool>(instance.nNodes, false);
            vector<int> myComponent;
            myComponent.push_back(n);
            inComponent[n] = true;
            myComponentHelper.sumPrize += instance.myPrizes[n];
            labels[n] = myCurrentLabel;

            queue<int> myQueue;
            myQueue.push(n);
            while (!myQueue.empty()) {
                int k = myQueue.front();

                myQueue.pop();
                for (int li : instance.adjList[k]) {
                    if (currentSolution[li] == 0) {
                        {
                            if (!componentBoundary[li]) {
                                componentBoundaryIndexed.push_back(li);
                                componentBoundary[li] = true;
                            }
                        }
                    } else if (labels[li] != myCurrentLabel) {
                        myComponent.push_back(li);
                        labels[li] = myCurrentLabel;
                        inComponent[li] = true;
                        myQueue.push(li);
                        if (instance.myPrizes[li] > 0) {
                            myComponentHelper.sumPrize += instance.myPrizes[li];
                        }
                    }
                }
            }
            myComponentHelper.boundaryIndexed = componentBoundaryIndexed;
            myComponentHelper.boundary = componentBoundary;
            myComponentHelper.components = myComponent;
            myComponents.push_back(myComponentHelper);

            myCurrentLabel++;
        }
    }

    numberOfComponents = myCurrentLabel - 1;
    if (numberOfComponents <= 1)
        return 0;

    std::sort(myComponents.begin(), myComponents.end(),
              std::greater<CompStruct>());

    for (unsigned i = 0; i < myComponents.size() - 1; ++i) {
        unsigned other = i + 1;
        if (other >= myComponents.size())
            other = 0;
        {

            cut myCut;

            myCut.lambda = 0;
            if (params.sepIterFreeze == 0)
                myCut.frozen = false;

            double prizeMin = -99999999;
            double rhsPrize = -99999999;
            myCut.rhs = vector<nodevaluepair>();

            nodevaluepair n;
            n.value = 1.0;
            n.id = myComponents[i].components[0];

            for (unsigned j = 0; j < myComponents[i].components.size(); ++j) {
                if (instance.trueTerminals[myComponents[i].components[j]]) {
                    n.id = myComponents[i].components[j];
                    break;
                }

                if (!instance.realTerminals[myComponents[i].components[j]])
                    continue;
                rhsPrize = realPrizes[myComponents[i].components[j]] /
                           (max(instance.myBudgetCost[myComponents[i].components[j]], 0.00001));

                if (rhsPrize > prizeMin) {
                    n.id = myComponents[i].components[j];
                    prizeMin = rhsPrize;
                }
            }
            myCut.rhs.push_back(n);
            if (myComponents[i].sumPrize + epsInt >= incumbentObj ||
                params.separation == 0) {
                nodevaluepair m;
                m.value = 1.0;
                m.id = myComponents[other].components[0];

                prizeMin = -99999999;
                rhsPrize = -99999999;

                for (unsigned j = 0; j < myComponents[other].components.size();
                     ++j) {
                    if (instance.trueTerminals[myComponents[other].components[j]]) {
                        m.id = myComponents[other].components[j];
                        break;
                    }

                    if (!instance
                            .realTerminals[myComponents[other].components[j]])
                        continue;
                    rhsPrize = realPrizes[myComponents[other].components[j]] /
                               (max(instance.myBudgetCost[myComponents[other].components[j]], 0.00001));

                    if (rhsPrize > prizeMin) {
                        m.id = myComponents[other].components[j];
                        prizeMin = rhsPrize;
                    }
                }
                myCut.rhsConst = 1;
                myCut.rhs.push_back(m);
                if (myCut.rhs[0].id > myCut.rhs[1].id) {
                    int helper = myCut.rhs[0].id;
                    myCut.rhs[0].id = myCut.rhs[1].id;
                    myCut.rhs[1].id = helper;
                }
            }

            vector<int> myBoundary;

            if (params.separation == 0) {

                int otherBoundarySize =
                        myComponents[other].boundaryIndexed.size();
                int otherBoundaryFound = 0;
                queue<int> myQueue;

                vector<int> inComponent = vector<int>(instance.nNodes, 0);
                for (unsigned j = 0; j < myComponents[i].components.size();
                     ++j) {
                    inComponent[myComponents[i].components[j]] = true;
                }

                for (unsigned n = 0; n < instance.nNodes; n++) {
                    if (inComponent[n])
                        continue;

                    if (myComponents[i].boundary[n]) {
                        inComponent[n] = true;
                        if (myComponents[other].boundary[n]) {
                            myBoundary.push_back(n);
                            otherBoundaryFound++;
                            if (otherBoundaryFound == otherBoundarySize) {
                                break;
                            }
                        } else {
                            myQueue.push(n);
                        }
                    }
                }

                while (!myQueue.empty() &&
                       otherBoundaryFound < otherBoundarySize) {
                    int m = myQueue.front();
                    myQueue.pop();

                    for (int k : instance.adjList[m]) {
                        if (fixedToZero[k])
                            continue;

                        if (!inComponent[k]) {
                            inComponent[k] = true;

                            if (myComponents[other].boundary[k]) {

                                myBoundary.push_back(k);
                                otherBoundaryFound++;
                                if (otherBoundaryFound == otherBoundarySize)
                                    break;
                            } else {
                                myQueue.push(k);
                            }
                        }
                    }
                }
            } else {
                for (unsigned j = 0; j < myComponents[i].boundaryIndexed.size();
                     ++j) {
                    if (fixedToZero[myComponents[i].boundaryIndexed[j]])
                        continue;
                    myBoundary.push_back(myComponents[i].boundaryIndexed[j]);
                }
            }

            myCut.lhs = vector<nodevaluepair>();

            sort(myBoundary.begin(), myBoundary.end());

            for (int m : myBoundary) {
                nodevaluepair nv;
                nv.id = m;
                nv.value = 1.0;
                myCut.lhs.push_back(nv);
            }

            myCut.myHash = myCut.myHasher(myCut.lhs, myCut.rhs, instance.nNodes);

            if (myCutHash.find(myCut.myHash) != myCutHash.end())
                continue;

            if (myCut.rhs.size() == 2) {
                vector<nodevaluepair> dummy;
                nodevaluepair n;
                n.id = myCut.rhs[0].id;
                dummy.push_back(n);
                long helper1 =
                        myCut.myHasher(myCut.lhs, dummy, instance.nNodes);
                if (myCutHash.find(helper1) != myCutHash.end())
                    continue;

                vector<nodevaluepair> dummy2;
                nodevaluepair n2;
                n2.id = myCut.rhs[1].id;
                dummy2.push_back(n2);
                long helper2 =
                        myCut.myHasher(myCut.lhs, dummy2, instance.nNodes);
                if (myCutHash.find(helper2) != myCutHash.end())
                    continue;
            }

            myCutHash.insert(myCut.myHash);

            myCut.violated = true;
            myCut.direction = 0;
            myCut.directionPrevious = 0;
            myCut.age = 0;
            if (iterations % sepIterFreeze == 0) {
                myCut.frozen = false;
            }
            myNewCuts.push_back(myCut);

            numberOfCuts++;
        }
    }

    return numberOfCuts;
}
