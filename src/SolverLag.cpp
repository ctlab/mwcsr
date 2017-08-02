/*
 * SolverLag.cpp
 *
 *  Created on: Mar 29, 2015
 *      Author: markus
 */

#include "utility/ProgramOptions.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <solverLag/SolverLag.h>
#include <stack>
//#include <time.h>
#include <boost/filesystem.hpp>
#include <chrono>
#include <iomanip> // std::setprecision

bool cutToRemove(SolverLag::cut x) { return x.toRemove; }

SolverLag::SolverLag(Instance &_instance, int _maxIterations)
    : instance(_instance), myCuts{list<cut>()}, myNewCuts{list<cut>()},
      realPrizes{vector<double>(instance.nNodes, 0)},
      currentSolution{vector<double>(instance.nNodes, 0)},
      previousSolution{vector<double>(instance.nNodes, 0)},
      sumSolution{vector<int>(instance.nNodes, 0)}, incumbent{vector<bool>(
                                                        instance.nNodes, 0)},
      dualIncumbent{vector<int>(instance.nNodes, 0)}, labels{vector<int>(
                                                          instance.nNodes, 0)},
      fixedToZero{vector<int>(instance.nNodes, 0)}, fixedToOne{vector<int>(
                                                        instance.nNodes, 0)},
      incumbentObj{0.0}, subgradientSquared{1}, subgradientNorm{0},
      directionPrevSquared{0}, alpha{params.beta}, noImprov{0},
      numberOfComponents{0}, bestBound{std::numeric_limits<double>::max()},
      currentBound{std::numeric_limits<double>::max()},
      previousBound{std::numeric_limits<double>::max()},
      bestBoundCFT{std::numeric_limits<double>::max()},
      worstBoundCFT{std::numeric_limits<double>::lowest()}, counterCFT{0},
      maxIterations{_maxIterations}, iterations{0}, sepIter(params.sepIter),
      sepIterFreeze(params.sepIterFreeze), inRins(false), savedObj(0.0),
      runtime(0.0) {}

SolverLag::~SolverLag() {}

int SolverLag::solve() {

    solveSubgradient(maxIterations);
    if (params.solver == 1) {
        // writeCutsToInstance();
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
                // cout<<i.id<<" ";
            }
            // cout<<endl;
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

    for (int i = 0; i < instance.nNodes; ++i) {
        instance.fixedToOne[i] = fixedToOne[i];
        instance.fixedToZero[i] = fixedToZero[i];
    }

    return 1;
}

int SolverLag::writeSolutionToInstance() {
    instance.incumbent = vector<bool>(instance.nNodes, false);

    for (int i = 0; i < instance.nNodes; ++i) {
        instance.incumbent[i] = incumbent[i];
        instance.incumbentFound = true;
    }

    return 0;
}

int SolverLag::solveSubgradient(int maxIterations) {
    // time_t startT;
    // time(&startT);

    chrono::time_point<std::chrono::system_clock> startTime =
        chrono::system_clock::now();
    iterations = 0;
    // if(params.outputlag)
    //	cout<<"fixed "<<costBasedFixing()<<" variables to zero due to cost in
    //component"<<endl;

    double eps = 1e-10;

    int nFixed = 0;

    if (params.integer) {
        eps = 0;
    }

    bool boundImprov = false;
    addInitCuts();

    /*
    for(int i:instance.myTrueTerminals)
    {
            fixedToOne[i]=true;
            //cout<<i<<endl;
    }*/

    while (iterations < maxIterations && sqrt(subgradientSquared) > epsOpt) {
        boundImprov = false;
        // cout<<sqrt(subgradientSquared)<<endl;
        subgradientSquared = 0.0;
        currentBound = calculateCurrentSolution(true);
        // cout<<currentBound<<endl;
        // cout<<"incumbentObj"<<incumbentObj<<endl;
        // cout<<currentBound<<" "<<previousBound<<"
        // "<<(currentBound>previousBound)<<endl;

        // double boundPCSTP=0.0;

        if (currentBound < bestBound) {
            // cout<<currentBound<<" "<<bestBound<<"
            // "<<(currentBound>bestBound)<<endl;
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

                // if(params.subgradient==1)
                {
                    for (int i = 0; i < instance.nNodes; ++i) {
                        dualIncumbent[i] = currentSolution[i];
                    }
                }
            } else
                noImprov++;
        }

        // previousBound=currentBound;

        if (params.subgradient == 2) {
            if (boundImprov) {
                bestBound = currentBound;
                for (cut &c : myCuts) {
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

        // if(iterations%50==0 && iterations>0) RINS(100);

        if (iterations % params.heurIter == 0) {
            primalHeuristic();
        }

        nFixed = 0;
        if ((boundImprov) && params.pegging > 0 && !inRins) {
            // calculateCurrentSolution(false);
            for (unsigned c = 0; c < instance.components.size(); ++c) {
                if (instance.maxRevenueInComponent[c] < incumbentObj &&
                    !instance.componentFixed[c]) {
                    // cout<<c<<" "<<instance.maxRevenueInComponent[c]<<endl;
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
            if (nFixed) {
                cout << "fixed " << nFixed << " variables" << endl;
            }
        }

        // cout<<"bestBoundCheck "<<bestBoundCheck<<" "<<incumbentObj<<endl;
        // cout<<bestBound<<" "<<currentBound<<" "<<incumbentObj<<"
        // "<<bestBoundCheck<<endl;

        // cout<<bestBoundCheck+eps<<" "<<incumbentObj<<endl;

        if (bestBoundCheck <= incumbentObj + eps) {
            cout << bestBound << endl;
            cout << "Bound check optimality" << endl;
            break;
        }

        if (params.outputlag) {
            if (inRins) {
                cout << "RINS ";
            }
            cout << std::setprecision(9) << "iteration: \t" << iterations
                 << "\t lagrangian bound: \t"
                 << instance.transformInternalValue(bestBound)
                 << "\t current bound: \t "
                 << instance.transformInternalValue(currentBound)
                 << "\t incumbent: \t "
                 << instance.transformInternalValue(incumbentObj)
                 << "\t number of active cuts: \t" << numberOfCuts << endl;
        }
        if (nFixed)
            myCuts.erase(
                std::remove_if(myCuts.begin(), myCuts.end(), cutToRemove),
                myCuts.end());
        upgradeMultipliers();
        for (int i = 0; i < instance.nNodes; ++i) {
            previousSolution[i] = currentSolution[i];
        }

        iterations++;
        // cout<<alpha<<" "<< sqrt(subgradientSquared)<<endl;
    }

    if (params.outputlag) {
        cout << "finished" << endl;
        cout << std::setprecision(9) << "iteration: \t" << iterations
             << "\t lagrangian bound: \t"
             << instance.transformInternalValue(bestBound)
             << "\t incumbent: \t "
             << instance.transformInternalValue(incumbentObj) << endl;
    }

    chrono::time_point<std::chrono::system_clock> endTime =
        chrono::system_clock::now();
    chrono::duration<double> elapsedSeconds = endTime - startTime;

    // time_t endT;
    // time(&endT);
    // runtime=difftime(endT,startT);
    runtime = elapsedSeconds.count();
    writeStatistics();

    // cout<<statf.c_str()<<endl;
    return 1;
}

double SolverLag::calculateReducedCosts() {
    double obj = 0.0;
    for (int n = 0; n < instance.nNodes; ++n) {
        realPrizes[n] = instance.myPrizes[n];
        // if(realPrizes[n]>0 ||instance.myPrizes[n]>0)
        // cout<<realPrizes[n]<<" "<<instance.myPrizes[n]<<endl;
    }

    for (cut &c : myCuts) {
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
        // cout<<c.rhsConst<<" "<<c.lambda<<" "<<obj<<endl;
    }

    /*
    for(int n=0;n<instance.nNodes;++n)
    {
            //if(realPrizes[n]!=instance.myPrizes[n])
            //cout<<realPrizes[n]<<" "<<instance.myPrizes[n]<<endl;
    }*/

    /*
    for(int n=0;n<instance.nNodes;++n)
    {
            //if(realPrizes[n]!=instance.myPrizes[n])
            cout<<realPrizes[n]<<" y("<<n<<")+ ";
    }
    cout<<endl;*/

    return obj;
}

void SolverLag::writeStatistics() {

    // std::string
    // filename(boost::filesystem::path(params.file).stem().string());

    instance.bestBoundLag = bestBound;
    instance.incumbentObjLag = incumbentObj;
    instance.iterationsLag = iterations;
    instance.runtimeLag = runtime;

    instance.incumbent = vector<bool>(instance.nTrueNodes, false);
    instance.solSize = 0;
    for (int i = 0; i < instance.nNodes; ++i) {
        if (incumbent[i]) {
            instance.solSize++;
            instance.incumbent[instance.map[i]] = true;
        }
    }

    instance.gapLag = (instance.transformInternalValue(bestBound) -
                       instance.transformInternalValue(incumbentObj)) /
                      instance.transformInternalValue(incumbentObj) * 100;
    if (params.inputformat == 0) {
        instance.gapLag = (instance.transformInternalValue(incumbentObj) -
                           instance.transformInternalValue(bestBound)) /
                          instance.transformInternalValue(bestBound) * 100;
    }

    if (instance.gapLag < epsOpt)
        instance.gapLag = 0;
}

int SolverLag::createCuts(int iter) {
    int numberOfCuts = 0;
    int cnt = 0;

    if (iter % sepIter == 0 && iter > 0)
        cnt += separateCuts();
    // cout<<cnt<<endl;
    // numberOfCuts+=cnt;
    if (iter % sepIterFreeze == 0) {
        // cout<<"here"<<endl;
        for (cut &c : myCuts) {
            c.frozen = false;
        }
        numberOfCuts += cnt;
        // alpha=alpha*2;
    }

    numberOfCuts += checkPreviousCuts(true);

    return numberOfCuts;
}

void SolverLag::updateMultipliersSherali() {
    if (noImprov > params.betaIter) {
        noImprov = 0;
        alpha /= 2;
        if (params.outputlag)
            cout << "no improvement for" << params.betaIter
                 << " iterations, new alpha is " << alpha << endl;
        currentBound = bestBound;
        for (int n = 0; n < instance.nNodes; ++n) {
            currentSolution[n] = dualIncumbent[n];
        }
        subgradientSquared = 0;
        checkPreviousCuts(false);
        for (cut &c : myCuts) {
            if (c.frozen)
                continue;
            c.directionPrevious = 0;
        }
    }
    directionPrevSquared = 0.0;
    for (cut &c : myCuts) {
        if (c.frozen)
            continue;
        directionPrevSquared += (c.directionPrevious * c.directionPrevious);
    }

    double sigma = 0;
    // cout<<subgradientSquared<<endl;
    // cout<<directionPrevSquared<<endl;

    if (directionPrevSquared > epsOpt) {
        sigma = sqrt(subgradientSquared) / sqrt(directionPrevSquared);
    }

    // cout<<"sigma "<<sigma<<endl;
    double directionSquared = 0.0;
    // cout<<"sherali "<<myCuts.size()<<endl;
    for (cut &c : myCuts) {
        if (c.frozen)
            continue;
        // cerr<<c.direction<<endl;
        c.direction = c.subgradient + sigma * c.directionPrevious;
        c.directionPrevious = c.direction;
        directionSquared += (c.direction * c.direction);
    }

    if (directionSquared < epsOpt) {
        // cout<<"directionSquared "<<directionSquared<<endl;
        directionSquared = subgradientSquared;
        for (cut &c : myCuts) {
            if (c.frozen)
                continue;
            c.direction = c.subgradient;
        }
    }

    double theta = alpha * (currentBound - incumbentObj) / (directionSquared);
    // cout<<subgradientSquared<<" "<<directionSquared<<endl;

    // cout<<theta<<" "<<alpha<<" "<<incumbentObj<<" "<<directionSquared<<"
    // "<<subgradientSquared<<endl;

    for (cut &c : myCuts) {
        if (c.frozen)
            continue;
        // cerr<<c.lambda<<" "<<c.lambda-theta*c.direction<<"
        // "<<c.direction<<endl;
        c.lambda = std::max(0.0, c.lambda - theta * c.direction);
        // c.directionPrevious=c.direction;
        // if(c.lambda<epsOpt)
        //	c.lambda=0;
        // cout<<c.lambda<<endl;
        // if(c.lambda>0)
        //	cerr<<c.lambda<<" "<<c.subgradient<<endl;
    }

    // directionPrevSquared=directionSquared;
}

void SolverLag::updateMultipliersLucena() {
    if (noImprov > params.betaIter) {
        noImprov = 0;
        alpha /= 2;
        if (params.outputlag)
            cout << "no improvement for" << params.betaIter
                 << " iterations, new alpha is " << alpha << endl;
    }

    double theta = alpha * (currentBound - incumbentObj) / (subgradientSquared);
    // cout<<subgradientSquared<<endl;
    // cout<<theta<<" "<<incumbentObj<<" "<<alpha<<" "<<currentBound<<"
    // "<<subgradientSquared<<endl;
    // incumbentObj=7;

    for (cut &c : myCuts) {
        if (c.frozen)
            continue;
        // cerr<<"cut "<<c.lambda<<" "<<c.lambda-theta*c.subgradient<<"
        // "<<c.subgradient<<endl;
        c.lambda = std::max(0.0, c.lambda - theta * c.subgradient);
        // cerr<<"cut "<<c.lambda<<endl;
    }
}

void SolverLag::updateMultipliersCFT() {
    // cout<<noImprov<<endl;
    if (noImprov >= params.betaIter) {
        noImprov = 0;
        alpha /= 2;
        for (cut &c : myCuts) {
            c.lambda = c.bestLambda;
        }
    }
#ifdef COMMENT
    else {
        cftBounds.push_back(currentBound);
        cftBoundsBest.push_back(bestBound);

        if (cftBounds.size() > 100) {
            cftBounds.pop_front();
            double oldBest = cftBoundsBest.front();
            cftBoundsBest.pop_front();
            double delta = bestBound - incumbentObj;
            if (oldBest - bestBound < 0.01 * delta) {
                // cout<<oldBest<<" "<<bestBound<<endl;
                double boundSum = 0.0;

                for (double bound : cftBounds) {
                    boundSum += bound;
                }

                double avgBound = boundSum / (double)cftBounds.size();
                // cout<<avgBound<<" "<<bestBound<<" "<<delta<<"
                // "<<0.001*delta<<endl;
                // cout<<"alpha"<<alpha<<endl;
                if (avgBound - bestBound < 0.001 * delta) {
                    // cout<<"delta10"<<endl;
                    alpha *= 10;
                } else if (avgBound - bestBound < 0.01 * delta) {
                    // cout<<"delta2"<<endl;
                    alpha *= 2;
                } else /*if(alpha/2>epsOpt)*/
                {
                    alpha /= 2;
                }
                /*
                else
                {

                        for(cut& c: myCuts)
                        {
                                c.lambda=c.bestLambda;
                        }
                        alpha=params.beta;
                }*/
                // cout<<"alpha"<<alpha<<endl;
            }
        }
    }
#endif

    double theta = alpha * (currentBound - incumbentObj) / (subgradientSquared);
    // cout<<theta<<endl;

    for (cut &c : myCuts) {
        if (c.frozen)
            continue;
        // cerr<<c.lambda<<" "<<c.lambda-theta*c.subgradient<<"
        // "<<c.subgradient<<endl;
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

    // int notViolated=0;
    // cout<<"myCuts.size() "<<myCuts.size()<<endl;
    for (cut &c : myCuts) {
        if (c.frozen)
            continue;

        c.violated = true;
        c.subgradient = calculateSubgradientCuts(c);
        // cout<<c.subgradient<<endl;

        double lhs = 0.0;
        unsigned nFixedToZero = 0;
        for (nodevaluepair n : c.lhs) {
            lhs += n.value * currentSolution[n.id];
            if (fixedToZero[n.id])
                nFixedToZero++;
        }

        if (nFixedToZero == c.lhs.size() && c.rhs.size() == 1 && addCuts) {
            // cout<<"here"<<endl;
            for (nodevaluepair n : c.rhs) {
                // cout<<n.id<<" "<<fixedToZero.size()<<endl;
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

                // if(incumbent[n.id])
                //	cout<<n.id<<" "<<incumbent[n.id]<<" "<<c.lhs.size()<<"
                //"<<instance.adjList[n.id].size()<<endl;
            }
            c.toRemove = true;
            c.subgradient = 0;

            // continue;
        }

        double rhs = -c.rhsConst;
        nFixedToZero = 0;

        for (nodevaluepair n : c.rhs) {
            rhs += n.value * currentSolution[n.id];
            if (fixedToZero[n.id])
                nFixedToZero++;
        }

        if (nFixedToZero == c.rhs.size() + c.rhsConst && addCuts) {
            // cout<<c.rhsConst<<endl;
            c.subgradient = 0;
            c.toRemove = true;
        }

        // cout<<lhs<<" "<<rhs<<endl;

        if (lhs < rhs) {
            numberOfCuts++;
            c.age = 0;

            /*
            for(nodevaluepair l:c.lhs)
                    cout<<l.id<<" ";
            cout<<endl;
            for(nodevaluepair r:c.rhs)
                    cout<<r.id<<" ";
            cout<<endl;*/

            if (addCuts && false) {
                for (cut c2 : myNewCuts) {
                    if (c2.rhs[0].id == c.rhs[0].id && c.myHash == c2.myHash) {
                        if (c.rhs[0].value >= c2.rhs[0].value) {
                            c2.toRemove = true;
                            // c.subgradient=0;
                            // numberOfCuts--;
                            // cout<<"remove"<<endl;
                            break;
                        }
                        /*
                        else
                        {
                                c2.toRemove=true;
                                break;
                        }*/
                    }
                }
            }
            // numberOfCuts++;

            // cerr<<c.subgradient<<endl;

            /*
            for(nodevaluepair n: c.lhs)
            {
                    cout<<n.id<<" ";
            }
            cout<<endl;

            for(nodevaluepair n: c.rhs)
            {
                    cout<<n.id<<" ";
            }
            cout<<endl;
            cout<<c.lambda<<endl;*/

        } else {
            c.violated = false;
            c.age++;

            if (c.lambda == 0 && c.subgradient > 0 && c.age > params.maxAge) {
                c.subgradient = 0;
            }

            if (c.lambda == 0 && c.subgradient == 0 &&
                c.directionPrevious > 0 && c.age > params.maxAge) {
                // directionPrevSquared-=(c.directionPrevious*c.directionPrevious);
                // c.directionPrevious=0.0;
            }
        }
        subgradientSquared += (c.subgradient * c.subgradient);

        // cout<<"prev "<<c.subgradient<<" "<<c.violated<<" "<<lhs<<"
        // "<<subgradientSquared<<" "<<c.age<<endl;
    }

    // numberOfCuts+=myNewCuts.size();
    if (addCuts) {
        // myNewCuts.erase(std::remove_if(myNewCuts.begin(), myNewCuts.end(),
        // cutToRemove),	myNewCuts.end());

        for (cut &c : myNewCuts) {
            if (!c.frozen) {
                c.subgradient = calculateSubgradientCuts(c);
                c.directionPrevious = calculateSubgradientCutsPrevious(c);
                // directionPrevSquared+=(c.subgradient*c.subgradient);
                subgradientSquared += (c.subgradient * c.subgradient);
            }
            // cout<<"push"<<" "<<myNewCuts.size()<<endl;
            // cout<<c.rhs.size()<<" "<<c.lhs.size()<<" "<<c.rhsConst<<endl;
            myCuts.push_back(c);
        }
    }

    // cout<<"subgradientSquared "<<subgradientSquared<<endl;
    // cout<<myCuts.size()<<endl;

    // cout<<subgradientSquared<<" "<<numberOfCuts<<endl;

    // cerr<<myCuts.size()<<endl;
    // myCuts.remove_if([&](const cut& c){ return (c.lambda==0.0 &&
    // c.subgradient==0 ); } );
    // cerr<<myCuts.size()<<endl;
    // cerr<<notViolated<<endl;

    return numberOfCuts;
}

double SolverLag::calculateSubgradientCuts(const cut &myCut) {
    double subg = myCut.rhsConst;
    for (nodevaluepair n : myCut.lhs) {
        // cout<<n<<" "<<currentSolution[n]<<" "<<subg<<endl;
        subg += (n.value * currentSolution[n.id]);
    }

    for (nodevaluepair n : myCut.rhs) {
        // cout<<n<<" "<<currentSolution[n]<<" "<<subg<<endl;
        subg -= (n.value * currentSolution[n.id]);
    }
    // cout<<myCut.rhs1<<" "<<myCut.rhs2<<endl;

    // cout<<"subg "<<subg<<endl;

    return subg;
}

double SolverLag::calculateSubgradientCutsPrevious(const cut &myCut) {
    double subg = myCut.rhsConst;
    for (nodevaluepair n : myCut.lhs) {
        // cout<<n<<" "<<currentSolution[n]<<" "<<subg<<endl;
        subg += (n.value * previousSolution[n.id]);
    }

    for (nodevaluepair n : myCut.rhs) {
        // cout<<n<<" "<<currentSolution[n]<<" "<<subg<<endl;
        subg -= (n.value * previousSolution[n.id]);
    }
    // cout<<myCut.rhs1<<" "<<myCut.rhs2<<endl;

    // cout<<"subg "<<subg<<endl;

    return subg;
}

int SolverLag::setVariableFixing(const vector<int> &toZero,
                                 const vector<int> &toOne) {
    int numberFixed = toZero.size() + toOne.size();

    for (unsigned i = 0; i < toZero.size(); ++i) {
        // cout<<toZero[i]<<" ";
        fixedToZero[toZero[i]] = true;
    }
    // cout<<endl;

    for (unsigned i = 0; i < toOne.size(); ++i) {
        // cout<<toOne[i]<<" ";
        fixedToOne[toOne[i]] = true;
    }
    // cout<<endl;

    return numberFixed;
}

void SolverLag::initCuts(list<SolverLag::cut> &cuts) {
    myCuts = cuts;
    for (cut &c : myCuts) {
        // cerr<<c.direction<<endl;
        c.subgradient = 0;
        c.direction = 0;
        c.directionPrevious = 0;
        c.age = 0;
        // c.lambda=0;
    }
}

int SolverLag::separateCuts() {
    int numberOfCuts = 0;
    myNewCuts.clear();
    // int n=-1;

    int myCurrentLabel = 1;
    std::fill(labels.begin(), labels.end(), 0);
    myComponents.clear();

    for (int n = 0; n < instance.nNodes; ++n) {
        // cerr<<n<<" "<<currentSolution[n]<<endl;
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
            // cout<<endl;
            myComponentHelper.sumPrize += instance.myPrizes[n];
            // cout<<instance.myPrizes[n]<<endl;
            labels[n] = myCurrentLabel;

            queue<int> myQueue;
            myQueue.push(n);
            // cerr<<"component "<<n<<" "<<myCurrentLabel<<endl;
            while (!myQueue.empty()) {
                // node k=myQueue.front();
                int k = myQueue.front();

                myQueue.pop();
                // forall_adj_edges(e,k)
                for (int li : instance.adjList[k]) {
                    if (currentSolution[li] == 0) {
                        // componentBoundary[l->index()]=true;
                        // if(!fixedToZero[li])
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
                            // cout<<li<<" "<<myComponentHelper.sumPrize<<"
                            // "<<instance.myPrizes[li]<<endl;
                            myComponentHelper.sumPrize += instance.myPrizes[li];
                        }
                    }
                }
            }
            myComponentHelper.boundaryIndexed = componentBoundaryIndexed;
            myComponentHelper.boundary = componentBoundary;
            myComponentHelper.components = myComponent;
            myComponents.push_back(myComponentHelper);
            // componentsNested.push_back(myComponent);

            // cout<<"myComponent: "<<myComponent.size()<<" myCurrentLabel:
            // "<<myCurrentLabel<<endl;

            myCurrentLabel++;
        }
    }

    numberOfComponents = myCurrentLabel - 1;
    cout << "numberOfComponents " << numberOfComponents << endl;
    if (numberOfComponents <= 1)
        return 0;

    std::sort(myComponents.begin(), myComponents.end(),
              std::greater<CompStruct>());

    for (unsigned i = 0; i < myComponents.size() - 1; ++i) {
        // cout<<myComponents[i].sumPrize<<endl;
        // if(myComponents[i].sumPrize<incumbentObj/10)
        //	break;
        unsigned other = i + 1;
        // if(other>=i)
        //	continue;
        if (other >= myComponents.size())
            other = 0;
        // cout<<myComponents[i].sumPrize<<"
        // "<<myComponents[other].sumPrize<<endl;

        // cout<<i<<" "<<other<<endl;
        // for(unsigned other=i+1;other<myComponents.size();++other)
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

            // cout<<endl;
            for (unsigned j = 0; j < myComponents[i].components.size(); ++j) {
                // if(myComponents[i].sumPrize<=incumbentObj)
                //	cout<<myComponents[i].components[j]<<endl;
                if (instance.trueTerminals[myComponents[i].components[j]]) {
                    n.id = myComponents[i].components[j];
                    break;
                }

                // rhsPrize=instance.myPrizes[myComponents[i].components[j]];
                if (!instance.realTerminals[myComponents[i].components[j]])
                    continue;
                rhsPrize =
                    realPrizes[myComponents[i].components[j]] /
                    (max(instance.myBudgetCost[myComponents[i].components[j]],
                         0.00001));

                if (rhsPrize > prizeMin) {
                    n.id = myComponents[i].components[j];
                    prizeMin = rhsPrize;
                }
            }
            myCut.rhs.push_back(n);
            // cout<<myComponents[i].sumPrize<<" "<<incumbentObj<<endl;
            if (myComponents[i].sumPrize + epsInt >= incumbentObj ||
                params.separation == 0) {
                nodevaluepair m;
                m.value = 1.0;
                m.id = myComponents[other].components[0];

                // cout<<n<<" "<<instance.trueTerminals[n]<<endl;
                prizeMin = -99999999;
                rhsPrize = -99999999;

                for (unsigned j = 0; j < myComponents[other].components.size();
                     ++j) {
                    if (instance
                            .trueTerminals[myComponents[other].components[j]]) {
                        m.id = myComponents[other].components[j];
                        break;
                    }

                    if (!instance
                             .realTerminals[myComponents[other].components[j]])
                        continue;
                    // cout<<"j "<<j<<" "<<components[other][j]<<endl;
                    // rhsPrize=instance.myPrizes[myComponents[other].components[j]];
                    rhsPrize = realPrizes[myComponents[other].components[j]] /
                               (max(instance.myBudgetCost[myComponents[other]
                                                              .components[j]],
                                    0.00001));

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

            // cout<<n<<" "<<instance.trueTerminals[n]<<endl;
            // cout<<endl;

            vector<int> myBoundary;

            if (params.separation == 0) {

                int otherBoundarySize =
                    myComponents[other].boundaryIndexed.size();
                int otherBoundaryFound = 0;
                // queue<node> myQueue=queue<node>();
                queue<int> myQueue;

                vector<int> inComponent = vector<int>(instance.nNodes, 0);
                for (unsigned j = 0; j < myComponents[i].components.size();
                     ++j) {
                    // myQueue.push(components[i][j]);
                    inComponent[myComponents[i].components[j]] = true;
                }

                for (int n = 0; n < instance.nNodes; n++) {
                    // if(fixedToZero[n])
                    //	continue;

                    if (inComponent[n])
                        continue;

                    if (myComponents[i].boundary[n]) {
                        inComponent[n] = true;
                        if (myComponents[other].boundary[n]) {
                            myBoundary.push_back(n);
                            otherBoundaryFound++;
                            //	cout<<otherBoundaryFound<<"
                            //"<<otherBoundarySize<<endl;
                            if (otherBoundaryFound == otherBoundarySize) {
                                break;
                            }
                        } else {
                            myQueue.push(n);
                        }
                    }
                }

                // cerr<<"component "<<n<<" "<<myCurrentLabel<<endl;

                while (!myQueue.empty() &&
                       otherBoundaryFound < otherBoundarySize) {
                    int m = myQueue.front();
                    // int m=myQueue.top();

                    // cerr<<m<<endl;
                    myQueue.pop();

                    // forall_adj_edges(e,m)
                    for (int k : instance.adjList[m]) {
                        if (fixedToZero[k])
                            continue;

                        // if(k>=instance.nNodes)

                        if (!inComponent[k]) {
                            inComponent[k] = true;
                            // cout<<k<<" "<<instance.nNodes<<endl;

                            if (myComponents[other].boundary[k]) {

                                myBoundary.push_back(k);
                                otherBoundaryFound++;
                                // cout<<otherBoundaryFound<<"
                                // "<<otherBoundarySize<<endl;
                                if (otherBoundaryFound == otherBoundarySize)
                                    break;
                            } else {
                                myQueue.push(k);
                            }
                        }
                    }
                }

                // cout<<otherBoundaryFound<<" "<<otherBoundarySize<<endl;
            } else {
                for (unsigned j = 0; j < myComponents[i].boundaryIndexed.size();
                     ++j) {
                    if (fixedToZero[myComponents[i].boundaryIndexed[j]])
                        continue;
                    myBoundary.push_back(myComponents[i].boundaryIndexed[j]);
                }
                // myBoundary=boundariesNested[i];
            }

            myCut.lhs = vector<nodevaluepair>();

            sort(myBoundary.begin(), myBoundary.end());

            for (int m : myBoundary) {
                nodevaluepair nv;
                // cout<<m<<" ";
                nv.id = m;
                nv.value = 1.0;
                myCut.lhs.push_back(nv);
            }

            // myCut.lambda=myCut.lhs.size();
            myCut.myHash =
                myCut.myHasher(myCut.lhs, myCut.rhs, instance.nNodes);

            // cout<<"cut"<<endl;
            if (myCutHash.find(myCut.myHash) != myCutHash.end())
                continue;
            // cout<<"cut1"<<endl;

            if (myCut.rhs.size() == 2) {
                vector<nodevaluepair> dummy;
                nodevaluepair n;
                n.id = myCut.rhs[0].id;
                dummy.push_back(n);
                long helper1 =
                    myCut.myHasher(myCut.lhs, dummy, instance.nNodes);
                if (myCutHash.find(helper1) != myCutHash.end())
                    continue;
                // cout<<"cut2"<<endl;

                vector<nodevaluepair> dummy2;
                nodevaluepair n2;
                n2.id = myCut.rhs[1].id;
                dummy2.push_back(n2);
                long helper2 =
                    myCut.myHasher(myCut.lhs, dummy2, instance.nNodes);
                if (myCutHash.find(helper2) != myCutHash.end())
                    continue;
                // cout<<"cut3"<<endl;
            }

            // cout<<"cut4"<<endl;

            myCutHash.insert(myCut.myHash);
            /*cout<<"insert"<<endl;
            for(nodevaluepair l:myCut.lhs)
                    cout<<l.id<<" ";
            cout<<endl;
            for(nodevaluepair r:myCut.rhs)
                    cout<<r.id<<" ";
            cout<<endl;
            cout<<"insertend"<<endl;*/

            // cout<<myCut.myHash<<endl;
            myCut.violated = true;
            myCut.direction = 0;
            myCut.directionPrevious = 0;
            myCut.age = 0;
            // myCut.lambda
            if (iterations % sepIterFreeze == 0) {
                myCut.frozen = false;
                //	myCut.subgradient=calculateSubgradientCuts(myCut);
                //	subgradientSquared+=(myCut.subgradient*myCut.subgradient);
            }
            // myCuts.push_back(myCut);
            // cout<<myCut.rhs.size()<<" "<<myCut.lhs.size()<<"
            // "<<myCut.rhsConst<<endl;
            myNewCuts.push_back(myCut);

            numberOfCuts++;
        }
    }

    // cout<<myNewCuts.size()<<endl;

    return numberOfCuts;
}
