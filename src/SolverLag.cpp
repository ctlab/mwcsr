/*
 * SolverLag.cpp
 *
 *  Created on: Mar 29, 2015
 *      Author: markus
 */

#include "include/Parameters.h"
#include <queue>
#include "include/SolverLag.h"
#include <stack>
#include <chrono>

using std::vector;
using std::list;
using std::queue;
using std::max;

using Rcpp::Rcout;

namespace chrono = std::chrono;

bool cutToRemove(SolverLag::cut x) { return x.toRemove; }

SolverLag::SolverLag(Instance& instance, Parameters& params)
        : instance(instance), params(params), myCuts{list<cut>()}, myNewCuts{list<cut>()},
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
          maxIterations{params.maxIter}, iterations{0}, sepIter(params.sepIter),
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
                // Rcout<<i.id<<" ";
            }
            // Rcout<<"\n";
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

    chrono::time_point <std::chrono::system_clock> startTime =
            chrono::system_clock::now();
    iterations = 0;
    // if(params.outputlag)
    //	Rcout<<"fixed "<<costBasedFixing()<<" variables to zero due to cost in
    //component"<<"\n";

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
            //Rcout<<i<<"\n";
    }*/

    while (iterations < maxIterations && sqrt(subgradientSquared) > epsOpt) {
        boundImprov = false;
        // Rcout<<sqrt(subgradientSquared)<<"\n";
        subgradientSquared = 0.0;
        currentBound = calculateCurrentSolution(true);
        // Rcout<<currentBound<<"\n";
        // Rcout<<"incumbentObj"<<incumbentObj<<"\n";
        // Rcout<<currentBound<<" "<<previousBound<<"
        // "<<(currentBound>previousBound)<<"\n";

        // double boundPCSTP=0.0;

        if (currentBound < bestBound) {
            // Rcout<<currentBound<<" "<<bestBound<<"
            // "<<(currentBound>bestBound)<<"\n";
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
                    // Rcout<<c<<" "<<instance.maxRevenueInComponent[c]<<"\n";
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

        // Rcout<<"bestBoundCheck "<<bestBoundCheck<<" "<<incumbentObj<<"\n";
        // Rcout<<bestBound<<" "<<currentBound<<" "<<incumbentObj<<"
        // "<<bestBoundCheck<<"\n";

        // Rcout<<bestBoundCheck+eps<<" "<<incumbentObj<<"\n";

        if (bestBoundCheck <= incumbentObj + eps) {
            break;
        }

        if (params.outputlag) {
            if (inRins) {
                Rcout << "RINS ";
            }
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
        for (int i = 0; i < instance.nNodes; ++i) {
            previousSolution[i] = currentSolution[i];
        }

        iterations++;
        // Rcout<<alpha<<" "<< sqrt(subgradientSquared)<<"\n";
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

    // time_t endT;
    // time(&endT);
    // runtime=difftime(endT,startT);
    runtime = elapsedSeconds.count();
    writeStatistics();

    // Rcout<<statf.c_str()<<"\n";
    return 1;
}

double SolverLag::calculateReducedCosts() {
    double obj = 0.0;
    for (int n = 0; n < instance.nNodes; ++n) {
        realPrizes[n] = instance.myPrizes[n];
        // if(realPrizes[n]>0 ||instance.myPrizes[n]>0)
        // Rcout<<realPrizes[n]<<" "<<instance.myPrizes[n]<<"\n";
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
        // Rcout<<c.rhsConst<<" "<<c.lambda<<" "<<obj<<"\n";
    }

    /*
    for(int n=0;n<instance.nNodes;++n)
    {
            //if(realPrizes[n]!=instance.myPrizes[n])
            //Rcout<<realPrizes[n]<<" "<<instance.myPrizes[n]<<"\n";
    }*/

    /*
    for(int n=0;n<instance.nNodes;++n)
    {
            //if(realPrizes[n]!=instance.myPrizes[n])
            Rcout<<realPrizes[n]<<" y("<<n<<")+ ";
    }
    Rcout<<"\n";*/

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
/*
    instance.gapLag = (instance.transformInternalValue(bestBound) -
                       instance.transformInternalValue(incumbentObj)) /
                      instance.transformInternalValue(incumbentObj) * 100;
    if (params.inputformat == 0) {
        instance.gapLag = (instance.transformInternalValue(incumbentObj) -
                           instance.transformInternalValue(bestBound)) /
                          instance.transformInternalValue(bestBound) * 100;
    }*/

    if (instance.gapLag < epsOpt)
        instance.gapLag = 0;
}

int SolverLag::createCuts(int iter) {
    int numberOfCuts = 0;
    int cnt = 0;

    if (iter % sepIter == 0 && iter > 0)
        cnt += separateCuts();
    // Rcout<<cnt<<"\n";
    // numberOfCuts+=cnt;
    if (iter % sepIterFreeze == 0) {
        // Rcout<<"here"<<"\n";
        for (cut& c : myCuts) {
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
        currentBound = bestBound;
        for (int n = 0; n < instance.nNodes; ++n) {
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
    // Rcout<<subgradientSquared<<"\n";
    // Rcout<<directionPrevSquared<<"\n";

    if (directionPrevSquared > epsOpt) {
        sigma = sqrt(subgradientSquared) / sqrt(directionPrevSquared);
    }

    // Rcout<<"sigma "<<sigma<<"\n";
    double directionSquared = 0.0;
    // Rcout<<"sherali "<<myCuts.size()<<"\n";
    for (cut& c : myCuts) {
        if (c.frozen)
            continue;
        // cerr<<c.direction<<"\n";
        c.direction = c.subgradient + sigma * c.directionPrevious;
        c.directionPrevious = c.direction;
        directionSquared += (c.direction * c.direction);
    }

    if (directionSquared < epsOpt) {
        // Rcout<<"directionSquared "<<directionSquared<<"\n";
        directionSquared = subgradientSquared;
        for (cut& c : myCuts) {
            if (c.frozen)
                continue;
            c.direction = c.subgradient;
        }
    }

    double theta = alpha * (currentBound - incumbentObj) / (directionSquared);
    // Rcout<<subgradientSquared<<" "<<directionSquared<<"\n";

    // Rcout<<theta<<" "<<alpha<<" "<<incumbentObj<<" "<<directionSquared<<"
    // "<<subgradientSquared<<"\n";

    for (cut& c : myCuts) {
        if (c.frozen)
            continue;
        // cerr<<c.lambda<<" "<<c.lambda-theta*c.direction<<"
        // "<<c.direction<<"\n";
        c.lambda = std::max(0.0, c.lambda - theta * c.direction);
        // c.directionPrevious=c.direction;
        // if(c.lambda<epsOpt)
        //	c.lambda=0;
        // Rcout<<c.lambda<<"\n";
        // if(c.lambda>0)
        //	cerr<<c.lambda<<" "<<c.subgradient<<"\n";
    }

    // directionPrevSquared=directionSquared;
}

void SolverLag::updateMultipliersLucena() {
    if (noImprov > params.betaIter) {
        noImprov = 0;
        alpha /= 2;
    }

    double theta = alpha * (currentBound - incumbentObj) / (subgradientSquared);
    // Rcout<<subgradientSquared<<"\n";
    // Rcout<<theta<<" "<<incumbentObj<<" "<<alpha<<" "<<currentBound<<"
    // "<<subgradientSquared<<"\n";
    // incumbentObj=7;

    for (cut& c : myCuts) {
        if (c.frozen)
            continue;
        // cerr<<"cut "<<c.lambda<<" "<<c.lambda-theta*c.subgradient<<"
        // "<<c.subgradient<<"\n";
        c.lambda = std::max(0.0, c.lambda - theta * c.subgradient);
        // cerr<<"cut "<<c.lambda<<"\n";
    }
}

void SolverLag::updateMultipliersCFT() {
    // Rcout<<noImprov<<"\n";
    if (noImprov >= params.betaIter) {
        noImprov = 0;
        alpha /= 2;
        for (cut& c : myCuts) {
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
                // Rcout<<oldBest<<" "<<bestBound<<"\n";
                double boundSum = 0.0;

                for (double bound : cftBounds) {
                    boundSum += bound;
                }

                double avgBound = boundSum / (double)cftBounds.size();
                // Rcout<<avgBound<<" "<<bestBound<<" "<<delta<<"
                // "<<0.001*delta<<"\n";
                // Rcout<<"alpha"<<alpha<<"\n";
                if (avgBound - bestBound < 0.001 * delta) {
                    // Rcout<<"delta10"<<"\n";
                    alpha *= 10;
                } else if (avgBound - bestBound < 0.01 * delta) {
                    // Rcout<<"delta2"<<"\n";
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
                // Rcout<<"alpha"<<alpha<<"\n";
            }
        }
    }
#endif

    double theta = alpha * (currentBound - incumbentObj) / (subgradientSquared);
    // Rcout<<theta<<"\n";

    for (cut& c : myCuts) {
        if (c.frozen)
            continue;
        // cerr<<c.lambda<<" "<<c.lambda-theta*c.subgradient<<"
        // "<<c.subgradient<<"\n";
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
    // Rcout<<"myCuts.size() "<<myCuts.size()<<"\n";
    for (cut& c : myCuts) {
        if (c.frozen)
            continue;

        c.violated = true;
        c.subgradient = calculateSubgradientCuts(c);
        // Rcout<<c.subgradient<<"\n";

        double lhs = 0.0;
        unsigned nFixedToZero = 0;
        for (nodevaluepair n : c.lhs) {
            lhs += n.value * currentSolution[n.id];
            if (fixedToZero[n.id])
                nFixedToZero++;
        }

        if (nFixedToZero == c.lhs.size() && c.rhs.size() == 1 && addCuts) {
            // Rcout<<"here"<<"\n";
            for (nodevaluepair n : c.rhs) {
                // Rcout<<n.id<<" "<<fixedToZero.size()<<"\n";
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
                //	Rcout<<n.id<<" "<<incumbent[n.id]<<" "<<c.lhs.size()<<"
                //"<<instance.adjList[n.id].size()<<"\n";
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
            // Rcout<<c.rhsConst<<"\n";
            c.subgradient = 0;
            c.toRemove = true;
        }

        // Rcout<<lhs<<" "<<rhs<<"\n";

        if (lhs < rhs) {
            numberOfCuts++;
            c.age = 0;

            /*
            for(nodevaluepair l:c.lhs)
                    Rcout<<l.id<<" ";
            Rcout<<"\n";
            for(nodevaluepair r:c.rhs)
                    Rcout<<r.id<<" ";
            Rcout<<"\n";*/

            if (addCuts && false) {
                for (cut c2 : myNewCuts) {
                    if (c2.rhs[0].id == c.rhs[0].id && c.myHash == c2.myHash) {
                        if (c.rhs[0].value >= c2.rhs[0].value) {
                            c2.toRemove = true;
                            // c.subgradient=0;
                            // numberOfCuts--;
                            // Rcout<<"remove"<<"\n";
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

            // cerr<<c.subgradient<<"\n";

            /*
            for(nodevaluepair n: c.lhs)
            {
                    Rcout<<n.id<<" ";
            }
            Rcout<<"\n";

            for(nodevaluepair n: c.rhs)
            {
                    Rcout<<n.id<<" ";
            }
            Rcout<<"\n";
            Rcout<<c.lambda<<"\n";*/

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

        // Rcout<<"prev "<<c.subgradient<<" "<<c.violated<<" "<<lhs<<"
        // "<<subgradientSquared<<" "<<c.age<<"\n";
    }

    // numberOfCuts+=myNewCuts.size();
    if (addCuts) {
        // myNewCuts.erase(std::remove_if(myNewCuts.begin(), myNewCuts.end(),
        // cutToRemove),	myNewCuts.end());

        for (cut& c : myNewCuts) {
            if (!c.frozen) {
                c.subgradient = calculateSubgradientCuts(c);
                c.directionPrevious = calculateSubgradientCutsPrevious(c);
                // directionPrevSquared+=(c.subgradient*c.subgradient);
                subgradientSquared += (c.subgradient * c.subgradient);
            }
            // Rcout<<"push"<<" "<<myNewCuts.size()<<"\n";
            // Rcout<<c.rhs.size()<<" "<<c.lhs.size()<<" "<<c.rhsConst<<"\n";
            myCuts.push_back(c);
        }
    }

    // Rcout<<"subgradientSquared "<<subgradientSquared<<"\n";
    // Rcout<<myCuts.size()<<"\n";

    // Rcout<<subgradientSquared<<" "<<numberOfCuts<<"\n";

    // cerr<<myCuts.size()<<"\n";
    // myCuts.remove_if([&](const cut& c){ return (c.lambda==0.0 &&
    // c.subgradient==0 ); } );
    // cerr<<myCuts.size()<<"\n";
    // cerr<<notViolated<<"\n";

    return numberOfCuts;
}

double SolverLag::calculateSubgradientCuts(const cut& myCut) {
    double subg = myCut.rhsConst;
    for (nodevaluepair n : myCut.lhs) {
        // Rcout<<n<<" "<<currentSolution[n]<<" "<<subg<<"\n";
        subg += (n.value * currentSolution[n.id]);
    }

    for (nodevaluepair n : myCut.rhs) {
        // Rcout<<n<<" "<<currentSolution[n]<<" "<<subg<<"\n";
        subg -= (n.value * currentSolution[n.id]);
    }
    // Rcout<<myCut.rhs1<<" "<<myCut.rhs2<<"\n";

    // Rcout<<"subg "<<subg<<"\n";

    return subg;
}

double SolverLag::calculateSubgradientCutsPrevious(const cut& myCut) {
    double subg = myCut.rhsConst;
    for (nodevaluepair n : myCut.lhs) {
        // Rcout<<n<<" "<<currentSolution[n]<<" "<<subg<<"\n";
        subg += (n.value * previousSolution[n.id]);
    }

    for (nodevaluepair n : myCut.rhs) {
        // Rcout<<n<<" "<<currentSolution[n]<<" "<<subg<<"\n";
        subg -= (n.value * previousSolution[n.id]);
    }
    // Rcout<<myCut.rhs1<<" "<<myCut.rhs2<<"\n";

    // Rcout<<"subg "<<subg<<"\n";

    return subg;
}

int SolverLag::setVariableFixing(const vector<int>& toZero,
                                 const vector<int>& toOne) {
    int numberFixed = toZero.size() + toOne.size();

    for (unsigned i = 0; i < toZero.size(); ++i) {
        // Rcout<<toZero[i]<<" ";
        fixedToZero[toZero[i]] = true;
    }
    // Rcout<<"\n";

    for (unsigned i = 0; i < toOne.size(); ++i) {
        // Rcout<<toOne[i]<<" ";
        fixedToOne[toOne[i]] = true;
    }
    // Rcout<<"\n";

    return numberFixed;
}

void SolverLag::initCuts(list<SolverLag::cut>& cuts) {
    myCuts = cuts;
    for (cut& c : myCuts) {
        // cerr<<c.direction<<"\n";
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
        // cerr<<n<<" "<<currentSolution[n]<<"\n";
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
            // Rcout<<"\n";
            myComponentHelper.sumPrize += instance.myPrizes[n];
            // Rcout<<instance.myPrizes[n]<<"\n";
            labels[n] = myCurrentLabel;

            queue<int> myQueue;
            myQueue.push(n);
            // cerr<<"component "<<n<<" "<<myCurrentLabel<<"\n";
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
                            // Rcout<<li<<" "<<myComponentHelper.sumPrize<<"
                            // "<<instance.myPrizes[li]<<"\n";
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

            // Rcout<<"myComponent: "<<myComponent.size()<<" myCurrentLabel:
            // "<<myCurrentLabel<<"\n";

            myCurrentLabel++;
        }
    }

    numberOfComponents = myCurrentLabel - 1;
    if (numberOfComponents <= 1)
        return 0;

    std::sort(myComponents.begin(), myComponents.end(),
              std::greater<CompStruct>());

    for (unsigned i = 0; i < myComponents.size() - 1; ++i) {
        // Rcout<<myComponents[i].sumPrize<<"\n";
        // if(myComponents[i].sumPrize<incumbentObj/10)
        //	break;
        unsigned other = i + 1;
        // if(other>=i)
        //	continue;
        if (other >= myComponents.size())
            other = 0;
        // Rcout<<myComponents[i].sumPrize<<"
        // "<<myComponents[other].sumPrize<<"\n";

        // Rcout<<i<<" "<<other<<"\n";
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

            // Rcout<<"\n";
            for (unsigned j = 0; j < myComponents[i].components.size(); ++j) {
                // if(myComponents[i].sumPrize<=incumbentObj)
                //	Rcout<<myComponents[i].components[j]<<"\n";
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
            // Rcout<<myComponents[i].sumPrize<<" "<<incumbentObj<<"\n";
            if (myComponents[i].sumPrize + epsInt >= incumbentObj ||
                params.separation == 0) {
                nodevaluepair m;
                m.value = 1.0;
                m.id = myComponents[other].components[0];

                // Rcout<<n<<" "<<instance.trueTerminals[n]<<"\n";
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
                    // Rcout<<"j "<<j<<" "<<components[other][j]<<"\n";
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

            // Rcout<<n<<" "<<instance.trueTerminals[n]<<"\n";
            // Rcout<<"\n";

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
                            //	Rcout<<otherBoundaryFound<<"
                            //"<<otherBoundarySize<<"\n";
                            if (otherBoundaryFound == otherBoundarySize) {
                                break;
                            }
                        } else {
                            myQueue.push(n);
                        }
                    }
                }

                // cerr<<"component "<<n<<" "<<myCurrentLabel<<"\n";

                while (!myQueue.empty() &&
                       otherBoundaryFound < otherBoundarySize) {
                    int m = myQueue.front();
                    // int m=myQueue.top();

                    // cerr<<m<<"\n";
                    myQueue.pop();

                    // forall_adj_edges(e,m)
                    for (int k : instance.adjList[m]) {
                        if (fixedToZero[k])
                            continue;

                        // if(k>=instance.nNodes)

                        if (!inComponent[k]) {
                            inComponent[k] = true;
                            // Rcout<<k<<" "<<instance.nNodes<<"\n";

                            if (myComponents[other].boundary[k]) {

                                myBoundary.push_back(k);
                                otherBoundaryFound++;
                                // Rcout<<otherBoundaryFound<<"
                                // "<<otherBoundarySize<<"\n";
                                if (otherBoundaryFound == otherBoundarySize)
                                    break;
                            } else {
                                myQueue.push(k);
                            }
                        }
                    }
                }

                // Rcout<<otherBoundaryFound<<" "<<otherBoundarySize<<"\n";
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
                // Rcout<<m<<" ";
                nv.id = m;
                nv.value = 1.0;
                myCut.lhs.push_back(nv);
            }

            // myCut.lambda=myCut.lhs.size();
            myCut.myHash =
                    myCut.myHasher(myCut.lhs, myCut.rhs, instance.nNodes);

            // Rcout<<"cut"<<"\n";
            if (myCutHash.find(myCut.myHash) != myCutHash.end())
                continue;
            // Rcout<<"cut1"<<"\n";

            if (myCut.rhs.size() == 2) {
                vector<nodevaluepair> dummy;
                nodevaluepair n;
                n.id = myCut.rhs[0].id;
                dummy.push_back(n);
                long helper1 =
                        myCut.myHasher(myCut.lhs, dummy, instance.nNodes);
                if (myCutHash.find(helper1) != myCutHash.end())
                    continue;
                // Rcout<<"cut2"<<"\n";

                vector<nodevaluepair> dummy2;
                nodevaluepair n2;
                n2.id = myCut.rhs[1].id;
                dummy2.push_back(n2);
                long helper2 =
                        myCut.myHasher(myCut.lhs, dummy2, instance.nNodes);
                if (myCutHash.find(helper2) != myCutHash.end())
                    continue;
                // Rcout<<"cut3"<<"\n";
            }

            // Rcout<<"cut4"<<"\n";

            myCutHash.insert(myCut.myHash);
            /*Rcout<<"insert"<<"\n";
            for(nodevaluepair l:myCut.lhs)
                    Rcout<<l.id<<" ";
            Rcout<<"\n";
            for(nodevaluepair r:myCut.rhs)
                    Rcout<<r.id<<" ";
            Rcout<<"\n";
            Rcout<<"insertend"<<"\n";*/

            // Rcout<<myCut.myHash<<"\n";
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
            // Rcout<<myCut.rhs.size()<<" "<<myCut.lhs.size()<<"
            // "<<myCut.rhsConst<<"\n";
            myNewCuts.push_back(myCut);

            numberOfCuts++;
        }
    }

    // Rcout<<myNewCuts.size()<<"\n";

    return numberOfCuts;
}
