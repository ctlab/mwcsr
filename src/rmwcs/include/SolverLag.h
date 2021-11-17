/*
 * SolverLag.h
 *
 *  Created on: Mar 29, 2015
 *      Author: markus
 */

#ifndef SOLVERLAG_H_
#define SOLVERLAG_H_

#include <list>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "Instance.h"
#include "../../include/monitor.h"

class SolverLag {
protected:
    Instance& instance;
    Parameters& params;
    mwcsr::monitor int_monitor;

public:
    struct nodevaluepair {
        int id = -1;
        double value = 1.0;

        bool operator<(const nodevaluepair& other) const {
            return value < other.value;
        }

        bool operator>(const nodevaluepair& other) const {
            return (value > other.value);
        }
    };

    struct CompStruct {
        // int id=-1;
        double sumPrize = 0.0;
        std::vector<bool> boundary;
        std::vector<int> boundaryIndexed;
        std::vector<int> components;
        // std::vector<vector<int>> componentsNested;
        std::vector<int> boundaryIndexedNested;

        bool operator<(const CompStruct& other) const {
            return sumPrize < other.sumPrize;
        }

        bool operator>(const CompStruct& other) const {
            return (sumPrize > other.sumPrize);
        }
    };

    struct cut {
        bool violated = true;
        double rhsConst = 0;
        std::vector<nodevaluepair> lhs;
        std::vector<nodevaluepair> rhs;

        int age = 0;
        double lambda = 0.0;
        double bestLambda = 0.0;
        double subgradient = 0.0;
        double direction = 0.0;
        double directionPrevious = 0.0;
        bool frozen = 1;
        bool toRemove = false;

        std::size_t myHash = -1;

        std::size_t myHasher(std::vector<nodevaluepair> const& vec,
                             std::vector<nodevaluepair> const& vec2,
                             int offset) const {
            std::size_t seed = 0;
            for (auto& i : vec) {
                seed ^= ((std::size_t) i.id) + 0x9e3779b9 + (seed << 6) +
                        (seed >> 2);
            }
            for (auto& i : vec2) {
                seed ^= ((std::size_t) i.id + offset) + 0x9e3779b9 +
                        (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

protected:
    std::vector<CompStruct> myComponents;
    std::list<double> cftBounds;
    std::list<double> cftBoundsBest;
    std::list<cut> myCuts;
    std::list<cut> myNewCuts;

    std::unordered_set<long> myCutHash;

    std::vector<double> realPrizes;
    std::vector<double> currentSolution;
    std::vector<double> previousSolution;
    std::vector<int> sumSolution;
    std::vector<bool> incumbent;
    std::vector<int> dualIncumbent;
    std::vector<int> labels;

    std::vector<int> fixedToZero;
    std::vector<int> fixedToOne;

    double incumbentObj;

    double subgradientSquared;
    double subgradientNorm;
    double directionPrevSquared;

    int solveSubgradient(int maxIterations);

    virtual double calculateCurrentSolution(bool save) = 0;

    virtual bool primalHeuristic() = 0;

    virtual int lagrangianPegging() = 0;

    virtual int separateCuts();

    virtual int addInitCuts() = 0;

    double calculateReducedCosts();

    void upgradeMultipliers();

    double calculateSubgradientCuts(const cut& myCut);

    double calculateSubgradientCutsPrevious(const cut& myCut);

    int checkPreviousCuts(bool addCuts);

    int createCuts(int iter);

    void writeStatistics();

    double alpha;
    int noImprov;

    int numberOfComponents;

    double bestBound;
    double currentBound;
    double previousBound;

    double bestBoundCFT;
    double worstBoundCFT;
    int counterCFT;

    int maxIterations;
    int iterations;

    int sepIter;
    int sepIterFreeze;

    double savedObj = 0.0;

    double runtime;

    void updateMultipliersCFT();

    void updateMultipliersLucena();

    void updateMultipliersSherali();

    // inline bool

    struct myclass {
        bool operator()(int i, int j) { return (i < j); }
    } nodeComparator;

    static constexpr double epsInt = 1e-3;
    static constexpr double epsOpt = 1e-6;

    int setVariableFixing(const std::vector<int>& toZero, const std::vector<int>& toOne);

    int initCuts();

    // bool initCuts;

public:
    SolverLag(Instance& instance, Parameters& params, mwcsr::monitor int_monitor);

    virtual ~SolverLag();

    std::string statf;

    int writeCutsToInstance();

    int writeFixingToInstance();

    int writeSolutionToInstance();

    inline double getBestBound() const { return bestBound; }

    inline double getIncumbentObj() const { return incumbentObj; }

    inline std::vector<bool> getIncumbent() const { return incumbent; }

    inline std::list<SolverLag::cut> getCuts() const { return myCuts; }

    inline double getDirectionPrevSquared() const {
        return directionPrevSquared;
    }

    inline void setDirectionPrevSquared(double dps) {
        directionPrevSquared = dps;
    }

    void initCuts(std::list<SolverLag::cut>& cuts);

    std::string getStatistics() { return statf; }

    int solve();
};

#endif /* SOLVERLAG_H_ */
