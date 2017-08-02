/*
 * SolverLag.h
 *
 *  Created on: Mar 29, 2015
 *      Author: markus
 */

#ifndef SOLVERLAG_H_
#define SOLVERLAG_H_

#include "instance/Instance.h"
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

class SolverLag {
  protected:
    Instance &instance;

  public:
    struct nodevaluepair {
        int id = -1;
        double value = 1.0;

        bool operator<(const nodevaluepair &other) const {
            // cout<<"sdsfdsf"<<endl;
            return value < other.value;
        }

        bool operator>(const nodevaluepair &other) const {
            return (value > other.value);
        }
    };

    struct CompStruct {
        // int id=-1;
        double sumPrize = 0.0;
        vector<bool> boundary;
        vector<int> boundaryIndexed;
        vector<int> components;
        // vector<vector<int>> componentsNested;
        vector<int> boundaryIndexedNested;

        bool operator<(const CompStruct &other) const {
            // cout<<"sdsfdsf"<<endl;
            return sumPrize < other.sumPrize;
        }

        bool operator>(const CompStruct &other) const {
            return (sumPrize > other.sumPrize);
        }
    };

    struct cut {
        bool violated = true;
        double rhsConst = 0;
        vector<nodevaluepair> lhs;
        vector<nodevaluepair> rhs;

        int age = 0;
        double lambda = 0.0;
        double bestLambda = 0.0;
        double subgradient = 0.0;
        double direction = 0.0;
        double directionPrevious = 0.0;
        bool frozen = 1;
        bool toRemove = false;

        std::size_t myHash = -1;

        std::size_t myHasher(std::vector<nodevaluepair> const &vec,
                             std::vector<nodevaluepair> const &vec2,
                             int offset) const {
            std::size_t seed = 0;
            for (auto &i : vec) {
                seed ^= ((std::size_t)i.id) + 0x9e3779b9 + (seed << 6) +
                        (seed >> 2);
            }
            for (auto &i : vec2) {
                seed ^= ((std::size_t)i.id + offset) + 0x9e3779b9 +
                        (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

  protected:
    vector<CompStruct> myComponents;
    list<double> cftBounds;
    list<double> cftBoundsBest;
    list<cut> myCuts;
    list<cut> myNewCuts;

    unordered_set<long> myCutHash;

    vector<double> realPrizes;
    vector<double> currentSolution;
    vector<double> previousSolution;
    vector<int> sumSolution;
    vector<bool> incumbent;
    vector<int> dualIncumbent;
    vector<int> labels;

    vector<int> fixedToZero;
    vector<int> fixedToOne;

    double incumbentObj;

    double subgradientSquared;
    double subgradientNorm;
    double directionPrevSquared;

    int solveSubgradient(int maxIterations);
    int RINS(int maxIterations);

    virtual double calculateCurrentSolution(bool save) = 0;
    virtual bool primalHeuristic() = 0;
    virtual int lagrangianPegging() = 0;
    virtual int separateCuts();
    virtual int addInitCuts() = 0;

    double calculateReducedCosts();
    void upgradeMultipliers();
    double calculateSubgradientCuts(const cut &myCut);
    double calculateSubgradientCutsPrevious(const cut &myCut);
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

    bool inRins;
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

    int setVariableFixing(const vector<int> &toZero, const vector<int> &toOne);
    int initCuts();

    // bool initCuts;

  public:
    SolverLag(Instance &_instance, int _maxIterations);
    virtual ~SolverLag();
    string statf;

    int writeCutsToInstance();
    int writeFixingToInstance();
    int writeSolutionToInstance();

    inline double getBestBound() const { return bestBound; }

    inline double getIncumbentObj() const { return incumbentObj; }

    inline vector<bool> getIncumbent() const { return incumbent; }

    inline list<SolverLag::cut> getCuts() const { return myCuts; }

    inline double getDirectionPrevSquared() const {
        return directionPrevSquared;
    }

    inline void setDirectionPrevSquared(double dps) {
        directionPrevSquared = dps;
    }

    void initCuts(list<SolverLag::cut> &cuts);

    string getStatistics() { return statf; }

    int solve();
};

#endif /* SOLVERLAG_H_ */
