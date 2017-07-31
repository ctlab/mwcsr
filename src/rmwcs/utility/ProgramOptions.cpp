/*
 * ProgramOptions.cpp
 *
 *  Created on: Mar 28, 2015
 *      Author: markus
 */

#include "ProgramOptions.h"
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

ProgramOptions::Parameters params;

ProgramOptions::ProgramOptions(int &argc, char ** &argv)
{
	po::options_description general_options("General options");
	general_options.add_options()
			("help,h", "produce help message")
			("file,f", po::value<string>(&params.file), "instance file to process")
			//("statfile", po::value<string>(&params.statfile), "output stat file")
			("inputformat,i", po::value<int>(&params.inputformat), "instance format: 0-PCSTP, 1-MWCS")
			("problem", po::value<int>(&params.problem), "problem: 0-class, 1-cardinality, 2-budget")
			("budget", po::value<int>(&params.budget), "budget as integer")
			("budgetMultiplier", po::value<double>(&params.budgetMultiplier), "budgetMultiplier: set budget as fraction of complete weight; parameter budget will be ignored")
			("cardMultiplier", po::value<double>(&params.cardMultiplier), "cardMultiplier: set cardinality as fraction of number of nodes; parameter cardCons will be ignored")
			("cardCons", po::value<int>(&params.cardCons), "cardinality constraint as integer")
			("integer", po::value<int>(&params.integer), "input is integer")
			("outputlag", po::value<int>(&params.outputlag), "output of lagrangian")
			("startcons", po::value<int>(&params.startcons), "add flow-conservation/degree cons at start?")
			("subgradient", po::value<int>(&params.subgradient), "subgradient: 0-classic, 1-average direction, 2-cft")
			("separation", po::value<int>(&params.separation), "separation: 0-strong, 1-fast")
			("maxAge", po::value<int>(&params.maxAge), "extending the life of non-violated ineqs")
			("beta", po::value<double>(&params.beta), "beta for subgradient")
			("pegging", po::value<int>(&params.pegging), "pegging (variable fixing)")
			("sepIterFreeze", po::value<int>(&params.sepIterFreeze), "after how many iterations we are checking added ineqs (i.e., >0 gives delayed relax-and-cut) ")
			("heurIter", po::value<int>(&params.heurIter), "after how many iterations we are doing heuristics")
			("sepIter", po::value<int>(&params.sepIter), "after how many iterations we are separating (i.e., >0 gives delayed relax-and-cut) ")
			("betaIterations", po::value<int>(&params.betaIter), "number of nonimproving iterations until beta is halfed")
			("maxIterations,m", po::value<int>(&params.maxIter), "maximum number of subgradient iterations")
			("timelimit", po::value<int>(&params.timelimit), "timelimit")
			("setting", po::value<int>(&params.setting), "to quickly set parameters to settings from the paper, see ProgramOptions.cpp for the details of the available settings; e.g., setting 23 should work stable for most problems (it is setting FAC from the paper with 1000 iterations)");


	po::options_description all("Allowed parameters");
	all.add(general_options);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
	    cout << general_options << endl;
	    exit(0);
	}
	else if (!vm.count("file"))
	{
		cout << "No input file given!" << endl;
		exit(0);
	}


	if(params.setting==1)
	{
		params.solver=0;
		params.maxIter=400;
		params.betaIter=20;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=0;
		params.separation=1;
		params.maxAge=3;
		params.beta=2;
		params.pegging=1;
	}
	if(params.setting==2)
	{
		params.solver=0;
		params.maxIter=400;
		params.betaIter=20;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=0;
		params.separation=1;
		params.maxAge=1000;
		params.beta=2;
		params.pegging=1;
	}
	if(params.setting==3)
	{
		params.solver=0;
		params.maxIter=400;
		params.betaIter=20;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=0;
		params.separation=0;
		params.maxAge=3;
		params.beta=2;
		params.pegging=1;
	}

	if(params.setting==11)
	{
		params.solver=0;
		params.maxIter=400;
		params.betaIter=20;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=1;
		params.separation=1;
		params.maxAge=3;
		params.beta=2;
		params.pegging=1;
	}

	if(params.setting==12)
	{
		params.solver=0;
		params.maxIter=400;
		params.betaIter=20;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=1;
		params.separation=1;
		params.maxAge=1000;
		params.beta=2;
		params.pegging=1;
	}
	if(params.setting==13)
	{
		params.solver=0;
		params.maxIter=400;
		params.betaIter=20;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=1;
		params.separation=0;
		params.maxAge=3;
		params.beta=2;
		params.pegging=1;
	}


	if(params.setting==21)
	{
		params.solver=0;
		params.maxIter=1000;
		params.betaIter=50;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=0;
		params.separation=1;
		params.maxAge=3;
		params.beta=2;
		params.pegging=1;
	}
	if(params.setting==22)
	{
		params.solver=0;
		params.maxIter=1000;
		params.betaIter=50;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=0;
		params.separation=1;
		params.maxAge=1000;
		params.beta=2;
		params.pegging=1;
	}
	if(params.setting==23)
	{
		params.solver=0;
		params.maxIter=1000;
		params.betaIter=50;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=0;
		params.separation=0;
		params.maxAge=3;
		params.beta=2;
		params.pegging=1;
	}

	if(params.setting==31)
	{
		params.solver=0;
		params.maxIter=1000;
		params.betaIter=50;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=1;
		params.separation=1;
		params.maxAge=3;
		params.beta=2;
		params.pegging=1;
	}

	if(params.setting==32)
	{
		params.solver=0;
		params.maxIter=1000;
		params.betaIter=50;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=1;
		params.separation=1;
		params.maxAge=1000;
		params.beta=2;
		params.pegging=1;
	}
	if(params.setting==33)
	{
		params.solver=0;
		params.maxIter=1000;
		params.betaIter=50;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=1;
		params.separation=0;
		params.maxAge=3;
		params.beta=2;
		params.pegging=1;
	}

	if(params.setting==81)
	{
		params.solver=1;
		params.maxIter=400;
		params.betaIter=20;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=0;
		params.separation=1;
		params.maxAge=3;
		params.beta=2;
		params.pegging=1;
	}

	if(params.setting==82)
	{
		params.solver=1;
		params.maxIter=400;
		params.betaIter=20;
		params.sepIter=1;
		params.sepIterFreeze=1;
		params.heurIter=1;
		params.subgradient=0;
		params.separation=1;
		params.maxAge=3;
		params.beta=2;
		params.pegging=1;
	}


}


ProgramOptions::~ProgramOptions()
{

}

