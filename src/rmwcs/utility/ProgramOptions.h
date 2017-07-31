/*
 * ProgramOptions.h
 *
 * class for handling the program options
 *
 *  Created on: Mar 28, 2015
 *      Author: markus
 */

#ifndef PROGRAMOPTIONS_H_
#define PROGRAMOPTIONS_H_

#include <iostream>

class ProgramOptions
{

public:
	struct Parameters
	{
		/// Input file (instance)
		std::string file;
		std::string statfile;

		int problem=0;
		int budget=0;
		int solver=0;
		double budgetMultiplier=-1;
		int maxIter=1000;
		int memlimit=2800;
		int inputformat=0;
		int betaIter=5;
		int integer=false;
		int cardCons=0;
		int separation=0;
		int maxAge=10;
		int outputlag=1;
		int startcons=1;
		int pegging=1;
		int sepIterFreeze=50;
		int sepIter=10;
		int heurIter=10;
		int sepIterFreezeRins=10;
		int sepIterRins=10;
		int subgradient=0;
		double cardMultiplier=0;
		double beta=2.0;
		int timelimit=1800;
		int setting=0;
	};

	ProgramOptions(int &argc, char ** &argv);
	virtual ~ProgramOptions();

};

extern ProgramOptions::Parameters params;

#endif /* PROGRAMOPTIONS_H_ */
