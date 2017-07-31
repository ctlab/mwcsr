/*
 * mwcs.cpp
 *
 *  Created on: Mar 28, 2015
 *      Author: markus
 */

#include "utility/ProgramOptions.h"
#include "instance/Instance.h"
#include "solverLag/SolverLag.h"
#include "solverLag/SolverClassic.h"
#include "solverLag/SolverCardinality.h"
#include "solverLag/SolverBudget.h"
#include <boost/filesystem.hpp>
#include <fstream>
#include <iomanip>

int main(int argc, char* argv[])
{
	ProgramOptions po(argc, argv);

	string stat;
	Instance myInstance=Instance();


		if(params.problem==0)
		{
			SolverClassic mySolver=SolverClassic(myInstance,params.maxIter);
			mySolver.solve();
			stat=mySolver.getStatistics();
		}
		else if(params.problem==1)
		{
			SolverBudget mySolver=SolverBudget(myInstance,params.maxIter);
			mySolver.solve();
			stat=mySolver.getStatistics();
		}
		else if(params.problem==2)
		{
			SolverCardinality mySolver=SolverCardinality(myInstance,params.maxIter);
			mySolver.solve();
			stat=mySolver.getStatistics();
		}


	std::string filename(boost::filesystem::path(params.file).stem().string());


	cout<<"STAH,fullname,filename,nodes,edges,iterationsLag,runtimeLag,bestBoundLag,incumbentObjLag,gapLag,nFixedZero,nFixedOne"
			",problem,budget,budgetMultiplier,maxIter,inputformat,betaIter,integer"
			",cardCons,cardMultiplier,separation,maxAge,outputlag,startcons,pegging,sepIterFreeze,sepIter"
			",heurIter,subgradient,timelimit,setting"<<endl;

	cout<<"STAT,"
	<<params.file<<","
	<<filename<<","
	<<myInstance.nTrueNodes<<","
	<<myInstance.nTrueEdges<<","
	<<myInstance.iterationsLag<<","
	<<myInstance.runtimeLag<<","
	<<std::setprecision(9)<<myInstance.transformInternalValue(myInstance.bestBoundLag)<<","
	<<std::setprecision(9)<<myInstance.transformInternalValue(myInstance.incumbentObjLag)<<","
	<<myInstance.gapLag<<","
	<<myInstance.nFixedZero<<","
	<<myInstance.nFixedOne<<","
	<<params.problem<<","
	<<params.budget<<","
	<<params.budgetMultiplier<<","
	<<params.maxIter<<","
	<<params.inputformat<<","
	<<params.betaIter<<","
	<<params.integer<<","
	<<params.cardCons<<","
	<<params.cardMultiplier<<","
	<<params.separation<<","
	<<params.maxAge<<","
	<<params.outputlag<<","
	<<params.startcons<<","
	<<params.pegging<<","
	<<params.sepIterFreeze<<","
	<<params.sepIter<<","
	<<params.heurIter<<","
	<<params.subgradient<<","
	<<params.timelimit<<","
	<<params.setting<<endl;
	
	cout<<"SOL,";
	cout<<myInstance.transformInternalValue(myInstance.incumbentObjLag)<<",";
	cout<<myInstance.solSize<<endl;
	for(int i=0;i<myInstance.nTrueNodes;++i)
	{
		if(myInstance.incumbent[i])
		{
			cout<<"N,"<<i<<endl;
		}
	}
	
	stringstream solfile;
	solfile<<filename<<"_"<<params.problem<<"_"<<params.budget<<"_"<<params.budgetMultiplier<<"_"<<params.cardCons<<"_"<<params.cardMultiplier<<".sol";
	ofstream os(solfile.str());
	os<<"SOL,";
	os<<myInstance.transformInternalValue(myInstance.incumbentObjLag)<<",";
	os<<myInstance.solSize<<endl;
	for(int i=0;i<myInstance.nTrueNodes;++i)
	{
		if(myInstance.incumbent[i])
		{
			os<<"N,"<<i<<endl;
		}
	}
	os.close();	


	if (!params.statfile.empty())
	{
		ofstream os(params.statfile);
		os<<stat<<endl;
		os.close();
	}


	return 0;
}


