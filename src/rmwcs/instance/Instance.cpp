/*
 * Instance.cpp
 *
 *  Created on: Mar 29, 2015
 *      Author: markus
 */

#include "Instance.h"
#include "utility/ProgramOptions.h"

#include <sstream>
#include <string>
#include <queue>
#include <limits>
#include <fstream>
#include <algorithm>
#include <cmath>

Instance::Instance():
components {vector<vector<int>>()},
maxRevenueInComponent{ vector<double>() },
nComponents{0},
maxPrize {0},
minWeight{std::numeric_limits<double>::max()},
sumPrizes {0},
budget(params.budget),
nRealTerminals {0}

{


	if(params.inputformat==0)
		readPCSTP();
	else if(params.inputformat==1)
		readMWCS();
	else if(params.inputformat==2)
		readDilkina();
	else if(params.inputformat==3)
		readSTP();
	else if(params.inputformat==4)
		readMWCSBudget();
	cout<<"Instance reading finished"<<endl;
	cout<<"Nodes: \t \t \t"<<nNodes<<endl;
	cout<<"Edges: \t \t \t"<<nEdges<<endl;
	cout<<"TerminalsRead: \t \t"<<nTerminals<<endl;

	/*
	myTerminals.clear();


	for(int i=0;i<nNodes;++i)
	{
		if(realTerminals[i])
		{
			//cout<<i<<" "<<myPrizes[i]<<endl;
			bool remainTerminal=false;
			for(int k:adjList[i])
			{
				if(myPrizes[k]+myPrizes[i]>0)
				{
					remainTerminal=true;
					break;
				}
			}
			if(!remainTerminal)
			{
				realTerminals[i]=false;
			}
			else
			{
				myTerminals.push_back(i);
			}
		}
	}
	nTerminals=myTerminals.size();
	cout<<"TerminalsReal: \t \t"<<nTerminals<<endl;*/



	nTrueNodes=nNodes;
	nTrueEdges=nEdges;

	nComponents=calculateComponents();
	cout<<"Components:  \t \t "<<nComponents<<endl;
	if(params.problem==1)
		cout<<"Budget:  \t \t "<<budget<<endl;
	if(params.problem==2)
		cout<<"Cardinality:  \t \t "<<params.cardCons<<endl;
	//exit(0);
	if(preprocessing()>0)
	{
		rebuildDatastructures();
		cout<<"Size after preprocessing"<<endl;
		cout<<"Nodes: \t \t \t"<<nNodes<<endl;
		cout<<"Edges: \t \t \t"<<nEdges<<endl;
		cout<<"Terminals:  \t \t"<<nTerminals<<endl;
		nComponents=calculateComponents();
		cout<<"Components:  \t \t "<<nComponents<<endl;
	}



	fixedToOne=vector<int>(nNodes,0);
	fixedToZero=vector<int>(nNodes,0);

	nFlowNodes=2*nNodes;
	nFlowArcs=2*nEdges+nNodes;

	flowArcs=new std::pair<int,int>[nFlowArcs];
	int counter=0;
	for(int i=0;i<nNodes;++i)
	{
		flowArcs[counter]=std::make_pair(i,i+nNodes);
		counter++;
	}

	for(int i=0;i<nNodes;++i)
	{
		for(int j:adjList[i])
		{
			flowArcs[counter]=std::make_pair(i+nNodes,j);
			counter++;
		}
	}

}


void Instance::readSTP()
{
	std::ifstream infile(params.file);
	std::string line;

	std::string magicNumber;
	std::getline(infile, line);
	std::istringstream iss(line);
	iss>>magicNumber;

	if(magicNumber!="33D32945")
	{
		cout<<"Not in STP Format"<<endl;
		exit(0);
	}

	int srcID, dstID;
	double cost;
	int edgeCounter=0;
	int termID;
	int count=0;
	while(std::getline(infile, line))
	{
		if(line.empty())
			continue;

		iss.clear();
		iss.str(line);
		std::string token;

		iss>>token;
		//cerr<<count<<" "<<token<<" "<<line<<" "<<iss.str()<<endl;
		count++;
		if(token=="Nodes")
		{
			iss>>nNodes;
			adjList=vector<vector<int>>(nNodes);
			//cerr<<nNodes<<endl;
			//myNodes=vector<node>(nNodes,nullptr);
			realTerminals=vector<bool>(nNodes,false);
			trueTerminals=vector<bool>(nNodes,false);
			myPrizes=vector<double>(nNodes);
			myBudgetCost=vector<double>(nNodes);

			nodesToRemove=vector<bool>(nNodes);
			for(int i=0;i<nNodes;i++)
			{
				//myNodes[i]=G.newNode(i);
				myPrizes[i]=-1;
				myBudgetCost[i]=1;
				adjList[i]=vector<int>();
			}
		}
		else if(token=="Edges")
		{
			iss>>nEdges;
			//cerr<<nEdges<<endl;

			//myEdges=vector<edge>(nEdges,nullptr);
		}
		else if(token=="E")
		{
			iss>>srcID>>dstID>>cost;
			if(cost<minWeight)
				minWeight=cost;
			//cerr<<edgeCounter<<endl;
			//myEdges[edgeCounter]=G.newEdge(myNodes[srcID-1],myNodes[dstID-1]);
			//myCost.push_back(cost);

			adjList[srcID-1].push_back(dstID-1);
			adjList[dstID-1].push_back(srcID-1);

			edgeCounter++;
		}
		else if(token=="Terminals")
		{
			iss>>nTerminals;
			//myTerminals=vector<int>(nTerminals,0);
		}
		else if(token=="T")
		{
			iss>>termID;
			trueTerminals[termID-1]=true;
			myTrueTerminals.push_back(termID-1);
			realTerminals[termID-1]=true;
			myTerminals.push_back(termID-1);

		}
	}
	//cout<<myTrueTerminals.size()<<endl;
}

void Instance::readDilkina()
{
	std::ifstream infile(params.file);
	std::string line;

	std::getline(infile, line);
	std::getline(infile, line);
	std::istringstream iss(line);
	iss>>nNodes;
	myBudgetCost=vector<double>(nNodes);
	myPrizes=vector<double>(nNodes);
	realTerminals=vector<bool>(nNodes,false);
	trueTerminals=vector<bool>(nNodes,false);
	adjList=vector<vector<int>>(nNodes);
	nEdges=0;

	std::getline(infile, line);
	std::getline(infile, line);

	double prizes, cost, dummy;

	for(int i=0;i<nNodes;++i)
	{
		std::getline(infile, line);
		//cerr<<i<<" "<<line<<endl;
		std::istringstream iss(line);

		iss.clear();
		iss.str(line);
		iss>>dummy>>dummy>>dummy>>prizes>>cost;
		myPrizes[i]=prizes;
		myBudgetCost[i]=cost;
		adjList[i]=vector<int>();
	}

	std::getline(infile, line);
	std::getline(infile, line);

	int state=0;

	int n1, n2;
	//, termID;
	int totalCost=0;
	nTerminals=0;
	while(std::getline(infile, line))
	{
		iss.clear();
		iss.str(line);
		std::string token;

		iss>>token;
		//cerr<<token<<endl;
		if(token=="Reserves")
		{
			state=1;
			continue;
		}
		else if(token=="TotalCost")
		{
			state=2;
			continue;
		}
		else if(token=="SteinerNodes")
		{
			state=3;
			continue;
		}


		if(state==0)
		{
			iss>>n1>>n2;
			nEdges++;
			//cerr<<n1<<" "<<n2<<endl;
			adjList[n1].push_back(n2);
			adjList[n2].push_back(n1);
		}
		else if(state==1)
		{
			//cerr<<line<<endl;
			//iss>>termID;
			//termID=atoi(token.c_str());
			//nRealTerminals++;
			//cout<<termID<<" ";
			//nTerminals++;
			//realTerminals[termID]=true;
			//trueTerminals[termID]=true;
			//myTerminals.push_back(termID);
			//myTrueTerminals.push_back(termID);
		}
		else if(state==2)
		{
			//totalCost=atoi(token.c_str());;
			//cout<<iss.str().c_str()<<endl;
		}
	}

	totalCost=0;
	for(int i=0;i<nNodes;++i)
	{
		totalCost+=myBudgetCost[i];
		if(myPrizes[i]>0)
		{
			realTerminals[i]=true;
			myTerminals.push_back(i);
		}
	}
	nTerminals=myTerminals.size();

	if(params.budgetMultiplier>=1)
	{
		budget=floor(totalCost*params.budgetMultiplier/100.0);
	}
	//cout<<nRealTerminals<<endl;
}

void Instance::readPCSTP()
{

	std::ifstream infile(params.file);
	std::string line;

	std::string magicNumber;
	std::getline(infile, line);
	std::istringstream iss(line);
	iss>>magicNumber;

	if(magicNumber!="33D32945")
	{
		cout<<"Not in STP Format"<<endl;
		exit(0);
	}

	int srcID, dstID;
	double cost;
	int edgeCounter=0;
	int termID;
	double prize;
	int termCounter=0;
	int count=0;
	while(std::getline(infile, line))
	{
		if(line.empty())
			continue;

		iss.clear();
		iss.str(line);
		std::string token;

		iss>>token;
		//cerr<<count<<" "<<token<<" "<<line<<" "<<iss.str()<<endl;
		count++;
		if(token=="Nodes")
		{
			iss>>nNodes;
			adjList=vector<vector<int>>(nNodes);
			//cerr<<nNodes<<endl;
			//myNodes=vector<node>(nNodes,nullptr);
			realTerminals=vector<bool>(nNodes,false);
			trueTerminals=vector<bool>(nNodes,false);
			myPrizes=vector<double>(nNodes);
			myBudgetCost=vector<double>(nNodes,1);
			nodesToRemove=vector<bool>(nNodes);
			for(int i=0;i<nNodes;i++)
			{
				//myNodes[i]=G.newNode(i);
				adjList[i]=vector<int>();
			}
		}
		else if(token=="Edges")
		{
			iss>>nEdges;
			//cerr<<nEdges<<endl;

			//myEdges=vector<edge>(nEdges,nullptr);
		}
		else if(token=="E")
		{
			iss>>srcID>>dstID>>cost;
			if(cost<minWeight)
				minWeight=cost;
			//cerr<<edgeCounter<<endl;
			//myEdges[edgeCounter]=G.newEdge(myNodes[srcID-1],myNodes[dstID-1]);
			//myCost.push_back(cost);

			adjList[srcID-1].push_back(dstID-1);
			adjList[dstID-1].push_back(srcID-1);

			edgeCounter++;
		}
		else if(token=="Terminals")
		{
			iss>>nTerminals;
			myTerminals=vector<int>(nTerminals,0);
		}
		else if(token=="TP")
		{
			iss>>termID>>prize;
			myTerminals[termCounter]=termID-1;
			myPrizes[myTerminals[termCounter]]=prize;
			sumPrizes+=prize;
			termCounter++;
		}
	}
	nTerminals=0;
	for(int i=0;i<nNodes;i++)
	{
		myPrizes[i]-=minWeight;
		//myPrizes[i]*=-1;
		if(myPrizes[i]>0)
		{
			nTerminals++;
			realTerminals[i]=true;
		}
		if(myPrizes[i]>maxPrize)
		{
			maxPrize=myPrizes[i];
		}
	}
}

void Instance::readMWCS()
{

	std::ifstream infile(params.file);
	std::string line;

	std::string magicNumber;
	std::getline(infile, line);
	std::istringstream iss(line);
	iss.precision(16);

	iss>>magicNumber;

	if(magicNumber!="33D32945")
	{
		cout<<"Not in STP Format"<<endl;
		exit(0);
	}

	minWeight=0;
	int srcID, dstID;
	int edgeCounter=0;
	int termID;
	double prize;
	int termCounter=0;
	int count=0;
	int countRealTerminals=0;
	while(std::getline(infile, line))
	{
		if(line.empty())
			continue;

		iss.clear();
		iss.str(line);
		std::string token;

		iss>>token;
		//cerr<<count<<" "<<token<<" "<<line<<" "<<iss.str()<<endl;
		count++;
		if(token=="Nodes")
		{

			iss>>nNodes;
			adjList=vector<vector<int>>(nNodes);
			//cerr<<nNodes<<endl;
			adjList=vector<vector<int>>(nNodes);
			realTerminals=vector<bool>(nNodes,false);
			trueTerminals=vector<bool>(nNodes,false);
			myPrizes=vector<double>(nNodes);
			myBudgetCost=vector<double>(nNodes,1);
			nodesToRemove=vector<bool>(nNodes);
			for(int i=0;i<nNodes;i++)
			{
				adjList[i]=vector<int>();
			}

		}
		else if(token=="Edges")
		{
			iss>>nEdges;
			//cerr<<nEdges<<endl;

		}
		else if(token=="E")
		{
			iss>>srcID>>dstID;
			//cerr<<edgeCounter<<endl;
			adjList[srcID-1].push_back(dstID-1);
			adjList[dstID-1].push_back(srcID-1);
			//myCost[edgeCounter]=0;
			edgeCounter++;
			edgeCounter++;
		}
		else if(token=="Terminals")
		{
			iss>>nTerminals;
			myTerminals=vector<int>(nTerminals,-1);
		}
		else if(token=="T")
		{
			iss>>termID>>prize;
			myTerminals[termCounter]=termID-1;
			myPrizes[termID-1]=prize;
			//cout.precision(16);
			//cout<<(termID-1)<<" "<<prize<<endl;


			if(myPrizes[termID-1]>0)
			{
				realTerminals[termID-1]=true;
				countRealTerminals++;
			}

			if(myPrizes[termID-1]>maxPrize)
			{
				maxPrize=myPrizes[termID-1];
			}

			termCounter++;

		}
	}



	if(params.problem==2 && params.cardMultiplier>=1)
	{
		int nPositive=0;
		for(unsigned i=0;i<myPrizes.size();++i)
		{
			if(myPrizes[i]>0)
			{
				nPositive++;
			}
		}
		params.cardCons=floor(nPositive*params.cardMultiplier/100.0);

	}

	//cout<<"Real Terminals \t"<<countRealTerminals<<endl;

}


void Instance::readMWCSBudget()
{

	std::ifstream infile(params.file);
	std::string line;

	std::string magicNumber;
	std::getline(infile, line);
	std::istringstream iss(line);
	iss.precision(16);

	iss>>magicNumber;

	if(magicNumber!="33D32945")
	{
		cout<<"Not in STP Format"<<endl;
		exit(0);
	}

	minWeight=0;
	int srcID, dstID;
	int edgeCounter=0;
	int termID;
	double prize;
	int termCounter=0;
	int count=0;
	int countRealTerminals=0;
	int myMaxPrize=0;
	while(std::getline(infile, line))
	{
		if(line.empty())
			continue;

		iss.clear();
		iss.str(line);
		std::string token;

		iss>>token;
		//cerr<<count<<" "<<token<<" "<<line<<" "<<iss.str()<<endl;
		count++;
		if(token=="Nodes")
		{

			iss>>nNodes;
			adjList=vector<vector<int>>(nNodes);
			//cerr<<nNodes<<endl;
			adjList=vector<vector<int>>(nNodes);
			realTerminals=vector<bool>(nNodes,false);
			trueTerminals=vector<bool>(nNodes,false);
			myPrizes=vector<double>(nNodes);
			myBudgetCost=vector<double>(nNodes);
			nodesToRemove=vector<bool>(nNodes);
			for(int i=0;i<nNodes;i++)
			{
				adjList[i]=vector<int>();
			}

		}
		else if(token=="Edges")
		{
			iss>>nEdges;
			//cerr<<nEdges<<endl;

		}
		else if(token=="E")
		{
			iss>>srcID>>dstID;
			//cerr<<edgeCounter<<endl;
			adjList[srcID-1].push_back(dstID-1);
			adjList[dstID-1].push_back(srcID-1);
			//myCost[edgeCounter]=0;
			edgeCounter++;
			edgeCounter++;
		}
		else if(token=="Terminals")
		{
			iss>>nTerminals;
			myTerminals=vector<int>(nTerminals,-1);
		}
		else if(token=="T")
		{
			iss>>termID>>prize>>myMaxPrize;
			myTerminals[termCounter]=termID-1;
			myPrizes[termID-1]=prize;
			myBudgetCost[termID-1]=myMaxPrize;
			//cout.precision(16);
			//cout<<(termID-1)<<" "<<prize<<endl;
			if(myPrizes[termID-1]>0)
			{
				realTerminals[termID-1]=true;
				countRealTerminals++;
			}

			if(myPrizes[termID-1]>maxPrize)
			{
				maxPrize=myPrizes[termID-1];
			}

			termCounter++;

		}
	}


	if(params.problem==2 && params.cardMultiplier>=1)
	{
		int nPositive=0;
		for(unsigned i=0;i<myPrizes.size();++i)
		{
			if(myPrizes[i]>0)
			{
				nPositive++;
			}
		}
		params.cardCons=floor(nPositive*params.cardMultiplier/100.0);

	}


	if(params.problem==1 && params.budgetMultiplier>=1)
	{
		double totalCost=0;
		for(int i=0;i<nNodes;++i)
		{
			totalCost+=myBudgetCost[i];
			//cout<<myBudgetCost[i]<<" "<<myPrizes[i]<<endl;
		}

		budget=floor(totalCost*params.budgetMultiplier/100.0);
	}

	cout<<"Real Terminals \t"<<countRealTerminals<<endl;

}



void Instance::rebuildDatastructures()
{


	vector<double> myPrizes2;
	vector<double> myBudgetCost2;
	vector<bool> realTerminals2;
	vector<int> myTerminals2;
	vector<vector<int>> adjList2;

	
	vector<int> backMap=vector<int>(nNodes);

	int newNNodes=0;

	for(int i=0;i<nNodes;++i)
	{
		if(!nodesToRemove[i])
		{
			backMap[i]=newNNodes;
			map.push_back(i);
			adjList2.push_back(vector<int>());
			myPrizes2.push_back(myPrizes[i]);
			myBudgetCost2.push_back(myBudgetCost[i]);
			realTerminals2.push_back(realTerminals[i]);
			if(realTerminals[i])
				myTerminals2.push_back(newNNodes);
			newNNodes++;
		}
	}

	int newNEdges=0;
	for(int i=0;i<nNodes;++i)
	{
		if(!nodesToRemove[i])
		{
			for(unsigned j=0;j<adjList[i].size();++j)
			{
				if(!nodesToRemove[adjList[i][j]])
				{
					adjList2[backMap[i]].push_back(backMap[adjList[i][j]]);
					newNEdges++;
				}
			}
		}
	}

	newNEdges/=2;


	myPrizes=myPrizes2;
	myBudgetCost=myBudgetCost2;

	realTerminals=realTerminals2;
	myTerminals=myTerminals2;
	adjList=adjList2;

	nNodes=newNNodes;
	nEdges=newNEdges;
	nTerminals=myTerminals.size();

	trueTerminals=vector<bool>(nNodes,false);


	//componentArray.init(G);

}

int Instance::preprocessing()
{
	//return 0;

	int numberRemoved=0;
	numberRemoved+=uselessComponentsTest();
	numberRemoved+=degreeOneTest();
	numberRemoved+=degreeZeroTest();

	return numberRemoved;

}


int Instance::uselessComponentsTest()
{
	int numberRemoved=0;

	for(int i=0;i<nComponents;++i)
	{
		//cerr<<maxRevenueInComponent[i]<<" "<<maxPrize<<endl;
		if(maxRevenueInComponent[i]<=maxPrize)
		{
			numberRemoved+=components[i].size();
			for(unsigned j=0;j<components[i].size();++j)
			{
				nodesToRemove[components[i][j]]=true;
			}
		}

	}


	return numberRemoved;

}

int Instance::degreeOneTest()
{
	int numberRemoved=0;

	vector<int> toRemove;
	do
	{
		toRemove.clear();

		for(int n=0;n<nNodes;++n)
		{
			if((adjList[n].size()==1 && !realTerminals[n]))
			{
				toRemove.push_back(n);
			}
			/*if((adjList[n].size()==1 && realTerminals[n]) && params.problem==0)
			{
				toRemove.push_back(n);
				int adjacentNode=adjList[n][0];
				myPrizes[adjacentNode]+=myPrizes[n];
				if(myPrizes[adjacentNode]>0)
					realTerminals[adjacentNode]=true;
			}*/

		}

		for(unsigned i=0;i<toRemove.size();++i)
		{

			int node=toRemove[i];
			int adjacentNode=adjList[node][0];
			adjList[node].clear();

			if(adjList[adjacentNode].size()==0)
				continue;

			unsigned j=0;
			for(;j<adjList[adjacentNode].size();++j)
			{
				if(adjList[adjacentNode][j]==node)
				{
					break;
				}
			}
			//cout<<j<<" "<<adjList[adjacentNode].size()<<endl;
			//cout<<node<<" "<<adjList[adjacentNode][j]<<endl;
			adjList[adjacentNode].erase(adjList[adjacentNode].begin()+j);
		}

		numberRemoved+=toRemove.size();
	}
	while(toRemove.size()>0);

	return numberRemoved;

}

int Instance::degreeZeroTest()
{
	int numberRemoved=0;

	for(int n=0;n<nNodes;++n)
	{
		if(adjList[n].size()==0 && !nodesToRemove[n])
		{
			nodesToRemove[n]=true;
			numberRemoved++;
		}
	}
	return numberRemoved;
}


int Instance::calculateComponents()
{
	int numberOfComponents=0;

	componentArray=vector<int>(nNodes,0);
	components.clear();
	maxRevenueInComponent.clear();

	for(int n=0;n<nNodes;n++)
	{
		if(componentArray[n]==0)
		{
			//cerr<<n<<endl;
			vector<int> componentHelper;
			double revInComp=0.0;
			numberOfComponents++;
			componentArray[n]=numberOfComponents;
			//cout<<n<<endl;
			//cout<<myPrizes[n]<<endl;
			if(myPrizes[n]>0)
				revInComp+=myPrizes[n];
			componentHelper.push_back(n);
			queue<int> myQueue;
			myQueue.push(n);
			//cerr<<"component "<<n<<" "<<myCurrentLabel<<endl;
			while(!myQueue.empty())
			{
				int m=myQueue.front();
				myQueue.pop();
				for(int k:adjList[m])
				{
					if(componentArray[k]==0)
					{
						componentArray[k]=numberOfComponents;
						componentHelper.push_back(k);
						if(myPrizes[k]>0)
							revInComp+=myPrizes[k];
						myQueue.push(k);
					}
				}
			}
			components.push_back(componentHelper);
			maxRevenueInComponent.push_back(revInComp);
			componentFixed.push_back(0);
		}
	}

	/*cout<<"number of components "<<numberOfComponents<<endl;
	for(int i=0;i<numberOfComponents;++i)
	{
		cout<<"size: \t"<<components[i].size()<<"\t "<<maxRevenueInComponent[i]<<endl;
	}*/


	return numberOfComponents;
}

double Instance::transformInternalValue(double value) const
{
	//cout<<"value "<<value<<" "<<params.inputformat<<" "<<nRealTerminals<<" "<<bigM<<endl;
	if(params.inputformat==0)
		return (-1)*(value-sumPrizes+minWeight);
	else if(params.inputformat==1)
		return value;
	else if(params.inputformat==2)
		return value;
	else if(params.inputformat==3)
		return (-1)*(value+1);
	else if(params.inputformat==4)
		return value;
	return -1;
}


Instance::~Instance() {
	// TODO Auto-generated destructor stub
}

