#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <string>

using namespace std;
#ifndef INITIALIZATION_H_
#define INITIALIZATION_H_

class Initialization {
	private:
//	vector<vector<int>> D;
	public:
		Initialization();
		virtual ~Initialization();

		////////////////////CPLEX needed
	    IloModel *lp;
	    IloEnv env;
	    double gap_p;
	    double  objval;
	    IloCplex *cplexPtr = nullptr;
	    IloBool status;


	    ////////////////Sets
	    int t=0;
	    int s=0;
	    int d=0;
	    //////////////////////Parameters
//	    double C1=1;					GET VALUE IN MACRO.h
//	    double C2=1.5;
//	    double C3=0.5;
//	    double C4=10;
	    int P;
//	    int X1max=10, X2max=10;
	    int numScenarios, numTimePeriods;


	    // Declare variables
	    IloArray<IloNumVarArray> *X1;
	    IloArray<IloNumVarArray> *X2;
	    IloArray<IloNumVarArray> *I;
	    IloArray<IloNumVarArray> *S;
	    IloArray<IloNumVarArray> *Z1;
	    IloArray<IloNumVarArray> *Z2;
	    IloArray<IloNumVarArray> *Y;


	    //////pointer to another class member
	    vector<vector<int>> D;
	    //vector<vector<int>> X1max;
	    //vector<vector<int>> X2max;
	    int X1max;
	    int X2max;
		////////////Functions
		void Init();
		void CreateProblem();
		void AllocateMemory();
		void ObjectiveFunction();
		void DefineConstraints();
		void Solve();
		void FreeProblem();
		void PrepareCplexParams();

//const vector<vector<int>>& D
};





#endif /* INITIALIZATION_H_ */
