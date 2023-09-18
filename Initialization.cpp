#include "Initialization.h"
#include <ilcplex/ilocplex.h>
#include <time.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "Macro.h"

using namespace std;

Initialization::Initialization() {

}

Initialization::~Initialization() {
}


/////////////////////////////////Let's start doing something

void Initialization::Init()
{
	Initialization inite;
	cout << " Initializing has started " << endl;
	CreateProblem();
	AllocateMemory();
	ObjectiveFunction();
	this->DefineConstraints();

	//inite.DefineConstraints(D);


}
void Initialization::CreateProblem()
{
	cout << " CreateProblem has started " << endl;

    lp = new IloModel(env);
    if (lp == NULL)
    {
        cout << "Cannot create problem!!" << endl;
        exit(502);
    }
}
void Initialization::AllocateMemory()
{

	cout << " AllocateMemory has started " << endl;

	////initialization of parameter D



    /////Initialize variables
    X1 = new IloArray<IloNumVarArray> (env, numScenarios);
	{
    	for(s=0; s< numScenarios; s++)
    	{
        	(*X1)[s] = IloNumVarArray(env, numTimePeriods);
        	for(t=0; t< numTimePeriods; t++)
        	{
                string valname = "X1("+to_string(s)+")("+to_string(t)+")";
                (*X1)[s][t]= IloNumVar(env, 0, M, IloNumVar::Float,valname.c_str());
        	}
    	}
	}

////////////////////////////////////////////////////////////////

    X2 = new IloArray<IloNumVarArray> (env, numScenarios);
	{
    	for(s=0; s< numScenarios; s++)
    	{
        	(*X2)[s] = IloNumVarArray(env, numTimePeriods);
        	for(t=0; t< numTimePeriods; t++)
        	{
                string valname = "X2("+to_string(s)+")("+to_string(t)+")";
                (*X2)[s][t]= IloNumVar(env, 0, M, IloNumVar::Float,valname.c_str());
        	}

    	}
	}
///////////////////////////////////////////////////////////////
    I = new IloArray<IloNumVarArray> (env, numScenarios);
	{
    	for(s=0; s< numScenarios; s++)
    	{
        	(*I)[s] = IloNumVarArray(env, numTimePeriods);
        	for(t=0; t< numTimePeriods; t++)
        	{
                string valname = "I("+to_string(s)+")("+to_string(t)+")";
                (*I)[s][t]= IloNumVar(env, 0, M, IloNumVar::Float,valname.c_str());
        	}

    	}
	}
/////////////////////////////////////////////////////////////////
    S = new IloArray<IloNumVarArray> (env, numScenarios);
	{
    	for(s=0; s< numScenarios; s++)
    	{
        	(*S)[s] = IloNumVarArray(env, numTimePeriods);
        	for(t=0; t< numTimePeriods; t++)
        	{
                string valname = "S("+to_string(s)+")("+to_string(t)+")";
                (*S)[s][t]= IloNumVar(env, 0, M, IloNumVar::Float,valname.c_str());
        	}
    	}
	}
/////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
    Y = new IloArray<IloNumVarArray> (env, numScenarios);
	{
    	for(s=0; s< numScenarios; s++)
    	{
        	(*Y)[s] = IloNumVarArray(env, numTimePeriods);
        	for(t=0; t< numTimePeriods; t++)
        	{
                string valname = "Y("+to_string(s)+")("+to_string(t)+")";
                (*Y)[s][t]= IloNumVar(env, 0, 1, IloNumVar::Bool,valname.c_str());
        	}
    	}
	}

}
void Initialization::ObjectiveFunction()
{
	cout << " DefineObjectiveFunction has started " << endl;

	/////////define objective function
    IloExpr objective(env);

    //////////first component
    for(s=0; s<numScenarios; s++)
    {
    	for(t=0; t<numTimePeriods; t++)
    	{
    		objective+=(C1* (*X1)[s][t]);
    	}
    }
    //////////second component
    for(s=0; s<numScenarios; s++)
    {
    	for(t=0; t<numTimePeriods; t++)
    	{
    		objective.operator+=(C2* (*X2)[s][t]);
    	}
    }
    ////////third component
    for(s=0; s<numScenarios; s++)
    {
    	for(t=0; t<numTimePeriods; t++)
    	{
    		objective.operator+=(C3* (*I)[s][t]);
    	}
    }
    ///////fourth component
    for(s=0; s<numScenarios; s++)
    {
    	for(t=0; t<numTimePeriods; t++)
    	{
    		objective.operator+=(C4* (*S)[s][t]);
    	}
    }

    IloObjective obj(env, objective, IloObjective::Minimize);

    lp->add(obj);
    objective.end();

}

void Initialization::DefineConstraints()
{


	cout << "DefineConstraints has started" << endl;

    for (int s = 0; s < numScenarios; s++)
    {
          for (int t = 0; t < numTimePeriods; t++)
          {
              cout << "Scenario " << s << ", Period " << t << ", Demand: " << D[s][t] << std::endl;
          }
      }

	cout << "DefineConstraints has started 1" << endl;

	////////////////////////C2
    for (int s = 0; s < numScenarios; s++)
    {
        for (int t = 1; t < numTimePeriods; t++)
        {
            // Add constraints to ensure Y[s][t] is 1 if D[s][t] <= (*I)[s][t-1] + (*X1)[s][t-1] + (*X2)[s][t-1]
            lp->add((*Y)[s][t] * M >= (*I)[s][t-1] + (*X1)[s][t-1] + (*X2)[s][t-1] - D[s][t]);
        }
    }
	cout << "DefineConstraints has started 2" << endl;

////////////////////////C3
/*    for (int s = 0; s < numScenarios; s++)
    {
        for (int t = 1; t < numTimePeriods; t++)
        {
            // Add constraints to ensure Y[s][t] is 1 if D[s][t] <= (*I)[s][t-1] + (*X1)[s][t-1] + (*X2)[s][t-1]
            lp->add((1 - (*Y)[s][t]) * M >= - D[s][t] + ((*I)[s][t-1] + (*X1)[s][t-1] + (*X2)[s][t-1]));
        }
    }
	cout << "DefineConstraints has started 3" << endl;
*/
//////////////////////////C4
    for (int s = 0; s < numScenarios; s++)
    {												//// line 220 is commented for my version / for Martin's version line 221 should be activated
        for (int t = 1; t < numTimePeriods; t++)
        {
            // Add constraint: (*I)[s][t] == (*I)[s][t-1] + (*X1)[s][t-1] + (*X2)[s][t-1] - D[s][t] if Y[s][t] == 1
       //     lp->add((*I)[s][t] - ((*I)[s][t-1] - (*X1)[s][t-1] - (*X2)[s][t-1]) + D[s][t] <= M * (1 - (*Y)[s][t]));
        	lp->add( (*X1)[s][t] + (*X2)[s][t] + (*I)[s][t-1] + (*S)[s][t] >= D[s][t] + (*I)[s][t]   );
        }
    }
////////////////////////C5
/*    for (int s = 0; s < numScenarios; s++)
    {
        for (int t = 1; t < numTimePeriods; t++)
        {
            // Add constraint: (*I)[s][t] == 0 if Y[s][t] == 0
            lp->add((*I)[s][t] <= M * (*Y)[s][t]);
        }
    }
*/
	cout << "DefineConstraints has started 4" << endl;

	////////////////////////C6 & C7
    for (int s = 0; s < numScenarios; s++) {
        for (int t = 1; t < numTimePeriods; t++) {
            // Add constraint: (*S)[s][t] >= D[s][t] - (*I)[s][t-1] - (*X1)[s][t-1] - (*X2)[s][t-1] if Y[s][t] == 0
            lp->add((*S)[s][t] - (D[s][t] - (*I)[s][t-1] - (*X1)[s][t-1] - (*X2)[s][t-1]) >= -M * (*Y)[s][t]);

            // Add constraint: (*S)[s][t] == 0 if Y[s][t] == 1
            lp->add((*S)[s][t] <= M * (1 - (*Y)[s][t]));
        }
    }
	cout << "DefineConstraints has started 5" << endl;

    for (int s = 0; s < numScenarios - 1; s++)
    {
    	for(int t=0; t < numTimePeriods-1; t++)
    	{
    		if(t==0)
    			lp->add((*I)[s][t] == (*I)[s+1][t]);
    	}

    }
	cout << "DefineConstraints has started 6" << endl;

/*    for (int s = 0; s < numScenarios; s++)
    {
    	for(int s_prime=0; s_prime < numScenarios; s_prime++ )
    	{
    		if(s!=s_prime)
    		{
    			for (int t = 0; t < numTimePeriods; t++)
    			{
    				if (t==1)
    				{
    					lp->add( (*I)[s][t] == (*I)[s_prime][t]);
    				}
    			}
    		}
    	}
     }*/
	////////////////////////C9
    for(int s=0; s< numScenarios; s++)
    {
    	for(int t=1; t< numTimePeriods; t++)
    	{
    		lp->add((*X1)[s][t] <= X1max/*[s][t]*/);
    	}
    }
	cout << "DefineConstraints has started 7" << endl;

	////////////////////////C10
    for(int s=0; s< numScenarios; s++)
    {
    	for(int t=1; t< numTimePeriods; t++)
    	{
    		lp->add((*X2)[s][t] <= X2max/*[s][t]*/);
    	}
    }


	cout << "DefineConstraints has finished" << endl;

	PrepareCplexParams();
}

void Initialization::PrepareCplexParams()
{
	cout << "preparing Cplex Parameters has started" << endl;

    if(cplexPtr == nullptr)
    {
        cplexPtr = new IloCplex(*lp);
        std::string ModelName = std::string("Model-") + ".lp";
        cplexPtr->exportModel(ModelName.c_str());
        double tL = 3600;
        cplexPtr->setParam(IloCplex::TiLim, tL);		        //Run Time = 1 Hours
        cplexPtr->setParam(IloCplex::EpGap, 0.0);
    }
}


void Initialization::Solve()
{
	cout << "solving has started" << endl;

    IloEnv env;
    IloModel model(env);
	cout << "solving has started 2" << endl;

	   try{
	        double t1 = clock();
	        status = cplexPtr->solve();						//cout<<"Solving . . .";
	        objval = cplexPtr->getObjValue();
	        double LB = cplexPtr->getBestObjValue();
	        gap_p = (objval-LB)/objval;
	    } catch(const IloException& e) {
	        std::cerr << "\n\nCPLEX Raised an exception:\n";
	        std::cerr << e << "\n";
	    } catch(const std::exception& e) {
	        std::cerr << "\n\nStandard exception:\n";
	        std::cerr << e.what() << "\n";
	    } catch(...) {
	        std::cerr << "\n\nAn exception of unknown type was thrown.\n";
	    }

		cout << "solving has started 3" << endl;


	    if (!status)
	    {
	        std::cerr << "\n\nCplex error!\n";
	        std::cerr << "\tStatus: " << cplexPtr->getStatus() << "\n";
	        std::cerr << "\tSolver status: " << cplexPtr->getCplexStatus() << "\n";
	    }
	    else
	    {
	        std::cout << "\n\nCplex success!\n";
	        std::cout << "\tStatus: " << cplexPtr->getStatus() << "\n";
	        std::cout << "\tObjective value: " << cplexPtr->getObjValue() << "\n";
	    }

	    cout << " solving has finished " << endl;


	    int objective1=0;
	    int objective2=0;
	    int objective3=0;
	    int objective4=0;

	    try{
	    for (int s = 0; s < numScenarios; s++)
	    {
	        for (int t = 0; t < numTimePeriods; t++)
	        {
	                    objective1 += (C1*cplexPtr->getValue((*X1)[s][t]));
	        }
	    }
	    }catch(const IloException &e)
	    {
	        std::cerr << __FUNCTION__ << ":X1 Parameters :: " << e << "\n";
	    }
	    cerr << __FUNCTION__ << ": objective 1 :: " << objective1 << endl;


	    try{
		    for (int s = 0; s < numScenarios; s++)
		    {
		        for (int t = 0; t < numTimePeriods; t++)
		        {
		        	objective2 += (C2*cplexPtr->getValue((*X2)[s][t]));
		        }
		    }
		    }catch(const IloException &e)
		    {
		        std::cerr << __FUNCTION__ << ":X2 Parameters :: " << e << "\n";
		    }
		    cerr << __FUNCTION__ << ": objective 2 :: " << objective2 << endl;


		    try{
			    for (int s = 0; s < numScenarios; s++)
			    {
			        for (int t = 0; t < numTimePeriods; t++)
			        {
			        	objective3 += (C3*cplexPtr->getValue((*I)[s][t]));
			        }
			    }
			    }catch(const IloException &e)
			    {
			        std::cerr << __FUNCTION__ << ":I Parameters :: " << e << "\n";
			    }
			    cerr << __FUNCTION__ << ": objective 3 :: " << objective3 << endl;


			    try{
			 	    for (int s = 0; s < numScenarios; s++)
			 			    {
			 		        for (int t = 0; t < numTimePeriods; t++)
			 			        {
			 			        	objective4 += (C3*cplexPtr->getValue((*S)[s][t]));
			 			        }
			 			    }
			 			    }catch(const IloException &e)
			 			    {
			 			        std::cerr << __FUNCTION__ << ":S Parameters :: " << e << "\n";
			 			    }
			 			    cerr << __FUNCTION__ << ": objective 4 :: " << objective4 << endl;

}
void Initialization::FreeProblem()
{
	cout << "FreeProblem has started" << endl;

    lp->end();
    env.end();
    if (lp != NULL)
    {
        delete lp;
    }

}

