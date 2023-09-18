#include "Writer.h"
#include <time.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <string>
#include <map>
#include <cmath>
#include <math.h>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <numeric>
#include "fstream"
#include "iomanip"
#include <sstream>
#include "unordered_map"
#include <unordered_set>
#include "Macro.h"

using namespace std;

Writer::Writer() {
}

Writer::~Writer() {
}

void Writer::WritingData(const Initialization *initialization,const string &outputname)
{


	cout << " Writing Data has started " << endl;

//    const char* OutName = new char;

    std::string StrInfName = "OUTPUT/" + outputname + "-" + ".txt";

    std::ofstream OutFile(StrInfName, std::ios::out); // Open the file for writing

    if (!OutFile.is_open()) {
        std::cerr << "Failed to open the file!" << std::endl;
        return; // or handle the error in a way that makes sense for your application
    }


    ////////////////////////////////////////////////

//    fprintf(Out, "Obj:	%f\n", initialization->objval);


    //////////////start printing


/*	    for (int s = 0; s < numScenarios; s++) {
	        for (int t = 1; t < numTimePeriods; t++) {
	            IloNum x = initialization->cplexPtr->getValue((*initialization->X1)[s][t]);
	            OutFile << "X1[" << s << "][" << t << "]= " << x << "\n";
	        }
	    }
*/
	    try {
	        for (int s = 0; s < numScenarios; s++) {
	            for (int t = 0; t < numTimePeriods; t++) {
	                double value = initialization->cplexPtr->getValue((*initialization->X1)[s][t]);
	                OutFile << "X1[" << s << "][" << t << "] = " << value << endl;
	            }
	        }
	    } catch(const IloException &e) {
	        cerr << "Error getting X1 values: " << e << "\n";
	    }


	///////////////////////////////Print X2
/*	try {
	    for (int s = 0; s < numScenarios; s++) {
	        for (int t = 1; t < numTimePeriods; t++) {
	         //   IloNum value = initialization->cplexPtr->getValue((*initialization->X2)[s][t]);
	        //    OutFile << "X2[" << s << "][" << t << "]= " << value << "\n";
                OutFile << "X2[" << s << "][" << t << "]" << "=" << initialization->cplexPtr->getValue((*initialization->X2)[s][t]) << endl;

	        }
	    }
	} catch(const IloException& e) {
	    std::cerr << __FUNCTION__ << ":X2 Parameters :: " << e << "\n";
	}
*/
    try {
        for (int s = 0; s < numScenarios; s++) {
            for (int t = 0; t < numTimePeriods; t++) {
                double value = initialization->cplexPtr->getValue((*initialization->X2)[s][t]);
                OutFile << "X2[" << s << "][" << t << "] = " << value << endl;
            }
        }
    } catch(const IloException &e) {
        cerr << "Error getting X2 values: " << e << "\n";
    }

    ////////////////////////////////Print I
    try{
        for (int s = 0; s < numScenarios; s++)
        {
            for (int t = 0; t < numTimePeriods; t++)
            {
                OutFile << "I[" << s << "][" << t << "]" << "=" << initialization->cplexPtr->getValue((*initialization->I)[s][t]) << endl;
//                fprintf(Out, "I[%d][%d]= %f\n", s, t, initialization->cplexPtr->getValue((*initialization->I)[s][t]));
            }
        }

    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":I Parameters :: " << e << "\n";
    }
    ///////////////////////////////PRINT S
    try{
        for (int s = 0; s < numScenarios; s++)
        {
            for (int t = 1; t < numTimePeriods; t++)
            {
//                fprintf(Out, "S[%d][%d]= %f\n", s, t, initialization->cplexPtr->getValue((*initialization->S)[s][t]));
                OutFile << "S[" << s << "][" << t << "]" << "=" << initialization->cplexPtr->getValue((*initialization->S)[s][t]) << endl;
            }
        }

    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":S Parameters :: " << e << "\n";
    }
    ///////////////////////////////Printing Y
    try{
        for (int s = 0; s < numScenarios; s++)
        {
            for (int t = 1; t < numTimePeriods; t++)
            {
//                fprintf(Out, "S[%d][%d]= %f\n", s, t, initialization->cplexPtr->getValue((*initialization->S)[s][t]));
                OutFile << "Y[" << s << "][" << t << "]" << "=" << initialization->cplexPtr->getValue((*initialization->Y)[s][t]) << endl;
            }
        }

    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":Y Parameters :: " << e << "\n";
    }



///////////////////////////////////Finish printing


    cout << "numScenarios: " << numScenarios << ", numTimePeriods: " << numTimePeriods << endl;
/*
    try {
        for (int s = 0; s < numScenarios; s++) {
            for (int t = 0; t < numTimePeriods; t++) {
                double value = initialization->cplexPtr->getValue((*initialization->X1)[s][t]);
                cout << "X1[" << s << "][" << t << "] = " << value << endl;
            }
        }
    } catch(const IloException &e) {
        cerr << "Error getting X1 values: " << e << "\n";
    }
*/

//    fclose(Out);
    cout << "objval :: " << initialization->objval << endl;

    cout << "writing data has finished " << endl;
}
