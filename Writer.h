#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <string>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "Initialization.h"
using namespace std;

#ifndef WRITER_H_
#define WRITER_H_

class Writer {
public:
	Writer();
	virtual ~Writer();
    vector<vector<int>> D;
    vector<vector<int>> X1max;
    vector<vector<int>> X2max;

    void WritingData(const Initialization *initialization,const string& outputname);


    FILE* Out;
    string AddOutput = "OUTPUT/";

};

#endif /* WRITER_H_ */
