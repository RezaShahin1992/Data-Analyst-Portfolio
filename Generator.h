#include <time.h>
#include "Initialization.h"
#include "Writer.h"
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
using namespace std;
#ifndef GENERATOR_H_
#define GENERATOR_H_

class Generator {
public:
	Generator();
	virtual ~Generator();
	vector<vector<int>>  DemandGenerator(int from, int to , int time, int scenario, Initialization *Ob1, Writer *wr);
	vector<vector<int>>  X1ProductionGenerator(int from, int to , int time, int scenario, Initialization *Ob1, Writer *wr);
	vector<vector<int>>  X2ProductionGenerator(int from, int to , int time, int scenario, Initialization *Ob1, Writer *wr);

};


#endif /* GENERATOR_H_ */
