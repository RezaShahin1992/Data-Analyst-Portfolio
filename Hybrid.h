#ifndef HYBRID_H
#define HYBRID_H
#include "Macro.h"
#include <map>
#include <cmath>
#include <math.h>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <set>


using namespace std;

class Hybrid
{
public:
    Hybrid();

    typedef unsigned int Index;
    typedef pair<double,double> Coord;

    typedef set<pair<double,double>> CoordSet;

public:
    Coord StartPoint;
    Coord EndPoint;
    Coord CurrentPoint;
    int nextTrip = 0;
    int currentTrip = 0;
    int Pick = -1;
    int Drop = -1;
    int R_Type = -1;
    double R_Tau = -1;
    bool isCP = false;
};

#endif // HYBRID_H
