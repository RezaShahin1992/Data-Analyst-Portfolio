#ifndef TRIPSINFO_H
#define TRIPSINFO_H

#include <set>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

class TripsInfo
{
public:
    TripsInfo();
    ~TripsInfo();

    typedef unsigned int Index;
    typedef pair<double,double> Coord;
    typedef set<pair<double,double>> CoordSet;

    int id = 0;
    int RequestNo=0, NodesNo=0, VehicleIndex=0,C=0;
    int VehicleCapacity = 0;
    int SuspendTime = 0;
    int SourcePoint = -1;
    int DestinationPoint = -1;
    int MiddleNode = 0;

    bool isNew = false;

    double StartTime = -1;
    double RealStartTime = 0;
    double RealArrivalTime = 0;
    double TripComputedTime = 0.0;

    set<int> Nodes;
    set<int> CP_Index;
    set<int> NonCP_Index;
    set<Coord> TripPickDropSetNonCP;
    set<pair<double,double>> CPCoords; ///< Set Of Checkpoint Nodes
    set<pair<double,double>> NonCPCoord; ///< Set Of Non Checkpoint Nodes
    set<pair<double,double>> NodesCoord; ///< A set of nodes

    vector<int> q;
    vector<int> NodesList;
    vector <int> StopsType;
    vector <double> StopsTeta;
    vector <double> HiddenTeta;
    vector <double> StopsB;
    vector <int> RequestsType;
    vector <int> RequestsPickStopNo;//ps(k)
    vector <int> RequestsDropStopNo;//ds(k)
    vector <int> HeuristicRequestsPickStopNo;//H_ps(k)
    vector <int> HeuristicRequestsDropStopNo;//H_ds(k)
    vector <double> RequestsTau;
    vector <int> RequestsXX1;
    vector <int> RequestsYY1;
    vector <int> RequestsXX2;
    vector <int> RequestsYY2;

    vector <double> W;
    vector<vector<int>> ArcsMatrix; ///< A matrix of input arcs
    vector<double> PickStopXList;
    vector<double> PickStopYList;
    vector<double> DropStopXList;
    vector<double> DropStopYList;

    vector<int> m_PossibleTripsForRequests;

     void ClearTripsInfo();
};

#endif // TRIPSINFO_H
