#ifndef HEURISTIC_H
#define HEURISTIC_H

#include "Solver.h"
#include "SolverBase.h"
#include "Utilities.h"

using namespace std;

class Heuristic
{
public:
    Heuristic();
    Heuristic(Solver* s,SolverBase *sb);
    virtual ~Heuristic();

    typedef set<pair<double,double>> CoordSet;
    typedef vector<pair<double,double>> CoordList;
    typedef pair<double,double> Coord;
    typedef unsigned int Index;

    void doHeuristicSolve();
    void ApplyHeuToCplex();
    void CleanUnusedArcsMatrix();
    void LoadTrips();
    void Compute_q();
    void Compute_q_Trip(TripsInfo &ti);
    int GlobalToLocal(const int& globalIdx);
    void GetNearestCPToCurrentPos(const SolverBase* sb,const Coord& rootPos,vector<int>& cpIndices,CoordList& cpPosition);
    void CreateRawGraph();
    void SaveSolverTrips();
    void FindWarmStartParams();
    void CheckThetaCondition(TripsInfo& ti);
    void CheckCapacityCondition(TripsInfo& ti,const int& qMax);
    void HandlePDTrap(TripsInfo& ti,const int trapIdx);
    void UpdateParameterOnRemove(TripsInfo& ti,const int& ReqIdx,const int& RemovedIdx);
    void UpdateParameterOnRemoveNPND(TripsInfo& ti,const int& ReqIdx,const Coord& RemovedPosPick,const Coord& RemovedPosDrop);
    void SaveHybrid(TripsInfo& ti,const int& ReqIdx,const int& Pick,const int& Drop,const REQUESTTYPE& type);
    void SaveCpacityHybrid(TripsInfo& ti,const int& ReqIdx,const int& Pick,const int& Drop,const REQUESTTYPE& type);
    Coord ComputeRootPosInTrip(const Coord& rootPos,const int& lastTripIdx,const int& nextTripIdx);
    void RefreshRequests(TripsInfo &ti);
    void RefreshRequestsTotal();
    void RefreshCheckpointsIndices();
    void ApplyThetaTransmission(TripsInfo& ti);
    void ApplyCapacityTransmission(TripsInfo& ti,const int& qMax);
    double doFixedPos(const double& p);
    bool RemoveUnusedNodes(TripsInfo& ti,const Coord& rootPos);
    void UpdateGraph();
    void ComputeGlobalTheta();
    void PrepareGlobalParams();
    int GetNextCPIndexFromRoot(const Coord& rootPos);
    int GetNextCPIndexFromRoot(const int& Idx);
    int GetIndexFromList(const vector<Coord>& input,const Coord& coord);
    inline Coord GetCoordinationFromIndex(const int& idx);
    bool isEqual(const Coord& a,const Coord& b);
    int GetTripNumberOfNode(const Coord& root);

    const map<int,int>& GetZ(){return Zvar;}
    const map<int,int>& GetX(){return Xvar;}

    int maxTripId = -1;

    Utilities* utility = nullptr;

private:
    Solver* solver;
    SolverBase* base;

    Index TotalReq = 0;
    vector <double> StopsTeta;
    vector <double> HiddenTeta;
    vector <double> HiddenStopTeta;
    vector <double> StopsB;
    vector<int> TotalReqType;
    vector<int> Nodes;
    vector <int> RPS;
    vector <int> RDS;
    vector <int> HRPS;
    vector <int> HRDS;
    vector <double> RequestsTau;
    vector <int> StopsType;
    vector<Hybrid> m_HybridTransferList;
    vector<Hybrid> m_CapacityTransferList;
    CoordList GraphNodesList;
    CoordList RawGraphNodesList;
    CoordSet RawGraphNodes;
    vector<vector<int>> RejectedArcsMatrix; ///< A matrix of input arcs
    set<int> CPNodes;
    set<Coord> GraphCPCoords;
    set<Coord> RawGraphCPCoords;
    vector<TripsInfo> SolverTrips;
    map<int,int> Zvar;
    map<int,int> Xvar;
    map<int,int> Qvar;
    set<int> insertedRequests;
};

#endif // HEURISTIC_H
