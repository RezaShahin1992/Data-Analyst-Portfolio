#ifndef UTILITIES_H
#define UTILITIES_H

#include "Macro.h"
#include "SolverBase.h"
#include "TripsInfo.h"

class Utilities
{
public:
    Utilities();
    ~Utilities(){}

    typedef set<pair<double,double>> CoordSet;
    typedef vector<pair<double,double>> CoordList;
    typedef pair<double,double> Coord;
    typedef unsigned int Index;


    void ComputeTripTime(SolverBase *sb);
    int GetNbReqFromType(const SolverBase* sb,const REQUESTTYPE &type);
    int GetNbReqFromRPS(const SolverBase* sb,const int ps,const REQUESTTYPE &type);
    int GetNbReqFromRDS(const SolverBase* sb, const int ds,const REQUESTTYPE &type);
    int GetNbReqFromPC(const vector<vector<vector<int> > > &S_PND, const SolverBase* sb, const int pick);
    int GetNbReqFromDC(const vector <vector<vector<int>>>& S_NPD,const SolverBase* sb, const int drop);
    int GetTripNumberOfRequest(const SolverBase* sb,const int& reqID);
    Coord GetCoordinationFromIndex(const SolverBase* sb,const int& idx);
    int GetCurrentRequestNo(const SolverBase* sb,const REQUESTTYPE& rt);
    vector<Index> GetRequestsInNode(const vector <vector<vector<int>>>& S_PND,const SolverBase* sb,const int& idx);
    /**
     * @brief Compute the minimum amount of time needed to go from pickup to pickup point
     * @param sb
     * @param a
     * @param b
     * @return minimum delta from Pick to Drop without any deviation
     */
    double ComputeDeltaMinReqFromPickToDrop(const set<int>& CPNodes,const SolverBase* sb,const int& PickUp,const int& DropOff);
    //
    string ConvertToRealTime(const double& time);
    double CalculateDistance(const Coord& a,const Coord& b);
    vector<int> GetSameCPIndexFromRoot(const SolverBase* sb,const CoordList& GraphNodesList,const Coord& rootPos);
    CoordList GetSimilarPositionOfNode(const SolverBase* sb,const Coord& rootPos);
    Coord GetSimilarPositionOfNodeInTrip(const SolverBase* sb,const int& TripIdx,const Coord& rootPos);
    vector<int> GetNearestCPToCurrentPos(const CoordList& GraphCPGlobalList, const CoordSet& GraphNodes, const Coord& rootPos);
    Index GetVehicleNumberOfCurrentNode(const SolverBase* sb,const Index& idx);
    int GetNextCPIndexFromRoot(const set<int> CPNodes,const CoordList& GraphNodesList,const Coord& rootPos);
    /**
     * @brief map source coordination to destination trip and compute distance from des to src
     * @param sb
     * @param src
     * @param des
     * @return distance
     */
    double GetDistanceFromSrcToDes(const SolverBase* sb,const int& src,const int& des);
    int GetIndexFromSet(const CoordSet& input,const Coord& coord);
    int GetIndexFromList(const vector<Coord>& input,const Coord& coord);
    bool isEqual(const double& a,const double& b);
    bool isEqual(const Coord& a,const Coord& b);
    double doFixedPos(const double& p);

};

#endif // UTILITIES_H
