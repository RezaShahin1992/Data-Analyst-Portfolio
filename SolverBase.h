#ifndef SOLVERBASE_H
#define SOLVERBASE_H

#include "Macro.h"
#include "TripsInfo.h"
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


class SolverBase
{
public:
    SolverBase();

    typedef unsigned int Index;
    typedef pair<double,double> Coord;
    typedef vector<pair<double,double>> CoordList;

    typedef set<pair<double,double>> CoordSet;

    struct Hybrid
    {
    public:
        Coord StartPoint;
        Coord EndPoint;
        Coord CurrentPoint;
        int Pick = -1;
        int Drop = -1;
        int R_Type = -1;
        double R_Tau = -1;
        bool isCP = false;
    };

    class VehicleContainer
    {
    public:
        //Default Constructor
        VehicleContainer(){}

        VehicleContainer(bool inUse){InUse = inUse;}

        VehicleContainer(bool inUse,int idx){InUse = inUse;V_Index = idx;}

        bool InUse = false;
        int V_Index = 0; // Vehicle Index

    };

    //Setter
    void Set_SA_X(const int& sx){Input_ServiceAreaX = sx;}
    void Set_SA_Y(const int& sy){Input_ServiceAreaY = sy;}
    void Set_CP_Number(const int& cpno){Input_CheckpointsNumber = cpno;}
    void Set_VehicleCap(const int& vc){Input_VehiclesCapacity = vc;}
    void Set_Interval(const double& interval){Input_StartTimeInterval = interval;}
    void Set_PD_Per(const int& pd){Input_PD = pd;}
    void Set_PND_Per(const int& pnd){Input_PND = pnd;}
    void Set_NPD_Per(const int& npd){Input_NPD = npd;}
    void Set_NPND_Per(const int& npnd){Input_NPND = npnd;}
    void Set_StopB(const double& b){Input_StopB = b;}
    void Set_SlackTime(const double& slack){Input_SlackTime = slack;}
    void Set_W0(const double& w0){Input_W0 = w0;}
    void Set_W1(const double& w1){Input_W1 = w1;}
    void Set_W2(const double& w2){Input_W2 = w2;}
    void Set_XWS(const vector<vector<vector<double> > >& xws){X_WS = xws;}
    void Set_ZWS(const vector<vector<vector<double> > >& zws){Z_WS = zws;}

    void SetTripsInfoList(vector<TripsInfo>& tripsinfo){TripsInfoList = tripsinfo;}

    void execute();

    double CalculateDistance(const Coord &a, const Coord &b);
    double doFixedPos(const double& p);

    //Private Function
private :
    void GenerateCheckPoints(TripsInfo& ti,const TripsInfo& ti_pre);
    void ComputeNonCPNumber(TripsInfo& ti);
    void GenerateNonCheckPoints(TripsInfo& ti,const TripsInfo& ti_pre);
    void InitNodeLoad(const TripsInfo& ti);
    void ComputeTripTeta(TripsInfo& ti,const double& slacktime);
    void ComputeGlobalTeta(const double& slacktime);
    void ComputeHiddenTeta(TripsInfo& ti);
    void GenerateRequests(TripsInfo& ti,const TripsInfo& ti_pre,bool needtau = true);
    void GeneratePDRequests(TripsInfo& ti,const TripsInfo& ti_pre,bool needtau = true);
    void GenerateNPDRequests(TripsInfo& ti,const TripsInfo& ti_pre,bool needtau = true);
    void GeneratePNDRequests(TripsInfo& ti,const TripsInfo& ti_pre,bool needtau = true);
    void GenerateNPNDRequests(TripsInfo& ti,const TripsInfo& ti_pre,bool needtau = true);
    void TripClearPreProcess(TripsInfo &ti);
    void TripPreProcess(TripsInfo& ti);
    void TripCompute_q(TripsInfo& ti);
    void ClearPreProcess();
    bool Handle_PND_NPD_Condition(TripsInfo &ti);
    void doVehiclesGeneration(const int& nbTrips);
    void GenerateVehicles(TripsInfo& ti);
    std::vector<int> ShuffleGeneration(int from,int to,int no);
    std::vector<int> ShuffleGenerationUnSorted(int from,int to,int no);
    std::vector<double> ShuffleGenerationReal(double from,double to,int no);
    std::vector<double> ShuffleGenerationUnSortedReal(double from,double to,int no);
    void ComputeTau(TripsInfo& ti,const TripsInfo& ti_pre,const double& start_x,const double& start_y);
    bool isAReturnTrip(const int& tripID);
    int GetIndexFromSet(const CoordSet& input,const Coord& coord);

    Coord GetCoordinationFromIndex(const int& idx);
    /**
     * @brief isTripSourceCoord
     * @param p
     * @return if p as a source coord of a trip we return current trip idx
     */
    int isTripSourceCoord(const Coord& p);
public:
//    vector<TripsInfo>

    int Input_ServiceAreaX = 0;
    int Input_ServiceAreaY = 0;
    int Input_CheckpointsNumber = 0;
    double TotalStartTime = 0.0;
    int Input_VehiclesCapacity = 0;
    double Input_StartTimeInterval = 0.0;
    int Input_PD = 0;
    int Input_PND = 0;
    int Input_NPD = 0;
    int Input_NPND = 0;
    double Input_StopB = 0;
    double Input_SlackTime = 0.0;
    double Input_W0;
    double Input_W1;
    double Input_W2;
    int VehiclesNo=0;
    int nbPoints = 0;
    int TC0,TC;
    vector<vector<double>> delta;
    double m_TripsGenerationTime = 0.0;
    double m_SlackTime = 0.0;
    double m_StopB = 0.0;

    int noPD_Req = 0;
    int noNPD_Req = 0;
    int noPND_Req = 0;
    int noNPND_Req = 0;
    int NonCPNo = 0;

    int TripsNo = 0;

    std::map<Coord,int> Xmap;
    std::map<int,double> Tmap;
    std::map<int,double> TBarmap;
    std::map<int,double> Pmap;
    std::map<int,double> Dmap;
    std::map<Coord,int> Zmap;
    std::map<int,int> qmap;
    std::map<Coord,int> Qmap;

    vector<Hybrid> m_InvalidPos;
    vector<TripsInfo> TripsInfoList;
    set<int> VehiclesSet;
    std::map<int,double> VehicleFreeTime;
    std::vector<VehicleContainer> VehicleContainerList;
    vector <pair<int,int>> RequestsTypeIndex;
    vector<vector<vector<double>>> X_WS; // Warm start list for X variable
    vector<vector<vector<double>>> Z_WS; // Warm start list for Z variable

    std::set<int> VisitedIndices_Hybrid;
    vector <double> StopsTeta;
    vector <double> StopsB;

    vector<std::pair<Coord,Coord>> PD_RequestCoordinations;
    vector<std::pair<Coord,Coord>> PND_RequestCoordinations;
    vector<std::pair<Coord,Coord>> NPD_RequestCoordinations;
    vector<std::pair<Coord,Coord>> NPND_RequestCoordinations;

    vector<pair<int,pair<Coord,Coord> > > TripsHeadTailList;
};

#endif // SOLVERBASE_H
