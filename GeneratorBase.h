#ifndef GENERATORBASE_H
#define GENERATORBASE_H
#include <fstream>
#include <cmath>
#include <math.h>
#include <map>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <random>
#include <algorithm>
#include <numeric>
#include "TripsInfo.h"
#include "Macro.h"


using namespace std;

class GeneratorBase
{
public:

    typedef unsigned int Index;
    typedef set<pair<double,double>> CoordSet;
    typedef pair<double,double> Coord;

    enum Teta_Method
    {
        MaximumOptimal = 1,
        Optimal = 2
    };

    GeneratorBase();

    void ComputeNonCPNumber(TripsInfo& ti);

    bool isAReturnTrip(const int& tripID);

    void InitNodeLoad();

    void GenerateRequests(TripsInfo& ti,const TripsInfo& ti_pre,bool needtau = true);
    void GeneratePDRequests(TripsInfo& ti,const TripsInfo& ti_pre,bool needtau = true);
    void GenerateNPDRequests(TripsInfo& ti,const TripsInfo& ti_pre,bool needtau = true);
    void GeneratePNDRequests(TripsInfo& ti,const TripsInfo& ti_pre,bool needtau = true);
    void GenerateNPNDRequests(TripsInfo& ti,const TripsInfo& ti_pre,bool needtau = true);

    bool IsGeneratedCoordExist(const Coord& pos);

    void ComputeTau(TripsInfo& ti,vector<double>& reqtau,const double& start_x,const double& start_y);

    void GenerateCheckPoints(TripsInfo& ti,const TripsInfo& ti_pre);

    void GenerateNonCheckPoints(TripsInfo& ti);

    std::vector<int> ShuffleGeneration(int from,int to,int no);
    std::vector<int> ShuffleGenerationUnSorted(int from,int to,int no);
    std::vector<double> ShuffleGenerationReal(double from,double to,int no);
    std::vector<double> ShuffleGenerationUnSortedReal(double from,double to,int no);

    void PrintRawDataToFile();

    double doFixedPos(const double& p);

    void StartGenerationOnTimeHorizons();

    void doGenerationTrips();

    int ComputeNumberOfTrips(const double &timehorizon,const double& slackTime);

    void ComputeTripTeta(TripsInfo& ti,const double& slacktime);

public:

    int noPD_Req = 0;
    int noNPD_Req = 0;
    int noPND_Req = 0;
    int noNPND_Req = 0;

    int addedPD_Req = 0;
    int addedNPD_Req = 0;
    int addedPND_Req = 0;
    int addedNPND_Req = 0;

    int Input_CntInput = 0;

    int NonCPNo = 0;

    double TotalStartTime = 0.0;

    std::vector<Coord> VisitedCoord_Hybrid;
    std::vector<Coord> VisitedCoord_Hybrid_PND;

    vector<std::pair<Coord,Coord>> PD_RequestCoordinations;
    vector<std::pair<Coord,Coord>> PND_RequestCoordinations;
    vector<std::pair<Coord,Coord>> NPD_RequestCoordinations;
    vector<std::pair<Coord,Coord>> NPND_RequestCoordinations;

    vector <double> PD_RequestsTau;
    vector <double> PND_RequestsTau;
    vector <double> NPD_RequestsTau;
    vector <double> NPND_RequestsTau;

    //Setter
    void SetTripsNumber(const int& tno){Input_TripsNumber = tno;}
    void SetCheckpointsNumber(const int& cpno){Input_CheckpointsNumber = (2*cpno-1);}
    void SetStartTimeInterval(const double& interval){Input_StartTimeInterval = (interval);}
    void SetVehiclesCapacity(const int& c){Input_VehiclesCapacity = c;}
    void SetStopB(const double& b){Input_StopB = b;}
    void SetRequestNumber(const int& rno){Input_RequestNumber = rno;}
    void SetSuspendTime(const int& st){Input_SuspendTime = (st);}
    void SetTetaMethod(Teta_Method tm){Input_TetaMethod = tm;}
    void SetTetaIncrementalValue(const int& iv){Input_TetaIncrementalValue = iv;}
    void SetServiceAreaX(const int& x){Input_ServiceAreaX = x;}
    void SetServiceAreaY(const int& y){Input_ServiceAreaY = y;}
    void SetRequestPercentagePD(const int& pd){Input_PD = pd;}
    void SetRequestPercentagePND(const int& pnd){Input_PND = pnd;}
    void SetRequestPercentageNPD(const int& npd){Input_NPD = npd;}
    void SetRequestPercentageNPND(const int& npnd){Input_NPND = npnd;}
    void SetW0(const int& w0){Input_W0 = w0;}
    void SetW1(const int& w1){Input_W1 = w1;}
    void SetW2(const int& w2){Input_W2 = w2;}
    void SetSlackTime(const double& st){Input_SlackTime = st;}
    void SetTimeHorizons(const double& th){Input_TimeHorizon = th;}
    void SetNormal(const int& n)
    {
        if(n == 1)
            Input_IsNormalDistributaion = true;
        else
            Input_IsNormalDistributaion = false;
    }

private:
    int Input_TripsNumber = 0;
    int Input_CheckpointsNumber = 0;
    int Input_VehiclesCapacity = 0;
    int Input_RequestNumber = 0;
    int Input_SuspendTime = 0;
    Teta_Method Input_TetaMethod = MaximumOptimal;
    int Input_TetaIncrementalValue = 0;
    int Input_ServiceAreaX = 0;
    int Input_ServiceAreaY = 0;
    int Input_PD = 0;
    int Input_PND = 0;
    int Input_NPD = 0;
    int Input_NPND = 0;
    int Input_W0 = 0;
    int Input_W1 = 0;
    int Input_W2 = 0;

    int m_nbTrips = 0;
    int m_nbRequests = 0;

    double m_TripsGenerationTime = 0.0;

    double Input_StartTimeInterval = 0.0;
    double Input_StopB = 0;
    double Input_SlackTime = 0.0;
    double Input_TimeHorizon = 0.0;

    bool Input_IsNormalDistributaion = false;


private:

    vector<TripsInfo> GeneratorTrips;
};

#endif // GENERATORBASE_H
