#ifndef READERHELPER_H
#define READERHELPER_H

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include "SolverBase.h"
#include "unordered_map"
#include <unordered_set>


using namespace std;



class ReaderHelper
{
public:
    ReaderHelper();
    ~ReaderHelper()
    {
    }

    typedef unsigned int Index;
    typedef pair<double,double> Coord;

    typedef set<pair<double,double>> CoordSet;

    int m_ReservedTripHeu = 5;


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

    std::vector<VehicleContainer> VehicleContainerList;

    vector<TripsInfo> TripsInfoList;
    CoordSet m_ImpossibleCoord;

public:
    string AddInput = "INPUT/";
    string AddOutput = "OUTPUT/";
    string Address="";
    vector<vector<vector<double>>> X_WS; // Warm start list for X variable
    vector<vector<vector<double>>> Z_WS; // Warm start list for Z variable

    char* trash = new char;

    FILE *Input;
    FILE* Infrastructure;
    FILE* Instance;

    int TripsNo = 0;
    int LocalTripsNo = 0;
    int LocalTotalReqNo = 0;
    int nbPoints = 0;
    int VehiclesNo=0;
    int TC0,TC;
    int MaxTripsNo = 0;
    int MaxVehiclesNo = 0;
    vector<vector<double>> delta;
    vector <pair<int,int>> RequestsTypeIndex;
    set<int> VehiclesSet;
    vector<int> VehicleCapList;
    vector<int> GammaList;
    double GammaImpact = 0.0;

    double m_SlackTime = 0.0;
    double m_StopB = 0.0;
    double MinStopTime = 0.0;
    double RMax = 0.0;
    double WMax = 0.0;
    double HMax = 0.0; // W(max) variable in article
    int m_SAX = 0;
    int m_SAY = 0;
    double m_Interval = 0.0;
    int m_PD = 0;
    int m_PND = 0;
    int m_NPD = 0;
    int m_NPND = 0;
    double m_TripsGenerationTime = 0.0;


    //Base Data
    int m_CPno = 0;
    int m_VehicleCap = 0;
    double m_W0 = 0;
    double m_W1 = 0;
    double m_W2 = 0;

    std::map<Coord,int> Xmap;
    std::map<int,double> Tmap;
    std::map<int,double> TBarmap;
    std::map<int,double> Pmap;
    std::map<int,double> Dmap;
    std::map<Coord,int> Zmap;
    std::map<int,int> qmap;
    std::map<Coord,int> Qmap;

public:
    void ReadData(int index0,SolverBase* base = nullptr);
    void ReadOutputFile(int idx);
    void SendOutputDataToBase(SolverBase* base = nullptr);
    double doFixedPos(const double& p);
    int Read_Address(int& it);
    int Read_Address_Inf(int& it);
    int Read_Address_Ins(int& it);
    void ClearData();
    void split(std::string str, std::string splitBy, std::vector<std::string>& tokens);
    bool ImportNonCPToGraph(TripsInfo &ti,const Coord& noncpcoord);

    inline bool IsExists (const std::string& name) {
        std::ifstream infile(name);
        return infile.good();
    }


private:
    void ReadInfrastructure(const int &index0);
    void ReadInstance(const int &index0);
    void IncreaseTrips(const int& currentTripsNo,const int& targetTripsNo);
    void AddReservedTrips(const int& currentTripsNo,const int& targetTripsNo);
    void PreProcessInstanceFile(const int &index0);
    void NodeCreation(TripsInfo& ti);
    void RequestCreation(TripsInfo& ti);
    void GetNonCPFromInsFile(TripsInfo& ti);
};

#endif // READERHELPER_H
