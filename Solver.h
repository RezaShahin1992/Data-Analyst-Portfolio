#ifndef SOLVER_H
#define SOLVER_H

//#define _CRT_SECURE_NO_WARNINGS
//#define _USE_MATH_DEFINES

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include "ReaderHelper.h"
#include "fstream"
#include "iomanip"
#include "Macro.h"
#include "SolverBase.h"
#include "Hybrid.h"
#include "Utilities.h"

class Solver
{
public:
    Solver();
    ~Solver(){}

    typedef set<pair<double,double>> CoordSet;
    typedef vector<pair<double,double>> CoordList;
    typedef pair<double,double> Coord;
    typedef unsigned int Index;

    IloBool status; // integer... for what??
    int cntInputFile = 0;

    double gap_p;
    double  objval;

    double TotalTravelTime = 0.0;
    double TotalRequestTime = 0.0;
    double TotalRequestWaitingTime = 0.0;
    double SumRej = 0.0;

    bool IsSecondRun = false;

    IloEnv env;
    IloModel *lp;
    IloCplex *cplexPtr = nullptr;
    Utilities *utility = nullptr;
    
    IloArray<IloArray<IloNumVarArray>> * Xval;
    IloNumVarArray * Tval;
    IloNumVarArray * Tbarval;
    IloNumVarArray * Pval;
    IloNumVarArray * Dval;
    IloNumVarArray * qval;
    IloArray<IloNumVarArray> * Qval;
    IloNumVarArray * Eval;
    IloArray<IloArray<IloNumVarArray>> * Zval;

    set<int> CPNodes;
    set<int> NonCPNodes;
    CoordSet GraphNodes;
    CoordList GraphNodesList;
    CoordSet GraphCPGlobal;
    CoordList GraphCPGlobalList;
    vector<vector<int>> ArcsMatrix; ///< A matrix of input arcs
    double minDisCP = 0.0;

    vector <vector<vector<int>>> S_PND;///PND[k][r][v]=a; pickup stop number of request k in trip r by Vehicle v is equal to a //pc(k,r,v)
    vector <vector<vector<int>>> S_NPD;///NPD[k][r][v]=a; dropoff stop number of request k in trip r by Vehicle v is equal to a //dc(k,r,v)

    vector<pair<int,vector<int>>> SameCheckpointIndicesList;
    set<int> m_VisitedCP;
    set<int> m_Sources;
    set<int> m_Destinations;
    vector<double> m_DeltaMinReq;
    vector<vector<vector<int>>> LC;

public:
    void Init(const SolverBase* rh);
    void PreProcess(SolverBase *sb);
    void HandleXZQParams(const SolverBase *sb, const map<int,int>& Z, const map<int,int>& X, const map<int,int>& Q);
    void Create_Problem();
    void AllocateMemory(const SolverBase* rh);
    void DefineConstraints(const SolverBase* rh);
    /**
     * @brief Compute the minimum amount of time needed to go from pickup to pickup point between all trips
     * @default ONLY NPND Requests
     * @param sb :: A pointer to SolverBase Class
     */
    void ComputeDeltaMinRequest(const SolverBase* sb);
    /**
     * @brief ComputeLC
     * @param sb
     */
    void ComputeLC(const SolverBase* sb);

    void Free_Problem();

    void AddDataToWarmStart(const int& ins,const SolverBase *rh);
    void doSolve(const int& ins,const SolverBase *rh);
    void PrepareCplexParams(const int& ins);
    void CplexSolve();

    void CreateArc(SolverBase *sb, const int& NodesNo, const int& TripsNo);

};

#endif // SOLVER_H
