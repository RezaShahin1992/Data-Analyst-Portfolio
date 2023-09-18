#include "Solver.h"
#include <time.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>

Solver::Solver()
{
    utility = new Utilities();
}

void Solver::Init(const SolverBase *sb)
{
    Create_Problem();
    AllocateMemory(sb);
    DefineConstraints(sb);
}

void Solver::PreProcess(SolverBase *rh)
{
    cerr << __FUNCTION__ << ": Started ... " << endl;
    double t1 = clock();
    int TripsNo = rh->TripsNo;							//etelaat marboot be safar haro migirim
    vector<int> RequestsType ;							//noe darkhast haro migirim
    int TotalReq = 0;
    CoordSet TotalCP;
    CoordSet TotalNonCP;
    for(int i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo& ti = rh->TripsInfoList[i];		//miad tripsinfo object ro ba indice i az tripsinfolist migire k too rh zakhire shode
        const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
        const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination
        for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)		//istgah ha va gheyre istgah haro miad migire
        {																							// va mirize tooye set graphnodes.
            GraphNodes.insert(*cit);
            TotalCP.insert(*cit);
            GraphCPGlobal.insert(*cit);
        }
        for(CoordSet::const_iterator cit = NonCPCoord.begin() ; cit != NonCPCoord.end() ; cit++)
        {
            GraphNodes.insert(*cit);
            TotalNonCP.insert(*cit);
        }
        TotalReq+= ti.RequestNo;
    }
    for(CoordSet::const_iterator cit = GraphNodes.begin() ; cit != GraphNodes.end() ; cit++)
    {
        GraphNodesList.push_back(*cit);
    }
    for(CoordSet::const_iterator cit = GraphCPGlobal.begin() ; cit != GraphCPGlobal.end() ; cit++)
    {
        GraphCPGlobalList.push_back(*cit);
    }
    rh->nbPoints = GraphNodes.size();				//enja mifahmim ke size set graphnodes tedade noghtehash ro mide behemoon bar oon asas

    CoordSet::const_iterator cpCit = TotalCP.begin();
    double xcp0 = cpCit->first;
    double xcp1 = (++cpCit)->first;
    minDisCP = xcp1 - xcp0;

    fstream req;
    req.open("req.txt",ios_base::out | ios_base::trunc); //enja file ro baz mikonim baraye en 4ta file
    fstream pcdc;
    pcdc.open("PC_DC.txt",ios_base::out | ios_base::trunc);
    fstream hyb;
    hyb.open("hybridtrips.txt",ios_base::out | ios_base::trunc);
    fstream graphfile;
    graphfile.open("graphfile.txt",ios_base::out | ios_base::trunc);
    int ReqID = 0;

    for(CoordSet::const_iterator cit = GraphNodes.begin() ; cit != GraphNodes.end() ; cit++)
    {
        const Coord& _c = *cit;
        if(TotalCP.find(_c) != TotalCP.end())          //file graph ro be ezaye har istgah ye setare
        {												// be ezaye har gheyre istgah ye x
            graphfile << "*" << "\t";
        }
        if(TotalNonCP.find(_c) != TotalNonCP.end())
        {
            graphfile << "x" << "\t";
        }
    }
    graphfile.close();

    S_PND.resize(TotalReq);
    S_NPD.resize(TotalReq);
    LC.resize(TotalReq);
    for(int i=0 ; i < TotalReq ; i++)
    {
        S_PND[i].resize(TripsNo);
        S_NPD[i].resize(TripsNo);
        LC[i].resize(TripsNo);
        for(int j=0 ; j < TripsNo ; j++)
        {
            S_PND[i][j].resize(rh->VehiclesNo); // enja matrix haye marboot be pc va dc por mishe
            S_NPD[i][j].resize(rh->VehiclesNo);
            LC[i][j].resize(rh->VehiclesNo);
        }
    }

    for(Index i=0 ; i < S_PND.size() ; i++)
    {
        for(Index j=0 ; j < S_PND[i].size() ; j++)
        {
            for(Index k=0 ; k < S_PND[i][j].size() ; k++)
            {
                S_PND[i][j][k] = -1;
            }
        }
    }

    for(Index i=0 ; i < S_NPD.size() ; i++)
    {
        for(Index j=0 ; j < S_NPD[i].size() ; j++)
        {
            for(Index k=0 ; k < S_NPD[i][j].size() ; k++)
            {
                S_NPD[i][j][k] = -1;
            }
        }
    }

    for(Index i=0 ; i < LC.size() ; i++)
    {
        for(Index j=0 ; j < LC[i].size() ; j++)
        {
            for(Index k=0 ; k < LC[i][j].size() ; k++)
            {
                LC[i][j][k] = -1;
            }
        }
    }

    for(int it=0 ; it < TripsNo ; it++)				//hala bia bbin chanta safar darim har safar ro bia barresi kon
    {
        TripsInfo& ti = rh->TripsInfoList[it];
//        int VID = ti.VehicleIndex;
        ti.RequestsType.resize(ti.RequestNo);			//tedade darkhast haro miaym enja resize mikonim listeshoon ro
        ti.RequestsPickStopNo.resize(ti.RequestNo);
        ti.RequestsDropStopNo.resize(ti.RequestNo);
        ti.HeuristicRequestsPickStopNo.resize(ti.RequestNo);
        ti.HeuristicRequestsDropStopNo.resize(ti.RequestNo);

        const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
        const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination
        const vector<double>& PickStopXList = ti.PickStopXList;							//bbinim k moshtari ha too che noghte e savar
        const vector<double>& PickStopYList = ti.PickStopYList;							// va piade mishan enjaha
        const vector<double>& DropStopXList = ti.DropStopXList;
        const vector<double>& DropStopYList = ti.DropStopYList;

        CoordSet LocalGraphNodes;
        CoordSet LocalGraphCPNodes;
        for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)
        {
            LocalGraphNodes.insert(*cit);
            LocalGraphCPNodes.insert(*cit);
        }
        for(CoordSet::const_iterator cit = NonCPCoord.begin() ; cit != NonCPCoord.end() ; cit++)
        {
            LocalGraphNodes.insert(*cit);
        }

        for(size_t i=0 ; i < PickStopXList.size() ; i++)			//ye halghe baraye tamame noghate moshtari ha
        {
            const double& px = PickStopXList[i];
            const double& py = PickStopYList[i];
            const double& dx = DropStopXList[i];
            const double& dy = DropStopYList[i];
            Coord coord_p = make_pair(px,py);						//joft kardan noghate savar o piade shodan moshtari ha
            Coord coord_d = make_pair(dx,dy);
            if(CPCoords.find(coord_p) != CPCoords.end() && CPCoords.find(coord_d) != CPCoords.end())			//noe moshtari ha moshakhas shode enja
            {
                ti.RequestsType[i] = static_cast<int>(PD);
            }
            else if (CPCoords.find(coord_p) != CPCoords.end() && CPCoords.find(coord_d) == CPCoords.end())
            {
                ti.RequestsType[i] = static_cast<int>(PND);
            }
            else if (CPCoords.find(coord_p) == CPCoords.end() && CPCoords.find(coord_d) != CPCoords.end())
            {
                ti.RequestsType[i] = static_cast<int>(NPD);
            }
            else if (CPCoords.find(coord_p) == CPCoords.end() && CPCoords.find(coord_d) == CPCoords.end())
            {
                ti.RequestsType[i] = static_cast<int>(NPND);
            }
            int idxP = utility->GetIndexFromSet(GraphNodes,coord_p);					//moadele andisi ro be dast avordim, yani masalan oon mokhtasate
            int idxD = utility->GetIndexFromSet(GraphNodes,coord_d);					//baraye savar shodan tooye gerafemoon moadele che andisi shode

            pair<int,int> typeIndex = make_pair(it,rh->RequestsTypeIndex.size());		//be har darkhasti ye andisi bede
            rh->RequestsTypeIndex.push_back(typeIndex);									//va oono tooye ye set darkhasta zakhire kon

            RequestsType.push_back(ti.RequestsType[i]);

            if(ti.RequestsType[i] == static_cast<int>(PND))						//hala migim age darkhast az jense PND bood,
            {
                ti.RequestsPickStopNo[i] = -1;									//savar shodanesh ro paen tar meghdar midim
                ti.RequestsDropStopNo[i] = idxD;								//piade shodanesh ro ke darim meghdar midim

                ti.HeuristicRequestsPickStopNo[i] = idxP;
                ti.HeuristicRequestsDropStopNo[i] = idxD;

                vector<int> sameCP;
                for(std::vector<int>::const_iterator cit = ti.m_PossibleTripsForRequests.begin() ; cit != ti.m_PossibleTripsForRequests.end() ; cit++)
                {
                    int LocalidxP_CP = utility->GetIndexFromSet(LocalGraphCPNodes,coord_p);
                    int CPID = 0;
                    int CP_IDX = 0;
                    for(CoordSet::const_iterator a = rh->TripsInfoList[*cit].CPCoords.begin() ;
                        a != rh->TripsInfoList[*cit].CPCoords.end() ; a++)
                    {
                        if(CPID == LocalidxP_CP)
                        {
                            CP_IDX = utility->GetIndexFromSet(GraphNodes,*a);
                            break;
                        }
                        CPID++;
                    }
                    for(int p=0 ; p < rh->VehiclesNo ; p++)
                    {
                        if(cit != ti.m_PossibleTripsForRequests.begin())
                            S_PND[ReqID][*cit][p] = CP_IDX;
                        else
                            S_PND[ReqID][*cit][p] = idxP;
                    }
                    if(CP_IDX != idxP)
                    {
                        sameCP.push_back(CP_IDX);
                    }
                }																	//dar paen miaym hamin kare bala ro baraye mahali ham anjam midim
                if(m_VisitedCP.find(idxP) == m_VisitedCP.end())
                {
                    SameCheckpointIndicesList.push_back(make_pair(idxP,sameCP));
                    m_VisitedCP.insert(idxP);
                }
            }
            else if(ti.RequestsType[i] == static_cast<int>(NPD))					//hamin karo miaym baraye npd ha ham anjam midim
            {
                int LocalidxD_CP = utility->GetIndexFromSet(LocalGraphCPNodes,coord_d);

                ti.RequestsPickStopNo[i] = idxP;
                ti.RequestsDropStopNo[i] = -1;

                ti.HeuristicRequestsPickStopNo[i] = idxP;
                ti.HeuristicRequestsDropStopNo[i] = idxD;

                vector<int> sameCP;
                for(std::vector<int>::const_iterator cit = ti.m_PossibleTripsForRequests.begin() ; cit != ti.m_PossibleTripsForRequests.end() ; cit++)
                {
                    int CPID = 0;
                    int CP_IDX = 0;
                    for(CoordSet::const_iterator a = rh->TripsInfoList[*cit].CPCoords.begin() ;
                        a != rh->TripsInfoList[*cit].CPCoords.end() ; a++)
                    {
                        if(CPID == LocalidxD_CP)
                        {
                            CP_IDX = utility->GetIndexFromSet(GraphNodes,*a);
                            break;
                        }
                        CPID++;
                    }
                    for(int p=0 ; p < rh->VehiclesNo ; p++)
                    {
                        if(cit != ti.m_PossibleTripsForRequests.begin())
                            S_NPD[ReqID][*cit][p] = CP_IDX;
                        else
                            S_NPD[ReqID][*cit][p] = idxD;
                    }
                    if(CP_IDX != idxD)
                    {
                        sameCP.push_back(CP_IDX);
                    }
                }
                if(m_VisitedCP.find(idxD) == m_VisitedCP.end())
                {
                    SameCheckpointIndicesList.push_back(make_pair(idxD,sameCP));
                    m_VisitedCP.insert(idxD);
                }
            }
            else																	//hala migim age az 2ta noe dige boodan moshtarihamoon
            {																		//varede en else besho va meghdar dehi kon savar o piade shodan ro
                ti.RequestsPickStopNo[i] = idxP;
                ti.RequestsDropStopNo[i] = idxD;
                ti.HeuristicRequestsPickStopNo[i] = idxP;
                ti.HeuristicRequestsDropStopNo[i] = idxD;
            }
            ReqID++;
        }

        for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++) // ye halghe baraye istgah ha
        {
            const Coord& cpcoord = *cit;
            int idxcp = utility->GetIndexFromSet(GraphNodes,cpcoord);   // ye tabee getindex for farakhani mikonim  ba parameter
            CPNodes.insert(idxcp);							// graphnode va mokhtasate istgah ha va gharar midim toye
            ti.Nodes.insert(idxcp);							//idxcp ke andiseshoon hast. 3 ta jaye mokhtalef en ro vared
            ti.CP_Index.insert(idxcp);						// mikone.
            if(cit == CPCoords.begin())			//agar be avaling element mokhtasate istgah ha eshare mikone
                m_Sources.insert(idxcp);		// mirizimesh tooye msource
            if(cit == --CPCoords.end())			//agar be akharin element eshare mikone,
                m_Destinations.insert(idxcp);		//mirizimesh tooye mdestination

            if(cit == CPCoords.begin())				//age be avaling element mokhtasate istgah ha eshare mikone
            {
                ti.SourcePoint = idxcp;				//mirizimesh tooye sourcepoint
            }
            if(cit == (--CPCoords.end()))			// age be akharin element eshare mikne
            {
                ti.DestinationPoint = idxcp;		// mirizimesh toye destinationpoint
            }
        }

        for(CoordSet::const_iterator cit = NonCPCoord.begin() ; cit != NonCPCoord.end() ; cit++) //halgheye gheyre istgah ha
        {
            const Coord& noncpcoord = *cit;
            int idxnoncp = utility->GetIndexFromSet(GraphNodes,noncpcoord); //andis haro migire az function .. parameter graph va gheyre istgahi
            NonCPNodes.insert(idxnoncp);			//insert mikonimesh tooye nodehaye gheyre istgahi
            ti.Nodes.insert(idxnoncp);				//tooye en
            ti.NonCP_Index.insert(idxnoncp);		// va tooye en ke badan too gheyda azashun mikhaym estefade konim
        }
    }
    req.close();

    int idxGlobalSource = *m_Sources.begin();
    int idxGlobalDestination = *--m_Destinations.end();
    m_Sources.clear();
    m_Destinations.clear();
    m_Sources.insert(idxGlobalSource);
    m_Destinations.insert(idxGlobalDestination);

    vector<int> TotalReqType;   // ye vectori miaym az darkhast ha misazim
    for(int i=0 ; i <TripsNo ; i++)				// ye halghe baraye safar ha mirim
    {
        for(int j=0 ; j < rh->TripsInfoList[i].RequestNo ; j++)   // ye halghe dige mirim baraye darkhast haye ith safar
        {
            TotalReqType.push_back(rh->TripsInfoList[i].RequestsType[j]);  // va entori miaym oon vectoremoon ro por mikonim
        }
    }
    for(Index i=0 ; i < S_PND.size() ; i++)
    {
        for(Index j=0 ; j < S_PND[i].size() ; j++)
        {
            for(Index k=0 ; k < S_PND[i][j].size() ; k++)
            {
                if(TotalReqType[i] == PND)			// enja migirm migim age oon darkhast az jense pnd bod hala khorooji bede
                {
                    pcdc << "PND["<<i<<"]["<<j<<"]["<<k<<"] = " << S_PND[i][j][k] << endl;
                }
            }
        }
    }

    for(Index i=0 ; i < S_NPD.size() ; i++)
    {
        for(Index j=0 ; j < S_NPD[i].size() ; j++)
        {
            for(Index k=0 ; k < S_NPD[i][j].size() ; k++)
            {
                if(TotalReqType[i] == NPD)					// hamin karo baraye npd ha mikonim ta matris haro namayesh bedim
                {
                    pcdc << "NPD["<<i<<"]["<<j<<"]["<<k<<"] = " << S_NPD[i][j][k] << endl;
                }
            }
        }
    }

    for(int i=0 ; i < TripsNo ; i++)  // enja etelaat marboot be enke oon HYBR (k) chi hast ro namayesh midim
    {
        TripsInfo& ti = rh->TripsInfoList[i];
        hyb << "Trips["<<ti.id << "] :: ";
        for(std::vector<int>::const_iterator cit = ti.m_PossibleTripsForRequests.begin() ; cit != ti.m_PossibleTripsForRequests.end() ; cit++)
        {
            hyb << *cit << "," ;
        }
        hyb << endl;
    }

    pcdc.close();
    hyb.close();

    for(int i=0 ; i < TripsNo ; i++)												     	//rajebe q kar mikonm enja
    {
        TripsInfo& ti = rh->TripsInfoList[i];
        for(set<int>::const_iterator cit = ti.Nodes.begin() ; cit != ti.Nodes.end() ;cit++)	//yani masalan safare aval, q esh ro be tedade node haye
        {																					//oon safar resize mikonim
            ti.NodesList.push_back(*cit);
        }
        ti.q.resize(ti.NodesList.size());
    }

    rh->TC0 = ((int)CPNodes.size() - 1) * TripsNo + 1;
    rh->TC = rh->TC0 * rh->VehiclesNo;
    ////////////////////////////////// Arcs & Path & Delta ///////////////////////////////

    int NodesNo = static_cast<int>(GraphNodes.size());			//tedade node haro mishmorim va ekhtesas midim be nodesno
    vector<Coord> GraphList;									//ye vector dorost mikonim
    for(CoordSet::const_iterator cit = GraphNodes.begin() ; cit != GraphNodes.end() ; cit++) //nodehaye graph ro miaym mirizim too oon vector
    {
        GraphList.push_back(*cit);
    }
    rh->delta.resize(NodesNo);								//enja miaym mohasebate marboot be delta ro namayesh midim
    for (int i = 0; i < NodesNo; i++)
    {
        rh->delta[i].resize(NodesNo);
    }
    for (int i = 0; i < NodesNo; i++) {
        for (int j = 0; j < NodesNo; j++) {
            if(i == j)
            {
                rh->delta[i][j] = M;									//bara mohasebe delta age nodehaye barabar bashe m mizarim
            }
            else
            {
                if(rh->delta[j][i] != 0.0)
                    rh->delta[i][j] = rh->delta[j][i]/*M*/;				//miaym bala mosalasi haro too matris delta m mizarim
                else
                    rh->delta[i][j] = utility->CalculateDistance(GraphList[i],GraphList[j]);				//baraye baghiye jaha zaman residan behesh ro mohasebe mikonim

            }
        }
    }

    CreateArc(rh,NodesNo,TripsNo);


    fstream arc;
    arc.open("arc.txt",ios_base::out | ios_base::trunc);			//enja arc haro namayesh midim
    arc.clear();
    for(size_t i=0 ; i < ArcsMatrix.size() ; i++)
    {
        for(size_t j=0 ; j < ArcsMatrix.size() ; j++)
        {
            arc << ArcsMatrix[i][j] << "\t";
        }
        arc << "\n";
    }
    arc << "\n";
    arc.close();


    fstream deltafile;
    deltafile.open("delta.txt",ios_base::out | ios_base::trunc);

    vector <int> RDS;
    vector <int> RPS;
    for(int i=0 ; i < TripsNo ; i++)
    {
        for(Index j=0 ; j < rh->TripsInfoList[i].RequestsDropStopNo.size() ; j++)
            RDS.push_back(rh->TripsInfoList[i].RequestsDropStopNo[j]);
        for(Index j=0 ; j < rh->TripsInfoList[i].RequestsPickStopNo.size() ; j++)
            RPS.push_back(rh->TripsInfoList[i].RequestsPickStopNo[j]);
    }

    for(Index i=0 ; i < S_PND.size() ; i++)
    {
        int dIdx = RDS[i];
        int detectedTrip = utility->GetTripNumberOfRequest(rh,i);
        if(rh->TripsInfoList[detectedTrip].m_PossibleTripsForRequests.empty())
            continue;
        int firstTripOfHybrid = rh->TripsInfoList[detectedTrip].m_PossibleTripsForRequests[0];
        for(Index j=0 ; j < S_PND[i].size() ; j++)
        {
            for(Index k=0 ; k < S_PND[i][j].size() ; k++)
            {
                if(TotalReqType[i] == PND)			// enja migirm migim age oon darkhast az jense pnd bod hala khorooji bede
                {
                    int pndIdx = S_PND[i][j][k];
                    if(pndIdx == -1)
                        continue;
                    if(ArcsMatrix[pndIdx][dIdx] == 1)
                    {
                        int originalIdx = S_PND[i][firstTripOfHybrid][k];
                        const Coord& original_Coord = utility->GetCoordinationFromIndex(rh,originalIdx);
                        const Coord& drop_Coord     = utility->GetCoordinationFromIndex(rh,dIdx);
                        const Coord& pnd_Coord      = utility->GetCoordinationFromIndex(rh,pndIdx);
                        double xdiff = std::fabs(original_Coord.first - drop_Coord.first);
                        double realX = pnd_Coord.first + xdiff;
                        Coord realCoord;
                        realCoord.first  = realX;
                        realCoord.second = drop_Coord.second;

//                        rh->delta[pndIdx][dIdx] = rh->delta[originalIdx][dIdx];
                        for(int q = pndIdx ; q < NodesNo ; q++)
                        {
                            const Coord& q_Coord = utility->GetCoordinationFromIndex(rh,q);
                            if(TotalCP.find(q_Coord) != TotalCP.end())
                            {
                                if(q_Coord.first < realCoord.first)
                                {
                                    rh->delta[q][dIdx] = utility->CalculateDistance(q_Coord,realCoord);
                                }
                                else if(q_Coord.first > realCoord.first)
                                {
                                    rh->delta[dIdx][q] = utility->CalculateDistance(realCoord,q_Coord);
                                    break;
                                }
                            }
                            else
                            {
                                if(q_Coord.first > realCoord.first)
                                {
                                    rh->delta[dIdx][q] = utility->CalculateDistance(realCoord,q_Coord);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for(Index i=0 ; i < S_NPD.size() ; i++)
    {
        int pIdx = RPS[i];
        int detectedTrip = utility->GetTripNumberOfRequest(rh,i);
        if(rh->TripsInfoList[detectedTrip].m_PossibleTripsForRequests.empty())
            continue;
        int firstTripOfHybrid = rh->TripsInfoList[detectedTrip].m_PossibleTripsForRequests[0];
        for(Index j=0 ; j < S_NPD[i].size() ; j++)
        {
            for(Index k=0 ; k < S_NPD[i][j].size() ; k++)
            {
                if(TotalReqType[i] == NPD)			// enja migirm migim age oon darkhast az jense pnd bod hala khorooji bede
                {
                    int npdIdx = S_NPD[i][j][k];
                    if(npdIdx == -1)
                        continue;
                    if(ArcsMatrix[pIdx][npdIdx] == 1)
                    {
                        int originalIdx = S_NPD[i][firstTripOfHybrid][k];
                        const Coord& original_Coord = utility->GetCoordinationFromIndex(rh,originalIdx);
                        const Coord& pick_Coord     = utility->GetCoordinationFromIndex(rh,pIdx);
                        const Coord& npd_Coord      = utility->GetCoordinationFromIndex(rh,npdIdx);
                        double xdiff = std::fabs(original_Coord.first - pick_Coord.first);
                        double realX = npd_Coord.first - xdiff;
                        Coord realCoord;
                        realCoord.first  = realX;
                        realCoord.second = pick_Coord.second;

                        rh->delta[pIdx][npdIdx] = utility->CalculateDistance(realCoord,npd_Coord);
                    }
                }
            }
        }
    }

    for(int i=0 ; i < rh->delta.size() ; i++)
    {
        for(int j=0 ; j < rh->delta[i].size() ; j++)
        {
            stringstream tmp;
            tmp << setprecision(3) << fixed << rh->delta[i][j];
            double new_val = stod(tmp.str());
            deltafile << new_val << "\t";
        }
        deltafile << "\n";
    }
    deltafile.close();

    utility->ComputeTripTime(rh);
    ComputeDeltaMinRequest(rh);
    ComputeLC(rh);

    double t2 =clock();
    double diff = (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
    cerr << __FUNCTION__ << ": Finished ... dt = " << diff << endl;
}

void Solver::Create_Problem()
{
    lp = new IloModel(env);
    if (lp == NULL)
    {
        cout << "Cannot create problem!!" << endl;
        exit(502);
    }
}

void Solver::AllocateMemory(const SolverBase *sb)
{
    //////////Number of Decision Variables

    double t1 = clock();
    int NodesNo = 0;
    int VehiclesNo = sb->VehiclesNo;
    int RequestNo = 0;
    double VehicleCap = sb->Input_VehiclesCapacity;
    int TripsNo = sb->TripsNo;
    vector <double> Weights;
    vector <double> ReqTau;
    for(int i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo &ti = sb->TripsInfoList[i];
        cerr << "VehicleIndex : " << sb->TripsInfoList[i].VehicleIndex << endl;
        RequestNo += ti.RequestNo;
        Weights = ti.W;
        for(int j=0 ; j < ti.RequestsTau.size() ; j++)
            ReqTau.push_back(ti.RequestsTau[j]);
    }

    NodesNo = GraphNodes.size();

    Index TotalReq = 0;
    vector<int> TotalReqType;

    for(int i=0 ; i <TripsNo ; i++)
    {
        TotalReq += sb->TripsInfoList[i].RequestNo;
        for(int j=0 ; j < sb->TripsInfoList[i].RequestNo ; j++)
        {
            TotalReqType.push_back(sb->TripsInfoList[i].RequestsType[j]);
        }
    }

    //===========================================================================================
    Xval = new IloArray<IloArray<IloNumVarArray>>(env, NodesNo);
    for (int i = 0; i < NodesNo; i++)
    {
        (*Xval)[i] = IloArray<IloNumVarArray>(env, NodesNo);
        for (int j = 0; j < NodesNo; j++)
        {
            (*Xval)[i][j] = IloNumVarArray(env, VehiclesNo);
            for (int v = 0; v < VehiclesNo; v++)
            {
                int detectedTrip = 0;
                for(int q=0 ; q < TripsNo ; q++)
                {
                    if(sb->TripsInfoList[q].Nodes.find(i) != sb->TripsInfoList[q].Nodes.end()
                            && sb->TripsInfoList[q].Nodes.find(j) != sb->TripsInfoList[q].Nodes.end())
                    {
                        detectedTrip = q;
                        break;
                    }
                }
                string valname = "X("+to_string(i)+")("+to_string(j)+")("+to_string(v)+")";
                if(!IsSecondRun)
                    (*Xval)[i][j][v] = IloNumVar(env, 0, 1, IloNumVar::Float,valname.c_str());
                else
                    (*Xval)[i][j][v] = IloNumVar(env, 0, 1, IloNumVar::Bool,valname.c_str());

            }
        }
    }

    //////////////
    Tval = new IloNumVarArray(env, NodesNo);
    for(int i=0 ;  i < NodesNo ; i++)
    {
        string valname = "T("+to_string(i)+")";
        (*Tval)[i] = IloNumVar(env, 0, M , IloNumVar::Float,valname.c_str());
    }
    //////////////
    Tbarval = new IloNumVarArray(env, NodesNo);
    for(int i=0 ;  i < NodesNo ; i++)
    {
        string valname = "Tbar("+to_string(i)+")";
        (*Tbarval)[i] = IloNumVar(env, 0, M , IloNumVar::Float,valname.c_str());
    }
    /////////////////
    Pval = new IloNumVarArray(env, RequestNo);
    for(int k=0 ;  k < RequestNo ; k++)
    {
        string valname = "P("+to_string(k)+")";
        (*Pval)[k] = IloNumVar(env, 0, M , IloNumVar::Float,valname.c_str());
    }

    //////////////////
    Dval = new IloNumVarArray(env, RequestNo);
    for(int k=0 ;  k < RequestNo ; k++)
    {
        string valname = "D("+to_string(k)+")";
        (*Dval)[k] = IloNumVar(env, 0, M , IloNumVar::Float,valname.c_str());
    }
    //////////////////
    Zval = new IloArray<IloArray<IloNumVarArray>>(env, RequestNo);
    for (int k = 0; k < RequestNo; k++)
    {
        (*Zval)[k] = IloArray<IloNumVarArray>(env, TripsNo);
        for (int r = 0; r < TripsNo; r++)
        {
            (*Zval)[k][r] = IloNumVarArray(env, VehiclesNo);
            for (int v = 0; v < VehiclesNo; v++)
            {
                string valname = "Z("+to_string(k)+")("+to_string(r)+")("+to_string(v)+")";
                (*Zval)[k][r][v] = IloNumVar(env, 0, 1, IloNumVar::Bool,valname.c_str());
            }
        }
    }

    /////////////////
    Qval = new IloArray<IloNumVarArray>(env,NodesNo);
    for(int i=0 ; i < NodesNo ; i++)
    {
        (*Qval)[i] = IloNumVarArray(env,VehiclesNo);
        for(int v=0 ; v < VehiclesNo ; v++)
        {
            string valname = "Q("+to_string(i)+")("+to_string(v)+")";
            (*Qval)[i][v] = IloNumVar(env, 0.0, M , IloNumVar::Float,valname.c_str());
        }
    }

    /////////////////
    qval = new IloNumVarArray(env, NodesNo);
    for(int i=0 ;  i < NodesNo ; i++)
    {
        string valname = "q("+to_string(i)+")";
        (*qval)[i] = IloNumVar(env, -M, M , IloNumVar::Float,valname.c_str());
    }

    /////////////////


    // Total shuttle travel time
    IloExpr objective(env);
    for(int i=0; i<NodesNo; i++)
    {
        for(int j=0; j<NodesNo; j++)
        {
            for(int v=0 ; v < VehiclesNo ; v++)
            {
                if(ArcsMatrix[i][j] == 1)
                    objective+=(Weights[0]*sb->delta[i][j]*(*Xval)[i][j][v]);
            }
        }
    }

    for(int k=0; k<RequestNo; k++)
        objective.operator+=(Weights[1]*((*Dval)[k] - (*Pval)[k]));

    for(int k=0; k<RequestNo; k++)
        objective.operator+=(Weights[2]*((*Pval)[k] - ReqTau[k]));




    IloObjective obj(env, objective, IloObjective::Minimize);

    lp->add(obj);

    objective.end();

    double t2 =clock();
    double diff = (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
}

void Solver::DefineConstraints(const SolverBase *rh)
{
    const vector<vector<double>>& delta = rh->delta;
    const int& VehiclesNo = rh->VehiclesNo;
    int TripsNo = rh->TripsNo;

    Index TotalReq = 0;
    vector<int> TotalReqType;
    vector<int> Nodes;
    vector <int> RPS;
    vector <int> RDS;
    vector <double> StopsB = rh->StopsB;
    vector <double> RequestsTau;
    vector <int> StopsType;
    vector <double> StopsTeta = rh->StopsTeta;
    Index nbNodes = rh->nbPoints;
    Index VehicleIndex = 0;

    fstream ConsoleOutput;
    ConsoleOutput.open("OUTPUT/Console"+to_string(cntInputFile)+".txt",ios_base::out | ios_base::trunc);
    for(int tno = 0 ; tno < rh->TripsNo ; tno++)
    {
        ConsoleOutput << "Trip["<<rh->TripsInfoList[tno].id<<"] . Start Time : " << utility->ConvertToRealTime(rh->TripsInfoList[tno].RealStartTime) << std::endl;
        ConsoleOutput << "Trip["<<rh->TripsInfoList[tno].id<<"] . Arrival Time : " << utility->ConvertToRealTime(rh->TripsInfoList[tno].RealArrivalTime) << std::endl;
    }

    for(int i=0 ; i <TripsNo ; i++)
    {
        TotalReq += rh->TripsInfoList[i].RequestNo;
        for(int j=0 ; j < rh->TripsInfoList[i].RequestNo ; j++)
        {
            TotalReqType.push_back(rh->TripsInfoList[i].RequestsType[j]);
        }
    }
    for(int i=0 ; i < TripsNo ; i++)
    {
        for(Index j=0 ; j < rh->TripsInfoList[i].NodesList.size() ; j++)
            Nodes.push_back(rh->TripsInfoList[i].NodesList[j]);

        for(Index j=0 ; j < rh->TripsInfoList[i].RequestsPickStopNo.size() ; j++)
            RPS.push_back(rh->TripsInfoList[i].RequestsPickStopNo[j]);

        for(Index j=0 ; j < rh->TripsInfoList[i].RequestsDropStopNo.size() ; j++)
            RDS.push_back(rh->TripsInfoList[i].RequestsDropStopNo[j]);

        for(Index j=0 ; j < rh->TripsInfoList[i].RequestsTau.size() ; j++)
            RequestsTau.push_back(rh->TripsInfoList[i].RequestsTau[j]);

        for(Index j=0 ; j < rh->TripsInfoList[i].StopsType.size() ; j++)
            StopsType.push_back(rh->TripsInfoList[i].StopsType[j]);
    }


    /**
     *  Capacity Constraints
     *  **/

    /**
      * @brief :: Constraint 44
      * **/
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::NPD || TotalReqType[k] == REQUESTTYPE::PND)
        {
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
//            for(int v=0 ; v < VehiclesNo ; v++)
            {
                for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                     r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
                {
                    int sourceIdx           = rh->TripsInfoList[DetectedTrip].SourcePoint;
                    int destinationIdx      = rh->TripsInfoList[DetectedTrip].DestinationPoint;
                    double ThetaSource      = StopsTeta[sourceIdx];
                    double ThetaDestination = StopsTeta[destinationIdx];
                    double diffTheta        = ThetaDestination - ThetaSource;
                    int v = rh->TripsInfoList[*r].VehicleIndex;
                    lp->add((*Dval)[k] - (*Pval)[k] <= diffTheta + M*(1-(*Zval)[k][*r][v]));
                }
            }
        }
    }
    /**
      * @brief :: Constraint 46
      * **/
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::NPD || TotalReqType[k] == REQUESTTYPE::NPND)
        {
            int nbReq = utility->GetNbReqFromRPS(rh,RPS[k],static_cast<REQUESTTYPE>(TotalReqType[k]));
            lp->add((*qval)[RPS[k]] == nbReq);
        }
    }
    /**
      * @brief :: Constraint 47
      * **/
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::PND || TotalReqType[k] == REQUESTTYPE::NPND)
        {
            int nbReq = utility->GetNbReqFromRDS(rh,RDS[k],static_cast<REQUESTTYPE>(TotalReqType[k]));
            lp->add((*qval)[RDS[k]] == -nbReq);
        }
    }
    /**
      * @brief :: Constraint 45
      * **/
    for(Index i=0 ; i< nbNodes ; i++)
    {
        int idx_i = Nodes[i];
        if (StopsType[i] == 1) // As a checkpoints node
        {
            IloExpr sumPND(env);
            IloExpr sumNPD(env);
            int nb_pd_ps = 0;
            int nb_pd_ds = 0;
            for (Index k = 0; k < TotalReq; k++)
            {
                if (TotalReqType[k] == REQUESTTYPE::PD)
                {
                    if(RPS[k] == idx_i)
                    {
                        nb_pd_ps += utility->GetNbReqFromRPS(rh,RPS[k],static_cast<REQUESTTYPE>(TotalReqType[k]));
                    }
                }

                if (TotalReqType[k] == REQUESTTYPE::PD)
                {
                    if(RDS[k] == idx_i)
                    {
                        nb_pd_ds += utility->GetNbReqFromRDS(rh,RDS[k],static_cast<REQUESTTYPE>(TotalReqType[k]));
                    }
                }

                for (int r=0 ; r < TripsNo ; r++)
                {
                    if (TotalReqType[k] == REQUESTTYPE::PND)
                    {
                        int nbReq_PND = utility->GetNbReqFromPC(S_PND,rh,idx_i);
//                        for (int v = 0; v < VehiclesNo; v++)
                        {
                            int v = rh->TripsInfoList[r].VehicleIndex;
                            if(S_PND[k][r][v] == idx_i)
                            {
                                sumPND+=((*Zval)[k][r][v] * nbReq_PND);
                            }
                        }
                    }

                    if (TotalReqType[k] == REQUESTTYPE::NPD)
                    {
                        int nbReq_NPD = utility->GetNbReqFromDC(S_NPD,rh,idx_i);
//                        for (int v = 0; v < VehiclesNo; v++)
                        {
                            int v = rh->TripsInfoList[r].VehicleIndex;
                            if(S_NPD[k][r][v] == idx_i)
                            {
                                sumNPD+=((*Zval)[k][r][v] * nbReq_NPD);
                            }
                        }
                    }
                }
            }
            lp->add((*qval)[idx_i] == (nb_pd_ps + sumPND - sumNPD - nb_pd_ds));
        }
    }
    /**
      * @brief :: Constraint 43
      * **/
//    for (int v = 0; v < VehiclesNo; v++)
    {
        int v = utility->GetVehicleNumberOfCurrentNode(rh,0);
        lp->add((*Qval)[0][v] >= (*qval)[0]);
    }
    for (Index i = 0; i < nbNodes; i++)
    {
        int idx_i = Nodes[i];
        for (Index j = 0; j < nbNodes; j++)
        {
            int idx_j = Nodes[j];
//            if(ArcsMatrix[idx_i][idx_j] == 0)
//            {
//                continue;
//            }
//            for (int r=0 ; r < TripsNo ; r++)
            {
//                for (int v = 0; v < VehiclesNo; v++)
                {
                    int v = utility->GetVehicleNumberOfCurrentNode(rh,idx_i);
                    if(v >= rh->VehiclesNo)
                        continue;
                    lp->add((*Qval)[idx_j][v] >= ((*Qval)[idx_i][v] + (*qval)[idx_j]) - M*(1-(*Xval)[idx_i][idx_j][v]));
                }
            }
        }
    }

    /**
     * End Capacity Constraints
     * **/

    /////////////////////////////////C2//////////////////////
    for (Index j = 0; j < nbNodes; j++)
    {
        int idx_j = Nodes[j];
        IloExpr RejExp(env);
        vector<Index> ReqInNode = utility->GetRequestsInNode(S_PND,rh,j);

        if(m_Sources.find(idx_j) == m_Sources.end())
        {
            IloExpr sumX(env);
            bool isAssign = false;
//            for (int v = 0; v < VehiclesNo; v++)
            {
                for (Index i = 0; i < nbNodes; i++)
                {
                    int idx_i = Nodes[i];
                    if(ArcsMatrix[i][j] == 1)
                    {
                        int vi = utility->GetVehicleNumberOfCurrentNode(rh,i);
                        int vj = utility->GetVehicleNumberOfCurrentNode(rh,j);
                        int v = 0;
                        if(vi != vj)
                        {
                            if(i > j)
                                v = vi;
                            else
                                v = vj;
                        }
                        else
                        {
                            v = vi;
                        }
                        sumX+=((*Xval)[idx_i][idx_j][v]);
                        isAssign = true;
                    }
                }
            }
            if(isAssign)
                lp->add(sumX == 1 );
        }
    }
    /////////////////////////////////C3//////////////////////
    for (Index i = 0; i < nbNodes; i++)
    {
        int idx_i = Nodes[i];
        IloExpr RejExp(env);
        vector<Index> ReqInNode = utility->GetRequestsInNode(S_PND,rh,i);
        if (m_Destinations.find(idx_i) == m_Destinations.end())
        {
            IloExpr sumX(env);
            bool isAssign = false;
//            for (int v = 0; v < VehiclesNo; v++)
            {
                for (Index j = 0; j < nbNodes; j++)
                {
                    int idx_j = Nodes[j];
                    if(ArcsMatrix[i][j] == 1)
                    {
                        int vi = utility->GetVehicleNumberOfCurrentNode(rh,i);
                        int vj = utility->GetVehicleNumberOfCurrentNode(rh,j);
                        int v = 0;
                        if(vi != vj)
                        {
                            if(i > j)
                                v = vi;
                            else
                                v = vj;
                        }
                        else
                        {
                            v = vi;
                        }
                        sumX+=((*Xval)[idx_i][idx_j][v]);
                        isAssign = true;
                    }
                }
            }
            if(isAssign)
                lp->add(sumX == 1);
        }
    }
    /////////////////////////////////C4//////////////////////
    for (Index j = 0; j < nbNodes; j++)
    {
        int idx_j = Nodes[j];
        if(m_Sources.find(idx_j) == m_Sources.end() && m_Destinations.find(idx_j) == m_Destinations.end())
        {
            IloExpr sumX(env);
            bool isAssign = false;
//            for (int v = 0; v < VehiclesNo; v++)
            {
                for (Index i = 0; i < nbNodes; i++)
                {
                    int idx_i = Nodes[i];
                    if(ArcsMatrix[i][j] == 1)
                    {
                        int vi = utility->GetVehicleNumberOfCurrentNode(rh,i);
                        int vj = utility->GetVehicleNumberOfCurrentNode(rh,j);
                        int v = 0;
                        if(vi != vj)
                        {
                            if(i > j)
                                v = vi;
                            else
                                v = vj;
                        }
                        else
                        {
                            v = vi;
                        }
                        sumX+=((*Xval)[idx_i][idx_j][v]);
                    }
                }
            }
//            for (int v = 0; v < VehiclesNo; v++)
            {
                for (Index i = 0; i < nbNodes; i++)
                {
                    int idx_i = Nodes[i];
                    if(ArcsMatrix[j][i] == 1)
                    {
                        int vi = utility->GetVehicleNumberOfCurrentNode(rh,i);
                        int vj = utility->GetVehicleNumberOfCurrentNode(rh,j);
                        int v = 0;
                        if(vi != vj)
                        {
                            if(i > j)
                                v = vi;
                            else
                                v = vj;
                        }
                        else
                        {
                            v = vi;
                        }
                        sumX-=((*Xval)[idx_j][idx_i][v]);
                    }
                }
            }
            lp->add(sumX == 0);
        }
    }
    /////////////////////////////////C5//////////////////////
    for(Index i=0 ; i< nbNodes ; i++)
    {
        int idx_i = Nodes[i];
        if (StopsType[i] == 1)
        {
            lp->add((*Tval)[idx_i] == StopsTeta[i]);
        }
    }
    /////////////////////////////////C6//////////////////////
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] != 2)
        {
            lp->add((*Pval)[k] == (*Tval)[RPS[k]]);
        }
    }
    /////////////////////////////////C7//////////////////////
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] != REQUESTTYPE::NPD)
        {
            lp->add((*Dval)[k] == (*Tbarval)[RDS[k]]);
        }
    }
    /////////////////////////////////C8//////////////////////
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::NPD || TotalReqType[k] == REQUESTTYPE::PND)
        {
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            IloExpr sumZ(env);
//            for (int v = 0;v < VehiclesNo; v++)
            {
                for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                     r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
                {
                    int v = rh->TripsInfoList[*r].VehicleIndex;
                    sumZ+=((*Zval)[k][*r][v]);
                }
            }
            lp->add(sumZ == 1 );
        }
    }
//    lp->add((*Zval)[0][0][0] == 1);
//    lp->add((*Zval)[1][1][1] == 1);
//    lp->add((*Zval)[0][1][0] == 1);
//    lp->add((*Zval)[1][0][0] == 0);
    /////////////////////////////////C9//////////////////////
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::PND) {
            for (int r=0 ; r < TripsNo ; r++)
            {
//                for (int v = 0; v < VehiclesNo; v++)
                {
                    int v = rh->TripsInfoList[r].VehicleIndex;
                    if(S_PND[k][r][v] == -1)
                        continue;
                    lp->add((*Pval)[k] >= ((*Tval)[S_PND[k][r][v]]-M*(1-((*Zval)[k][r][v]))));
                }
            }
        }
    }
    /////////////////////////////////C10//////////////////////
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::PND) {
            for (int r=0 ; r < TripsNo ; r++)
            {
//                for (int v = 0; v < VehiclesNo; v++)
                {
                    int v = rh->TripsInfoList[r].VehicleIndex;
                    if(S_PND[k][r][v] == -1)
                        continue;
                    lp->add((*Pval)[k]<= ((*Tval)[S_PND[k][r][v]]+M*(1-((*Zval)[k][r][v]))));
                }
            }
        }
    }
    /////////////////////////////////C11//////////////////////
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::NPD) {
            for (int r=0 ; r < TripsNo ; r++)
            {
//                for (int v = 0; v < VehiclesNo; v++)
                {
                    int v = rh->TripsInfoList[r].VehicleIndex;
                    if(S_NPD[k][r][v] == -1)
                        continue;
                    lp->add((*Dval)[k] >= ((*Tbarval)[S_NPD[k][r][v]]-M*(1-((*Zval)[k][r][v]))));
                }
            }
        }
    }
    /////////////////////////////////C12//////////////////////
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::NPD) {
            for (int r=0 ; r < TripsNo ; r++)
            {
//                for (int v = 0; v < VehiclesNo; v++)
                {
                    int v = rh->TripsInfoList[r].VehicleIndex;
                    if(S_NPD[k][r][v] == -1)
                        continue;
                    lp->add((*Dval)[k] <= ((*Tbarval)[S_NPD[k][r][v]]+M*(1-((*Zval)[k][r][v]))));
                }
            }
        }
    }
    /////////////////////////////////C13//////////////////////
    for (Index k = 0; k < TotalReq; k++)
    {
        lp->add(((*Pval)[k] >= RequestsTau[k]));
    }
    /////////////////////////////////C14//////////////////////
    for (Index k = 0; k < TotalReq; k++)
    {
        lp->add((*Dval)[k] >= (*Pval)[k]);
    }
    /////////////////////////////////C15//////////////////////
    for (Index i = 0; i < nbNodes; i++)
    {
        int idx_i = Nodes[i];
        for (Index j = 0; j < nbNodes; j++)
        {
            int idx_j = Nodes[j];
            if(ArcsMatrix[idx_i][idx_j] == 0)
            {
                continue;
            }
            IloExpr sumXdelta(env), sumX(env), sumRHS(env);
//            for (int v = 0; v < VehiclesNo; v++)
            {
//                int v = GetVehicleNumberOfCurrentNode(rh,i);
                int vi = utility->GetVehicleNumberOfCurrentNode(rh,i);
                int vj = utility->GetVehicleNumberOfCurrentNode(rh,j);
                int v = 0;
                if(vi != vj)
                {
                    if(i > j)
                        v = vi;
                    else
                        v = vj;
                }
                else
                {
                    v = vi;
                }
                sumXdelta+=((delta[idx_i][idx_j]));
                sumX+=((*Xval)[idx_i][idx_j][v]);


            }
            lp->add((*Tbarval)[idx_j] >= (*Tval)[idx_i]-M+(M+sumXdelta)*sumX);
  //          lp->add(-M*-sumX);
        }
    }
    /////////////////////////////////C16//////////////////////
    for (Index i = 0; i < nbNodes; i++)
    {
        int idx_i = Nodes[i];
        if(m_Sources.find(idx_i) == m_Sources.end() && m_Destinations.find(idx_i) == m_Destinations.end())
        {
            lp->add((*Tval)[idx_i] >= (*Tbarval)[idx_i]+StopsB[i]);
        }
    }
    /////////////////////////////////C17//////////////////////
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k]==REQUESTTYPE::PD || TotalReqType[k]==REQUESTTYPE::NPND){
//            for (int v = 0; v < VehiclesNo; v++)
            {
                IloExpr sum(env);
                bool isAssign = false;
                for (Index j = 0; j < nbNodes; j++)
                {
                    int idx_j = Nodes[j];

                    int vi = utility->GetVehicleNumberOfCurrentNode(rh,RPS[k]);
                    int vj = utility->GetVehicleNumberOfCurrentNode(rh,j);
                    int v = 0;
                    if(vi != vj)
                    {
                        if(RPS[k] > j)
                            v = vi;
                        else
                            v = vj;
                    }
                    else
                    {
                        v = vi;
                    }
                    sum+=((*Xval)[RPS[k]][idx_j][v]);

                    vi = utility->GetVehicleNumberOfCurrentNode(rh,j);
                    vj = utility->GetVehicleNumberOfCurrentNode(rh,RDS[k]);

                    if(vi != vj)
                    {
                        if(RDS[k] > j)
                            v = vi;
                        else
                            v = vj;
                    }
                    else
                    {
                        v = vi;
                    }
                    sum-=((*Xval)[idx_j][RDS[k]][v]);
                    isAssign = true;

                }
                if(isAssign)
                    lp->add(sum == 0);
            }
        }
    }
    ///////////////////////////////C18//////////////////////
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::PND)
        {
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            int Drop = RDS[k];
//            for (int v = 0; v < VehiclesNo; v++)
            {
                IloExpr LHS(env);

                bool rhsAssign = false;
                bool lhsAssign = false;
                for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin(); r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
                {
                    for (Index j = 0; j < nbNodes; j++)
                    {
                        int idx_j = Nodes[j];
                        int v = rh->TripsInfoList[*r].VehicleIndex;
                        int li = S_PND[k][*r][v];
                        if(li == -1)
                            continue;
                        int vi = utility->GetVehicleNumberOfCurrentNode(rh,li);
                        int vj = utility->GetVehicleNumberOfCurrentNode(rh,j);
                        v = 0;
                        if(vi != vj)
                        {
                            if(li > j)
                                v = vi;
                            else
                                v = vj;
                        }
                        else
                        {
                            v = vi;
                        }
                        LHS+=((*Xval)[li][j][v]);
                    }
                }
                IloExpr RHS(env);
                for (int j = 0; j < nbNodes; j++)
                {
                    int vi = utility->GetVehicleNumberOfCurrentNode(rh,j);
                    int vj = utility->GetVehicleNumberOfCurrentNode(rh,Drop);
                    int v = 0;
                    if(vi != vj)
                    {
                        if(j > Drop)
                            v = vi;
                        else
                            v = vj;
                    }
                    else
                    {
                        v = vi;
                    }
                    RHS+=((*Xval)[j][Drop][v]);
                }
                lp->add((LHS == RHS));
            }
        }
    }
    /////////////////////////////////C19//////////////////////
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::NPD)
        {
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            for (int v = 0; v < VehiclesNo; v++)
            {
                IloExpr sumXRHS(env);
                IloExpr sumXLHS(env);
                bool rhsAssign = false;
                bool lhsAssign = false;
                for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                     r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
                {
                    for (Index j = 0; j < nbNodes; j++)
                    {
                        int idx_j = Nodes[j];
                        int v = rh->TripsInfoList[*r].VehicleIndex;
                        if(S_NPD[k][*r][v] == -1)
                            continue;

                        int vi = utility->GetVehicleNumberOfCurrentNode(rh,j);
                        int vj = utility->GetVehicleNumberOfCurrentNode(rh,S_NPD[k][*r][v]);
                        v = 0;
                        if(vi != vj)
                        {
                            if(j > S_NPD[k][*r][v])
                                v = vi;
                            else
                                v = vj;
                        }
                        else
                        {
                            v = vi;
                        }
                        sumXRHS+=((*Xval)[idx_j][S_NPD[k][*r][v]][v]);
                    }
                }
                for (Index j = 0; j < nbNodes; j++)
                {
                    int idx_j = Nodes[j];

                    int vi = utility->GetVehicleNumberOfCurrentNode(rh,RPS[k]);
                    int vj = utility->GetVehicleNumberOfCurrentNode(rh,j);
                    int v = 0;
                    if(vi != vj)
                    {
                        if(RPS[k] > j)
                            v = vi;
                        else
                            v = vj;
                    }
                    else
                    {
                        v = vi;
                    }
                    sumXLHS+=((*Xval)[RPS[k]][idx_j][v]);
                }
                lp->add(sumXLHS == sumXRHS);
            }
        }
    }

#ifdef GROUP_1
        ////////////////////////////// C-22 (Constraint-16 In Second Paper)
    for(Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::PND) // PND
        {
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                 r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
            {
                if(*r != rh->TripsNo-1)
                {
                    int m_NextTripIdx = ((*r) + 1);
                    int v = rh->TripsInfoList[m_NextTripIdx].VehicleIndex; // Vehicle Index For Next Trip
                    int vc = rh->TripsInfoList[*r].VehicleIndex; //  Vehicle Index For Current Trip
                    int pc = S_PND[k][m_NextTripIdx][v];
                    if (pc != -1)
                    {
                        lp->add((*Tval)[RDS[k]] <= (*Zval)[k][*r][vc] * StopsTeta[pc] + M*(1-(*Zval)[k][*r][vc]));
                    }
                }
            }
        }
    }
    //////////////////////////// C-23 (Constraint-17 In Second Paper)
    for(Index k=0 ; k < TotalReq ; k++)
    {
        if(TotalReqType[k] == REQUESTTYPE::NPD)
        {
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                 r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
            {
                if(*r != 0)
                {
                    int m_PreviousTripIdx = ((*r) - 1);
                    int vPre = rh->TripsInfoList[m_PreviousTripIdx].VehicleIndex; // Vehicle Index For Previous Trip
                    int vc   = rh->TripsInfoList[*r].VehicleIndex; //  Vehicle Index For Current Trip
                    int dc   = S_NPD[k][m_PreviousTripIdx][vPre];
                    if (dc != -1)
                    {
                        lp->add((*Tval)[RPS[k]] >= ((*Zval)[k][*r][vc]) * StopsTeta[dc] - M * (1-(*Zval)[k][*r][vc]));
                    }
                }
            }
        }
    }
    ////////////////////////////// C-24 (Constraint-18 In Second Paper)
    for(Index k=0 ; k < TotalReq ; k++)
    {
        if(TotalReqType[k] == REQUESTTYPE::PND)
        {
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            int i = RDS[k];
            for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                 r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
            {
                int cTrip = *r;
                int nTrip = ((*r) + 1);
                if(cTrip != rh->TripsNo-1)
                {
                    int v  = rh->TripsInfoList[cTrip].VehicleIndex;
                    int vn = rh->TripsInfoList[nTrip].VehicleIndex;
                    int pc = S_PND[k][cTrip][v];
                    int pc_next = S_PND[k][nTrip][vn];

                    if(pc == -1 || pc_next == -1)
                        continue;

                    for(Index j=0 ; j < nbNodes ; j++)
                    {
                        // Check 1 : j variable condition
                        // Check 2 : i variable condition
                        if(j > pc && j <= pc_next
                                && CPNodes.find(j) != CPNodes.end()
                                && ArcsMatrix[i][j] == 1)
                        {
                            int vi = utility->GetVehicleNumberOfCurrentNode(rh,i);
                            int vj = utility->GetVehicleNumberOfCurrentNode(rh,j);
                            int vnode = 0;
                            if(vi != vj)
                            {
                                if(i > j)
                                    vnode = vi;
                                else
                                    vnode = vj;
                            }
                            else
                            {
                                vnode = vi;
                            }
                            lp->add((*Xval)[i][j][vnode] <= (*Zval)[k][cTrip][v]);
                        }
                    }
                }
            }
        }
    }
    ////////////////////////////// C-25 (Constraint-19 In Second Paper)
    for(Index k=0 ; k < TotalReq ; k++)
    {
        if(TotalReqType[k] == REQUESTTYPE::PND)
        {
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            int j = RDS[k];
            for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                 r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
            {
                int cTrip = *r;
                int nTrip = ((*r) + 1);
                if(cTrip != rh->TripsNo-1)
                {
                    int v  = rh->TripsInfoList[cTrip].VehicleIndex;
                    int vn = rh->TripsInfoList[nTrip].VehicleIndex;

                    int pc = S_PND[k][cTrip][v];
                    int pc_next = S_PND[k][nTrip][vn];

                    if(pc == -1 || pc_next == -1)
                        continue;

                    for(Index i=0 ; i < nbNodes ; i++)
                    {
                        // Check 1 : j variable condition
                        // Check 2 : i variable condition
                        if(i >= pc && i < pc_next
                                && CPNodes.find(i) != CPNodes.end()
                                && ArcsMatrix[i][j] == 1)
                        {
                            int vi = utility->GetVehicleNumberOfCurrentNode(rh,i);
                            int vj = utility->GetVehicleNumberOfCurrentNode(rh,j);
                            int vnode = 0;
                            if(vi != vj)
                            {
                                if(i > j)
                                    vnode = vi;
                                else
                                    vnode = vj;
                            }
                            else
                            {
                                vnode = vi;
                            }
                            lp->add((*Xval)[i][j][vnode] <= (*Zval)[k][cTrip][v]);
                        }
                    }
                }
            }
        }
    }
    ////////////////////////////// C-26 (Constraint-20 In Second Paper)
    for(Index k=0 ; k < TotalReq ; k++)
    {
        if(TotalReqType[k] == REQUESTTYPE::NPD)
        {
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            int j = RPS[k];
            for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                 r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
            {
                int cTrip = *r;
                int pTrip = ((*r) - 1);
                if(cTrip != 0)
                {
                    int v  = rh->TripsInfoList[cTrip].VehicleIndex;
                    int vp = rh->TripsInfoList[pTrip].VehicleIndex;

                    int dc = S_NPD[k][cTrip][v];
                    int dc_previous = S_NPD[k][pTrip][vp];

                    if(dc == -1 || dc_previous == -1)
                        continue;

                    for(Index i=0 ; i < nbNodes ; i++)
                    {
                        // Check 1 : j variable condition
                        // Check 2 : i variable condition
                        if(i >= dc_previous && i < dc
                                && CPNodes.find(i) != CPNodes.end()
                                && ArcsMatrix[i][j] == 1)
                        {
                            int vi = utility->GetVehicleNumberOfCurrentNode(rh,i);
                            int vj = utility->GetVehicleNumberOfCurrentNode(rh,j);
                            int vnode = 0;
                            if(vi != vj)
                            {
                                if(i > j)
                                    vnode = vi;
                                else
                                    vnode = vj;
                            }
                            else
                            {
                                vnode = vi;
                            }
                            lp->add((*Xval)[i][j][vnode] <= (*Zval)[k][cTrip][v]);
                        }
                    }
                }
            }
        }
    }
    ////////////////////////////// C-27 (Constraint-21 In Second Paper)
    for(Index k=0 ; k < TotalReq ; k++)
    {
        if(TotalReqType[k] == REQUESTTYPE::NPD)
        {
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            int i = RPS[k];
            for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                 r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
            {
                int cTrip = *r;
                int pTrip = ((*r) - 1);
                if(cTrip != 0)
                {
                    int v  = rh->TripsInfoList[cTrip].VehicleIndex;
                    int vp = rh->TripsInfoList[pTrip].VehicleIndex;

                    int dc = S_NPD[k][cTrip][v];
                    int dc_previous = S_NPD[k][pTrip][vp];

                    if(dc == -1 || dc_previous == -1)
                        continue;

                    for(Index j=0 ; j < nbNodes ; j++)
                    {
                        // Check 1 : j variable condition
                        // Check 2 : i variable condition
                        if(j > dc_previous && j <= dc
                                && CPNodes.find(j) != CPNodes.end()
                                && ArcsMatrix[i][j] == 1)
                        {
                            int vi = utility->GetVehicleNumberOfCurrentNode(rh,i);
                            int vj = utility->GetVehicleNumberOfCurrentNode(rh,j);
                            int vnode = 0;
                            if(vi != vj)
                            {
                                if(i > j)
                                    vnode = vi;
                                else
                                    vnode = vj;
                            }
                            else
                            {
                                vnode = vi;
                            }
                            lp->add((*Xval)[i][j][vnode] <= (*Zval)[k][cTrip][v]);
                        }
                    }
                }
            }
        }
    }
#endif // End Group_1

#ifndef GROUP_2

    // (C-56)
//    for (Index k = 0; k < TotalReq; k++)
//    {
//        if (TotalReqType[k] == REQUESTTYPE::NPD) {
//            int DetectedTrip = GetTripNumberOfRequest(rh,k);

//            for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
//                 r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
//            {
//                int v = rh->TripsInfoList[*r].VehicleIndex;
//                int dc = S_NPD[k][*r][v];
//                if(dc == -1)
//                    continue;
//                IloExpr sumZ(env);
//                // Modification to strengthen the constraint
//                for (int rp=*r ; rp < TripsNo ; rp++)
//                {
//                    int vl = rh->TripsInfoList[rp].VehicleIndex;
//                    sumZ.operator+=((*Zval)[k][rp][vl]);
//                }
//                lp->add((*Dval)[k] >=((*Tbarval)[dc]-M*(1-((*Zval)[k][*r][v]))));
//                lp->add((*Dval)[k] >=((*Tbarval)[dc]-M*(1-sumZ)));
//            }
//        }
//    }
    // New constraint for the pickup of PND customers: (C-57)
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::PND) {
            IloExpr sumZtheta(env);
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                 r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
            {
                int v = rh->TripsInfoList[*r].VehicleIndex;
                int pc = S_PND[k][*r][v];
                if(pc != -1) // does this define HYBR(k)
                    sumZtheta.operator+=((*Zval)[k][*r][v]*StopsTeta[pc]);
            }
            lp->add((*Pval)[k] == sumZtheta);
        }
    }
    // New valid inequality for the dropoff of NPD customers: (C-58)
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::NPD) {
            IloExpr sumZtheta(env);
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                 r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
            {
                int v = rh->TripsInfoList[*r].VehicleIndex;
                int dc = S_NPD[k][*r][v];
                if(dc != -1) // does this define HYBR(k)
                    sumZtheta.operator+=((*Zval)[k][*r][v]*(StopsTeta[dc]-rh->Input_SlackTime));
            }
            lp->add((*Dval)[k] >= sumZtheta);
        }
    }
    // (C-59)
    fstream deltamin;
    deltamin.open("OUTPUT/deltamin"+to_string(cntInputFile)+".txt",ios_base::out | ios_base::trunc);
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::NPND)
        {
            deltamin << "Minimum Delta For NPND Request ["<<k<<"] = " << m_DeltaMinReq[k] << endl;
            lp->add((*Dval)[k] >= (*Pval)[k] + m_DeltaMinReq[k]);
        }
    }

    // New valid inequality for the checkpoints: (C-60)
    for(Index i=0 ; i< nbNodes ; i++)
    {
        int idx_i = Nodes[i];
        if (StopsType[i] == 1)
        {
            lp->add((*Tbarval)[idx_i] >= (*Tval)[idx_i]-rh->Input_SlackTime);
        }
    }
    // (C-61)
    fstream lcfile;
    lcfile.open("OUTPUT/lc"+to_string(cntInputFile)+".txt",ios_base::out | ios_base::trunc);
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::PND)
        {
            IloExpr sumZtheta(env);
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                 r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
            {
                int v = rh->TripsInfoList[*r].VehicleIndex;
                int lc = LC[k][*r][v];

                if(lc == -1)
                    continue;
                lcfile << "LC["<<k<<"]["<<*r<<"]["<<v<<"] = " << lc << endl;
                double t = utility->ComputeDeltaMinReqFromPickToDrop(CPNodes,rh,lc,RDS[k]);
                deltamin << "Minimum Delta For Request ["<<k<<"] From " << "LC["<<k<<"]["<<*r<<"]["<<v<<"] = " << lc << " To DS["<<k<<"] = " << RDS[k] << " == " << t<< endl;
                sumZtheta+=(StopsTeta[lc]*(*Zval)[k][*r][v] + t);
            }

            lp->add((*Dval)[k] >= sumZtheta);
        }
    }
    // (C-62)
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::NPD)
        {
            IloExpr sumZtheta(env);
            int DetectedTrip = utility->GetTripNumberOfRequest(rh,k);
            for (vector<int>::const_iterator r = rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.begin();
                 r != rh->TripsInfoList[DetectedTrip].m_PossibleTripsForRequests.end(); r++)
            {
                int v = rh->TripsInfoList[*r].VehicleIndex;
                int lc = LC[k][*r][v];
                if(lc == -1)
                    continue;
                lcfile << "LC["<<k<<"]["<<*r<<"]["<<v<<"] = " << lc << endl;
                double t = utility->ComputeDeltaMinReqFromPickToDrop(CPNodes,rh,lc,RPS[k]);
                deltamin << "Minimum Delta For Request ["<<k<<"] From " << "LC["<<k<<"]["<<*r<<"]["<<v<<"] = " << lc << " To PS["<<k<<"] = " << RPS[k] << " == " << t<< endl;
                sumZtheta+=(StopsTeta[lc]*(*Zval)[k][*r][v] + t);
            }

            lp->add((*Pval)[k] >= sumZtheta);
        }
    }
    deltamin.close();
    lcfile.close();
#endif
    ConsoleOutput.close();
//    double t2 = clock();
//    double diff = (double)(t2 - t1) / (double)CLOCKS_PER_SEC;
}

void Solver::ComputeDeltaMinRequest(const SolverBase* sb)
{
    int TripsNo = sb->TripsNo;

    Index TotalReq = 0;
    vector<int> TotalReqType;
    vector <int> RPS;
    vector <int> RDS;
    vector <double> StopsB = sb->StopsB;
    vector <double> RequestsTau;
    vector <int> StopsType;
    vector <double> StopsTeta = sb->StopsTeta;
    Index nbNodes = sb->nbPoints;

    for(int i=0 ; i <TripsNo ; i++)
    {
        TotalReq += sb->TripsInfoList[i].RequestNo;
        for(int j=0 ; j < sb->TripsInfoList[i].RequestNo ; j++)
        {
            TotalReqType.push_back(sb->TripsInfoList[i].RequestsType[j]);
        }
    }
    for(int i=0 ; i < TripsNo ; i++)
    {
        for(Index j=0 ; j < sb->TripsInfoList[i].RequestsPickStopNo.size() ; j++)
            RPS.push_back(sb->TripsInfoList[i].RequestsPickStopNo[j]);

        for(Index j=0 ; j < sb->TripsInfoList[i].RequestsDropStopNo.size() ; j++)
            RDS.push_back(sb->TripsInfoList[i].RequestsDropStopNo[j]);

        for(Index j=0 ; j < sb->TripsInfoList[i].RequestsTau.size() ; j++)
            RequestsTau.push_back(sb->TripsInfoList[i].RequestsTau[j]);

        for(Index j=0 ; j < sb->TripsInfoList[i].StopsType.size() ; j++)
            StopsType.push_back(sb->TripsInfoList[i].StopsType[j]);
    }

    for(Index k=0 ; k < TotalReq ; k++)
    {
        if(TotalReqType[k] == NPND)
        {
            int PickUp  = RPS[k];
            int DropOff = RDS[k];
            int lastIdx = PickUp;
            double t = 0;
            for(int i = PickUp+1 ; i <= DropOff ; i++)
            {
                if(CPNodes.find(i) != CPNodes.end() || i == DropOff)
                {
                    Coord LastCoord = utility->GetCoordinationFromIndex(sb,lastIdx);
                    Coord CheckpointCoord = utility->GetCoordinationFromIndex(sb,i);
                    t += utility->CalculateDistance(LastCoord,CheckpointCoord);
                    lastIdx = i;
                }
            }

            m_DeltaMinReq.push_back(t);
        }
        else
        {
            m_DeltaMinReq.push_back(0.0);
        }
    }
}

void Solver::ComputeLC(const SolverBase *sb)
{
    int TripsNo = sb->TripsNo;

    Index TotalReq = 0;
    vector<int> TotalReqType;
    vector <int> RPS;
    vector <int> RDS;
    vector <double> StopsB = sb->StopsB;
    vector <double> RequestsTau;
    vector <int> StopsType;
    vector <double> StopsTeta = sb->StopsTeta;
    Index nbNodes = sb->nbPoints;

    for(int i=0 ; i <TripsNo ; i++)
    {
        TotalReq += sb->TripsInfoList[i].RequestNo;
        for(int j=0 ; j < sb->TripsInfoList[i].RequestNo ; j++)
        {
            TotalReqType.push_back(sb->TripsInfoList[i].RequestsType[j]);
        }
    }
    for(int i=0 ; i < TripsNo ; i++)
    {
        for(Index j=0 ; j < sb->TripsInfoList[i].RequestsPickStopNo.size() ; j++)
            RPS.push_back(sb->TripsInfoList[i].RequestsPickStopNo[j]);

        for(Index j=0 ; j < sb->TripsInfoList[i].RequestsDropStopNo.size() ; j++)
            RDS.push_back(sb->TripsInfoList[i].RequestsDropStopNo[j]);

        for(Index j=0 ; j < sb->TripsInfoList[i].RequestsTau.size() ; j++)
            RequestsTau.push_back(sb->TripsInfoList[i].RequestsTau[j]);

        for(Index j=0 ; j < sb->TripsInfoList[i].StopsType.size() ; j++)
            StopsType.push_back(sb->TripsInfoList[i].StopsType[j]);
    }

    for(Index k=0 ; k < TotalReq ; k++)
    {
        if(TotalReqType[k] == PND)
        {
            int DropOff = RDS[k];
            int r = utility->GetTripNumberOfRequest(sb,k);
            int v = sb->TripsInfoList[r].VehicleIndex;
            Coord DropCoord = utility->GetCoordinationFromIndex(sb,DropOff);
            vector<int> NearestCP = utility->GetNearestCPToCurrentPos(GraphCPGlobalList,GraphNodes,DropCoord);
            int BackCP = NearestCP[0];
            LC[k][r][v] = BackCP;

        }
        if(TotalReqType[k] == NPD)
        {
            int PickUp = RPS[k];
            int r = utility->GetTripNumberOfRequest(sb,k);
            int v = sb->TripsInfoList[r].VehicleIndex;
            Coord PickCoord = utility->GetCoordinationFromIndex(sb,PickUp);
            vector<int> NearestCP = utility->GetNearestCPToCurrentPos(GraphCPGlobalList,GraphNodes,PickCoord);
            int BackCP = NearestCP[0];
            LC[k][r][v] = BackCP;
        }
    }
}

void Solver::Free_Problem()
{
    lp->end();
    env.end();
    if (lp != NULL)
    {
        delete lp;
    }

    if(utility != nullptr)
    {
        delete utility;
    }

}

void Solver::HandleXZQParams(const SolverBase *sb, const map<int, int> &Z, const map<int, int> &X, const map<int, int> &Q)
{
    IloNumVarArray WarmStartVar_X(env);
    IloNumArray WarmStartVal_X(env);
    IloNumVarArray WarmStartVar_Z(env);
    IloNumArray WarmStartVal_Z(env);
    IloNumVarArray WarmStartVar_Q(env);
    IloNumArray WarmStartVal_Q(env);

    int NodesNo = sb->nbPoints;
    int VehiclesNo = sb->VehiclesNo;
    int RequestNo = 0;
    int TripsNo = sb->TripsNo;
    vector <int> ReqType;
    for(int i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo &ti = sb->TripsInfoList[i];
        RequestNo += ti.RequestNo;
        for(int j=0 ; j < ti.RequestsTau.size() ; j++)
            ReqType.push_back(ti.RequestsType[j]);
    }


    for(map<int,int>::const_iterator it = Z.begin() ; it != Z.end() ; it++)
    {
        int k = it->first;
        int r = it->second;
        int v = sb->TripsInfoList[r].VehicleIndex;

        if(ReqType[k] == PND)
        {
            WarmStartVar_Z.add((*Zval)[k][r][v]);
            WarmStartVal_Z.add(1.0);
        }

        if(ReqType[k] == NPD)
        {
            WarmStartVar_Z.add((*Zval)[k][r][v]);
            WarmStartVal_Z.add(1.0);
        }

        lp->add((*Zval)[k][r][v] == 1);
    }
    for(map<int,int>::const_iterator it = X.begin() ; it != X.end() ; it++)
    {
        int i = it->first;
        int j = it->second;
        int vi = utility->GetVehicleNumberOfCurrentNode(sb,i);
        int vj = utility->GetVehicleNumberOfCurrentNode(sb,j);
        int v = 0;
        if(vi != vj)
        {
            if(i > j)
                v = vi;
            else
                v = vj;
        }
        else
        {
            v = vi;
        }
        if(ArcsMatrix[i][j] == 0)
        {
            cerr << "No Arc Between " << i << "," << j << endl;
        }
        WarmStartVar_X.add((*Xval)[i][j][v]);
        WarmStartVal_X.add(1.0);

        lp->add((*Xval)[i][j][v] == 1);
    }
    for(map<int,int>::const_iterator it = Q.begin() ; it != Q.end() ; it++)
    {
        int i = it->first;
        int n = it->second;
        int v = utility->GetVehicleNumberOfCurrentNode(sb,i);
//        lp->add((*Qval)[i][v] == n);
        WarmStartVar_Q.add((*Qval)[i][v]);
        WarmStartVal_Q.add(n);
    }

    cplexPtr = new IloCplex(*lp);

    string ModelName = "Model-"+to_string(0)+".lp";
    cplexPtr->exportModel(ModelName.c_str());
    double tL = 18000;
    cplexPtr->setParam(IloCplex::TiLim, tL);		        //Run Time = 1 Hours
    cplexPtr->setParam(IloCplex::EpGap, 0.0);
//    cplexPtr->setParam(IloCplex::PreInd,false);



//    try{
//        cplexPtr->addMIPStart(WarmStartVar_Z,WarmStartVal_Z,IloCplex::MIPStartSolveMIP,"WarmStartVar_Z");
//    }catch(const IloException &e)
//    {
//        std::cerr << __FUNCTION__ << ":WarmStartVar_Z Start :: " << e << "\n";
//    }
//    try{
//        cplexPtr->addMIPStart(WarmStartVar_X,WarmStartVal_X,IloCplex::MIPStartAuto,"WarmStartVar_X");

//    }catch(const IloException &e)
//    {
//        std::cerr << __FUNCTION__ << ":WarmStartVar_X Start :: " << e << "\n";
//    }

//    cplexPtr->writeMIPStarts("WarmStart");

    WarmStartVal_X.end();
    WarmStartVar_X.end();
    WarmStartVal_Z.end();
    WarmStartVar_Z.end();
    WarmStartVal_Q.end();
    WarmStartVar_Q.end();

}

void Solver::AddDataToWarmStart(const int &ins, const SolverBase *sb)
{
    int NodesNo = sb->nbPoints;
    int VehiclesNo = sb->VehiclesNo;
    int RequestNo = 0;
    int TripsNo = sb->TripsNo;
    vector <int> ReqType;
    for(int i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo &ti = sb->TripsInfoList[i];
        RequestNo += ti.RequestNo;
        for(int j=0 ; j < ti.RequestsTau.size() ; j++)
            ReqType.push_back(ti.RequestsType[j]);
    }

    IloNumVarArray WarmStartVar_X(env);
    IloNumArray WarmStartVal_X(env);
    IloNumVarArray WarmStartVar_Z(env);
    IloNumArray WarmStartVal_Z(env);
    IloNumVarArray WarmStartVar_Q(env);
    IloNumArray WarmStartVal_Q(env);
    IloNumVarArray WarmStartVar_T(env);
    IloNumArray WarmStartVal_T(env);
    IloNumVarArray WarmStartVar_TBar(env);
    IloNumArray WarmStartVal_TBar(env);
    IloNumVarArray WarmStartVar_P(env);
    IloNumArray WarmStartVal_P(env);
    IloNumVarArray WarmStartVar_D(env);
    IloNumArray WarmStartVal_D(env);
    IloNumVarArray WarmStartVar_q(env);
    IloNumArray WarmStartVal_q(env);

    for(std::map<Coord,int>::const_iterator it = sb->Xmap.begin() ; it != sb->Xmap.end() ; it++)
    {
        int i = (int)it->first.first;
        int j = (int)it->first.second;
        int v = it->second;
        if(ArcsMatrix[i][j] == 0)
        {
            cerr << "No Arc Between " << i << "," << j << endl;
        }
        WarmStartVar_X.add((*Xval)[i][j][v]);
        WarmStartVal_X.add(1.0);
    }

    for(map<int,double>::const_iterator it = sb->Tmap.begin() ; it != sb->Tmap.end() ; it++)
    {
        int key = it->first;
        double val = it->second;
        WarmStartVar_T.add((*Tval)[key]);
        WarmStartVal_T.add(val);
    }

    for(map<int,double>::const_iterator it = sb->TBarmap.begin() ; it != sb->TBarmap.end() ; it++)
    {
        int key = it->first;
        double val = it->second;
        WarmStartVar_TBar.add((*Tbarval)[key]);
        WarmStartVal_TBar.add(val);
    }

    for(map<int,double>::const_iterator it = sb->Pmap.begin() ; it != sb->Pmap.end() ; it++)
    {
        int key = it->first;
        double val = it->second;
        WarmStartVar_P.add((*Pval)[key]);
        WarmStartVal_P.add(val);
    }

    for(map<int,double>::const_iterator it = sb->Dmap.begin() ; it != sb->Dmap.end() ; it++)
    {
        int key = it->first;
        double val = it->second;
        WarmStartVar_D.add((*Dval)[key]);
        WarmStartVal_D.add(val);
    }

    for(std::map<Coord,int>::const_iterator it = sb->Zmap.begin() ; it != sb->Zmap.end() ; it++)
    {
        int k = (int)it->first.first;
        int r = (int)it->first.second;
        int v = it->second;

        if(ReqType[k] == PND)
        {
            WarmStartVar_Z.add((*Zval)[k][r][v]);
            WarmStartVal_Z.add(1.0);
        }

        if(ReqType[k] == NPD)
        {
            WarmStartVar_Z.add((*Zval)[k][r][v]);
            WarmStartVal_Z.add(1.0);
        }

    }

    for(map<int,int>::const_iterator it = sb->qmap.begin() ; it != sb->qmap.end() ; it++)
    {
        int key = it->first;
        int val = it->second;
        WarmStartVar_q.add((*qval)[key]);
        WarmStartVal_q.add(val);
    }

    for(map<Coord,int>::const_iterator it = sb->Qmap.begin() ; it != sb->Qmap.end() ; it++)
    {
        int i = it->first.first;
        int v = it->first.second;
        int val = it->second;
        WarmStartVar_Q.add((*Qval)[i][v]);
        WarmStartVal_Q.add(val);
    }

    PrepareCplexParams(ins);

    cplexPtr->addMIPStart(WarmStartVar_X,WarmStartVal_X,IloCplex::MIPStartSolveMIP,"WarmStartVar_X");
    cplexPtr->addMIPStart(WarmStartVar_T,WarmStartVal_T,IloCplex::MIPStartSolveMIP,"WarmStartVar_T");
    cplexPtr->addMIPStart(WarmStartVar_TBar,WarmStartVal_TBar,IloCplex::MIPStartSolveMIP,"WarmStartVar_Tbar");
    cplexPtr->addMIPStart(WarmStartVar_P,WarmStartVal_P,IloCplex::MIPStartSolveMIP,"WarmStartVar_P");
    cplexPtr->addMIPStart(WarmStartVar_D,WarmStartVal_D,IloCplex::MIPStartSolveMIP,"WarmStartVar_D");
    cplexPtr->addMIPStart(WarmStartVar_Z,WarmStartVal_Z,IloCplex::MIPStartSolveMIP,"WarmStartVar_Z");
    cplexPtr->addMIPStart(WarmStartVar_q,WarmStartVal_q,IloCplex::MIPStartSolveMIP,"WarmStartVar_q");
    cplexPtr->addMIPStart(WarmStartVar_Q,WarmStartVal_Q,IloCplex::MIPStartSolveMIP,"WarmStartVar_Q");

    WarmStartVal_X.end();
    WarmStartVar_X.end();
    WarmStartVal_Z.end();
    WarmStartVar_Z.end();
    WarmStartVal_Q.end();
    WarmStartVar_Q.end();
    WarmStartVar_T.end();
    WarmStartVal_T.end();
    WarmStartVar_TBar.end();
    WarmStartVal_TBar.end();
    WarmStartVar_P.end();
    WarmStartVal_P.end();
    WarmStartVar_D.end();
    WarmStartVal_D.end();
    WarmStartVar_q.end();
    WarmStartVal_q.end();

}

void Solver::doSolve(const int &ins, const SolverBase *sb)
{
    cerr << __FUNCTION__ << ": Started ... " << endl;
    int NodesNo = sb->nbPoints;
    int VehiclesNo = sb->VehiclesNo;
    int RequestNo = 0;
    int TripsNo = sb->TripsNo;
    vector<vector<double>> m_GraphDelta = sb->delta;
    vector<double> ReqTau;
//    vector <int> ReqType;
    for(int i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo &ti = sb->TripsInfoList[i];
        RequestNo += ti.RequestNo;
        for(int j=0 ; j < ti.RequestsTau.size() ; j++)
        {
            ReqTau.push_back(ti.RequestsTau[j]);
        }
    }

    PrepareCplexParams(ins);

    CplexSolve();

    TotalTravelTime = 0.0;
    TotalRequestTime = 0.0;
    TotalRequestWaitingTime = 0.0;

    try{
    for (int i = 0; i < NodesNo; i++)
    {
        for (int j = 0; j < NodesNo; j++)
        {
            for (int v = 0; v < VehiclesNo; v++)
            {
                if(ArcsMatrix[i][j] == 0)
                    continue;
                if (cplexPtr->getValue((*Xval)[i][j][v]) > 0.5) {
                    TotalTravelTime += (m_GraphDelta[i][j]*cplexPtr->getValue((*Xval)[i][j][v]));
                }
            }
        }
    }
    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":X Parameters :: " << e << "\n";
    }
    TotalTravelTime *= sb->Input_W0;

    double SumP = 0.0;
    double SumD = 0.0;
    try{
        for (int k = 0; k < RequestNo; k++)
        {
            SumP = cplexPtr->getValue((*Pval)[k]);
            SumD = cplexPtr->getValue((*Dval)[k]);
            TotalRequestWaitingTime += (cplexPtr->getValue((*Pval)[k]) - ReqTau[k]);
            TotalRequestTime += (SumD - SumP < std::numeric_limits<double>::epsilon() ? 0 : SumD - SumP);
        }
    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":P Parameters :: " << e << "\n";
    }

    TotalRequestWaitingTime *= sb->Input_W2;
    TotalRequestTime *= sb->Input_W1;

    SumRej = 0.0;

    cerr << __FUNCTION__ << ": TotalTravelTime :: " << TotalTravelTime << endl;
    cerr << __FUNCTION__ << ": TotalRequestTime :: " << TotalRequestTime << endl;
    cerr << __FUNCTION__ << ": TotalRequestWaitingTime :: " << TotalRequestWaitingTime << endl;

    cerr << __FUNCTION__ << ": Objective Value :: " << objval << endl;
    cerr << __FUNCTION__ << ": GAP :: " << gap_p << endl;
}

void Solver::PrepareCplexParams(const int& ins)
{
    if(cplexPtr == nullptr)
    {
        cplexPtr = new IloCplex(*lp);
        string ModelName = "Model-"+to_string(ins)+".lp";
        cplexPtr->exportModel(ModelName.c_str());
        double tL = 18000;
        cplexPtr->setParam(IloCplex::TiLim, tL);		        //Run Time = 1 Hours
        cplexPtr->setParam(IloCplex::EpGap, 0.0);
    }
}

void Solver::CplexSolve()
{
    try{
        double t1 = clock();
        status = cplexPtr->solve();						//cout<<"Solving . . .";
        objval = cplexPtr->getObjValue();
        double LB = cplexPtr->getBestObjValue();
        gap_p = (objval-LB)/objval;
    }catch(const IloException& e)
    {
        std::cerr << "\n\nCPLEX Raised an exception:\n";
        std::cerr << e << "\n";
    }

    if (!status)
    {
        std::cerr << "\n\nCplex error!\n";
        std::cerr << "\tStatus: " << cplexPtr->getStatus() << "\n";
        std::cerr << "\tSolver status: " << cplexPtr->getCplexStatus() << "\n";
    }
    else
    {
        std::cout << "\n\nCplex success!\n";
        std::cout << "\tStatus: " << cplexPtr->getStatus() << "\n";
        std::cout << "\tObjective value: " << cplexPtr->getObjValue() << "\n";
    }
}

void Solver::CreateArc(SolverBase* sb,const int &NodesNo, const int &TripsNo)
{
    ArcsMatrix.resize(GraphNodes.size());						//arcsmatrix ro resize mikonim bar asase tedade nodeha,
    for(int i=0 ; i < NodesNo ; i++)
    {
        ArcsMatrix[i].resize(NodesNo); 							//har element az arcs matriz shamele yal ha va har row be andaze tedade node ha resize mishe
        for(int j = 0 ; j < NodesNo ; j++)
        {
            ArcsMatrix[i][j] = 0;								//va init mikonim be sefr dar ebtedaye kar ke hanooz khaliye
        }
    }

    CoordSet d_CPSet;
    set<int> d_CPIndices;
    for(int i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo& ti = sb->TripsInfoList[i];		//miad tripsinfo object ro ba indice i az tripsinfolist migire k too rh zakhire shode
        const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
        const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination
        for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)		//istgah ha va gheyre istgah haro miad migire
        {																							// va mirize tooye set graphnodes.
            d_CPSet.insert(*cit);
        }
    }

    for(int i=0 ; i < GraphNodesList.size() ; i++)
    {
        const Coord& i_coord = GraphNodesList[i]; // current node
        for(int j=i+1 ; j < GraphNodesList.size() ; j++)
        {
            const Coord& j_coord = GraphNodesList[j]; // next node
            if(d_CPSet.find(j_coord) != d_CPSet.end())
            {
                // Next node is an checkpoint
                ArcsMatrix[i][j] = 1;
                break;
            }
            else
            {
                //next node is an non-checkpoint
                ArcsMatrix[i][j] = 1;
            }
        }
    }

    int TotalReq = 0;
    vector<int> TotalReqType;
    vector <int> RPS;
    vector <int> RDS;

    for(int i=0 ; i <TripsNo ; i++)
    {
        TotalReq += sb->TripsInfoList[i].RequestNo;
        for(int j=0 ; j < sb->TripsInfoList[i].RequestNo ; j++)
        {
            TotalReqType.push_back(sb->TripsInfoList[i].RequestsType[j]);
        }
    }
    for(int i=0 ; i < TripsNo ; i++)
    {
        for(int j=0 ; j < sb->TripsInfoList[i].RequestsPickStopNo.size() ; j++)
            RPS.push_back(sb->TripsInfoList[i].RequestsPickStopNo[j]);

        for(int j=0 ; j < sb->TripsInfoList[i].RequestsDropStopNo.size() ; j++)
            RDS.push_back(sb->TripsInfoList[i].RequestsDropStopNo[j]);
    }

    /**
     * @brief :: Handle Arcs for pnd requests
    **/
    for(int k=0 ; k < S_PND.size() ; k++)
    {
        const int& DropIdx = RDS[k];
        int detctedTrip = utility->GetTripNumberOfRequest(sb,k);
        for(int r=0 ; r < S_PND[k].size() ; r++)
        {
            for(int v=0 ; v < S_PND[k][r].size() ; v++)
            {
                if(TotalReqType[k] == PND)
                {
                    if(S_PND[k][r][v] == -1)
                        continue;
                    int orgPnd = S_PND[k][detctedTrip][v];
                    const int & PickIdx = S_PND[k][r][v];
                    const Coord& PickCoord = utility->GetCoordinationFromIndex(sb,PickIdx);
                    const Coord& DropCoord = utility->GetCoordinationFromIndex(sb,DropIdx);
                    if(orgPnd == PickIdx)
                    {
                        if(std::fabs(DropCoord.first - PickCoord.first) > minDisCP)
                        {
                            for(int t=DropIdx ; t < GraphNodesList.size() ; t++)
                            {
                                if(CPNodes.find(t) != CPNodes.end())
                                {
                                    ArcsMatrix[DropIdx][t] = 1;
                                    double distance = utility->GetDistanceFromSrcToDes(sb,DropIdx,t);
                                    sb->delta[DropIdx][t] = distance;
                                    break;
                                }
                            }
                            for(int t=DropIdx ; t >=0 ; t--)
                            {
                                if(CPNodes.find(t) != CPNodes.end())
                                {
                                    ArcsMatrix[t][DropIdx] = 1;
                                    double distance = utility->GetDistanceFromSrcToDes(sb,DropIdx,t);
                                    sb->delta[t][DropIdx] = distance;
                                    break;
                                }
                            }
                            CoordList DropPosSimilar_Next = utility->GetSimilarPositionOfNode(sb,DropCoord);

                            for(int a=0 ; a < DropPosSimilar_Next.size() ; a++)
                            {
                                Coord SimilarDropPos = DropPosSimilar_Next[a];
                                vector<int> nearestCPToCurrentNode = utility->GetNearestCPToCurrentPos(GraphCPGlobalList,GraphNodes,SimilarDropPos);
                                for(int b = nearestCPToCurrentNode[0] ; b <= nearestCPToCurrentNode[1] ; b++)
                                {
                                    int d_cpIdxLocal = b;
                                    if(d_cpIdxLocal == DropIdx)
                                        continue;
                                    Coord d_cpIdxLocal_Pos = utility->GetCoordinationFromIndex(sb,d_cpIdxLocal);
                                    if(d_cpIdxLocal_Pos.first < SimilarDropPos.first)
                                    {
                                        ArcsMatrix[d_cpIdxLocal][DropIdx] = 1;
                                        double distance = utility->GetDistanceFromSrcToDes(sb,DropIdx,d_cpIdxLocal);
                                        sb->delta[d_cpIdxLocal][DropIdx] = distance;
                                    }
                                    else
                                    {
                                        ArcsMatrix[DropIdx][d_cpIdxLocal] = 1;
                                        double distance = utility->GetDistanceFromSrcToDes(sb,DropIdx,d_cpIdxLocal);
                                        sb->delta[DropIdx][d_cpIdxLocal] = distance;
                                    }
                                }
                            }
                            continue;
                        }
                    }
                    else
                    {
                        if(ArcsMatrix[orgPnd][DropIdx] == 0)
                            continue;
                    }

                    ArcsMatrix[PickIdx][DropIdx] = 1;

                    double distance = utility->GetDistanceFromSrcToDes(sb,DropIdx,PickIdx);
                    sb->delta[PickIdx][DropIdx] = distance;

                    double xdiff = std::fabs(DropCoord.first - PickCoord.first);
                    for(int i=0 ; i < SameCheckpointIndicesList.size() ; i++)
                    {
                        const pair<int,vector<int>>& idxPair = SameCheckpointIndicesList[i];
                        if(idxPair.first == PickIdx)
                        {
                            for(int j=0 ; j < idxPair.second.size() ; j++)
                            {
                                const int& similarCP_Index = idxPair.second[j];
                                const Coord& similarCP_Coord = utility->GetCoordinationFromIndex(sb,similarCP_Index);
                                Coord implicitCoord;
                                implicitCoord.first = similarCP_Coord.first + xdiff;
                                implicitCoord.second = similarCP_Coord.second;
                                for(int q=similarCP_Index ; q < GraphNodesList.size() ; q++)
                                {
                                    const Coord& qCoord = GraphNodesList[q];
                                    for(int e=q+1 ; e < GraphNodesList.size() ; e++)
                                    {
                                        const Coord& eCoord = GraphNodesList[e];
                                        if(implicitCoord.first >= qCoord.first
                                                && implicitCoord.first <= eCoord.first)
                                        {
                                            ArcsMatrix[DropIdx][e] = 1;
                                            double distance = utility->GetDistanceFromSrcToDes(sb,DropIdx,e);
                                            sb->delta[DropIdx][e] = distance;
                                            if(d_CPSet.find(eCoord)!= d_CPSet.end())
                                                break;
                                        }
                                        else if(implicitCoord.first > qCoord.first
                                                && implicitCoord.first > eCoord.first)
                                        {
                                            if(d_CPSet.find(eCoord)!= d_CPSet.end())
                                                break;
                                            ArcsMatrix[e][DropIdx] = 1;
                                            double distance = utility->GetDistanceFromSrcToDes(sb,DropIdx,e);
                                            sb->delta[e][DropIdx] = distance;

                                        }
                                    }
                                }

                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief :: Handle Arcs for npd requests
    **/
    for(int k=0 ; k < S_NPD.size() ; k++)
    {
        const int& PickIdx = RPS[k];
        int detctedTrip = utility->GetTripNumberOfRequest(sb,k);
        for(int r=0 ; r < S_NPD[k].size() ; r++)
        {
            for(int v=0 ; v < S_NPD[k][r].size() ; v++)
            {
                if(TotalReqType[k] == NPD)
                {
                    if(S_NPD[k][r][v] == -1)
                        continue;
                    const int & DropIdx = S_NPD[k][r][v];
                    int orgDrop = S_NPD[k][detctedTrip][v];
                    const Coord& PickCoord = utility->GetCoordinationFromIndex(sb,PickIdx);
                    const Coord& DropCoord = utility->GetCoordinationFromIndex(sb,DropIdx);

                    if(orgDrop == DropIdx)
                    {
                        if(std::fabs(DropCoord.first - PickCoord.first) > minDisCP)
                        {
                            int nearestCPnext = -1;
                            for(int t=PickIdx ; t < GraphNodesList.size() ; t++)
                            {
                                if(CPNodes.find(t) != CPNodes.end())
                                {
                                    ArcsMatrix[PickIdx][t] = 1;
                                    double distance = utility->GetDistanceFromSrcToDes(sb,PickIdx,t);
                                    sb->delta[PickIdx][t] = distance;
                                    nearestCPnext = t;
                                    break;
                                }
                            }
                            for(int t=PickIdx ; t >=0 ; t--)
                            {
                                if(CPNodes.find(t) != CPNodes.end())
                                {
                                    ArcsMatrix[t][PickIdx] = 1;
                                    double distance = utility->GetDistanceFromSrcToDes(sb,PickIdx,t);
                                    sb->delta[t][PickIdx] = distance;
                                    break;
                                }
                            }
//                            Coord rootPosnext = GetCoordinationFromIndex(sb,nearestCPnext);
//                            vector<int> resnext = GetSameCPIndexFromRoot(sb,rootPosnext,nearestCPnext);
//                            for(int u=0 ; u < resnext.size() ; u++)
//                            {
//                                ArcsMatrix[PickIdx][resnext[u]] = 1;
//                            }
                            CoordList PickPosSimilar_Next = utility->GetSimilarPositionOfNode(sb,PickCoord);

                            for(int a=0 ; a < PickPosSimilar_Next.size() ; a++)
                            {
                                Coord SimilarPickPos = PickPosSimilar_Next[a];
                                vector<int> nearestCPToCurrentNode = utility->GetNearestCPToCurrentPos(GraphCPGlobalList,GraphNodes,SimilarPickPos);
                                for(int b = nearestCPToCurrentNode[0] ; b <= nearestCPToCurrentNode[1] ; b++)
                                {
                                    int d_cpIdxLocal = b;
                                    if(d_cpIdxLocal == PickIdx)
                                        continue;
                                    Coord d_cpIdxLocal_Pos = utility->GetCoordinationFromIndex(sb,d_cpIdxLocal);
                                    if(d_cpIdxLocal_Pos.first < SimilarPickPos.first)
                                    {
                                        ArcsMatrix[d_cpIdxLocal][PickIdx] = 1;
                                        double distance = utility->GetDistanceFromSrcToDes(sb,PickIdx,d_cpIdxLocal);
                                        sb->delta[d_cpIdxLocal][PickIdx] = distance;
                                    }
                                    else
                                    {
                                        ArcsMatrix[PickIdx][d_cpIdxLocal] = 1;
                                        double distance = utility->GetDistanceFromSrcToDes(sb,PickIdx,d_cpIdxLocal);
                                        sb->delta[PickIdx][d_cpIdxLocal] = distance;
                                    }
                                }
                            }
                            continue;
                        }
                    }
                    else
                    {
                        if(ArcsMatrix[orgDrop][DropIdx] == 0)
                            continue;
                    }

                    ArcsMatrix[PickIdx][DropIdx] = 1;

                    double distance = utility->GetDistanceFromSrcToDes(sb,PickIdx,DropIdx);
                    sb->delta[PickIdx][DropIdx] = distance;

                    double xdiff = std::fabs(DropCoord.first - PickCoord.first);
                    for(int i=0 ; i < SameCheckpointIndicesList.size() ; i++)
                    {
                        const pair<int,vector<int>>& idxPair = SameCheckpointIndicesList[i];
                        if(idxPair.first == DropIdx)
                        {
                            for(int j=0 ; j < idxPair.second.size() ; j++)
                            {
                                const int& similarCP_Index = idxPair.second[j];
                                const Coord& similarCP_Coord = utility->GetCoordinationFromIndex(sb,similarCP_Index);
                                Coord implicitCoord;
                                implicitCoord.first = similarCP_Coord.first - xdiff;
                                implicitCoord.second = PickCoord.second;
                                for(int q=similarCP_Index ; q >= 0 ; q--)
                                {
                                    const Coord& qCoord = GraphNodesList[q];
                                    for(int e=q-1 ; e >=0 ; e--)
                                    {
                                        const Coord& eCoord = GraphNodesList[e];
                                        if(implicitCoord.first <= qCoord.first
                                                && implicitCoord.first >= eCoord.first)
                                        {
                                            ArcsMatrix[e][PickIdx] = 1;
                                            double distance = utility->GetDistanceFromSrcToDes(sb,e,PickIdx);
                                            sb->delta[e][PickIdx] = distance;
                                            if(d_CPSet.find(eCoord)!= d_CPSet.end())
                                                break;
                                        }
                                        else if(implicitCoord.first < qCoord.first
                                                && implicitCoord.first < eCoord.first)
                                        {
                                            if(d_CPSet.find(eCoord)!= d_CPSet.end())
                                                break;
                                            ArcsMatrix[PickIdx][e] = 1;
                                            double distance = utility->GetDistanceFromSrcToDes(sb,e,PickIdx);
                                            sb->delta[PickIdx][e] = distance;

                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief :: Handle Arcs for npnd requests
    **/
    for(Index i=0 ; i < RPS.size() ; i++)
    {
        int Pick = RPS[i];
        if(CPNodes.find(Pick) != CPNodes.end()) // if Pick idx has been a cp idx , we must ignore it because the request is PD
            continue;
        if(Pick != -1)
        {
            Coord PickCoord = utility->GetCoordinationFromIndex(sb,Pick);
            vector<int> CoupleCP = utility->GetNearestCPToCurrentPos(GraphCPGlobalList,GraphNodes,PickCoord);
            if(CoupleCP.size() == 2) // base on graph, we have exatcly 2 checkpoint between each non-checkpoints
            {
                int BackCPIdx = CoupleCP[0];
                int FrontCPIdx = CoupleCP[1];
                Coord BackCPCoord = utility->GetCoordinationFromIndex(sb,BackCPIdx);
                Coord FrontCPCoord = utility->GetCoordinationFromIndex(sb,FrontCPIdx);
                CoordList similarCP_Back  = utility->GetSimilarPositionOfNode(sb,BackCPCoord);
                CoordList similarCP_Front = utility->GetSimilarPositionOfNode(sb,FrontCPCoord);

                double diffX_Pick_Back = std::fabs(PickCoord.first - BackCPCoord.first);

                for(Index k=1 ; k < similarCP_Back.size() ; k++)
                {
                    Coord cpsim_back  = similarCP_Back[k];
                    Coord cpsim_front = similarCP_Front[k];
                    Coord newPosPick;
                    newPosPick.first = cpsim_back.first+diffX_Pick_Back;
                    newPosPick.second = PickCoord.second;
                    int cpsim_back_idx  = utility->GetIndexFromSet(GraphNodes,cpsim_back);
                    int cpsim_front_idx = utility->GetIndexFromSet(GraphNodes,cpsim_front);
                    for(int j = cpsim_back_idx ; j <= cpsim_front_idx ; j++)
                    {
                        Coord nodeCoord = utility->GetCoordinationFromIndex(sb,j);
                        if(nodeCoord.first >= newPosPick.first)
                        {
                            ArcsMatrix[Pick][j] = 1;
                            double distance = utility->GetDistanceFromSrcToDes(sb,Pick,j);
                            sb->delta[Pick][j] = distance;
                        }
                        else
                        {
                            ArcsMatrix[j][Pick] = 1;
                            double distance = utility->GetDistanceFromSrcToDes(sb,Pick,j);
                            sb->delta[j][Pick] = distance;
                        }
                    }
                }
            }
        }
    }
    for(Index i=0 ; i < RDS.size() ; i++)
    {
        int Drop = RDS[i];
        if(CPNodes.find(Drop) != CPNodes.end()) // if Drop idx has been a cp idx , we must ignore it because the request is PD
            continue;
        if(Drop != -1)
        {
            Coord DropCoord = utility->GetCoordinationFromIndex(sb,Drop);
            vector<int> CoupleCP = utility->GetNearestCPToCurrentPos(GraphCPGlobalList,GraphNodes,DropCoord);
            if(CoupleCP.size() == 2) // base on graph, we have exatcly 2 checkpoint between each non-checkpoints
            {
                int BackCPIdx = CoupleCP[0];
                int FrontCPIdx = CoupleCP[1];
                Coord BackCPCoord = utility->GetCoordinationFromIndex(sb,BackCPIdx);
                Coord FrontCPCoord = utility->GetCoordinationFromIndex(sb,FrontCPIdx);
                CoordList similarCP_Back  = utility->GetSimilarPositionOfNode(sb,BackCPCoord);
                CoordList similarCP_Front = utility->GetSimilarPositionOfNode(sb,FrontCPCoord);

                double diffX_Pick_Back = std::fabs(DropCoord.first - BackCPCoord.first);

                for(Index k=1 ; k < similarCP_Back.size() ; k++)
                {
                    Coord cpsim_back  = similarCP_Back[k];
                    Coord cpsim_front = similarCP_Front[k];
                    Coord newPosPick;
                    newPosPick.first = cpsim_back.first+diffX_Pick_Back;
                    newPosPick.second = DropCoord.second;
                    int cpsim_back_idx  = utility->GetIndexFromSet(GraphNodes,cpsim_back);
                    int cpsim_front_idx = utility->GetIndexFromSet(GraphNodes,cpsim_front);
                    for(int j = cpsim_back_idx ; j <= cpsim_front_idx ; j++)
                    {
                        Coord nodeCoord = utility->GetCoordinationFromIndex(sb,j);
                        if(nodeCoord.first >= newPosPick.first)
                        {
                            ArcsMatrix[Drop][j] = 1;
                            double distance = utility->GetDistanceFromSrcToDes(sb,Drop,j);
                            sb->delta[Drop][j] = distance;
                        }
                        else
                        {
                            ArcsMatrix[j][Drop] = 1;
                            double distance = utility->GetDistanceFromSrcToDes(sb,Drop,j);
                            sb->delta[j][Drop] = distance;
                        }
                    }
                }
            }
        }
    }
}
