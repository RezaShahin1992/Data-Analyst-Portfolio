#include "Heuristic.h"
#define GAP_HEU 0.0

Heuristic::Heuristic()
{

}

Heuristic::Heuristic(Solver *s,SolverBase *sb):
    solver(s),
    base(sb)
{
    utility = new Utilities();
}

Heuristic::~Heuristic()
{

}

void Heuristic::doHeuristicSolve()
{
    SaveSolverTrips();							//tamame safarharo mirizim tooye ye araye be esme solvertrips
    PrepareGlobalParams();						//amadesazi voroodihamoon
    ComputeGlobalTheta();					//mohasebe teta be soorat glo

    int TripsNo = base->TripsNo;

    bool isConditionPass = false;

    while(!isConditionPass)				//agar sharte ma pass nashod, bemoon tooye en halghe ta moshkel ro
    {									//bar taraf koni
        for(int i = 0 ; i < TripsNo ; i++)
        {
            TripsInfo& ti = base->TripsInfoList[i];			//hala bia tak tak safararo bardar va
            CheckThetaCondition(ti);			//oon safaro bar asase sharayete tetash barresi kon
            ApplyThetaTransmission(ti);		//hala oonae ke reject shodan bayad beran be safare badi ba en tabe
//            ComputeGlobalTheta();		//dobare teta haro mohasebe mikonim
//            RefreshCheckpointsIndices();	//baraye enke khialemoon rahat beshe moshkeli nist
        }										//va mire safare badi
//        UpdateGraph();						//hame ena ke tamum shod gerafe localo miam be rooz resani mikonim
        Compute_q();				//miam eno farakhani mikonim chon niaz darim behesh baraye mohasebe zarfiat
        const int& qMax = base->Input_VehiclesCapacity; //az too voroodi max zarfiat ro mikhunim
        bool isMoved = false;
        ///:: last trip is your last chance to handle request with capacity condition
        /// so we don't check capacity of last trip
        for(int i=0 ; i < TripsNo-1 ; i++)				//roo kole safar ha ta yeki monde be akhar
        {												//safare akharo barresi nmikonim chon safare badi e nist
            TripsInfo& ti = base->TripsInfoList[i];		//ke bekhaym tebghe zarfiat barash ja be jae lahaz konim
            int Trip_qMax = *max_element(ti.q.begin(),ti.q.end());	//dar safare masalan 0 max q ha masalan cheghadre, masalan shode 20, yani ye
            if(Trip_qMax > qMax)			//ye gere e vojood dashte ke az max zarfiat bozorg tare
            {								//hala ke bishtar mishe
                CheckCapacityCondition(ti,qMax);			//mibinim ke bayad shamele ja be jae bshim
                ApplyCapacityTransmission(ti,qMax);		//hala mire bara ja be jae bekhatere zarfiat
                isMoved = true;
            }

        }

        if(!isMoved)
        {
            isConditionPass = true;
        }
    }

    cerr << "Heuristic Computation Done ... Detecting Solver Params " << endl;

    FindWarmStartParams();			//enja mirim baraye dastan warm start

    for(map<int,int>::const_iterator it = Xvar.begin() ; it!= Xvar.end() ; it++)
    {
        cerr << "X["<<it->first<<"]["<<it->second<<"] = 1" << endl;
    }

    for(map<int,int>::const_iterator it = Zvar.begin() ; it!= Zvar.end() ; it++)
    {
        cerr << "Z["<<it->first<<"]["<<it->second<<"] = 1" << endl;
    }

    fstream graphfile;
    graphfile.open("HeuristicGraph.txt",ios_base::out | ios_base::trunc);

    for(Index i=0 ; i < GraphNodesList.size() ; i++)
    {
        const Coord& _c = GraphNodesList[i];
        if(GraphCPCoords.find(_c) != GraphCPCoords.end())          //file graph ro be ezaye har istgah ye setare
        {												// be ezaye har gheyre istgah ye x
            graphfile << "*" << "\t";
        }
        else
        {
            graphfile << "x" << "\t";
        }
    }
    graphfile.close();
}

void Heuristic::ApplyHeuToCplex()
{
    solver->HandleXZQParams(base,Zvar,Xvar,Qvar);
}

void Heuristic::CleanUnusedArcsMatrix()
{
    for(Index i=0 ; i < solver->ArcsMatrix.size() ; i++)
    {
        for(Index j=0 ; j < solver->ArcsMatrix[i].size() ; j++)
        {
            solver->ArcsMatrix[i][j] = 0;
        }
    }
    for(map<int,int>::const_iterator it = Xvar.begin() ; it != Xvar.end() ; it++)
    {
        int i = it->first;
        int j = it->second;
        solver->ArcsMatrix[i][j] = 1;
    }
}

void Heuristic::LoadTrips()
{
    for(Index i=0 ; i < SolverTrips.size() ; i++)
    {
        base->TripsInfoList[i] = SolverTrips[i];
    }
}

void Heuristic::Compute_q()
{
    vector<int> Nodes;
//    vector <int> RequestsType;
    vector <int> LocalRequestsPickStopNo;
    vector <int> LocalRequestsDropStopNo;
    vector <vector<vector<int>>> PND;
    vector <vector<vector<int>>> NPD;
    int VehicleIndex;
    int TripsNo = base->TripsNo;
    for(int r = 0 ; r < TripsNo ; r++)						//roo kole safar ha
    {
        TripsInfo& ti = base->TripsInfoList[r];
        int TripNodesNo = ti.CPCoords.size() + ti.NonCPCoord.size(); //tedade gere ha
        ti.q.resize(TripNodesNo,0);
        const int& ReqNo = ti.RequestNo;					//tedade darkhasta
        const vector<int>& ReqType = ti.RequestsType;				//noe moshtari chie
        const vector<int>& ReqHeuDS = ti.HeuristicRequestsDropStopNo;		//koja savaro piade mishe
        const vector<int>& ReqHeuPS = ti.HeuristicRequestsPickStopNo;
        int VehicleLoad = 0;
        for(int i=0 ; i < TripNodesNo ; i++)			//halghe rooye tedade gereha chon tak takeshono mikhaym hesab konim
        {
            int lastLoad = VehicleLoad;
            for(int k = 0 ; k < ReqNo ; k++)				//baraye tak take darkhast ha
            {
                if(ReqType[k] == REQUESTTYPE::PND)					//hala age darkhaste ma az noe pnd bashe
                {
                    int localIdxPS = GlobalToLocal(ReqHeuPS[k]);	//koja mikhad savar beshe koja mikhad piade bshe
                    int localIdxDS = GlobalToLocal(ReqHeuDS[k]);	//chon safar be safar darim barresi mikonim miam be soorate local hesab mikonim, chon
                    if(localIdxPS == i) 		//too ye safar az 0 ta 10 e dobare safare badi az 0 ta 7 e masalan ama too gerafe koli masalan 30 ta 40 e
                    {				//bara hamin local mikonimesh.
                        VehicleLoad += Positive;		//age too ye gere ye savar shodan etefagh oftade
                        ti.q[i] = VehicleLoad;			//ye load mosbat behesh ezafe kon
                    }
                    else if(localIdxDS == i)			//age too en gere ye piade shodan etefagh oftade
                    {
                        VehicleLoad += Negative;			//ye load kam kon
                        ti.q[i] = VehicleLoad;			//beriz baraye zarfiate oon noghte.
                    }
                }
                else if(ReqType[k] == REQUESTTYPE::NPD) //hamin dobare baraye npd
                {
                    int localIdxPS = GlobalToLocal(ReqHeuPS[k]);
                    int localIdxDS = GlobalToLocal(ReqHeuDS[k]);
                    if(localIdxPS == i)
                    {
                        VehicleLoad += Positive;
                        ti.q[i] = VehicleLoad;
                    }
                    else if(localIdxDS == i)
                    {
                        VehicleLoad += Negative;
                        ti.q[i] = VehicleLoad;
                    }
                }
                else										//baraye npnd va pd ha
                {
                    int localIdxPS = GlobalToLocal(ReqHeuPS[k]);
                    int localIdxDS = GlobalToLocal(ReqHeuDS[k]);
                    if(localIdxPS == i)
                    {
                        VehicleLoad += Positive;
                        ti.q[i] = VehicleLoad;
                    }
                    else if(localIdxDS == i)
                    {
                        VehicleLoad += Negative;
                        ti.q[i] = VehicleLoad;
                    }
                }
            }
            if(VehicleLoad == lastLoad) // age loadi ham nabod ke akharin load ro bezar ke mamoolan entori nist
            {
                ti.q[i] = VehicleLoad;
            }
        }
    }
}

void Heuristic::Compute_q_Trip(TripsInfo& ti)
{
    int TripNodesNo = ti.CPCoords.size() + ti.NonCPCoord.size(); //tedade gere ha
    ti.q.resize(TripNodesNo,0);
    const int& ReqNo = ti.RequestNo;					//tedade darkhasta
    const vector<int>& ReqType = ti.RequestsType;				//noe moshtari chie
    const vector<int>& ReqHeuDS = ti.HeuristicRequestsDropStopNo;		//koja savaro piade mishe
    const vector<int>& ReqHeuPS = ti.HeuristicRequestsPickStopNo;
    int VehicleLoad = 0;
    for(int i=0 ; i < TripNodesNo ; i++)			//halghe rooye tedade gereha chon tak takeshono mikhaym hesab konim
    {
        int lastLoad = VehicleLoad;
        for(int k = 0 ; k < ReqNo ; k++)				//baraye tak take darkhast ha
        {
            if(ReqType[k] == REQUESTTYPE::PND)					//hala age darkhaste ma az noe pnd bashe
            {
                int localIdxPS = GlobalToLocal(ReqHeuPS[k]);	//koja mikhad savar beshe koja mikhad piade bshe
                int localIdxDS = GlobalToLocal(ReqHeuDS[k]);	//chon safar be safar darim barresi mikonim miam be soorate local hesab mikonim, chon
                if(localIdxPS == i) 		//too ye safar az 0 ta 10 e dobare safare badi az 0 ta 7 e masalan ama too gerafe koli masalan 30 ta 40 e
                {				//bara hamin local mikonimesh.
                    VehicleLoad += Positive;		//age too ye gere ye savar shodan etefagh oftade
                    ti.q[i] = VehicleLoad;			//ye load mosbat behesh ezafe kon
                }
                else if(localIdxDS == i)			//age too en gere ye piade shodan etefagh oftade
                {
                    VehicleLoad += Negative;			//ye load kam kon
                    ti.q[i] = VehicleLoad;			//beriz baraye zarfiate oon noghte.
                }
            }
            else if(ReqType[k] == REQUESTTYPE::NPD) //hamin dobare baraye npd
            {
                int localIdxPS = GlobalToLocal(ReqHeuPS[k]);
                int localIdxDS = GlobalToLocal(ReqHeuDS[k]);
                if(localIdxPS == i)
                {
                    VehicleLoad += Positive;
                    ti.q[i] = VehicleLoad;
                }
                else if(localIdxDS == i)
                {
                    VehicleLoad += Negative;
                    ti.q[i] = VehicleLoad;
                }
            }
            else										//baraye npnd va pd ha
            {
                int localIdxPS = GlobalToLocal(ReqHeuPS[k]);
                int localIdxDS = GlobalToLocal(ReqHeuDS[k]);
                if(localIdxPS == i)
                {
                    VehicleLoad += Positive;
                    ti.q[i] = VehicleLoad;
                }
                else if(localIdxDS == i)
                {
                    VehicleLoad += Negative;
                    ti.q[i] = VehicleLoad;
                }
            }
        }
        if(VehicleLoad == lastLoad) // age loadi ham nabod ke akharin load ro bezar ke mamoolan entori nist
        {
            ti.q[i] = VehicleLoad;
        }
    }
}

int Heuristic::GlobalToLocal(const int &globalIdx)	//global ro miam be local tabdim mikonim
{
    Coord ReqHeuPSPos = GetCoordinationFromIndex(globalIdx);		//en mige bia to be man andis bede, man behet mokhtasat midam
    int TripId_PS = GetTripNumberOfNode(ReqHeuPSPos);			//hala behem bego ke mokhtasati ke dari bara che safarie
    Coord FirstPosTrip = *base->TripsInfoList[TripId_PS].CPCoords.begin();	//safaresho bedast miare, va bego k en safar ba che istgahi shoroo mishe
    int FirstIdxTrip = GetIndexFromList(GraphNodesList,FirstPosTrip);	//miad mige masaln folan istgah ke mokhtasatesh 12,0 e tooye kolle geraf che andisi dare
    int LocalIdx = globalIdx - FirstIdxTrip;				//masalan globalesh mishe gereye 10, voroodi ma ye andis gelobal bood ke midunestim baraye en safare masalan
    return LocalIdx;		//14 bude menhaye 10 baraye avalin istgah mikonim ke entori mishe andis local ro be dast avord
}

void Heuristic::GetNearestCPToCurrentPos(const SolverBase *sb, const Heuristic::Coord &rootPos, vector<int> &cpIndices, Heuristic::CoordList &cpPosition)
{
    UNUSED(sb);
    double minDisback = std::numeric_limits<double>::max();
    double minDisnext = std::numeric_limits<double>::max();
    Coord next ;
    Coord back ;
    for(int i=0 ; i < solver->GraphCPGlobalList.size() ; i++)
    {
        double diffmin = std::fabs(rootPos.first - solver->GraphCPGlobalList[i].first);
        if(rootPos.first > solver->GraphCPGlobalList[i].first)
        {
            if(diffmin < minDisback)
            {
                minDisback = diffmin;
                back = solver->GraphCPGlobalList[i];
            }
        }
        else
        {
            if(diffmin < minDisnext)
            {
                minDisnext = diffmin;
                next = solver->GraphCPGlobalList[i];
            }
        }
    }
    int backIdx = utility->GetIndexFromSet(solver->GraphNodes,back);
    int nextIdx = utility->GetIndexFromSet(solver->GraphNodes,next);

    cpIndices.push_back(backIdx);
    cpIndices.push_back(nextIdx);

    cpPosition.push_back(back);
    cpPosition.push_back(next);

}

void Heuristic::CreateRawGraph()
{
    int TripsNo = base->TripsNo;

    for(int i=0 ; i < TripsNo ; i++)
    {
        TripsInfo ti = base->TripsInfoList[i];
        const CoordSet& CPCoords = ti.CPCoords;
        for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)
        {
            const Coord& cpcoord = *cit;
            RawGraphNodesList.push_back(cpcoord);
            RawGraphNodes.insert(cpcoord);
            RawGraphCPCoords.insert(cpcoord);
        }
    }
}

void Heuristic::SaveSolverTrips()
{
    int TripsNo = base->TripsNo;				//c etellat ro dare migire baraye shomare safarha
    for(int i=0 ; i < TripsNo ; i++)
    {
        TripsInfo ti = base->TripsInfoList[i];
        SolverTrips.push_back(ti);
    }
}

void Heuristic::FindWarmStartParams()
{
    int TripsNo = base->TripsNo;					//kolle safar haro doone doone barresi kon
    for(int i=0 ; i < TripsNo ; i++)
    {
        TripsInfo ti = base->TripsInfoList[i];
        Coord FirstPosOfTrip = *ti.CPCoords.begin();
        for(Index j=0 ; j < ti.PickStopXList.size() ; j++)
        {
            Coord Pick,Drop;
            const double& px = ti.PickStopXList[j];				//darkhast haro doone doone bardar masalan
            const double& py = ti.PickStopYList[j];					//too safar 0
            const double& dx = ti.DropStopXList[j];		//savaro piade shodaneshoono
            const double& dy = ti.DropStopYList[j];
            Pick = make_pair(px,py);											//zoj kon
            Drop = make_pair(dx,dy);
            double diffX_Pick = std::fabs(Pick.first - FirstPosOfTrip.first);	//che noghte e savar o piade shode
            double diffX_Drop = std::fabs(Drop.first - FirstPosOfTrip.first);	//masalan 2 savar mishe 3 piade mishe
            diffX_Pick = doFixedPos(diffX_Pick);	//entori migim ke az mabda kam mishe, masalan mabda 0 e,
            diffX_Drop = doFixedPos(diffX_Drop);	//migim 2vahed bara savar shodan jolo bude 3 vahed bara piade shodan

            for(int r=0 ; r <= ti.id ; r++)       // too kolle safar hae ke ta ghable en darkhast hastan
            {									//masalan en darkhast too safare 3 bude, ma halghe ta safare 3 mizanim
                bool isFound = false;
                const TripsInfo& ti_save = SolverTrips[r]; //hala too gerafe solver asli bia oon safaro begir va istgahe avalesho bedast biar
                Coord FirstPosOfTripOld = *ti_save.CPCoords.begin();
                Coord CheckPos_Pick = make_pair(doFixedPos(FirstPosOfTripOld.first+diffX_Pick),Pick.second); //mokhtasate oon istgaho bia + kon ba 2vahed ekhtelaf too savar
                Coord CheckPos_Drop = make_pair(doFixedPos(FirstPosOfTripOld.first+diffX_Drop),Drop.second);	//va 3 vahed ekhtelaf tooye piade shodan
                for(Index k=0 ; k < ti_save.PickStopXList.size() ; k++) //tamume darkhast haye oon safaro peyda kon
                {									//koja piade shodan koja savar shodan
                    Coord Pick_Old,Drop_Old;
                    const double& px_old = ti_save.PickStopXList[k];	//enja hame darkhasthaye oon safaro peyda mikne
                    const double& py_old = ti_save.PickStopYList[k];	// koja ha savar mishan
                    const double& dx_old = ti_save.DropStopXList[k];	//kojaha piade mishan
                    const double& dy_old = ti_save.DropStopYList[k];
                    Pick_Old = make_pair(px_old,py_old);
                    Drop_Old = make_pair(dx_old,dy_old);
//                    if(CheckPos_Pick == Pick_Old && CheckPos_Drop == Drop_Old)
                    if(isEqual(CheckPos_Pick,Pick_Old) && isEqual(CheckPos_Drop,Drop_Old))	//age oon 2 vahed o 3vahed ro age too safar asli peyda kardish
                    {
                        int sumReq = 0;
                        for(int a=0 ; a < ti_save.id ; a++)
                            sumReq += SolverTrips[a].RequestNo;	//entori mifahmim ke folan darkhast masalan az safare 0 bude oomade be safare 3

                        if(TotalReqType[sumReq+k] != PD && TotalReqType[sumReq+k] != NPND)
                            Zvar[sumReq+k] = ti.id;		//age darkhasta aza noe aval o akhar nabudan, entori mifahmim
                        //en darkhaste 6 masalan male safar 1 bude umade ama toye safare 3 servis gerefte
                        isFound = true;	//zed masalan be ma mige darkhast 6 ma too safere 0 bude
                        break;
                    }
                }
                if(isFound)
                    break;
            }
        }
    }
    for(Index i=0 ; i < GraphNodesList.size()-1 ; i++)	//hala miaym x ro por mikonim,
    {
        const Coord& pos = GraphNodesList[i];				// masalan az 0 be 1, 0,0 va 1,0
        const Coord& pos_next = GraphNodesList[i+1];		//moghe eyato peyda mikoni ya na
        int posIdx = utility->GetIndexFromList(solver->GraphNodesList,pos);
        int pos_nextIdx = utility->GetIndexFromList(solver->GraphNodesList,pos_next);
        if(posIdx == -1)
        {
            int rootTripIdx = GetTripNumberOfNode(pos);			//en moghe eyati ke donbaleshi, male che safari bude,
            const TripsInfo& ti = base->TripsInfoList[rootTripIdx]; //masalan mire en baraye safare 5 bude
            Coord FirstPosOfTrip = *ti.CPCoords.begin();		//mige ok boro safare 5 node avalesho begir
            double diffX = std::fabs(pos.first - FirstPosOfTrip.first);	//va moghe eyate vahedesh ro begir, cheghad rooye oon safar rafte jelo
            diffX = doFixedPos(diffX);			///

            for(int r=0 ; r < TripsNo ; r++)			///kolle safararo barresi kon
            {
                const TripsInfo& ti_save = SolverTrips[r];	//hala shorooe oon safar + oon ekhtelafe koni va position jadid peyda koni
                Coord FirstPosOfTripOld = *ti_save.CPCoords.begin();	//hala begu oono mituni peyda koni too gerafe aslit ya na
                Coord newPos = make_pair(doFixedPos(FirstPosOfTripOld.first + diffX),pos.second); //
                int newPosIdx = utility->GetIndexFromList(solver->GraphNodesList,newPos); //entori label asli ro peyda mikonim tooye gerafe asli
                if(newPosIdx != -1)
                {
                    posIdx = newPosIdx;
                    break;
                }
            }

            if(posIdx == -1)
            {
                cerr << "WTF :: Idx not detected" << endl;
            }
        }
        if(pos_nextIdx == -1)			//enam baraye andis badiye agar peyda nakardi
        {
            int rootTripIdx = GetTripNumberOfNode(pos_next);	//man behet migam k en toye oon safar kojast
            const TripsInfo& ti = base->TripsInfoList[rootTripIdx];
            Coord FirstPosOfTrip = *ti.CPCoords.begin();
            double diffX = std::fabs(pos_next.first - FirstPosOfTrip.first);
            diffX = doFixedPos(diffX);

            for(int r=0 ; r < TripsNo ; r++)
            {
                const TripsInfo& ti_save = SolverTrips[r];
                Coord FirstPosOfTripOld = *ti_save.CPCoords.begin();
                Coord newPos = make_pair(doFixedPos(FirstPosOfTripOld.first + diffX),pos_next.second);
                int newPosIdx = utility->GetIndexFromList(solver->GraphNodesList,newPos);
                if(newPosIdx != -1)
                {
                    pos_nextIdx = newPosIdx;
                    break;
                }
            }

            if(pos_nextIdx == -1)
            {
                cerr << "WTF :: rootTripIdx :: " << rootTripIdx << endl;
                cerr << "WTF :: ti.id :: " << ti.id << endl;
                cerr << "WTF :: diffX :: " << diffX << endl;
                cerr << "WTF :: i :: " << i << endl;
                cerr << "WTF :: pos_next :: " << pos_next.first << "," << pos_next.second << endl;
                cerr << "WTF :: FirstPosOfTrip :: " << FirstPosOfTrip.first << "," << FirstPosOfTrip.second << endl;
                cerr << "WTF :: Idx Next not detected" << endl;
            }
        }
        if(posIdx != -1 && pos_nextIdx != -1)						//
        {
            Xvar[posIdx] = pos_nextIdx;				///masalan az 0 mikhad bere be 1 , entri x dar miad
            int TripIDHeu_Pos = GetTripNumberOfNode(pos);			//az enja be bad mohasebate q e ke be kar nemiad
            int TripIDHeu_PosNext = GetTripNumberOfNode(pos_next);
            int LocalPosIdx = GlobalToLocal(i);
            int LocalPosIdx_Next = GlobalToLocal(i+1);
            const TripsInfo& ti1 = base->TripsInfoList[TripIDHeu_Pos];
            const TripsInfo& ti2 = base->TripsInfoList[TripIDHeu_PosNext];
            Qvar[posIdx] = ti1.q[LocalPosIdx];
            Qvar[pos_nextIdx] = ti2.q[LocalPosIdx_Next];
        }
        else
        {
            cerr << "Neg Index , Out of bound" << endl;
            exit(0);
        }
    }
}

void Heuristic::CheckThetaCondition(TripsInfo &ti)
{
    bool isPass = false;
    double StopB = base->Input_StopB;
    while(!isPass)							//dobare ye halgheye binahayat mizanim
    {
        isPass = true;
        const int& ReqNo = ti.RequestNo;
        const vector<int>& ReqType = ti.RequestsType;
        const vector<int>& ReqHeuDS = ti.HeuristicRequestsDropStopNo;
        const vector<int>& ReqHeuPS = ti.HeuristicRequestsPickStopNo;
        for(int k=0 ; k < ReqNo ; k++)
        {
            if(ReqType[k] == PND)
            {
                int Drop = ReqHeuDS[k];				//koja savaro piade mishe
                int Pick = ReqHeuPS[k];				//oon tabe ma be ma har position gheyre istgahi ke bedim, behemoon mige istgahe badish chie
                int NextCP = GetNextCPIndexFromRoot(Drop);	//bia piade shodanesho bbin va bego istgahe badiye en gereye piade shodan chie
                if(HiddenStopTeta[NextCP] > (HiddenTeta[NextCP] - GAP_HEU)) //yani darkhast ghabeliate pasokhgooe nadare
                {
                    cerr << "Theta :: PND Request["<<k<<"] Rejected From Trip["<<ti.id<<"]" << endl;

                    SaveHybrid(ti,k,Pick,Drop,PND);			//en darkhast dar en sharaye reject mishe o zakhire mishe
                    UpdateParameterOnRemove(ti,k,Drop);		//hala oon darkhast bayad remove beshe

                    isPass = false;			//enja migim ok ma ye darkhast hazf kardim dobare en masiro az no bayad beri
                    			//dobare baayad baresi beshe bebinim kasi potansiel hazf shodan ro dare ya na
                    break;
                }
            }
            if(ReqType[k] == NPD)    //daghighan mesle bala eno pish mibarim
            {
                int Drop = ReqHeuDS[k];
                int Pick = ReqHeuPS[k];
                int NextCP = GetNextCPIndexFromRoot(Pick);
                if(HiddenStopTeta[NextCP] > (HiddenTeta[NextCP] - GAP_HEU))
                {
                    cerr << "Theta :: NPD Request["<<k<<"] Rejected From Trip["<<ti.id<<"]" << endl;

                    SaveHybrid(ti,k,Pick,Drop,NPD);
                    UpdateParameterOnRemove(ti,k,Pick);

                    isPass = false;

                    break;
                }
            }
            if(ReqType[k] == NPND)			//har 3ta savaro piade shodan ro lahaz kon, chon agar haro
            {									//ezafe ya kam beshan niaze ke barresi beshan
                int Drop = ReqHeuDS[k];
                int Pick = ReqHeuPS[k];
                int NextCP_Pick = GetNextCPIndexFromRoot(Pick);
                int NextCP_Drop = GetNextCPIndexFromRoot(Drop);
                if(HiddenStopTeta[NextCP_Pick] > (HiddenTeta[NextCP_Pick] - GAP_HEU)
                        || HiddenStopTeta[NextCP_Drop] > (HiddenTeta[NextCP_Drop] - GAP_HEU))
                {
                    cerr << "Theta :: NPND Request["<<k<<"] Rejected From Trip["<<ti.id<<"]" << endl;

                    SaveHybrid(ti,k,Pick,Drop,NPND);
                    Coord PickPos = GetCoordinationFromIndex(Pick);
                    Coord DropPos = GetCoordinationFromIndex(Drop);
                    UpdateParameterOnRemoveNPND(ti,k,PickPos,DropPos); //enja ham etelaat en moshtari ro hazf mikonim

                    isPass = false;

                    break;
                }
            }
        }
    }
}

void Heuristic::CheckCapacityCondition(TripsInfo &ti, const int &qMax)
{
    while(1)				//enghadr sharayeto barresi mikone ke pass beshe
    {
        bool isOverloadFound = false;
        int revertIdx = -1000;
        int VehicleIndex = ti.VehicleIndex;			//
        const int& ReqNo = ti.RequestNo;		//che tedad darkhast bude
        const vector<int>& ReqType = ti.RequestsType;		//noesh chi bude
        const vector<int>& ReqHeuDS = ti.HeuristicRequestsDropStopNo;	//koja savaro piade mishodan
        const vector<int>& ReqHeuPS = ti.HeuristicRequestsPickStopNo;
        for(int j=0 ; j < ti.q.size() ; j++)			// az 0 ta kolle gere haye oon safar boro jolo
        {
            if(ti.q[j] > qMax)						//agar folan gere az zarfiat bishtar bashe, niaz be ja be jae has
            {
                int candidateReq = -1;					//
                isOverloadFound = true;
                double maxTau = -1000000.0;
                for(Index k=0 ; k < ti.RequestNo ; k++)			//roo hame darkhasta
                {
                    if(ti.RequestsType[k] == static_cast<int>(PND))			// masalan darkhast pnd
                    {
                        int Drop = ReqHeuDS[k];
                        int Pick = ReqHeuPS[k];
                        int localPick = GlobalToLocal(Pick);
                        int localDrop = GlobalToLocal(Drop);

                        if(localPick == j)			//age savar shodanesh mosavi oon gere e shode k bishtar az qmax has
                        {
                            if(ti.RequestsTau[k] > maxTau) // hala oon darkhasti ke tau bishtari dare entekhab kon
                            {									//baraye ja be jae.
                                maxTau = ti.RequestsTau[k];
                                candidateReq = k;				//entori moshakhas mishe k ki bayad ja be ja beshe
                            }
                        }
                    }
                    else if(ti.RequestsType[k] == static_cast<int>(NPD))
                    {
                        int Drop = ReqHeuDS[k];
                        int Pick = ReqHeuPS[k];
                        int localPick = GlobalToLocal(Pick);
                        int localDrop = GlobalToLocal(Drop);

                        if(localPick == j)
                        {
                            if(ti.RequestsTau[k] > maxTau)
                            {
                                maxTau = ti.RequestsTau[k];
                                candidateReq = k;
                            }
                        }
                    }
                    else if(ti.RequestsType[k] == static_cast<int>(NPND))
                    {
                        int Drop = ReqHeuDS[k];
                        int Pick = ReqHeuPS[k];
                        int localPick = GlobalToLocal(Pick);
                        int localDrop = GlobalToLocal(Drop);

                        if(localPick == j)
                        {
                            if(ti.RequestsTau[k] > maxTau)
                            {
                                maxTau = ti.RequestsTau[k];
                                candidateReq = k;
                            }
                        }
                    }
                }

                if(candidateReq == -1)
                {
                    HandlePDTrap(ti,j);
                    break;
                }
                int i = candidateReq;

                bool isHybridFound = false;

                if(ti.RequestsType[i] == static_cast<int>(PND))
                {
                    int Drop = ReqHeuDS[i];
                    int Pick = ReqHeuPS[i];
                    int localPick = GlobalToLocal(Pick);
                    int localDrop = GlobalToLocal(Drop);

                    if(localPick == j)		//bia bego enja ye darkhast peyda kardim ke niaz be ja be jae dare
                    {

                        cerr << "Capacity :: PND Request["<<i<<"] Rejected From Trip["<<ti.id<<"]" << endl;
                        Coord PickCoord = GetCoordinationFromIndex(Pick);	///savaro piade shodano dar miarim
                        Coord DropCoord = GetCoordinationFromIndex(Drop);

                        SaveCpacityHybrid(ti,i,Pick,Drop,PND);			//zakhire mikonim etelaatesho
                        UpdateParameterOnRemove(ti,i,Drop);			//hala darkhast ja be ja shod, paramet ha be rooz beshe
                        isHybridFound = true;
                        Compute_q_Trip(ti);			//hala dobare zarfiat haro hesab kon

                        break;
                    }
                }
                else if(ti.RequestsType[i] == static_cast<int>(NPD))
                {
                    int Drop = ReqHeuDS[i];
                    int Pick = ReqHeuPS[i];
                    int localPick = GlobalToLocal(Pick);
                    int localDrop = GlobalToLocal(Drop);

                    if(localPick == j)
                    {
                        cerr << "Capacity :: NPD Request["<<i<<"] Rejected From Trip["<<ti.id<<"]" << endl;
                        Coord PickCoord = GetCoordinationFromIndex(Pick);
                        Coord DropCoord = GetCoordinationFromIndex(Drop);

                        SaveCpacityHybrid(ti,i,Pick,Drop,NPD);
                        UpdateParameterOnRemove(ti,i,Pick);
                        isHybridFound = true;
                        Compute_q_Trip(ti);

                        break;
                    }
                }
                else if(ti.RequestsType[i] == static_cast<int>(NPND))
                {
                    int Drop = ReqHeuDS[i];
                    int Pick = ReqHeuPS[i];
                    int localPick = GlobalToLocal(Pick);
                    int localDrop = GlobalToLocal(Drop);

                    if(localPick == j)
                    {
                        cerr << "Capacity :: NPND Request["<<i<<"] Rejected From Trip["<<ti.id<<"]" << endl;
                        Coord PickCoord = GetCoordinationFromIndex(Pick);
                        Coord DropCoord = GetCoordinationFromIndex(Drop);

                        SaveCpacityHybrid(ti,i,Pick,Drop,NPND);
                        UpdateParameterOnRemoveNPND(ti,i,PickCoord,DropCoord);
                        isHybridFound = true;
                        Compute_q_Trip(ti);

                        break;
                    }
                }

                if(isHybridFound)
                    break;
            }
        }
        if(!isOverloadFound)
            break;
    }
}

void Heuristic::HandlePDTrap(TripsInfo& ti, const int trapIdx)
{
    set<int> visitedIndices;
    bool isPass = false;
    while(!isPass)
    {
        int VehicleIndex = ti.VehicleIndex;			//
        const int& ReqNo = ti.RequestNo;		//che tedad darkhast bude
        const vector<int>& ReqType = ti.RequestsType;		//noesh chi bude
        const vector<int>& ReqHeuDS = ti.HeuristicRequestsDropStopNo;	//koja savaro piade mishodan
        const vector<int>& ReqHeuPS = ti.HeuristicRequestsPickStopNo;

        int maxq = std::numeric_limits<int>::min();
        int detectedIdx = -1;


        for(int j=0 ; j < trapIdx ; j++)
        {
            if(ti.q[j] > maxq && visitedIndices.find(j) == visitedIndices.end())
            {
                maxq = ti.q[j];
                detectedIdx = j;
            }
        }
        if(detectedIdx == -1)
        {
            cerr << __FUNCTION__ << ": index selection failed ... PD request trap in node [" << trapIdx <<"]" << endl;
            return;
        }
        visitedIndices.insert(detectedIdx);
        double maxTau = -1000000.0;
        int candidateReq = -1;
        for(Index k=0 ; k < ti.RequestNo ; k++)			//roo hame darkhasta
        {
            if(ti.RequestsType[k] == static_cast<int>(PND))			// masalan darkhast pnd
            {
                int Drop = ReqHeuDS[k];
                int Pick = ReqHeuPS[k];
                int localPick = GlobalToLocal(Pick);
                int localDrop = GlobalToLocal(Drop);

                if(localPick == detectedIdx)			//age savar shodanesh mosavi oon gere e shode k bishtar az qmax has
                {
                    if(ti.RequestsTau[k] > maxTau) // hala oon darkhasti ke tau bishtari dare entekhab kon
                    {									//baraye ja be jae.
                        maxTau = ti.RequestsTau[k];
                        candidateReq = k;				//entori moshakhas mishe k ki bayad ja be ja beshe
                    }
                }
            }
            else if(ti.RequestsType[k] == static_cast<int>(NPD))
            {
                int Drop = ReqHeuDS[k];
                int Pick = ReqHeuPS[k];
                int localPick = GlobalToLocal(Pick);
                int localDrop = GlobalToLocal(Drop);

                if(localPick == detectedIdx)
                {
                    if(ti.RequestsTau[k] > maxTau)
                    {
                        maxTau = ti.RequestsTau[k];
                        candidateReq = k;
                    }
                }
            }
            else if(ti.RequestsType[k] == static_cast<int>(NPND))
            {
                int Drop = ReqHeuDS[k];
                int Pick = ReqHeuPS[k];
                int localPick = GlobalToLocal(Pick);
                int localDrop = GlobalToLocal(Drop);

                if(localPick == detectedIdx)
                {
                    if(ti.RequestsTau[k] > maxTau)
                    {
                        maxTau = ti.RequestsTau[k];
                        candidateReq = k;
                    }
                }
            }
        }
        if(candidateReq != -1)
        {
            int i = candidateReq;
            bool isHybridFound = false;

            if(ti.RequestsType[i] == static_cast<int>(PND))
            {
                int Drop = ReqHeuDS[i];
                int Pick = ReqHeuPS[i];
                int localPick = GlobalToLocal(Pick);
                int localDrop = GlobalToLocal(Drop);

                if(localPick == detectedIdx)		//bia bego enja ye darkhast peyda kardim ke niaz be ja be jae dare
                {

                    cerr << "Capacity :: PND Request["<<i<<"] Rejected From Trip["<<ti.id<<"]" << endl;
                    Coord PickCoord = GetCoordinationFromIndex(Pick);	///savaro piade shodano dar miarim
                    Coord DropCoord = GetCoordinationFromIndex(Drop);

                    SaveCpacityHybrid(ti,i,Pick,Drop,PND);			//zakhire mikonim etelaatesho
                    UpdateParameterOnRemove(ti,i,Drop);			//hala darkhast ja be ja shod, paramet ha be rooz beshe
                    isHybridFound = true;
                    Compute_q_Trip(ti);			//hala dobare zarfiat haro hesab kon

                    break;
                }
            }
            else if(ti.RequestsType[i] == static_cast<int>(NPD))
            {
                int Drop = ReqHeuDS[i];
                int Pick = ReqHeuPS[i];
                int localPick = GlobalToLocal(Pick);
                int localDrop = GlobalToLocal(Drop);

                if(localPick == detectedIdx)
                {
                    cerr << "Capacity :: NPD Request["<<i<<"] Rejected From Trip["<<ti.id<<"]" << endl;
                    Coord PickCoord = GetCoordinationFromIndex(Pick);
                    Coord DropCoord = GetCoordinationFromIndex(Drop);

                    SaveCpacityHybrid(ti,i,Pick,Drop,NPD);
                    UpdateParameterOnRemove(ti,i,Pick);
                    isHybridFound = true;
                    Compute_q_Trip(ti);

                    break;
                }
            }
            else if(ti.RequestsType[i] == static_cast<int>(NPND))
            {
                int Drop = ReqHeuDS[i];
                int Pick = ReqHeuPS[i];
                int localPick = GlobalToLocal(Pick);
                int localDrop = GlobalToLocal(Drop);

                if(localPick == detectedIdx)
                {
                    cerr << "Capacity :: NPND Request["<<i<<"] Rejected From Trip["<<ti.id<<"]" << endl;
                    Coord PickCoord = GetCoordinationFromIndex(Pick);
                    Coord DropCoord = GetCoordinationFromIndex(Drop);

                    SaveCpacityHybrid(ti,i,Pick,Drop,NPND);
                    UpdateParameterOnRemoveNPND(ti,i,PickCoord,DropCoord);
                    isHybridFound = true;
                    Compute_q_Trip(ti);

                    break;
                }
            }
            if(isHybridFound)
                isPass = true;
        }
    }
}

void Heuristic::UpdateParameterOnRemove(TripsInfo& ti, const int &ReqIdx, const int &RemovedIdx)
{
    ti.RequestNo --;					//yedoone az tedade darkhastha kam kon chon bayad update konim

    ti.RequestsTau.erase(ti.RequestsTau.begin()+ReqIdx);	//che noe darkhasti bod, tau esh ro bardar
    ti.RequestsType.erase(ti.RequestsType.begin()+ReqIdx);	// noesh ro ham bardar

    ti.PickStopXList.erase(ti.PickStopXList.begin()+ReqIdx);		//jae ke savar o piade mishe ro ham bardar
    ti.PickStopYList.erase(ti.PickStopYList.begin()+ReqIdx);
    ti.DropStopXList.erase(ti.DropStopXList.begin()+ReqIdx);
    ti.DropStopYList.erase(ti.DropStopYList.begin()+ReqIdx);

    Coord RemovedCoord = GetCoordinationFromIndex(RemovedIdx);		//andise gheyr istgahi ke az gerafe
    bool isRemove = RemoveUnusedNodes(ti,RemovedCoord);			//ma mikhast hazf beshe be ma bede
    if(isRemove)												//en hamishe true khahad bud chon
    {						//hamishe tooye har 1doone gheyre istgahi hamishe ye darkhast hast
        GraphNodesList.erase(GraphNodesList.begin()+RemovedIdx);	//oon gere ro hazf mikonim chon dg darkhasi nist barash
        ti.NonCPCoord.erase(RemovedCoord);			//oon ghereye gheyre istgahi ro ham az safar hazf mikonim
    }

    RefreshRequests(ti);					//kolle darkhast haro ye bar refresh miknoim baraye oon safar
    RefreshCheckpointsIndices();			//istgah ha mian refresh mishan
    ComputeGlobalTheta();						//teta dobare hesab mishe - hala ke ye darkhast pak shode
}											//begoo hala ba tavajoh be en dobare hidden teta o enaro mohasebe kon
				//masalan be 3ta darkhast teta pass nmishod, sevomi ro pak kardim hala dobare mirim barresi mikonim ke ba 2 ta okeye ya na
void Heuristic::UpdateParameterOnRemoveNPND(TripsInfo &ti, const int &ReqIdx, const Coord &RemovedPosPick, const Coord &RemovedPosDrop)
{
    ti.RequestNo --;

    ti.RequestsTau.erase(ti.RequestsTau.begin()+ReqIdx);
    ti.RequestsType.erase(ti.RequestsType.begin()+ReqIdx);

    ti.PickStopXList.erase(ti.PickStopXList.begin()+ReqIdx);
    ti.PickStopYList.erase(ti.PickStopYList.begin()+ReqIdx);
    ti.DropStopXList.erase(ti.DropStopXList.begin()+ReqIdx);
    ti.DropStopYList.erase(ti.DropStopYList.begin()+ReqIdx);

    bool isRemovePick = RemoveUnusedNodes(ti,RemovedPosPick); //ham andis baraye savar shodan ro barresi mikonim
    if(isRemovePick)
    {
        int RemoveIdx = GetIndexFromList(GraphNodesList,RemovedPosPick);
        GraphNodesList.erase(GraphNodesList.begin()+RemoveIdx);
        ti.NonCPCoord.erase(RemovedPosPick);
    }
    bool isRemoveDrop = RemoveUnusedNodes(ti,RemovedPosDrop);	//ham andis baraye piade shodan ro barresi mikonim
    if(isRemoveDrop)
    {
        int RemoveIdx = GetIndexFromList(GraphNodesList,RemovedPosDrop);	//miam az rooye gerafemoon bar midarim
        GraphNodesList.erase(GraphNodesList.begin()+RemoveIdx);		//az rooye gheyre istgahi ha ham barmidarim
        ti.NonCPCoord.erase(RemovedPosDrop);
    }

    RefreshRequests(ti);				//dobare ye berooz resani rooye darkhast ha a baghie chiza mizanim
    RefreshCheckpointsIndices();
    ComputeGlobalTheta();
}

void Heuristic::SaveHybrid(TripsInfo &ti,const int& ReqIdx, const int &Pick, const int &Drop,const REQUESTTYPE& type)
{
    Coord PickCoord = GetCoordinationFromIndex(Pick); //koja savaro piade mishode
    Coord DropCoord = GetCoordinationFromIndex(Drop);
    const vector <double>& ReqTau = ti.RequestsTau;

    Hybrid h;
    h.Pick = Pick;
    h.Drop = Drop;		//moghe eyateshoon mesle mokhtasat
    h.StartPoint = PickCoord;		//andis hashoone
    h.EndPoint = DropCoord;
    h.currentTrip = ti.id;			//too che safari bood
    h.nextTrip = ti.id+1;		//too che safari mikhad bere
    h.R_Type = static_cast<int>(type);		//noe darkhast chi bud
    h.R_Tau = ReqTau[ReqIdx];		//tau esh ro ham zakhire kon

    if(h.nextTrip >= base->TripsNo)				//agar safare badi e ke mikhad bere az kolle safar
    {									//hae k mikhaym behesh assign konim bishtar bashe mige khate paeno mide
        cerr << "Theta Condition :: You reach the maximum of trips, increase your trips number" << endl;
    }

    m_HybridTransferList.push_back(h);				//agar moshkeli nabashe, oon darkhast tooye ye  list zakhire mishe
}

void Heuristic::SaveCpacityHybrid(TripsInfo &ti, const int &ReqIdx, const int &Pick, const int &Drop, const REQUESTTYPE &type)
{
    Coord PickCoord = GetCoordinationFromIndex(Pick);
    Coord DropCoord = GetCoordinationFromIndex(Drop);
    const vector <double>& ReqTau = ti.RequestsTau;

    Hybrid h;
    h.Pick = Pick;
    h.Drop = Drop;											//etellaat moshtari ro dar miarim
    h.StartPoint = PickCoord;
    h.EndPoint = DropCoord;
    h.currentTrip = ti.id;
    h.nextTrip = ti.id+1;
    h.R_Type = static_cast<int>(type);
    h.R_Tau = ReqTau[ReqIdx];

    if(h.nextTrip >= base->TripsNo)				//age safar badi bishtar az kolle safara bud,
    {
        cerr << "Capacity Condition :: You reach the maximum of trips, increase your trips number" << endl;
    }

    m_CapacityTransferList.push_back(h);
}

Heuristic::Coord Heuristic::ComputeRootPosInTrip(const Heuristic::Coord &rootPos, const int &lastTripIdx, const int &nextTripIdx)
{
    Coord res;
    Coord firstPosOfTrip = *base->TripsInfoList[lastTripIdx].CPCoords.begin();
    double diffX = std::fabs(rootPos.first - firstPosOfTrip.first);

    Coord firstPosOfCurrentTrip = *base->TripsInfoList[nextTripIdx].CPCoords.begin();
    Coord newPos;
    newPos.first = firstPosOfCurrentTrip.first + diffX;
    newPos.second = rootPos.second;
    res = (newPos);

    return res;
}

void Heuristic::RefreshRequests(TripsInfo &ti)
{
    ti.HeuristicRequestsPickStopNo.resize(ti.RequestNo);			//savaro piade shodan ro be tedade darkhastha
    ti.HeuristicRequestsDropStopNo.resize(ti.RequestNo);			//ye bar dige resize kon

    const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
    const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination
    const vector<double>& PickStopXList = ti.PickStopXList;
    const vector<double>& PickStopYList = ti.PickStopYList;	//savaro piade shodan haye oon safaro dobare lahaz kon
    const vector<double>& DropStopXList = ti.DropStopXList;
    const vector<double>& DropStopYList = ti.DropStopYList;

    CoordSet LocalGraphNodes;
    for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)
    {
        LocalGraphNodes.insert(*cit);
    }
    for(CoordSet::const_iterator cit = NonCPCoord.begin() ; cit != NonCPCoord.end() ; cit++)
    {
        LocalGraphNodes.insert(*cit);
    }
    for(size_t i=0 ; i < PickStopXList.size() ; i++)
    {
        const double& px = PickStopXList[i];
        const double& py = PickStopYList[i];
        const double& dx = DropStopXList[i];
        const double& dy = DropStopYList[i];
        Coord coord_p = make_pair(px,py);
        Coord coord_d = make_pair(dx,dy);

        int idxP = GetIndexFromList(GraphNodesList,coord_p);
        int idxD = GetIndexFromList(GraphNodesList,coord_d);

        if(ti.RequestsType[i] == static_cast<int>(PND))		//ba tavajoh be noe darkhast ha
        {													//savaro piade shodanashoono ok kn
            ti.HeuristicRequestsPickStopNo[i] = idxP;
            ti.HeuristicRequestsDropStopNo[i] = idxD;		//en karo lazeme bokonim chon ye gere pak shode
        }														//hamechiz bayad az no bashe
        else if(ti.RequestsType[i] == static_cast<int>(NPD))
        {
            ti.HeuristicRequestsPickStopNo[i] = idxP;
            ti.HeuristicRequestsDropStopNo[i] = idxD;
        }
        else
        {
            ti.HeuristicRequestsPickStopNo[i] = idxP;
            ti.HeuristicRequestsDropStopNo[i] = idxD;
        }
    }
}

void Heuristic::RefreshRequestsTotal()
{
    int TripsNo = base->TripsNo;
    for(int i=0 ; i < TripsNo ; i++)
    {
        TripsInfo& ti = base->TripsInfoList[i];
        RefreshRequests(ti);
    }
}

void Heuristic::RefreshCheckpointsIndices() //agar tooye hazf kardan ha az graph, yekseri node ha momkene
{												// hazf beshan agar ye nafar hazf beshe graph update mishe entori andis hash
    CPNodes.clear();
    GraphCPCoords.clear();
    int TripsNo = base->TripsInfoList.size();				//ye halghe roo hame safar ha
    for(int r=0 ; r < TripsNo ; r++)
    {
        const TripsInfo& ti = base->TripsInfoList[r];
        const CoordSet& CPCoords = ti.CPCoords;				//istgah haye har safaro gerefte
        for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++) // ye halghe baraye istgah ha
        {
            const Coord& cpcoord = *cit;
            int idxcp = GetIndexFromList(GraphNodesList,cpcoord);   // ye tabee getindex for farakhani mikonim  ba parameter
            CPNodes.insert(idxcp);								//andise oon istgah tooye kolle gerafo daravorde
            GraphCPCoords.insert(cpcoord);						//enja ham tooye graphcp ha zakhire mikonim ke badan
        }														//estefade konim azash
    }
}

void Heuristic::ApplyThetaTransmission(TripsInfo &ti)
{
    UNUSED(ti);
    vector<TripsInfo>& TripsList = base->TripsInfoList;
    int TripsNo = TripsList.size();
    for(Index i=0 ; i < m_HybridTransferList.size() ; i++)	//en list ro dashtim az kasae ke reject karidm.
    {
        const Hybrid& h = m_HybridTransferList[i];		//doone doone darkhast haro barresi kon
        int currentTripId = h.currentTrip;
        int nextTripId = h.nextTrip;
        int ReqType = h.R_Type;
        double ReqTau = h.R_Tau;
        Coord PickCoord = h.StartPoint;
        Coord DropCoord = h.EndPoint;			//ta enja etelaatesh ro dar miarim
        if(nextTripId < TripsNo)				//hala agar safari ke mikhad bere kamtar az safar haye reserve bod
        {
            TripsInfo& ti_next = base->TripsInfoList[nextTripId];			//shoroo kon be kare enteghal
            										//bego too che safari gharae bere (safar badi ro begu chie)
            if(ReqType == static_cast<int>(PND))			//noe moshtari chie
            {
                Coord NewDropPos = ComputeRootPosInTrip(DropCoord,currentTripId,nextTripId);//begu too safare jadid position esh chi mishe
                Coord NewPickPos = ComputeRootPosInTrip(PickCoord,currentTripId,nextTripId); //yani savar o piadash dar miare

                NewDropPos.first = doFixedPos(NewDropPos.first);		//khataye double ro az beyn mibarim
                NewDropPos.second = doFixedPos(NewDropPos.second);
                NewPickPos.first = doFixedPos(NewPickPos.first);
                NewPickPos.second = doFixedPos(NewPickPos.second);

                ti_next.RequestNo++;									//safare jadid ye darkhaste jadid dari
                ti_next.PickStopXList.push_back(NewPickPos.first);		//ena ham position esh hast
                ti_next.PickStopYList.push_back(NewPickPos.second);
                ti_next.DropStopXList.push_back(NewDropPos.first);
                ti_next.DropStopYList.push_back(NewDropPos.second);
                ti_next.NonCPCoord.insert(NewDropPos);				//chon pnd darim ye gheyre istgahi behet ezafe mishe k baraye piade shodane
                ti_next.RequestsTau.push_back(ReqTau);			// tau o type ham ene
                ti_next.RequestsType.push_back(ReqType);
                cerr << "Theta :: Added Request (PND) to Trip : " << ti_next.id << endl;//folan darkhast be folan safar ezafe shode
            }
            if(ReqType == static_cast<int>(NPD)) //hamin bara NPD va npnd ham etefagh miyofte
            {
                Coord NewDropPos = ComputeRootPosInTrip(DropCoord,currentTripId,nextTripId);
                Coord NewPickPos = ComputeRootPosInTrip(PickCoord,currentTripId,nextTripId);

                NewDropPos.first = doFixedPos(NewDropPos.first);
                NewDropPos.second = doFixedPos(NewDropPos.second);
                NewPickPos.first = doFixedPos(NewPickPos.first);
                NewPickPos.second = doFixedPos(NewPickPos.second);

                ti_next.RequestNo++;
                ti_next.PickStopXList.push_back(NewPickPos.first);
                ti_next.PickStopYList.push_back(NewPickPos.second);
                ti_next.DropStopXList.push_back(NewDropPos.first);
                ti_next.DropStopYList.push_back(NewDropPos.second);
                ti_next.NonCPCoord.insert(NewPickPos);
                ti_next.RequestsTau.push_back(ReqTau);
                ti_next.RequestsType.push_back(ReqType);
                cerr << "Theta :: Added Request (NPD) to Trip : " << ti_next.id << endl;
            }
            if(ReqType == static_cast<int>(NPND))
            {
                Coord NewDropPos = ComputeRootPosInTrip(DropCoord,currentTripId,nextTripId);
                Coord NewPickPos = ComputeRootPosInTrip(PickCoord,currentTripId,nextTripId);

                NewDropPos.first = doFixedPos(NewDropPos.first);
                NewDropPos.second = doFixedPos(NewDropPos.second);
                NewPickPos.first = doFixedPos(NewPickPos.first);
                NewPickPos.second = doFixedPos(NewPickPos.second);

                ti_next.RequestNo++;
                ti_next.PickStopXList.push_back(NewPickPos.first);
                ti_next.PickStopYList.push_back(NewPickPos.second);
                ti_next.DropStopXList.push_back(NewDropPos.first);
                ti_next.DropStopYList.push_back(NewDropPos.second);
                ti_next.NonCPCoord.insert(NewDropPos);		//enja 2ta gheyre istgahi behesh ezafe mishe
                ti_next.NonCPCoord.insert(NewPickPos);
                ti_next.RequestsTau.push_back(ReqTau);
                ti_next.RequestsType.push_back(ReqType);
                cerr << "Theta :: Added Request (NPND) to Trip : " << ti_next.id << endl;
            }
        }										//dobare ja be ja mishe be badi
    }
    UpdateGraph();						//hamechizo be rooz resani mikonim
    RefreshRequestsTotal();
    RefreshCheckpointsIndices();		//ye gheyre istgahi chon ezafe shode andis istgah ha be rooz mishe
    ComputeGlobalTheta();				//ke dobare ham bayad enja teta check beshe va agar baz pass nashe
    m_HybridTransferList.clear();				//enteghal ha ke anjam shod kolle listo pak miknoim
}

void Heuristic::ApplyCapacityTransmission(TripsInfo &ti, const int &qMax)
{
    UNUSED(qMax);
    int VehicleIndex = ti.VehicleIndex;
    vector<TripsInfo>& TripsList = base->TripsInfoList;
    int TripsNo = TripsList.size();

    for(Index i=0 ; i < m_CapacityTransferList.size() ; i++)			//che tedad darkhast por shodan too oon list
    {
        Hybrid &h = m_CapacityTransferList[i];
        int currentTripId = h.currentTrip;
        int nextTripId = h.nextTrip;					//etelaatesh ro biad
        int ReqType = h.R_Type;
        double ReqTau = h.R_Tau;
        Coord PickCoord = h.StartPoint;
        Coord DropCoord = h.EndPoint;

        if(h.nextTrip < TripsNo)										//
        {
            TripsInfo& nextTrip = TripsList[h.nextTrip];
            TripsInfo& currentTrip = TripsList[h.currentTrip];

            Coord NewDropPos = ComputeRootPosInTrip(DropCoord,currentTripId,nextTripId);
            Coord NewPickPos = ComputeRootPosInTrip(PickCoord,currentTripId,nextTripId);

            NewDropPos.first = doFixedPos(NewDropPos.first);			//position ro moshakhas kon too safar
            NewDropPos.second = doFixedPos(NewDropPos.second);		// ke mikhay beri
            NewPickPos.first = doFixedPos(NewPickPos.first);
            NewPickPos.second = doFixedPos(NewPickPos.second);

            nextTrip.RequestNo++;										//ye darkhast ezafe kon
            nextTrip.PickStopXList.push_back(NewPickPos.first);		//noghate moadel ro ejad kon
            nextTrip.PickStopYList.push_back(NewPickPos.second);
            nextTrip.DropStopXList.push_back(NewDropPos.first);
            nextTrip.DropStopYList.push_back(NewDropPos.second);

            if(h.R_Type == PND)
            {
                cerr << "Capacity :: Added Request (PND) to Trip : " << nextTrip.id << endl;
                nextTrip.NonCPCoord.insert(NewDropPos);
                nextTrip.RequestsTau.push_back(ReqTau);
                nextTrip.RequestsType.push_back(ReqType);
                nextTrip.q.push_back(0);						//ye zarfiat ham baraye gere jadid dar nazar begir

            }
            else if(h.R_Type == NPD)
            {
                cerr << "Capacity :: Added Request (NPD) to Trip : " << nextTrip.id << endl;
                nextTrip.NonCPCoord.insert(NewPickPos);
                nextTrip.RequestsTau.push_back(ReqTau);
                nextTrip.RequestsType.push_back(ReqType);
                nextTrip.q.push_back(0);
            }
            else if(h.R_Type == NPND)
            {
                cerr << "Capacity :: Added Request (NPND) to Trip : " << nextTrip.id << endl;
                nextTrip.NonCPCoord.insert(NewDropPos);
                nextTrip.NonCPCoord.insert(NewPickPos);
                nextTrip.RequestsTau.push_back(ReqTau);
                nextTrip.RequestsType.push_back(ReqType);
                nextTrip.q.push_back(0);
                nextTrip.q.push_back(0);
            }
        }
    }
    UpdateGraph();
    RefreshCheckpointsIndices();
    ComputeGlobalTheta();
    RefreshRequestsTotal();
    Compute_q();
    m_CapacityTransferList.clear();
}

double Heuristic::doFixedPos(const double &p)
{
    double a = /*std::ceil*/(p*100000.00000f)/100000.00000f;
    stringstream tmp;
    tmp << setprecision(5) << fixed << a;
    double new_val = stod(tmp.str());
    return new_val;
}

bool Heuristic::RemoveUnusedNodes(TripsInfo& ti,const Coord& rootPos)
{
    for(Index i=0 ; i < ti.PickStopXList.size() ; i++)
    {
        Coord p;
        Coord d;
        p.first = ti.PickStopXList[i];
        p.second = ti.PickStopYList[i];
        d.first = ti.DropStopXList[i];
        d.second = ti.DropStopYList[i];

        if(doFixedPos(p.first) == doFixedPos(rootPos.first) ||
                doFixedPos(d.first) == doFixedPos(rootPos.first))
        {
            return false;
        }
    }
    return true;
//    for(Index i=0 ; i < GraphNodesList.size() ; i++)
//    {
//        const Coord& pos = GraphNodesList[i];
//        int TripsNo = base->TripsNo;
//        bool isEmpty = true;
//        for(int j=0 ; j <TripsNo ; j++)
//        {
//            const TripsInfo& ti = base->TripsInfoList[j];
//            Coord p;
//            Coord d;
//            p.first = ti.PickStopXList[i];
//            p.second = ti.PickStopYList[i];
//            d.first = ti.DropStopXList[i];
//            d.second = ti.DropStopYList[i];

//            if(p == pos || d == pos)
//            {
//                isEmpty = false;
//                break;
//            }
//        }
//        if(isEmpty)
//        {
//            GraphNodesList.erase(GraphNodesList.begin()+i);
//        }
//    }
}

void Heuristic::UpdateGraph()			//geraf berooz resani mishe
{							//chon gheyre istgahi behesh ezafe shode, gerafi ke locale.
    GraphNodesList.clear();
    int TripsNo = base->TripsNo;
    CoordSet LocalGraphNodes;

    for(int i=0 ; i <TripsNo ; i++)
    {
        const TripsInfo& ti = base->TripsInfoList[i];
        const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
        const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination
        for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)
        {
            LocalGraphNodes.insert(*cit);
        }
        for(CoordSet::const_iterator cit = NonCPCoord.begin() ; cit != NonCPCoord.end() ; cit++)
        {
            LocalGraphNodes.insert(*cit);
        }
    }
    for(CoordSet::const_iterator cit = LocalGraphNodes.begin() ; cit != LocalGraphNodes.end() ; cit++)
    {
        GraphNodesList.push_back(*cit);
    }
}

void Heuristic::ComputeGlobalTheta()
{
    StopsTeta.clear();
    StopsB.clear();
    HiddenTeta.clear();
    HiddenStopTeta.clear();

    int LastNodeIndex = 0;
    int hiddenLastNodeIndex = 0;
    StopsTeta.resize(GraphNodesList.size(),0.0);
    for(size_t i=0 ; i < GraphNodesList.size()  ; i++)
    {
        if(StopsB.empty())													//avali ke hamishe sefer
            StopsB.push_back(0);
        else
            StopsB.push_back(base->Input_StopB);						//baghiye ro az voroodi migire
        const Coord& node = GraphNodesList[i];
        if(i == 0)										//noghteye aval gerafemoon teta haro hamaro sefr mikonim
        {
            StopsTeta[i] = 0.0;
            HiddenTeta.push_back(0.0);
            HiddenStopTeta.push_back(0.0);
        }
        else												//baraye baghiyeye gereha ke mabda nistan,
        {
            if(GraphCPCoords.find(node) != GraphCPCoords.end())		//age en gere jozve istgah haye ma bood, tooye refreshcheckpoint eno por karde bodim
            {
                const Coord& preNode = GraphNodesList[LastNodeIndex];		//bia istgah ghablish ro be ma bede

                double cTime = utility->CalculateDistance(preNode,node);			//va fasele en istgah ba istgah ghablish ro be ma bede
                double stopB_Add = cTime + base->Input_StopB;			// fasele be hamrahe hizhdah saniye
                double slack_Add = stopB_Add + base->Input_SlackTime;	//hame oona behamrahe slack
                double previousTheta = StopsTeta[LastNodeIndex];		//meghdare teta ro mirize toosh
                double ComputedTeta = previousTheta + slack_Add;	//mohasebeye teta
                ComputedTeta = doFixedPos(ComputedTeta);			// khataye double ro migirim ba en tabe ta 5 ragham
                StopsTeta[i] = ComputedTeta;						//mohasebeye teta ro mizarim tooye stops teta
                LastNodeIndex = i;								//akharin istgahi ke peyda kardi ro lahaz kon
            }				//ke too dore badi bedoonim akharin istgahe ma chi bude
            else						//agar gere istgah nabashe , tabiatan tetash ro sefr dar nazar migirim
            {
                StopsTeta[i] = 0.0;
            }
            /*!
              * Prepare hidden teta : hidden teta used in tau computation, this variable store teta value for all node of graph
              * checkpoints and non-checkpoints
              */
            double Distance = 0.0;
            double sumB = 0.0;
//            for(Index j=hiddenLastNodeIndex ; j <= i-1 ; j++)
            {
                const Coord& innerNode = GraphNodesList[hiddenLastNodeIndex];		//gere ghabli chi bude, va feli chie
                const Coord& innerNodeNext = GraphNodesList[i];					//beyne 2ta noghte mirim jolo
                double d = utility->CalculateDistance(innerNode,innerNodeNext);	// fasele
                Distance = d;														//mirizim too fasele
                sumB = base->Input_StopB;										//oon stopb ro mizanim
            }
            double computedTime = Distance ;							//fasele ro mirizim too time
            double computedTeta = 0.0;									//teta ro sefr mikonim
            if(GraphCPCoords.find(node) != GraphCPCoords.end())				//agar istgah bash, hiddenteta hamon stopteta hast
            {
                HiddenTeta.push_back(StopsTeta[i]);							//tetamakhfi haro dar miarim
                computedTeta = (computedTime+HiddenTeta[hiddenLastNodeIndex]+(sumB));	//zaman +makhtiteta + b
                computedTeta = doFixedPos(computedTeta);	//darsade khataro migirim
                if(computedTeta < StopsTeta[i])		//agar hiddenstopsteta, koochiktar az tetaye oon istgah bod
                    HiddenStopTeta.push_back(StopsTeta[i]); //pas mitone berese be teta istgahe badish va moshkeli nist
                else
                    HiddenStopTeta.push_back(computedTeta);	//agar be teta istgah badi naresim, mirizim tooye hidden
            }										//en komak mikone badan befahmim k teta pas nmishe
            else
            {
                computedTeta = (computedTime+HiddenTeta[hiddenLastNodeIndex]+(sumB)); //az ye gheyre gere
                computedTeta = doFixedPos(computedTeta);		//daram mirim be ye gheyre gere dige
                HiddenTeta.push_back(computedTeta);					// hiddenteta bezar meghdare mohasebe shode
                HiddenStopTeta.push_back(0.0);					//raftan be oon gere baes mishe teta pass beshe ya na
            }

            hiddenLastNodeIndex = i;
        }
    }
}

void Heuristic::PrepareGlobalParams()
{
    int TripsNo = base->TripsNo;				//yek seri init anjam midim enja ke etelaat ro farakhani konim
    int NodesNo = solver->GraphNodes.size();		// gerehamoon ro farakhani mikonim

    for(int i=0 ; i <TripsNo ; i++)											//be tedade safar ha mirim
    {
        TotalReq += base->TripsInfoList[i].RequestNo;						// darkhast haro miaym init mikonim
        for(int j=0 ; j < base->TripsInfoList[i].RequestNo ; j++)			// ye halghe baraye tamame darkhasta
        {
            TotalReqType.push_back(base->TripsInfoList[i].RequestsType[j]);		//mirizim tooye en moteghayer
        }
        base->TripsInfoList[i].q.resize(base->TripsInfoList[i].Nodes.size(),0);		//tedade q haro be andaze node ha resize mikonim
    }																	//q ro baraye har safar resize mikonim
    for(int i=0 ; i < TripsNo ; i++)												// be andaze ye safar mizanim
    {
        for(Index j=0 ; j < base->TripsInfoList[i].NodesList.size() ; j++)			//node haro migirim
            Nodes.push_back(base->TripsInfoList[i].NodesList[j]);

        for(Index j=0 ; j < base->TripsInfoList[i].RequestsPickStopNo.size() ; j++)
            RPS.push_back(base->TripsInfoList[i].RequestsPickStopNo[j]);			//rps haro migirim

        for(Index j=0 ; j < base->TripsInfoList[i].RequestsDropStopNo.size() ; j++)
            RDS.push_back(base->TripsInfoList[i].RequestsDropStopNo[j]);			// rds haro migirim

        for(Index j=0 ; j < base->TripsInfoList[i].HeuristicRequestsPickStopNo.size() ; j++)	//too solver por mishan khat 234
            HRPS.push_back(base->TripsInfoList[i].HeuristicRequestsPickStopNo[j]);			// hrps haro migirim
        																			//too heu niazi be pc va dc nadarim
        for(Index j=0 ; j < base->TripsInfoList[i].HeuristicRequestsDropStopNo.size() ; j++)
            HRDS.push_back(base->TripsInfoList[i].HeuristicRequestsDropStopNo[j]);			//hrds migiirm

        for(Index j=0 ; j < base->TripsInfoList[i].RequestsTau.size() ; j++)
            RequestsTau.push_back(base->TripsInfoList[i].RequestsTau[j]);				//tau haro migirim

        for(Index j=0 ; j < base->TripsInfoList[i].StopsType.size() ; j++)
            StopsType.push_back(base->TripsInfoList[i].StopsType[j]);					//stops type haro migirim
    }

    RejectedArcsMatrix.resize(NodesNo);						//arcsmatrix ro resize mikonim bar asase tedade nodeha,
    for(int i=0 ; i < NodesNo ; i++)
    {
        RejectedArcsMatrix[i].resize(NodesNo); 							//har element az arcs matriz shamele yal ha va har row be andaze tedade node ha resize mishe
        for(int j = 0 ; j < NodesNo ; j++)
        {
            RejectedArcsMatrix[i][j] = 0;								//va init mikonim be sefr dar ebtedaye kar ke hanooz khaliye
        }
    }

    for(Index i=0 ; i < solver->GraphNodesList.size() ; i++)
    {
        GraphNodesList.push_back(solver->GraphNodesList[i]);
    }
    for(set<int>::const_iterator it = solver->CPNodes.begin() ; it != solver->CPNodes.end() ; it++)
    {
        CPNodes.insert(*it);
    }
    RefreshCheckpointsIndices();
}

int Heuristic::GetNextCPIndexFromRoot(const Coord &rootPos)
{
    const int& rootIdx = GetIndexFromList(GraphNodesList,rootPos);
    for(Index i=rootIdx+1 ; i < GraphNodesList.size() ; i++)
    {
        if(CPNodes.find(i) != CPNodes.end())
            return i;
    }
    return -1;
}

int Heuristic::GetNextCPIndexFromRoot(const int &Idx)
{
    for(Index i=Idx+1 ; i < GraphNodesList.size() ; i++)
    {
        if(CPNodes.find(i) != CPNodes.end())
            return i;
    }
    return -1;
}

int Heuristic::GetIndexFromList(const vector<Heuristic::Coord> &input, const Heuristic::Coord &coord)
{
    for(Index i=0 ; i < input.size() ; i++)
    {
        const Coord& _c = input[i];
        if(doFixedPos(_c.first) == doFixedPos(coord.first))
        {
            return i;
        }
    }
    return -1;
}

Heuristic::Coord Heuristic::GetCoordinationFromIndex(const int &idx)
{
    return GraphNodesList[idx];
}

bool Heuristic::isEqual(const Coord &a, const Coord &b)
{
    double ep = std::numeric_limits<double>::epsilon();
//    if(fabs(a.first-b.first) <= epz && fabs(a.second-b.second) <= epz)
//    if(fabs(a.first-b.first) <= 0.000001 && fabs(a.second-b.second) <= 0.000001)
    if(a.first == b.first)
        return true;
    return false;
}

int Heuristic::GetTripNumberOfNode(const Coord &root)
{
    int TripsNo = base->TripsNo;
    for(int i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo& ti = base->TripsInfoList[i];
        if(ti.CPCoords.find(root) != ti.CPCoords.end() || ti.NonCPCoord.find(root) != ti.NonCPCoord.end())
            return ti.id;
    }
    return -1;
}
