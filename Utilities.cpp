#include "Utilities.h"

Utilities::Utilities()
{

}

void Utilities::ComputeTripTime(SolverBase *sb)
{
    int TripsNo = sb->TripsNo;
    for(int i=0 ; i < TripsNo ; i++)
    {
        TripsInfo& ti = sb->TripsInfoList[i];
        const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
        const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination
        CoordSet NodesCoordList;
        double Distance = 0;
        for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)
        {
            NodesCoordList.insert(*cit);
        }
        for(CoordSet::const_iterator cit = NonCPCoord.begin() ; cit != NonCPCoord.end() ; cit++)
        {
            NodesCoordList.insert(*cit);
        }
        std::vector<Coord> NodesList;
        for(CoordSet::const_iterator cit = NodesCoordList.begin() ; cit != NodesCoordList.end() ; cit++)
        {
            const Coord& node = *cit;
            NodesList.push_back(node);
        }
        for(size_t j=0 ; j < NodesList.size() -1 ; j++)
        {
            const Coord& node = NodesList[j];
            double x1 = node.first;
            double y1 = node.second;
            const Coord& nextNode = NodesList[j+1];
            double x2 = nextNode.first;
            double y2 = nextNode.second;
            double d = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
            Distance += d;
        }

        ti.TripComputedTime = Distance / M_Velocity;
        ti.RealStartTime = ti.StartTime;
        ti.RealArrivalTime = ti.RealStartTime + ti.StopsTeta[ti.StopsTeta.size()-1] - ti.SuspendTime;
        std::cerr << __FUNCTION__ << "Trip["<<i<<"] . Start Time : " <<   ti.RealStartTime << std::endl;
        std::cerr << __FUNCTION__ << "Trip["<<i<<"] . Arrival Time : " << ti.RealArrivalTime << std::endl;
    }
}

int Utilities::GetNbReqFromType(const SolverBase *sb, const REQUESTTYPE &type)
{
    int TripsNo = sb->TripsNo;
    int nbReq = 0;
    for(int i=0 ; i <TripsNo ; i++)
    {
        for(int j=0 ; j < sb->TripsInfoList[i].RequestNo ; j++)
        {
            if(static_cast<REQUESTTYPE>(sb->TripsInfoList[i].RequestsType[j]) == type) //?
            {
                nbReq++;
            }
        }
    }
    return nbReq;
}

int Utilities::GetNbReqFromRPS(const SolverBase *sb, const int ps, const REQUESTTYPE &type) // we fill the RPS and RDS
{
    int nbReq = 0;
    Index TotalReq = 0;
    vector<int> TotalReqType;
    vector<int> RPS;
    vector<int> RDS;
    int TripsNo = sb->TripsNo;
    for(int i=0 ; i <TripsNo ; i++)
    {
        TotalReq += sb->TripsInfoList[i].RequestNo;
        for(int j=0 ; j < sb->TripsInfoList[i].RequestNo ; j++)
        {
            TotalReqType.push_back(sb->TripsInfoList[i].RequestsType[j]);
        }
        for(Index j=0 ; j < sb->TripsInfoList[i].RequestsPickStopNo.size() ; j++)
            RPS.push_back(sb->TripsInfoList[i].RequestsPickStopNo[j]);

        for(Index j=0 ; j < sb->TripsInfoList[i].RequestsDropStopNo.size() ; j++)
            RDS.push_back(sb->TripsInfoList[i].RequestsDropStopNo[j]);
    }

    for (Index k = 0; k < TotalReq; k++)
    {
        if (static_cast<REQUESTTYPE>(TotalReqType[k]) == type)
        {
            if(RPS[k] == ps)
            {
                nbReq++;
            }
        }
    }
    return nbReq;
}

int Utilities::GetNbReqFromRDS(const SolverBase *sb, const int ds, const REQUESTTYPE &type)
{
    int nbReq = 0;
    Index TotalReq = 0;
    vector<int> TotalReqType;
    vector<int> RPS;
    vector<int> RDS;
    int TripsNo = sb->TripsNo;
    for(int i=0 ; i <TripsNo ; i++)
    {
        TotalReq += sb->TripsInfoList[i].RequestNo;
        for(int j=0 ; j < sb->TripsInfoList[i].RequestNo ; j++)
        {
            TotalReqType.push_back(sb->TripsInfoList[i].RequestsType[j]);
        }
        for(Index j=0 ; j < sb->TripsInfoList[i].RequestsPickStopNo.size() ; j++)
            RPS.push_back(sb->TripsInfoList[i].RequestsPickStopNo[j]);

        for(Index j=0 ; j < sb->TripsInfoList[i].RequestsDropStopNo.size() ; j++)
            RDS.push_back(sb->TripsInfoList[i].RequestsDropStopNo[j]);
    }

    for (Index k = 0; k < TotalReq; k++)
    {
        if (static_cast<REQUESTTYPE>(TotalReqType[k]) == type)
        {
            if(RDS[k] == ds)
            {
                nbReq++;
            }
        }
    }
    return nbReq;
}

int Utilities::GetNbReqFromPC(const vector <vector<vector<int>>>& S_PND, const SolverBase *sb, const int pick)
{
    int nbReq = 0;
    Index TotalReq = 0;
    vector<int> TotalReqType;
    int TripsNo = sb->TripsNo;
    int VehicleNo = sb->VehiclesNo;

    for(int i=0 ; i <TripsNo ; i++)
    {
        TotalReq += sb->TripsInfoList[i].RequestNo;
        for(int j=0 ; j < sb->TripsInfoList[i].RequestNo ; j++)
        {
            TotalReqType.push_back(sb->TripsInfoList[i].RequestsType[j]);
        }
    }

    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::PND)
        {
            for(int r=0 ; r < TripsNo ; r++)
            {
                for(int v=0 ; v < VehicleNo ; v++)
                {
                    if(S_PND[k][r][v] == pick)
                    {
                        nbReq++;
                    }
                }
            }
        }
    }
    return 1;
}

int Utilities::GetNbReqFromDC(const vector<vector<vector<int> > > &S_NPD, const SolverBase *sb, const int drop)
{
    int nbReq = 0;
    Index TotalReq = 0;
    vector<int> TotalReqType;
    int TripsNo = sb->TripsNo;
    int VehicleNo = sb->VehiclesNo;

    for(int i=0 ; i <TripsNo ; i++)
    {
        TotalReq += sb->TripsInfoList[i].RequestNo;
        for(int j=0 ; j < sb->TripsInfoList[i].RequestNo ; j++)
        {
            TotalReqType.push_back(sb->TripsInfoList[i].RequestsType[j]);
        }
    }

    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == REQUESTTYPE::NPD)
        {
            for(int r=0 ; r < TripsNo ; r++)
            {
                for(int v=0 ; v < VehicleNo ; v++)
                {
                    if(S_NPD[k][r][v] == drop)
                    {
                        nbReq++;
                    }
                }
            }
        }
    }
    return 1;
}

int Utilities::GetTripNumberOfRequest(const SolverBase *sb, const int &reqID)
{
    int DetectedTrip = -1;
    for(int i=0 ; i < sb->RequestsTypeIndex.size() ; i++)
    {
        int CurrentTrip = sb->RequestsTypeIndex[i].first;
        int CurrentReq  = sb->RequestsTypeIndex[i].second;
        if(CurrentReq == reqID)
        {
            DetectedTrip = CurrentTrip;
            break;
        }
    }
    return DetectedTrip;
}

Utilities::Coord Utilities::GetCoordinationFromIndex(const SolverBase* sb,const int &idx) /// save all nodes coordinations
{
    CoordSet d_GraphNodesCoord;
    CoordList d_GraphNodesCoordList;
    int TripsNo = sb->TripsInfoList.size();

    for(int i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo& ti = sb->TripsInfoList[i];		//miad tripsinfo object ro ba indice i az tripsinfolist migire k too rh zakhire shode
        const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
        const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination
        for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)		//istgah ha va gheyre istgah haro miad migire
        {																							// va mirize tooye set graphnodes.
            d_GraphNodesCoord.insert(*cit);
        }
        for(CoordSet::const_iterator cit = NonCPCoord.begin() ; cit != NonCPCoord.end() ; cit++)
        {
            d_GraphNodesCoord.insert(*cit);
        }
    }
    for(CoordSet::const_iterator cit = d_GraphNodesCoord.begin() ; cit != d_GraphNodesCoord.end() ; cit++)
    {
        d_GraphNodesCoordList.push_back(*cit);
    }

    return d_GraphNodesCoordList[idx];
}

int Utilities::GetCurrentRequestNo(const SolverBase *sb, const REQUESTTYPE &rt) //compute the type of each request
{
    int TripsNo = sb->TripsNo;
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

    int nbRes = 0;
    for (Index k = 0; k < TotalReq; k++)
    {
        if (TotalReqType[k] == rt)
        {
            nbRes++;
        }
    }
    return nbRes;
}

vector<Utilities::Index> Utilities::GetRequestsInNode(const vector <vector<vector<int>>>& S_PND,const SolverBase *sb, const int &idx)
{
    const int& VehiclesNo = sb->VehiclesNo;
    int TripsNo = sb->TripsNo;

    Index TotalReq = 0;
    vector<int> TotalReqType;
    vector<int> Nodes;
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
        for(Index j=0 ; j < sb->TripsInfoList[i].NodesList.size() ; j++)
            Nodes.push_back(sb->TripsInfoList[i].NodesList[j]);

        for(Index j=0 ; j < sb->TripsInfoList[i].RequestsPickStopNo.size() ; j++)
            RPS.push_back(sb->TripsInfoList[i].RequestsPickStopNo[j]);

        for(Index j=0 ; j < sb->TripsInfoList[i].RequestsDropStopNo.size() ; j++)
            RDS.push_back(sb->TripsInfoList[i].RequestsDropStopNo[j]);
    }

    set<Index> ReqRes;

    for(Index k=0 ; k < TotalReq ; k++)
    {
        if(TotalReqType[k] == REQUESTTYPE::PND)
        {
            for(Index r=0 ; r < TripsNo ; r++)
            {
                int v = sb->TripsInfoList[r].VehicleIndex;
                if(S_PND[k][r][v] == idx)
                {
                    ReqRes.insert(k);
                }
            }
        }

        if(TotalReqType[k] == REQUESTTYPE::NPND || TotalReqType[k] == REQUESTTYPE::PD || TotalReqType[k] == REQUESTTYPE::NPD)
        {
            if(RPS[k] == idx)
            {
                ReqRes.insert(k);
            }
        }
    }
    vector<Index> ReqList;
    for(set<Index>::const_iterator cit = ReqRes.begin() ; cit != ReqRes.end() ; cit++)
    {
        ReqList.push_back(*cit);
    }
    return ReqList;
}

double Utilities::ComputeDeltaMinReqFromPickToDrop(const set<int>& CPNodes, const SolverBase *sb, const int &PickUp, const int &DropOff)
{
    double t = 0;
    int lastIdx = PickUp;
    for(int i = PickUp+1 ; i <= DropOff ; i++)
    {
        if(CPNodes.find(i) != CPNodes.end() || i == DropOff)
        {
            Coord LastCoord = GetCoordinationFromIndex(sb,lastIdx);
            Coord CheckpointCoord = GetCoordinationFromIndex(sb,i);
            t += CalculateDistance(LastCoord,CheckpointCoord);
            lastIdx = i;
        }
    }
    return t;
}

string Utilities::ConvertToRealTime(const double &time)
{
    double Hour = 0.0;
    double lowerHour = 0.0;
    lowerHour = std::modf(time,&Hour);

    double min = lowerHour * 60.0;

    int i_Hour = Hour;
    int i_min = 0;

    double uppersec = 0.0;
    double lowersec = std::modf(min,&uppersec);

    i_min = uppersec;

    double sec = lowersec * 60.0;

    int i_sec = ceil(sec);


    if(i_sec == 60)
    {
        i_sec = 0;
        i_min++;
    }

    if(i_min == 60)
    {
        i_min = 0;
        i_Hour++;
    }


    string str_hour = to_string(i_Hour);
    string str_min = to_string(i_min);
    string str_sec = to_string(i_sec);

    return (str_hour + ":" + str_min + ":" + str_sec);
}

double Utilities::CalculateDistance(const Utilities::Coord &a, const Utilities::Coord &b)
{
//    double x1 = a.first;						//dar gheyre en soorat miaym mokhtasate 2ta noghte ro be dast miarim
//    double y1 = a.second;
//    double x2 = b.first;
//    double y2 = b.second;
//    double xDis = std::fabs(std::fabs(x2)-std::fabs(x1));
//    double yDis = std::fabs(std::fabs(y2)-std::fabs(y1));
//    double d = xDis + yDis;		// fasele 2ta noghte ro hesab mikonim
    double d = sqrt((a.first - b.first) * (a.first - b.first) + (a.second - b.second) * (a.second - b.second));
    double dis = d / (M_Velocity);

    dis = doFixedPos(dis);

    return dis;
}

vector<int> Utilities::GetSameCPIndexFromRoot(const SolverBase *sb, const CoordList &GraphNodesList, const Utilities::Coord &rootPos)
{
    vector<int> res;
    Index TripsNo = sb->TripsInfoList.size();
    CoordSet d_GraphNodesCoord;

    CoordSet d_CPSet;
    int detectedTripNo = 0;
    int globalIdx = 0;
    for(Index i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo& ti = sb->TripsInfoList[i];		//miad tripsinfo object ro ba indice i az tripsinfolist migire k too rh zakhire shode
        const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
        const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination

        int localCPIndex = 0;
        for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)		//istgah ha va gheyre istgah haro miad migire
        {																							// va mirize tooye set graphnodes.
            d_CPSet.insert(*cit);
            if(rootPos.first == cit->first)
            {
                globalIdx = localCPIndex;
                detectedTripNo = ti.id;
                break;
            }
            localCPIndex++;
        }
    }

    CoordList DetectedCPCoord;
    for(Index i = detectedTripNo+1 ; i < TripsNo ; i++)
    {
        const TripsInfo& ti = sb->TripsInfoList[i];		//miad tripsinfo object ro ba indice i az tripsinfolist migire k too rh zakhire shode
        const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
        int localCPIndex = 0;
        for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)		//istgah ha va gheyre istgah haro miad migire
        {																							// va mirize tooye set graphnodes.
            if(localCPIndex == globalIdx)
            {
                DetectedCPCoord.push_back(*cit);
                break;
            }
            localCPIndex++;
        }
    }
    for(int i=0 ; i < DetectedCPCoord.size() ; i++)
    {
        for(int j=0 ; j < GraphNodesList.size() ; j++)
        {
            if(DetectedCPCoord[i].first == GraphNodesList[j].first)
            {
                res.push_back(j);
            }
        }
    }
    return res;
}

Utilities::CoordList Utilities::GetSimilarPositionOfNode(const SolverBase *sb, const Utilities::Coord &rootPos)
{
    CoordList res;
    Index TripsNo = sb->TripsInfoList.size();

    CoordSet d_CPSet;
    int detectedTripNo = 0;
    int globalIdx = 0;
    for(Index i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo& ti = sb->TripsInfoList[i];		//miad tripsinfo object ro ba indice i az tripsinfolist migire k too rh zakhire shode
        const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
        const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination

        if(CPCoords.find(rootPos) != CPCoords.end() || NonCPCoord.find(rootPos) != NonCPCoord.end())
        {
            detectedTripNo = ti.id;
        }
    }
    Coord firstPosOfTrip = *sb->TripsInfoList[detectedTripNo].CPCoords.begin();
    double diffX = std::fabs(rootPos.first - firstPosOfTrip.first);
    for(int i=detectedTripNo ; i < TripsNo ; i++)
    {
        Coord firstPosOfCurrentTrip = *sb->TripsInfoList[i].CPCoords.begin();
        Coord newPos;
        newPos.first = firstPosOfCurrentTrip.first + diffX;
        newPos.second = rootPos.second;
        res.push_back(newPos);
    }
    return res;
}

Utilities::Coord Utilities::GetSimilarPositionOfNodeInTrip(const SolverBase *sb, const int &TripIdx, const Utilities::Coord &rootPos)
{
    Coord res;
    Index TripsNo = sb->TripsInfoList.size();

    int detectedTripNo = 0;
    for(Index i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo& ti = sb->TripsInfoList[i];		//miad tripsinfo object ro ba indice i az tripsinfolist migire k too rh zakhire shode
        const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
        const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination

        if(CPCoords.find(rootPos) != CPCoords.end() || NonCPCoord.find(rootPos) != NonCPCoord.end())
        {
            detectedTripNo = ti.id;
        }
    }
    Coord firstPosOfTrip = *sb->TripsInfoList[detectedTripNo].CPCoords.begin();
    double diffX = std::fabs(rootPos.first - firstPosOfTrip.first);

    Coord firstPosOfCurrentTrip = *sb->TripsInfoList[TripIdx].CPCoords.begin();
    Coord newPos;
    newPos.first = firstPosOfCurrentTrip.first + diffX;
    newPos.second = rootPos.second;
    res = (newPos);

    return res;
}

vector<int> Utilities::GetNearestCPToCurrentPos(const CoordList &GraphCPGlobalList, const CoordSet &GraphNodes, const Utilities::Coord &rootPos)
{
    double minDisback = std::numeric_limits<double>::max();
    double minDisnext = std::numeric_limits<double>::max();
    Coord next ;
    Coord back ;
    for(int i=0 ; i < GraphCPGlobalList.size() ; i++)
    {
        double diffmin = std::fabs(rootPos.first - GraphCPGlobalList[i].first);
        if(rootPos.first > GraphCPGlobalList[i].first)
        {
            if(diffmin < minDisback)
            {
                minDisback = diffmin;
                back = GraphCPGlobalList[i];
            }
        }
        else
        {
            if(diffmin < minDisnext)
            {
                minDisnext = diffmin;
                next = GraphCPGlobalList[i];
            }
        }
    }
    int backIdx = GetIndexFromSet(GraphNodes,back);
    int nextIdx = GetIndexFromSet(GraphNodes,next);

    vector<int> res;
    res.push_back(backIdx);
    res.push_back(nextIdx);

    return res;
}

Utilities::Index Utilities::GetVehicleNumberOfCurrentNode(const SolverBase *sb, const Utilities::Index &idx)
{
    Index TripNo = sb->TripsInfoList.size();
    for(Index i=0 ; i < TripNo ; i++)
    {
        const TripsInfo& ti = sb->TripsInfoList[i];
        if(ti.Nodes.find(idx) != ti.Nodes.end())
            return ti.VehicleIndex;
    }
    return 0;
}

int Utilities::GetNextCPIndexFromRoot(const set<int> CPNodes, const CoordList &GraphNodesList, const Utilities::Coord &rootPos)
{
    const int& rootIdx = GetIndexFromList(GraphNodesList,rootPos);
    for(Index i=rootIdx+1 ; i < GraphNodesList.size() ; i++)
    {
        if(CPNodes.find(i) != CPNodes.end())
            return i;
    }
    return -1;
}

double Utilities::GetDistanceFromSrcToDes(const SolverBase *sb, const int &src, const int &des)
{
    Coord srcCoord = GetCoordinationFromIndex(sb,src);
    Coord desCoord = GetCoordinationFromIndex(sb,des);
    Index TripsNo = sb->TripsInfoList.size();

    int srcTripIndex = -1;
    int desTripIndex = -1;
    for(Index i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo& ti = sb->TripsInfoList[i];		//miad tripsinfo object ro ba indice i az tripsinfolist migire k too rh zakhire shode
        const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
        const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination
        if(CPCoords.find(srcCoord) != CPCoords.end() || NonCPCoord.find(srcCoord) != NonCPCoord.end())
        {
            srcTripIndex = ti.id;
        }
        if(CPCoords.find(desCoord) != CPCoords.end() || NonCPCoord.find(desCoord) != NonCPCoord.end())
        {
            desTripIndex = ti.id;
        }
    }
    Coord firstPosOfSrc = *sb->TripsInfoList[srcTripIndex].CPCoords.begin();
    Coord firstPosOfDes = *sb->TripsInfoList[desTripIndex].CPCoords.begin();
    double diffX_Src = std::fabs(srcCoord.first - firstPosOfSrc.first);
    Coord newPos_Des;
    newPos_Des.first = firstPosOfDes.first + diffX_Src;
    newPos_Des.second = srcCoord.second;

    double dis = CalculateDistance(desCoord,newPos_Des);

    return dis;
}

int Utilities::GetIndexFromSet(const Utilities::CoordSet &input, const Utilities::Coord &coord)
{
    int cnt = 0;
    for(CoordSet::const_iterator cit = input.begin() ; cit != input.end() ; cit++)
    {
        const Coord& _c = *cit;
        if(doFixedPos(_c.first) == doFixedPos(coord.first))
        {
            return cnt;
        }
        cnt++;
    }
    return -1;
}

int Utilities::GetIndexFromList(const vector<Utilities::Coord> &input, const Utilities::Coord &coord)
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

bool Utilities::isEqual(const double &a, const double &b)
{
    double ep = std::numeric_limits<double>::epsilon();
    if(a == b)
        return true;
    return false;
}

bool Utilities::isEqual(const Utilities::Coord &a, const Utilities::Coord &b)
{
    double ep = std::numeric_limits<double>::epsilon();
    if(a.first == b.first && a.second == b.second)
        return true;
    return false;
}

double Utilities::doFixedPos(const double &p)
{
    double a = (p*100000.00000f)/100000.00000f;
    stringstream tmp;
    tmp << setprecision(5) << fixed << a;
    double new_val = stod(tmp.str());
    return new_val;
}
