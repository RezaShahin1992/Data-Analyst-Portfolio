#include "SolverBase.h"

SolverBase::SolverBase()
{

}

void SolverBase::execute()
{
    int nTrips = TripsInfoList.size();				//we assign the nTrips with the size of our TripsInfoList
    const double& slacktime = Input_SlackTime;		//read the slack time

    for(int i=0 ; i < nTrips ; i++)
    {
        int TripID = i;
        TripsInfo& ti = TripsInfoList[TripID];		//here we count the number of trips.

        ComputeTripTeta(ti,slacktime);
        ComputeGlobalTeta(slacktime);				//we compute the teta

        ti.RealArrivalTime = ti.RealStartTime + ti.StopsTeta[ti.StopsTeta.size()-1];

    }

    TripsNo = TripsInfoList.size();

    doVehiclesGeneration(TripsNo);					//we generate the vehicles

    VehiclesNo = VehiclesSet.size();
}

double SolverBase::CalculateDistance(const SolverBase::Coord &a, const SolverBase::Coord &b)
{
//    double x1 = a.first;						//dar gheyre en soorat miaym mokhtasate 2ta noghte ro be dast miarim
//    double y1 = a.second;
//    double x2 = b.first;
//    double y2 = b.second;
//    double xDis = std::fabs(x2-x1);
//    double yDis = std::fabs(y2-y1);
//    double d = xDis + yDis;		// fasele 2ta noghte ro hesab mikonim
//    double dis = d / (M_Velocity);
    double x1 = a.first;						//dar gheyre en soorat miaym mokhtasate 2ta noghte ro be dast miarim
    double y1 = a.second;
    double x2 = b.first;
    double y2 = b.second;
    double xDis = std::fabs(std::fabs(x2)-std::fabs(x1));
    double yDis = std::fabs(std::fabs(y2)-std::fabs(y1));
//    double d = xDis + yDis;		// fasele 2ta noghte ro hesab mikonim
    double d = sqrt((a.first - b.first) * (a.first - b.first) + (a.second - b.second) * (a.second - b.second));
    double dis = d / (M_Velocity);

    dis = doFixedPos(dis);

    return dis;
}

double SolverBase::doFixedPos(const double &p)
{
    double a = /*std::ceil*/(p*100000.00000f)/100000.00000f;
    stringstream tmp;
    tmp << setprecision(5) << fixed << a;
    double new_val = stod(tmp.str());
    return new_val;
}

void SolverBase::ComputeTripTeta(TripsInfo &ti, const double &slacktime)
{
    ti.StopsType.clear();
    ti.StopsB.clear();
    ti.StopsTeta.clear();

    const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
    const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination
    CoordSet GraphNodes;
    for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)
    {
        GraphNodes.insert(*cit);
    }
    for(CoordSet::const_iterator cit = NonCPCoord.begin() ; cit != NonCPCoord.end() ; cit++)
    {
        GraphNodes.insert(*cit);
    }

    std::vector<Coord> NodesList;
    for(CoordSet::const_iterator cit = GraphNodes.begin() ; cit != GraphNodes.end() ; cit++)
    {
        const Coord& node = *cit;
        NodesList.push_back(node);
    }

    for(size_t i=0 ; i < NodesList.size()  ; i++)
    {
        const Coord& node = NodesList[i];
        if(ti.CPCoords.find(node) != ti.CPCoords.end())
        {
            ti.StopsType.push_back(1);
        }
        else
        {
            ti.StopsType.push_back(2);
        }
        if(ti.StopsB.empty())
            ti.StopsB.push_back(0);
        else
            ti.StopsB.push_back(Input_StopB);
    }

    int LastNodeIndex = 0;
    int hiddenLastNodeIndex = 0;
    ti.StopsTeta.resize(NodesList.size(),0.0);
    for(size_t i=0 ; i < NodesList.size()  ; i++)
    {
        const Coord& node = NodesList[i];
        if(i== 0 && ti.id == 0) // First Node Of First Trip :: 0 value for teta
        {
            ti.StopsTeta[i] = 0.0;
            ti.HiddenTeta.push_back(0.0);
        }
        else
        {
            if(i == 0 && isAReturnTrip(ti.id)) // first node of retrun trip
            {
                ///@brief :: Use Zero Teta For Return Trip
                ti.StopsTeta[i] = 0.0;
                ti.HiddenTeta.push_back(0.0);
            }
            else
            {
                if(i == 0)
                {
                    ti.StopsTeta[i] = 0.0;
                    ti.HiddenTeta.push_back(0.0);
                    continue;
                }
                else
                {
                    if(ti.CPCoords.find(node) != ti.CPCoords.end())
                    {

                        double x2 = node.first;
                        double y2 = node.second;

                        const Coord& preNode = NodesList[LastNodeIndex];
                        double x1 = preNode.first;
                        double y1 = preNode.second;

                        double dis = sqrt((x2 - x1)*(x2-x1) + (y2-y1)*(y2-y1));
                        double cTime = dis / M_Velocity;
                        double ComputedTeta = cTime + ti.StopsTeta[LastNodeIndex] + slacktime + ti.StopsB[i];
                        ti.StopsTeta[i] = ComputedTeta;
                        LastNodeIndex = i;
                    }
                    else
                    {
                        ti.StopsTeta[i] = 0.0;
                    }
                    /*!
                      * Prepare hidden teta : hidden teta used in tau computation, this variable store teta value for all node of graph
                      * checkpoints and non-checkpoints
                      */
                    double Distance = 0.0;
                    double sumB = 0.0;
                    for(Index j=hiddenLastNodeIndex ; j <= i-1 ; j++)
                    {
                        const Coord& innerNode = NodesList[j];
                        const Coord& innerNodeNext = NodesList[j+1];
                        double x1 = innerNode.first;
                        double y1 = innerNode.second;
                        double x2 = innerNodeNext.first;
                        double y2 = innerNodeNext.second;
                        double d = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
                        Distance += d;
                        sumB+=ti.StopsB[j+1];
                    }
                    double computedTime = Distance / M_Velocity;
//                    computedTime = ceil(computedTime);
                    double computedTeta = -1.0;
                    if(ti.CPCoords.find(node) != ti.CPCoords.end())
                    {
                        double tmpHDN = (computedTime*PER_MINUTE+ti.HiddenTeta[hiddenLastNodeIndex]+(sumB));
                        if(tmpHDN < ti.StopsTeta[i])
                            computedTeta = ti.StopsTeta[i];
                        else
                            computedTeta = tmpHDN;
                    }
                    else
                    {
                        computedTeta = (computedTime*PER_MINUTE+ti.HiddenTeta[hiddenLastNodeIndex]+(sumB));
                    }
                    ti.HiddenTeta.push_back(computedTeta);
                    hiddenLastNodeIndex = i;
                }
            }
        }
    }
}

void SolverBase::ComputeGlobalTeta(const double &slacktime)
{
    StopsTeta.clear();
    StopsB.clear();
    TripsHeadTailList.clear();

    CoordSet GraphNodes;
    CoordSet GraphCP;
    for(int i=0 ; i < TripsInfoList.size() ; i++)
    {
        int TripID = i;
        TripsInfo& ti = TripsInfoList[TripID];
        const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
        const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination

        Coord TripSource;
        Coord TripDestination;
        for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)
        {
            if(cit == CPCoords.begin())
            {
                TripSource = *cit;
            }
            if(cit == --CPCoords.end())
            {
                TripDestination = *cit;
            }
            GraphNodes.insert(*cit);
            GraphCP.insert(*cit);
        }
        pair<int,pair<Coord,Coord> > SD_Trip = make_pair(i,make_pair(TripSource,TripDestination));
        TripsHeadTailList.push_back(SD_Trip);
        for(CoordSet::const_iterator cit = NonCPCoord.begin() ; cit != NonCPCoord.end() ; cit++)
        {
            GraphNodes.insert(*cit);
        }
    }
    std::vector<Coord> NodesList;
    for(CoordSet::const_iterator cit = GraphNodes.begin() ; cit != GraphNodes.end() ; cit++)
    {
        const Coord& node = *cit;
        NodesList.push_back(node);
    }

    int LastNodeIndex = 0;
    int hiddenLastNodeIndex = 0;
    StopsTeta.resize(NodesList.size(),0.0);
    for(size_t i=0 ; i < NodesList.size()  ; i++)
    {
        if(StopsB.empty())
            StopsB.push_back(0);
        else
            StopsB.push_back(Input_StopB);
        const Coord& node = NodesList[i];
        if(i == 0) // First Node Of First Trip :: 0 value for teta
        {
            StopsTeta[i] = 0.0;
        }
        else
        {
            if(GraphCP.find(node) != GraphCP.end())
            {
                int CurrentTripIdx = isTripSourceCoord(node);
                double x2 = node.first;
                double y2 = node.second;

                const Coord& preNode = NodesList[LastNodeIndex];
                double x1 = preNode.first;
                double y1 = preNode.second;

                double cTime = CalculateDistance(node,preNode);
                double stopB_Add = cTime + Input_StopB;
                double slack_Add = stopB_Add + slacktime;
                double previousTheta = StopsTeta[LastNodeIndex];
                double ComputedTeta = previousTheta + slack_Add;
                ComputedTeta = doFixedPos(ComputedTeta);
                if(CurrentTripIdx != -1)
                {
                    StopsTeta[i] = TripsInfoList[CurrentTripIdx].RealStartTime + ComputedTeta;
                }
                else
                    StopsTeta[i] = ComputedTeta;
                LastNodeIndex = i;
            }
            else
            {
                StopsTeta[i] = 0.0;
            }
        }
    }
}

void SolverBase::doVehiclesGeneration(const int &nbTrips)
{
    for(int i=0 ; i < nbTrips ; i++)
    {
        int TripID = i;
        TripsInfo& ti = TripsInfoList[TripID];
        std::vector<TripsInfo> tri_pre_list;
        GenerateVehicles(ti);
        VehiclesSet.insert(ti.VehicleIndex);
    }
}

void SolverBase::GenerateVehicles(TripsInfo &ti)
{
    bool isSameVehicle = false;
    for(std::map<int,double>::const_iterator cit = VehicleFreeTime.begin() ; cit != VehicleFreeTime.end() ; cit++)
    {
        if(ti.RealStartTime > cit->second)
        {
            ti.VehicleIndex = cit->first;
            VehicleFreeTime[cit->first] = ti.RealArrivalTime;
            isSameVehicle = true;
            break;
        }
    }
    if(!isSameVehicle)
    {
        ti.VehicleIndex = VehicleContainerList.size();
        VehicleFreeTime[ti.VehicleIndex] = ti.RealArrivalTime;
        VehicleContainer vc(true,ti.VehicleIndex);
        VehicleContainerList.push_back(vc);
    }
}

bool SolverBase::isAReturnTrip(const int &tripID)
{
    if(tripID % 2 != 0)
        return true;
    return false;
}

SolverBase::Coord SolverBase::GetCoordinationFromIndex(const int &idx)
{
    CoordSet d_GraphNodesCoord;
    CoordList d_GraphNodesCoordList;
    int TripsNo = TripsInfoList.size();

    for(int i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo& ti = TripsInfoList[i];		//miad tripsinfo object ro ba indice i az tripsinfolist migire k too rh zakhire shode
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

int SolverBase::isTripSourceCoord(const Coord &p)
{
    for(int i=0 ; i < TripsHeadTailList.size() ; i++)
    {
        int tripIdx = TripsHeadTailList[i].first;
        Coord Source = TripsHeadTailList[i].second.first;
        Coord Destination = TripsHeadTailList[i].second.second;

        if(p.first == Source.first)
        {
            return tripIdx;
        }
    }
    return -1;
}
