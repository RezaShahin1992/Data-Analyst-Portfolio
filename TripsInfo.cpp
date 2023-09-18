#include "TripsInfo.h"

TripsInfo::TripsInfo()
{

}

TripsInfo::~TripsInfo()
{

}

void TripsInfo::ClearTripsInfo()
{
    Nodes.clear();
    q.clear();
    NodesList.clear();
    Nodes.clear();
    StopsType.clear();
    StopsTeta.clear();
    HiddenTeta.clear();
    StopsB.clear();
    RequestsType.clear();
    RequestsPickStopNo.clear();
    RequestsDropStopNo.clear();
    RequestsTau.clear();
    RequestsXX1.clear();
    RequestsYY1.clear();
    RequestsXX2.clear();
    RequestsYY2.clear();
    CP_Index.clear();
    NonCP_Index.clear();
    CPCoords.clear();
    NonCPCoord.clear();
    NodesCoord.clear();
    ArcsMatrix.clear();
    PickStopXList.clear();
    PickStopYList.clear();
    DropStopXList.clear();
    DropStopYList.clear();
    W.clear();
    TripPickDropSetNonCP.clear();
    m_PossibleTripsForRequests.clear();

    id = 0;

    RequestNo=0, NodesNo=0, VehicleIndex=0,C=0;
    VehicleCapacity = 0;
    StartTime = -1;
    SuspendTime = 0;
    RealStartTime = 0;
    RealArrivalTime = 0;
    TripComputedTime = 0.0;
}
