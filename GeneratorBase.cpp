#include "GeneratorBase.h"
#include "unordered_map"

GeneratorBase::GeneratorBase()
{
}

void GeneratorBase::PrintRawDataToFile()
{

    int nbTrips = GeneratorTrips.size();

    fstream Inf;
    Inf.open("OUTPUT/Inf-"+to_string(Input_CntInput)+".txt",ios_base::out | ios_base::trunc);
    Inf << "TripsTotal:" << "\t" << nbTrips << "\n";
    Inf << "SlackTime:" << "\t" << Input_SlackTime << "\n";
    Inf << "StopB:" << "\t" << Input_StopB << "\n";
    Inf << "SAX:" << "\t" << Input_ServiceAreaX << "\n";
    Inf << "SAY:" << "\t" << Input_ServiceAreaY << "\n";
    Inf << "Interval:" << "\t" << Input_StartTimeInterval << "\n";
    Inf << "PD:" << "\t" << Input_PD << "\n";
    Inf << "PND:" << "\t" << Input_PND << "\n";
    Inf << "NPD:" << "\t" << Input_NPD << "\n";
    Inf << "NPND:" << "\t" << Input_NPND << "\n";
    Inf << "GenerationTime:" << "\t" << m_TripsGenerationTime << "\n";

    for(int i=0 ; i < nbTrips ; i++)
    {
        int TripID = i;
        TripsInfo& ti = GeneratorTrips[TripID];
        Inf << "New:" << "\t" << ti.isNew << "\n";
        Inf << "VehiclesIndex:" << "\t" << ti.VehicleIndex << "\n";
        Inf << "VehiclesCapacity:" << "\t" << ti.VehicleCapacity << "\n";
        Inf << "Time:" << "\t" << ti.StartTime << "\n";
        Inf << "SuspendTime:" << "\t" << ti.SuspendTime << "\n";
        Inf << "RequestNo:" << "\t" << ti.RequestNo << "\n";
        Inf << "NodesNo:" << "\t" << ti.NodesNo << "\n";
        Inf << "C:" << "\t" << ti.C << "\n";
        Inf << "\n";
        Inf << "Checkpoints:" << "\n";
        for(CoordSet::iterator it = ti.CPCoords.begin() ; it != ti.CPCoords.end() ; it++)
        {
            Inf << it->first << "\t" << it->second << "\n";
        }
        Inf << "\n";
//        Inf << "NonCheckpoints:" << "\n";
//        for(CoordSet::iterator it = ti.NonCPCoord.begin() ; it != ti.NonCPCoord.end() ; it++)
//        {
//            Inf << it->first << "\t" << it->second << "\n";
//        }
//        Inf << "\n";
        Inf << "stops.Type:" << "\t" ;
        for(int i=0 ; i < ti.StopsType.size() ; i++)
        {
            Inf << ti.StopsType[i] << "\t";
        }
        Inf << "\n";
        Inf << "stops.Teta:" << "\t" ;
        for(int i=0 ; i < ti.StopsTeta.size() ; i++)
        {
            Inf << ti.StopsTeta[i] << "\t";
        }
        Inf << "\n";
        Inf << "stops.b:" << "\t" ;
        for(int i=0 ; i < ti.StopsB.size() ; i++)
        {
            Inf << ti.StopsB[i] << "\t";
        }
        Inf << "\n";
        Inf << "\n";
    }
    Inf.close();
    fstream Ins;
    Ins.open("OUTPUT/Ins-"+to_string(Input_CntInput)+".txt",ios_base::out | ios_base::trunc);
    for(int i=0 ; i < nbTrips ; i++)
    {
        int TripID = i;
        TripsInfo& ti = GeneratorTrips[TripID];
        Ins << "Requests.PickStopX:" << "\t" ;
        for(int i=0 ; i < ti.PickStopXList.size() ; i++)
        {
            Ins << doFixedPos(ti.PickStopXList[i]) << "\t";
        }
        Ins << "\n";
        Ins << "Requests.PickStopY:" << "\t" ;
        for(int i=0 ; i < ti.PickStopYList.size() ; i++)
        {
            Ins << doFixedPos(ti.PickStopYList[i]) << "\t";
        }
        Ins << "\n";
        Ins << "Requests.DropStopX:" << "\t" ;
        for(int i=0 ; i < ti.DropStopXList.size() ; i++)
        {
            Ins << doFixedPos(ti.DropStopXList[i]) << "\t";
        }
        Ins << "\n";
        Ins << "Requests.DropStopY:" << "\t" ;
        for(int i=0 ; i < ti.DropStopYList.size() ; i++)
        {
            Ins << doFixedPos(ti.DropStopYList[i]) << "\t";
        }
        Ins << "\n";
        Ins << "Requests.tau:" << "\t" ;
        for(int i=0 ; i < ti.RequestsTau.size() ; i++)
        {
            Ins << doFixedPos(ti.RequestsTau[i]) << "\t";
        }
        Ins << "\n";
        Ins << "W:" << "\t" ;
        for(int i=0 ; i < ti.W.size() ; i++)
        {
            Ins << ti.W[i] << "\t";
        }
        Ins << "\n";
        Ins << "\n";
    }
    Ins.close();
}

double GeneratorBase::doFixedPos(const double &p)
{
    double a = (p*100000.00000f)/100000.00000f;
    stringstream tmp;
    tmp << setprecision(5) << fixed << a;
    double new_val = stod(tmp.str());
    return new_val;
}

void GeneratorBase::StartGenerationOnTimeHorizons()		//start generating data ...
{
    doGenerationTrips();								//call doGenerationtrips functions
}

void GeneratorBase::doGenerationTrips()					//here we generate trips
{
    double t1 = clock();
    //Input
    const double& timehorizon = Input_TimeHorizon;			//take the Input_TimeHorizon from input file
    const int& cpNo = Input_CheckpointsNumber;				//take the number of checkpoints
    const double& slacktime = Input_SlackTime;				//take the slack time
    const int& TotalReq = Input_RequestNumber;				//take the number of requests.


    //Generation
    int nTrips = ComputeNumberOfTrips(timehorizon,slacktime);			//call a function ComputeNumberOfTrips based on slack time and time horizon
    std::cerr << __FUNCTION__ << ":nTrips :: " << nTrips << std::endl;		//we print the pre-assigned trips
    std::cerr << __FUNCTION__ << ":TotalReq :: " << TotalReq << std::endl;	//and the number of customers
    m_nbTrips = nTrips;
    m_nbRequests = TotalReq;
    if(nTrips == 0)				//wrong input!
    {
        std::cerr << __FUNCTION__ << " :: Number of trips equal 0 , please update your time horizons OR other parameters in input file"  << std::endl;
        exit(0);
    }
    if(TotalReq < nTrips)
    {
        std::cerr << __FUNCTION__ << " :: Total Request are smaller than total trips"  << std::endl;
        exit(0);
    }
    int nTripReq = TotalReq / nTrips;					//total number of requests divided by number of trips = how many requests per trip
    int diffTripReq = TotalReq % nTrips;				//the remainder will be in diffTripReq

    GeneratorTrips.resize(nTrips);						//resize the generator trips to nTrips

    for(int i=0 ; i < nTrips ; i++)
    {
        GeneratorTrips[i].id = i;						//assign a id of a trip
        GeneratorTrips[i].RequestNo = nTripReq;			//assign the requests to that trips
    }
    int cntdiff=0;
    while(diffTripReq != 0)								//?
    {
        GeneratorTrips[cntdiff].RequestNo++;
        diffTripReq--;
        cntdiff++;
        if(cntdiff >= GeneratorTrips.size())
            cntdiff=0;
    }

    for(int i=0 ; i < nTrips ; i++)
    {
        int TripID = i;
        TripsInfo& ti = GeneratorTrips[TripID];
        ti.VehicleCapacity = Input_VehiclesCapacity;		//take the capacity of vehicles from the input

        ti.StartTime = TotalStartTime;						//what do we do here?
        ti.RealStartTime = ti.StartTime;
        ti.SuspendTime = Input_SuspendTime;

        ti.MiddleNode = (Input_CheckpointsNumber - 1) / 2;		//we compute the middle checkpoint

        ti.isNew = false;										//?

        TotalStartTime += Input_StartTimeInterval;				//?

        TripsInfo ti_Pre;
        if(ti.id > 0)
            ti_Pre = GeneratorTrips[ti.id-1];					//?
        else
            ti_Pre = ti;

        GenerateCheckPoints(ti,ti_Pre);					//why do we call it here again? we already called it in ComputeNumberOfTrips function

        ComputeNonCPNumber(ti);							//we compute the number of non-checkpoint stops

        GenerateNonCheckPoints(ti);						//generate non-checkpoint nodes

        InitNodeLoad();									//initialize coordinates of each request type

        ComputeTripTeta(ti,slacktime);				//why do we call it here again?

        GenerateRequests(ti,ti_Pre);					//generate requests

        ti.NodesNo = ti.CPCoords.size() + ti.NonCPCoord.size();		//NodesNo = CP + NonCP
        ti.C = ti.CPCoords.size();

        ti.W.resize(3);
        ti.W[0] = Input_W0;
        ti.W[1] = Input_W1;
        ti.W[2] = Input_W2;

    }

    PrintRawDataToFile();
}

int GeneratorBase::ComputeNumberOfTrips(const double &timehorizon, const double &slackTime) //Generate number of trips based on TimeHorizon and SlackTime
{
    TripsInfo first_ti;										//an object from class TripsInfo
    first_ti.id = 0;
    //Vehicle Capacity
    first_ti.VehicleCapacity = Input_VehiclesCapacity;		//take the vehicle capacity

    //Trip Start Time
    first_ti.StartTime = TotalStartTime;					//start time of a trip to be 0
    first_ti.RealStartTime = first_ti.StartTime;			//real start time to be zero
    first_ti.SuspendTime = Input_SuspendTime;				//not being used anymore

    TripsInfo *ti_Pre;										//a pointer to the class
    ti_Pre = &first_ti;										//an object address to be assigned to the pointer

    GenerateCheckPoints(first_ti,*ti_Pre);					//call a function to generate checkpoints

    ComputeTripTeta(first_ti,slackTime);					//do we need this function anymore?

    first_ti.RealArrivalTime = first_ti.RealStartTime + first_ti.StopsTeta[first_ti.StopsTeta.size()-1]; //the time we reach the last checkpoint of a trip

    double maxTripTime = first_ti.StopsTeta[first_ti.StopsTeta.size()-1]; // Teta value for last checkpoints of current trip

    double dt = timehorizon - maxTripTime; // we need to know when the current trip trips end in the time horizon

    double dtf = dt / Input_StartTimeInterval;		// we need to know how frequent we want a shuttle to start a trip

    int nTrips = dtf;								//thus we can compute the number of trips.
    if(nTrips <= 0)
        nTrips = 1;

    return nTrips;
}

void GeneratorBase::ComputeTripTeta(TripsInfo &ti, const double &slacktime) //to compute the trip theta
{
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

void GeneratorBase::ComputeNonCPNumber(TripsInfo& ti)			//we compute the number of non-checkpoint nodes
{
    noPD_Req = Input_PD * Input_RequestNumber / 100;			//making the percentage of each customer
    noNPD_Req = Input_NPD * Input_RequestNumber / 100;
    noPND_Req = Input_PND * Input_RequestNumber / 100;
    noNPND_Req = Input_NPND * Input_RequestNumber / 100;
    int allReqNo = noPD_Req+noNPD_Req+noPND_Req+noNPND_Req;		//putting all of the No of requests into one
    if(allReqNo < ti.RequestNo)									//if the AllReqNo is lower than RequestNo
    {
        int diff = ti.RequestNo - allReqNo;						//compute the difference
        noPND_Req+=diff;										//increase the difference to PND customers
    }

    int TripReqNo = ti.RequestNo;								//put all of the requests to TRipReqNo
    if(noPD_Req < m_nbTrips)									//if the number of PD requests are lower than
    {															//the number of trips
        if(addedPD_Req < noPD_Req)								//why do we need this?
        {
            noPD_Req = 1;
            addedPD_Req++;
        }
        else
        {
            noPD_Req = 0;
        }
    }
    else
    {
        int d = noPD_Req / m_nbTrips;
        int f = noPD_Req % m_nbTrips;
        if(addedPD_Req < noPD_Req)
        {
            noPD_Req = d;
            addedPD_Req += d;
        }
    }
    if(noNPD_Req < m_nbTrips)
    {
        if(addedNPD_Req < noNPD_Req)
        {
            noNPD_Req = 1;
            addedNPD_Req++;
        }
        else
        {
            noNPD_Req = 0;
        }
    }
    else
    {
        int d = noNPD_Req / m_nbTrips;
        int f = noNPD_Req % m_nbTrips;
        if(addedNPD_Req < noNPD_Req)
        {
            noNPD_Req = d;
            addedNPD_Req += d;
        }
    }
    if(noPND_Req < m_nbTrips)
    {
        if(addedPND_Req < noPND_Req)
        {
            noPND_Req = 1;
            addedPND_Req++;
        }
        else
        {
            noPND_Req = 0;
        }
    }
    else
    {
        int d = noPND_Req / m_nbTrips;
        int f = noPND_Req % m_nbTrips;
        if(addedPND_Req < noPND_Req)
        {
            noPND_Req = d;
            addedPND_Req += d;
        }
    }
    if(noNPND_Req < m_nbTrips)
    {
        if(addedNPND_Req < noNPND_Req)
        {
            noNPND_Req = 1;
            addedNPND_Req++;
        }
        else
        {
            noNPND_Req = 0;
        }
    }
    else
    {
        int d = noNPND_Req / m_nbTrips;
        int f = noNPND_Req % m_nbTrips;
        if(addedNPND_Req < noNPND_Req)
        {
            noNPND_Req = d;
            addedNPND_Req += d;
        }
    }

    ti.RequestNo = noNPD_Req+noPND_Req+noNPND_Req+noPD_Req;		//all of the request numbers

    NonCPNo = noNPD_Req+noPND_Req+(noNPND_Req*2);				//number of non-checkpoint nodes
}

bool GeneratorBase::isAReturnTrip(const int &tripID)
{
    if(tripID % 2 != 0)
        return true;
    return false;
}

void GeneratorBase::InitNodeLoad()
{
    PD_RequestCoordinations.resize(0);
    PND_RequestCoordinations.resize(0);
    NPD_RequestCoordinations.resize(0);
    NPND_RequestCoordinations.resize(0);
}

void GeneratorBase::GenerateRequests(TripsInfo &ti, const TripsInfo &ti_pre, bool needtau) //generate requests
{
    ti.PickStopXList.resize(0);							//initialize the x,y coordinates list and tau
    ti.PickStopYList.resize(0);
    ti.DropStopXList.resize(0);
    ti.DropStopYList.resize(0);
    ti.RequestsTau.resize(0);

    PD_RequestCoordinations.clear();
    PND_RequestCoordinations.clear();
    NPD_RequestCoordinations.clear();
    NPND_RequestCoordinations.clear();

    PD_RequestsTau.clear();
    PND_RequestsTau.clear();
    NPD_RequestsTau.clear();
    NPND_RequestsTau.clear();

    VisitedCoord_Hybrid.clear();
    VisitedCoord_Hybrid_PND.clear();

    GeneratePDRequests(ti,ti_pre,needtau);				//generate different types of customers.
    GeneratePNDRequests(ti,ti_pre,needtau);
    GenerateNPDRequests(ti,ti_pre,needtau);
    GenerateNPNDRequests(ti,ti_pre,needtau);

    for(int i=0 ; i < PD_RequestCoordinations.size() ; i++)
    {
        Coord start = PD_RequestCoordinations[i].first;			//initialize the start and finish points
        Coord finish = PD_RequestCoordinations[i].second;
        ti.PickStopXList.push_back(start.first);				//make a list of pick-up and drop-off points
        ti.PickStopYList.push_back(start.second);
        ti.DropStopXList.push_back(finish.first);
        ti.DropStopYList.push_back(finish.second);
        ti.RequestsTau.push_back(PD_RequestsTau[i]);
    }
    for(int i=0 ; i < PND_RequestCoordinations.size() ; i++)
    {
        Coord start = PND_RequestCoordinations[i].first;
        Coord finish = PND_RequestCoordinations[i].second;
        ti.PickStopXList.push_back(start.first);
        ti.PickStopYList.push_back(start.second);
        ti.DropStopXList.push_back(finish.first);
        ti.DropStopYList.push_back(finish.second);
        ti.RequestsTau.push_back(PND_RequestsTau[i]);
    }
    for(int i=0 ; i < NPD_RequestCoordinations.size() ; i++)
    {
        Coord start = NPD_RequestCoordinations[i].first;
        Coord finish = NPD_RequestCoordinations[i].second;
        ti.PickStopXList.push_back(start.first);
        ti.PickStopYList.push_back(start.second);
        ti.DropStopXList.push_back(finish.first);
        ti.DropStopYList.push_back(finish.second);
        ti.RequestsTau.push_back(NPD_RequestsTau[i]);
    }
    for(int i=0 ; i < NPND_RequestCoordinations.size() ; i++)
    {
        Coord start = NPND_RequestCoordinations[i].first;
        Coord finish = NPND_RequestCoordinations[i].second;

        ti.PickStopXList.push_back(start.first);
        ti.PickStopYList.push_back(start.second);
        ti.DropStopXList.push_back(finish.first);
        ti.DropStopYList.push_back(finish.second);
        ti.RequestsTau.push_back(NPND_RequestsTau[i]);
    }

    std::cerr << "Request Generation Finished For Trip -> " << ti.id << std::endl;
}

void GeneratorBase::GeneratePDRequests(TripsInfo &ti, const TripsInfo &ti_pre, bool needtau)
{
    ///<! Generate PD Requests
    if(noPD_Req == 0)
        return;
    const int range_from  = 0;
    const int range_to    = ti.CPCoords.size();						//range from 0 to CPCoords

    CoordSet m_GraphNodes;
    const CoordSet& m_CPCoords = ti.CPCoords; ///< Checkpoints Coordination
    for(CoordSet::const_iterator cit = m_CPCoords.begin() ; cit != m_CPCoords.end() ; cit++) //iterate over all checkpoints
    {
        m_GraphNodes.insert(*cit);									//insert them into GraphNodes
    }

    std::vector<Coord> m_NodesList;								//initialize a vector
    for(CoordSet::const_iterator cit = m_GraphNodes.begin() ; cit != m_GraphNodes.end() ; cit++) //iterate over GraphNodes
    {
        const Coord& node = *cit;
        m_NodesList.push_back(node);							//what is the logic?
    }
    Coord MiddleCood = m_NodesList[ti.MiddleNode];				//assign the middle checkpoint to MiddleCood

    int GenSize = 0;
    if(range_to < noPD_Req*2)									//what is the logic?
    {
        bool genTerminate = false;								//?
        while(PD_RequestCoordinations.size() < noPD_Req)
        {
            GenSize = range_to;
            std::vector<int> values = ShuffleGeneration(range_from,range_to,GenSize); //generate data

            CoordSet GraphNodes;
            const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
            for(int i=0 ; i < values.size() ; i++)						//?
            {
                for(int j=i+1 ; j < values.size() ; j++)
                {
                    Coord start = *std::next(CPCoords.begin(),values[i]);
                    Coord finish = *std::next(CPCoords.begin(),values[j]);

                    if((start.first < MiddleCood.first && finish.first <= MiddleCood.first)
                        ||(start.first >= MiddleCood.first && finish.first > MiddleCood.first))
                    {
                        if(needtau)
                            ComputeTau(ti,PD_RequestsTau,start.first,start.second);				//generate tau
                        PD_RequestCoordinations.push_back(make_pair(start,finish));		//make pair of the start and finish points
                    }
                    if(PD_RequestCoordinations.size() == noPD_Req)			//if this is true
                    {
                        genTerminate = true;					//finish generating PD customers
                        break;
                    }
                }
                if(genTerminate)
                    break;
            }
            if(genTerminate)
                break;
        }
    }											//koja migim PD ha beyne checkpoint ha generate beshan?
    else
    {
        while(PD_RequestCoordinations.size() < noPD_Req)				//why do we need the following?
        {
            GenSize = noPD_Req * 2;
            std::vector<int> values = ShuffleGeneration(range_from,range_to,GenSize);

            CoordSet GraphNodes;
            const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
            int step = values.size() / 2;
            for(int i=0 ; i < step ; i++)
            {
                Coord start = *std::next(CPCoords.begin(),values[i]);
                Coord finish = *std::next(CPCoords.begin(),values[i+step]);
                if((start.first < MiddleCood.first && finish.first <= MiddleCood.first)
                    ||(start.first >= MiddleCood.first && finish.first > MiddleCood.first))
                {
                    if(needtau)
                        ComputeTau(ti,PD_RequestsTau,start.first,start.second);
                    PD_RequestCoordinations.push_back(make_pair(start,finish));
                }
            }
        }
    }
}

void GeneratorBase::GenerateNPDRequests(TripsInfo &ti, const TripsInfo &ti_pre, bool needtau)
{
    if(noNPD_Req == 0)
        return;

    CoordSet m_GraphNodes;
    const CoordSet& m_CPCoords = ti.CPCoords; ///< Checkpoints Coordination
    for(CoordSet::const_iterator cit = m_CPCoords.begin() ; cit != m_CPCoords.end() ; cit++)
    {
        m_GraphNodes.insert(*cit);
    }

    std::vector<Coord> m_NodesList;
    for(CoordSet::const_iterator cit = m_GraphNodes.begin() ; cit != m_GraphNodes.end() ; cit++)
    {
        const Coord& node = *cit;
        m_NodesList.push_back(node);
    }
    Coord MiddleCood = m_NodesList[ti.MiddleNode];

    int CPGenSize = 0;
    if(noNPD_Req > ti.CPCoords.size())
        CPGenSize = ti.CPCoords.size();
    else
        CPGenSize = noNPD_Req;

    while(NPD_RequestCoordinations.size() < noNPD_Req)
    {
        NPD_RequestCoordinations.clear();
        NPD_RequestsTau.clear();
        VisitedCoord_Hybrid.clear();
        for(vector<Coord>::iterator it = VisitedCoord_Hybrid_PND.begin() ; it != VisitedCoord_Hybrid_PND.end() ; it++)
        {
            VisitedCoord_Hybrid.push_back(*it);
        }

        std::vector<int> nonCPGen = ShuffleGeneration(0,ti.NonCPCoord.size(),noNPD_Req);
        std::vector<int> CPGen = ShuffleGeneration(0,ti.CPCoords.size(),CPGenSize);

        bool isGenerationOk = false;
        int cntError = 0;
//        while(!isGenerationOk)
//        {
//            VisitedCoord_Hybrid.clear();
//            for(vector<Coord>::iterator it = VisitedCoord_Hybrid_PND.begin() ; it != VisitedCoord_Hybrid_PND.end() ; it++)
//            {
//                VisitedCoord_Hybrid.push_back(*it);
//            }

////            std::unordered_set<int> visitedIndices;
//            for(int i = nonCPGen.size() -1 ; i >= 0 ; i--)
//            {
//                Coord noncpCoord = *std::next(ti.NonCPCoord.begin(),nonCPGen[i]);
//                bool isFounded = false;
//                for(int j=CPGen.size()-1 ; j >= 0 ; j--)
//                {
//                    Coord cpCoord = *std::next(ti.CPCoords.begin(),CPGen[j]);
//                    if(cpCoord.first > noncpCoord.first /*&& visitedIndices.find(CPGen[j]) == visitedIndices.end()*/
//                        && /*VisitedCoord_Hybrid.find(noncpCoord) == VisitedCoord_Hybrid.end()*/
//                        (!IsGeneratedCoordExist(noncpCoord)))
//                    {
//                        VisitedCoord_Hybrid.push_back(noncpCoord);
////                        visitedIndices.insert(CPGen[j]);
//                        isFounded = true;
//                        break;
//                    }
//                }
//                if(!isFounded)
//                {
//                    CPGen = ShuffleGeneration(0,ti.CPCoords.size(),CPGenSize);
//                    nonCPGen = ShuffleGeneration(0,ti.NonCPCoord.size(),noNPD_Req);
//                    break;
//                }
//            }
//            if(visitedIndices.size() == CPGen.size())
//                isGenerationOk = true;
//            else
//            {
//                cntError++;
//                if(cntError > 1000)
//                    break;
//            }
//        }
        if(1)
        {

            VisitedCoord_Hybrid.clear();
            for(vector<Coord>::iterator it = VisitedCoord_Hybrid_PND.begin() ; it != VisitedCoord_Hybrid_PND.end() ; it++)
            {
                VisitedCoord_Hybrid.push_back(*it);
            }

            std::vector<int> nonCPGen = ShuffleGeneration(0,ti.NonCPCoord.size(),ti.NonCPCoord.size());
            std::vector<int> CPGen = ShuffleGeneration(0,ti.CPCoords.size(),ti.CPCoords.size());

            std::vector<int> nonCP_Array;
            std::vector<int> CP_Array;
            while(nonCP_Array.size() != noNPD_Req)
            {
                for(int i = 0 ; i < nonCPGen.size() ; i++)
                {
                    Coord noncpCoord = *std::next(ti.NonCPCoord.begin(),nonCPGen[i]);
                    bool isFounded = false;
                    for(int j=CPGen.size()-1 ; j >= 0 ; j--)
                    {
                        Coord cpCoord = *std::next(ti.CPCoords.begin(),CPGen[j]);
                        if((cpCoord.first > noncpCoord.first) &&
                            ((noncpCoord.first < MiddleCood.first && cpCoord.first <= MiddleCood.first)
                             ||(noncpCoord.first >= MiddleCood.first && cpCoord.first > MiddleCood.first)))
                        {
                            if(!IsGeneratedCoordExist(noncpCoord))
                            {
                                VisitedCoord_Hybrid.push_back(noncpCoord);
                                nonCP_Array.push_back(nonCPGen[i]);
                                CP_Array.push_back(CPGen[j]);
                                break;
                            }

                        }
                    }
                    if(nonCP_Array.size() == noNPD_Req)
                    {
                        break;
                    }
                }
                nonCPGen.clear();
                CPGen.clear();
                nonCPGen = ShuffleGeneration(0,ti.NonCPCoord.size(),ti.NonCPCoord.size());
                CPGen = ShuffleGeneration(0,ti.CPCoords.size(),ti.CPCoords.size());
                if(nonCP_Array.size() != noNPD_Req)
                {
                    VisitedCoord_Hybrid.clear();
                    for(vector<Coord>::iterator it = VisitedCoord_Hybrid_PND.begin() ; it != VisitedCoord_Hybrid_PND.end() ; it++)
                    {
                        VisitedCoord_Hybrid.push_back(*it);
                    }
                }
            }

            VisitedCoord_Hybrid.clear();
            for(vector<Coord>::iterator it = VisitedCoord_Hybrid_PND.begin() ; it != VisitedCoord_Hybrid_PND.end() ; it++)
            {
                VisitedCoord_Hybrid.push_back(*it);
            }

            for(size_t i = 0 ; i < nonCP_Array.size() ; i++)
            {
                Coord start = *std::next(ti.NonCPCoord.begin(),nonCP_Array[i]);
                Coord finish = *std::next(ti.CPCoords.begin(),CP_Array[i]);
                if((start.first < MiddleCood.first && finish.first <= MiddleCood.first)
                    ||(start.first >= MiddleCood.first && finish.first > MiddleCood.first))
                {
                    if(needtau)
                        ComputeTau(ti,NPD_RequestsTau,start.first,start.second);
                    NPD_RequestCoordinations.push_back(make_pair(start,finish));
                    VisitedCoord_Hybrid.push_back(start);
                }
            }

        }
        else
        {
            VisitedCoord_Hybrid.clear();
            for(vector<Coord>::iterator it = VisitedCoord_Hybrid_PND.begin() ; it != VisitedCoord_Hybrid_PND.end() ; it++)
            {
                VisitedCoord_Hybrid.push_back(*it);
            }

            for(size_t i = 0 ; i < nonCPGen.size() ; i++)
            {
                Coord start = *std::next(ti.NonCPCoord.begin(),nonCPGen[i]);
                Coord finish = *std::next(ti.CPCoords.begin(),CPGen[i]);
                if((start.first < MiddleCood.first && finish.first <= MiddleCood.first)
                    ||(start.first >= MiddleCood.first && finish.first > MiddleCood.first))
                {
                    if(needtau)
                        ComputeTau(ti,NPD_RequestsTau,start.first,start.second);
                    NPD_RequestCoordinations.push_back(make_pair(start,finish));
                    VisitedCoord_Hybrid.push_back(start);
                }
            }
        }
    }
}

void GeneratorBase::GeneratePNDRequests(TripsInfo &ti, const TripsInfo &ti_pre, bool needtau)
{
    if(noPND_Req == 0)
        return;
    int CPGenSize = 0;
    if(noPND_Req > ti.CPCoords.size())
        CPGenSize = ti.CPCoords.size();
    else
        CPGenSize = noPND_Req;

    CoordSet m_GraphNodes;
    const CoordSet& m_CPCoords = ti.CPCoords; ///< Checkpoints Coordination
    for(CoordSet::const_iterator cit = m_CPCoords.begin() ; cit != m_CPCoords.end() ; cit++)
    {
        m_GraphNodes.insert(*cit);
    }

    std::vector<Coord> m_NodesList;
    for(CoordSet::const_iterator cit = m_GraphNodes.begin() ; cit != m_GraphNodes.end() ; cit++)
    {
        const Coord& node = *cit;
        m_NodesList.push_back(node);
    }
    Coord MiddleCood = m_NodesList[ti.MiddleNode];

    while(PND_RequestCoordinations.size() < noPND_Req)
    {
        std::vector<int> CPGen = ShuffleGeneration(0,ti.CPCoords.size(),CPGenSize);
        std::vector<int> nonCPGen = ShuffleGeneration(0,ti.NonCPCoord.size(),noPND_Req);
        PND_RequestCoordinations.clear();
        PND_RequestsTau.clear();
        VisitedCoord_Hybrid.clear();
        VisitedCoord_Hybrid_PND.clear();

//        bool isGenerationOk = false;
//        int cntError = 0;
//        while(!isGenerationOk)
//        {
//            std::unordered_set<int> visitedIndices;
//            std::set<Coord> visitedCoord;
//            for(int i = CPGen.size() -1 ; i>= 0 ; i--)
//            {
//                Coord cpCoord = *std::next(ti.CPCoords.begin(),CPGen[i]);
//                bool isFounded = false;
//                for(int j=nonCPGen.size()-1 ; j>= 0 ; j--)
//                {
//                    Coord noncpCoord = *std::next(ti.NonCPCoord.begin(),nonCPGen[j]);
//                    if(noncpCoord.first > cpCoord.first && (!IsGeneratedCoordExist(noncpCoord)))
//                    {
//                        visitedCoord.insert(noncpCoord);
//                        visitedIndices.insert(nonCPGen[j]);
//                        VisitedCoord_Hybrid.push_back(noncpCoord);
//                        isFounded = true;
//                        break;
//                    }
//                }
//                if(!isFounded)
//                {
//                    CPGen = ShuffleGeneration(0,ti.CPCoords.size(),CPGenSize);
//                    nonCPGen = ShuffleGeneration(0,ti.NonCPCoord.size(),noPND_Req);
//                    break;
//                }
//            }
//            if(visitedIndices.size() == nonCPGen.size())
//                isGenerationOk = true;
//            cntError++;
//            if(cntError > 1000)
//            {
//                break;
//            }
//        }
        if(1)
        {

            VisitedCoord_Hybrid.clear();
            //TODO :: Set PND Requests With Another Approach
            CPGen = ShuffleGeneration(0,ti.CPCoords.size(),ti.CPCoords.size());
            nonCPGen = ShuffleGeneration(0,ti.NonCPCoord.size(),ti.NonCPCoord.size());

            std::vector<int> nonCP_Array;
            std::vector<int> CP_Array;
            while(CP_Array.size() != noPND_Req)
            {
                for(int i = 0 ; i< CPGen.size() ; i++)
                {
                    Coord cpCoord = *std::next(ti.CPCoords.begin(),CPGen[i]);
                    bool isFounded = false;
                    for(int j=nonCPGen.size()-1 ; j>= 0 ; j--)
                    {
                        Coord noncpCoord = *std::next(ti.NonCPCoord.begin(),nonCPGen[j]);
                        if(noncpCoord.first > cpCoord.first &&
                            ((cpCoord.first < MiddleCood.first && noncpCoord.first <= MiddleCood.first)
                             ||(cpCoord.first >= MiddleCood.first && noncpCoord.first > MiddleCood.first)))
                        {
                            if(!IsGeneratedCoordExist(noncpCoord))
                            {
                                VisitedCoord_Hybrid.push_back(noncpCoord);
                                nonCP_Array.push_back(nonCPGen[j]);
                                CP_Array.push_back(CPGen[i]);
                                isFounded = true;
                                break;
                            }
                        }
                    }
                    if(CP_Array.size() == noPND_Req)
                    {
                        break;
                    }
                }


            }

            VisitedCoord_Hybrid.clear();
            for(size_t i = 0 ; i < CP_Array.size() ; i++)
            {
                Coord start = *std::next(ti.CPCoords.begin(),CP_Array[i]);
                Coord finish = *std::next(ti.NonCPCoord.begin(),nonCP_Array[i]);
                if((start.first < MiddleCood.first && finish.first <= MiddleCood.first)
                    ||(start.first >= MiddleCood.first && finish.first > MiddleCood.first))
                {
                    if(needtau)
                        ComputeTau(ti,PND_RequestsTau,start.first,start.second);
                    PND_RequestCoordinations.push_back(make_pair(start,finish));
                    VisitedCoord_Hybrid.push_back(finish);
                    VisitedCoord_Hybrid_PND.push_back(finish);
                }
            }
        }
        else
        {
            VisitedCoord_Hybrid.clear();
            for(size_t i = 0 ; i < CPGen.size() ; i++)
            {
                Coord start = *std::next(ti.CPCoords.begin(),CPGen[i]);
                Coord finish = *std::next(ti.NonCPCoord.begin(),nonCPGen[i]);
                if((start.first < MiddleCood.first && finish.first <= MiddleCood.first)
                    ||(start.first >= MiddleCood.first && finish.first > MiddleCood.first))
                {
                    if(needtau)
                        ComputeTau(ti,PND_RequestsTau,start.first,start.second);
                    PND_RequestCoordinations.push_back(make_pair(start,finish));
                    VisitedCoord_Hybrid.push_back(finish);
                    VisitedCoord_Hybrid_PND.push_back(finish);
                }
            }
        }
    }
}

void GeneratorBase::GenerateNPNDRequests(TripsInfo &ti, const TripsInfo &ti_pre, bool needtau)
{
    if(noNPND_Req == 0)
        return;
    std::vector<Coord> SaveVisitedCoord_Hybrid;
    SaveVisitedCoord_Hybrid.clear();
    for(vector<Coord>::iterator it = VisitedCoord_Hybrid.begin() ; it != VisitedCoord_Hybrid.end() ; it++)
    {
        SaveVisitedCoord_Hybrid.push_back(*it);
    }
    CoordSet m_GraphNodes;
    const CoordSet& m_CPCoords = ti.CPCoords; ///< Checkpoints Coordination
    for(CoordSet::const_iterator cit = m_CPCoords.begin() ; cit != m_CPCoords.end() ; cit++)
    {
        m_GraphNodes.insert(*cit);
    }

    std::vector<Coord> m_NodesList;
    for(CoordSet::const_iterator cit = m_GraphNodes.begin() ; cit != m_GraphNodes.end() ; cit++)
    {
        const Coord& node = *cit;
        m_NodesList.push_back(node);
    }
    Coord MiddleCood = m_NodesList[ti.MiddleNode];

    int GenSize = 0;
    if(ti.NonCPCoord.size() < noNPND_Req*2)
    {
        GenSize = ti.NonCPCoord.size();
        bool genTerminate = false;
        while(NPND_RequestCoordinations.size() < noNPND_Req)
        {
            GenSize = ti.NonCPCoord.size();
            std::vector<int> values = ShuffleGeneration(0,ti.NonCPCoord.size(),GenSize);
            VisitedCoord_Hybrid.clear();
            for(vector<Coord>::iterator it = SaveVisitedCoord_Hybrid.begin() ; it != SaveVisitedCoord_Hybrid.end() ; it++)
            {
                VisitedCoord_Hybrid.push_back(*it);
            }
            CoordSet GraphNodes;
            const CoordSet& nonCPCoords = ti.NonCPCoord; ///< Non-Checkpoints Coordination
            for(int i=0 ; i < values.size() ; i++)
            {
                for(int j=i+1 ; j < values.size() ; j++)
                {
                    Coord start = *std::next(nonCPCoords.begin(),values[i]);
                    Coord finish = *std::next(nonCPCoords.begin(),values[j]);


                    if((start.first < finish.first)
                        && (!IsGeneratedCoordExist(start))
                        && (!IsGeneratedCoordExist(finish)))
                    {
                        NPND_RequestCoordinations.push_back(make_pair(start,finish));
                        if(needtau)
                            ComputeTau(ti,NPND_RequestsTau,start.first,start.second);
                        VisitedCoord_Hybrid.push_back(start);
                        VisitedCoord_Hybrid.push_back(finish);
                    }
                    if(NPND_RequestCoordinations.size() == noNPND_Req)
                    {
                        genTerminate = true;
                        break;
                    }
                }
                if(genTerminate)
                    break;
            }
            if(genTerminate)
                break;
        }
    }
    else
    {
        while(NPND_RequestCoordinations.size() < noNPND_Req)
        {
            VisitedCoord_Hybrid.clear();
            for(vector<Coord>::iterator it = SaveVisitedCoord_Hybrid.begin() ; it != SaveVisitedCoord_Hybrid.end() ; it++)
            {
                VisitedCoord_Hybrid.push_back(*it);
            }
            GenSize = noNPND_Req * 2;
            bool is_generationFinished = false;
            std::vector<int> values = ShuffleGeneration(0,ti.NonCPCoord.size(),/*noNPND_Req*2*/ti.NonCPCoord.size());
            CoordSet GraphNodes;
            const CoordSet& nonCPCoords = ti.NonCPCoord; ///< Non-Checkpoints Coordination
            int step = values.size() / 2;
            for(int i=0 ; i < values.size() ; i++)
            {
                for(int j=i+1 ; j < values.size() ; j++)
                {
                    int va = values[i];
                    int vb = values[j];
                    Coord start = *std::next(nonCPCoords.begin(),va);
                    Coord finish = *std::next(nonCPCoords.begin(),vb);


                    if((start.first < finish.first)
                        && (!IsGeneratedCoordExist(start))
                        && (!IsGeneratedCoordExist(finish)))
                    {
                        VisitedCoord_Hybrid.push_back(start);
                        VisitedCoord_Hybrid.push_back(finish);
                        NPND_RequestCoordinations.push_back(make_pair(start,finish));
                        if(needtau)
                            ComputeTau(ti,NPND_RequestsTau,start.first,start.second);
                        if(NPND_RequestCoordinations.size() >= noNPND_Req)
                        {
                            is_generationFinished = true;
                            break;
                        }

                    }
                }
                if(is_generationFinished)
                    break;
            }
            if(is_generationFinished)
                break;
        }
    }
}

bool GeneratorBase::IsGeneratedCoordExist(const Coord &pos)
{
    for(Index i=0 ; i < VisitedCoord_Hybrid.size() ; i++)
    {
        const Coord& vc = VisitedCoord_Hybrid[i];
        if((vc.first) == (pos.first)
                && (vc.second) == (pos.second))
        {
            return true;
        }
    }
    return false;
}

void GeneratorBase::ComputeTau(TripsInfo &ti, vector<double> &reqtau, const double &start_x, const double &start_y) //generate tau for each request
{
    const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
    const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination
    CoordSet NodesCoordList;
    double Distance = 0;
    for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++) //iterate over checkpoints coordinates
    {
        NodesCoordList.insert(*cit);				//insert into the nodesCoordList
    }
    for(CoordSet::const_iterator cit = NonCPCoord.begin() ; cit != NonCPCoord.end() ; cit++)//iterate over non-checkpoints coordinate
    {
        NodesCoordList.insert(*cit);				//insert into nodescoordlist
    }
    std::vector<Coord> NodesList;
    for(CoordSet::const_iterator cit = NodesCoordList.begin() ; cit != NodesCoordList.end() ; cit++) //iterate over all nodes
    {
        const Coord& node = *cit;
        NodesList.push_back(node);				//push back into nodesList - why ?
    }
    for(size_t i=0 ; i < NodesList.size()  ; i++)
    {
        const Coord& node = NodesList[i];
        double x1 = node.first;
        double y1 = node.second;
        if(x1 == start_x && y1 == start_y && i == 0)		//if the coordinates are the start position (origin depot)
        {
            reqtau.push_back(0);					//the tau is 0
            break;
        }

        if(x1 == start_x && y1 == start_y)
        {
            double range_from  = ti.StopsTeta[0];				//initialize
            double range_to = 0.0;								//initialize
            if(ti.CPCoords.find(node) == ti.CPCoords.end())		//checks whether the element node was not found in CPCoords
                range_to    = ti.HiddenTeta[i];					//if it was true, it means that the node was not in the set
            else
                range_to    = ti.StopsTeta[i];					//else, put the stopsteta

            std::vector<double> values = ShuffleGenerationUnSortedReal(range_from,range_to,1); //generate

            if(values.empty())						//if the vector values has no element
                reqtau.push_back(0);		//push back tau as 0
            else
                reqtau.push_back(values[0]);	//the first element of the vector value is assigned to Tau

            break;
        }
    }
}

void GeneratorBase::GenerateCheckPoints(TripsInfo &ti,const TripsInfo& ti_pre)		//to generate checkpoint nodes
{
    double sa_x = Input_ServiceAreaX;								//take the service area x
    double noCP = Input_CheckpointsNumber;							//take the number of checkpoints
    double step = sa_x / noCP;										//divide the two above

    double x = 0;
    double init_X = 0;

    if(ti_pre.CPCoords.empty())									//the first checkpoint to be zero in x coordinates
    {
        x = 0;
    }
    else
    {
        x = (*--ti_pre.CPCoords.end()).first + (sa_x / (noCP-1));		//assign x to be = size of service area/number of checkpoints-1
    }
    init_X = x;												//assign x to init_X

    int y = 0;

    vector<Coord> localCoord;								//local coordinates of checkpoints
    localCoord.resize(noCP);								//resize it to be as the size of NumberOfCheckpoints
    pair<double,double> firstPos = make_pair(x,y);			//first position is the 0,0
    pair<double,double> LastPos = make_pair(sa_x+init_X,y); // last position is as the size of service area x, y=0
    localCoord[0] = firstPos;								//assign it to localCoord vector
    localCoord[noCP-1] = LastPos;							//assign it to localCoord vector

    noCP -=1;												//remove one checkpoints from the set (to have manage over them) - removing the last checkpoint that we assigned
    step = (double)sa_x / (double)noCP;						//grow step for the coordination of checkpoints
    x+=step;												//step will be assigned to x
    ti.CPCoords.insert(firstPos);							//we insert the origin to the CPCoords
    for(int i=1 ; i < localCoord.size() -1 ; i++)			//starting from 1 beucase 0 assigned to the origin checkpoint
    {
        pair<double,double> pos = make_pair(x,y);			//make a pair of the x and y coordinates
        ti.CPCoords.insert(pos);							//insert it to CPCoords
        x+=step;											//add the step to x
    }
    ti.CPCoords.insert(LastPos);							//we insert the last position
}

void GeneratorBase::GenerateNonCheckPoints(TripsInfo &ti)	//to generate the non-checkpoint nodes
{
    int sa_x = (*--ti.CPCoords.end()).first;				//what are we doing here?
    int sa_y = Input_ServiceAreaY;							//take the y from input
    bool isFirstTrip = false;								//?
    int rangeFrom = 0;
    if(ti.id == 0)											//if the id=0
    {
        isFirstTrip = true;									//this is the first trip
        rangeFrom = 0;										//starting from 0
    }
    if(!isFirstTrip)										//if it is not a first trip then
    {
        rangeFrom = (*ti.CPCoords.begin()).first;			//the range starts from the first CPCoords
        rangeFrom++;										//then we increase the range
    }

    std::vector<double> X_NonCP;
    bool is_XnonCP_OK = false;								//it is initialize is_XNonCP to be false
    while(!is_XnonCP_OK)									//while it is not false
    {
        X_NonCP.clear();
        X_NonCP = ShuffleGenerationUnSortedReal(rangeFrom,sa_x,NonCPNo); //Go to this function to generate from,to,number
        bool isOK = true;
        for(int i=0 ; i < X_NonCP.size() ; i++)					//where did we fill X_NonCP?
        {
//            Coord searchCoord = make_pair(X_NonCP[i],0);
            for(CoordSet::const_iterator it = ti.CPCoords.begin() ; it != ti.CPCoords.end() ; it++) //iterate over CPCoords
            {
                const Coord& cp = *it;
                if(cp.first == X_NonCP[i])
                {
                    isOK = false;
                    break;
                }
            }
            if(!isOK)
                break;
        }
        if(isOK)
            is_XnonCP_OK = true;
    }


    std::sort(X_NonCP.begin(),X_NonCP.end());					//we sort the X_NonCP

    std::vector<double> Y_NonCP = ShuffleGenerationUnSortedReal(-(sa_y/2.0),sa_y/2.0,NonCPNo); //generate the y coordinate

    for(int i=0 ; i < Y_NonCP.size() ; i++)
    {
        pair<double,double> pos = make_pair(X_NonCP[i],Y_NonCP[i]);					//make a pair of x,y coordinates
        ti.NonCPCoord.insert(pos);													//insert into NonCPCoord
    }
}

std::vector<int> GeneratorBase::ShuffleGeneration(int from, int to, int no)
{
    const int range_from  = from;
    const int range_to    = to-1;

    std::vector<int> values;
    if(range_to < no-1)
        return values;
    std::srand(std::time(nullptr));

    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_int_distribution<int> distr(range_from, range_to);
    std::set<int> setGeneration;
    while(setGeneration.size() != no)
    {
        setGeneration.insert(distr(eng));
    }
    for(std::set<int>::const_iterator it = setGeneration.begin() ; it != setGeneration.end() ; it++)
    {
        values.push_back(*it);
    }
    //Sorting...
    std::sort(values.begin(),values.end());			//sorts the data
    return values;
}

std::vector<int> GeneratorBase::ShuffleGenerationUnSorted(int from, int to, int no)
{
    const int range_from  = from;
    const int range_to    = to-1;

    std::vector<int> values;
    if(range_to < no-1)
        return values;
    std::srand(std::time(nullptr));

    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_int_distribution<int> distr(range_from, range_to);
    for (int n = 0; n < no; ++n)
    {
        values.push_back(distr(eng));
    }

    return values;
}

std::vector<double> GeneratorBase::ShuffleGenerationReal(double from, double to, int no)
{
    const double range_from  = from;
    const double range_to    = to-1;

    std::vector<double> values;
    std::srand(std::time(nullptr));

    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(range_from, range_to);
    std::set<double> setGeneration;
    while(setGeneration.size() != no)
    {
        setGeneration.insert(distr(eng));
    }
    for(std::set<double>::const_iterator it = setGeneration.begin() ; it != setGeneration.end() ; it++)
    {
        values.push_back(*it);
    }
    //Sorting...
    std::sort(values.begin(),values.end());
    return values;
}

std::vector<double> GeneratorBase::ShuffleGenerationUnSortedReal(double from, double to, int no) //generate non-checkpoints
{
    const double range_from  = from;
    const double range_to    = to-std::numeric_limits<double>::epsilon();	//?

    std::vector<double> values;
    std::srand(std::time(nullptr));											//?

    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(range_from, range_to);
    std::set<double> setGeneration;
    while(setGeneration.size() != no)										//as long as generation is not equal to the number of non-checkpoints
    {
        setGeneration.insert(distr(eng));					//generate and insert them into setGeneration
    }
    for(std::set<double>::const_iterator it = setGeneration.begin() ; it != setGeneration.end() ; it++) //iterate all over the generated data
    {
        values.push_back(*it);			//push back into values. - ?1340 and 1344?
    }

    std::shuffle(values.begin(),values.end(),std::mt19937{});		//we shuffle over them

    return values;
}

