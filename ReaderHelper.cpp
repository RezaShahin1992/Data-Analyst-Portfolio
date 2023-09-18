#include "ReaderHelper.h"
#include "fstream"
#include <sstream>
#include <iomanip>
#define ReservedTripsNo 12

ReaderHelper::ReaderHelper()
{

}

void ReaderHelper::ReadData(int index0, SolverBase *base)    //
{															//
    PreProcessInstanceFile(index0);							// this function open the ins-0 file and with respect to the ins file would tell us how many trips do we have and each trip contain how many requests
    cout << __FUNCTION__ << ": Read Infrastructure Started " <<endl;
    ReadInfrastructure(index0);								//read the infrastructure file
    cout << __FUNCTION__ << ": Read Infrastructure Finished " <<endl;
    cout << __FUNCTION__ << ": Read Instance Started " <<endl;
    ReadInstance(index0);									//read the instance file
    cout << __FUNCTION__ << ": Read Instance Finished " <<endl;

    if(base != nullptr)													//
    {
        base->SetTripsInfoList(TripsInfoList);							//send some of the information that we read from input files to solverBase
        base->Set_CP_Number(m_CPno);									//
        base->Set_SA_X(m_SAX);
        base->Set_SA_Y(m_SAY);
        base->Set_VehicleCap(m_VehicleCap);
        base->Set_Interval(m_Interval);
        base->Set_PD_Per(m_PD);
        base->Set_PND_Per(m_PND);
        base->Set_NPD_Per(m_NPD);
        base->Set_NPND_Per(m_NPND);
        base->Set_StopB(m_StopB);
        base->Set_SlackTime(m_SlackTime);
        base->Set_W0(m_W0);
        base->Set_W1(m_W1);
        base->Set_W2(m_W2);
        base->Set_XWS(X_WS);
        base->Set_ZWS(Z_WS);
    }
}

void ReaderHelper::ReadOutputFile(int idx)
{
    string StrInfName = AddOutput +"OutputHeu-"+ to_string(idx) + ".txt";	// yek file e ro misaze ba peyvaste kardan yek string addinput ba meghdar index0 ke yek file txt besaze va dar stringinfname zakhire mishe
    fstream outFile;													//ye shey ro declare mikone
    outFile.open(StrInfName,ios_base::in | ios_base::out);				//file voroodi va khorooji ro baz mikone

    if(outFile.is_open())
    {
        int cnt = 0;
        while(!outFile.eof())
        {
            string line;
            getline(outFile,line);
            if(cnt < 3)
            {
                cnt++;
                continue;
            }

            std::vector<std::string> nameSplit;
            split(line,"=",nameSplit);
            std::vector<std::string> indicesList;
            split(nameSplit[0],"[",indicesList);
            if(indicesList[0] == "X")
            {
                for(int i=1 ; i < indicesList.size() ; i++)
                {
                    string &data = indicesList[i];
                    data.erase(--data.end());
                }
                Coord path = make_pair(stoi(indicesList[1]),stoi(indicesList[2]));
                int v = stoi(indicesList[3]);
                Xmap[path] = v;
            }
            if(indicesList[0] == "t")
            {
                string &data = indicesList[1];
                data.erase(--data.end());
                int idx = stoi(indicesList[1]);
                double value = stod(nameSplit[1]);
                Tmap[idx] = value;
            }

            if(indicesList[0] == "tbar")
            {
                string &data = indicesList[1];
                data.erase(--data.end());
                int idx = stoi(indicesList[1]);
                double value = stod(nameSplit[1]);
                TBarmap[idx] = value;
            }

            if(indicesList[0] == "p")
            {
                string &data = indicesList[1];
                data.erase(--data.end());
                int idx = stoi(indicesList[1]);
                double value = stod(nameSplit[1]);
                Pmap[idx] = value;
            }

            if(indicesList[0] == "d")
            {
                string &data = indicesList[1];
                data.erase(--data.end());
                int idx = stoi(indicesList[1]);
                double value = stod(nameSplit[1]);
                Dmap[idx] = value;
            }

            if(indicesList[0] == "Z")
            {
                for(int i=1 ; i < indicesList.size() ; i++)
                {
                    string &data = indicesList[i];
                    data.erase(--data.end());
                }
                Coord path = make_pair(stoi(indicesList[1]),stoi(indicesList[2]));
                int v = stoi(indicesList[3]);
                Zmap[path] = v;
            }

            if(indicesList[0] == "q")
            {
                string &data = indicesList[1];
                data.erase(--data.end());
                int idx = stoi(indicesList[1]);
                double value = stod(nameSplit[1]);
                qmap[idx] = value;
            }

            if(indicesList[0] == "Q")
            {
                for(int i=1 ; i < indicesList.size() ; i++)
                {
                    string &data = indicesList[i];
                    data.erase(--data.end());
                }
                Coord path = make_pair(stoi(indicesList[1]),stoi(indicesList[2]));
                int idx = stoi(indicesList[1]);
                double value = stod(nameSplit[1]);
                Qmap[path] = value;
            }

        }
    }
}

void ReaderHelper::SendOutputDataToBase(SolverBase *base)
{
    base->Xmap = Xmap;
    base->Tmap = Tmap;
    base->TBarmap = TBarmap;
    base->Pmap = Pmap;
    base->Dmap = Dmap;
    base->Zmap = Zmap;
    base->qmap = qmap;
    base->Qmap = Qmap;
}

double ReaderHelper::doFixedPos(const double &p)
{
    double a = /*std::ceil*/(p*100000.00000f)/100000.00000f;
    stringstream tmp;
    tmp << setprecision(5) << fixed << a;
    double new_val = stod(tmp.str());
    return new_val;
}

int ReaderHelper::Read_Address(int &it)
{
    AddInput = /*Address + "\\INPUT\\"*/"INPUT/";
    AddOutput = /*Address + "\\Outputs"*/"OUTPUT/";

    bool isOK_Inf = IsExists(AddInput+"Inf-"+to_string(it)+".txt");
    bool isOK_Ins = IsExists(AddInput+"Ins-"+to_string(it)+".txt");

    return (isOK_Inf && isOK_Ins);
}

int ReaderHelper::Read_Address_Inf(int &it)
{
    AddInput  = "INPUT/";
    AddOutput = "OUTPUT/";

    bool isOK_Inf = IsExists(AddInput+"Inf-"+to_string(it)+".txt");

    return (isOK_Inf);
}

int ReaderHelper::Read_Address_Ins(int &it)
{
    AddInput  = "INPUT/";
    AddOutput = "OUTPUT/";

    bool isOK_Ins = IsExists(AddInput+"Ins-"+to_string(it)+".txt");

    return (isOK_Ins);
}

void ReaderHelper::ClearData()
{
    TripsNo = 0;
    LocalTripsNo = 0;
    LocalTotalReqNo = 0;
    nbPoints = 0;
    VehiclesNo = 0;
    TC0 = TC = 0;
    delta.clear();
    RequestsTypeIndex.clear();
    VehiclesSet.clear();
    VehicleCapList.clear();
    GammaList.clear();
    TripsInfoList.resize(0);
}

void ReaderHelper::split(string str, string splitBy, std::vector<string> &tokens)
{
    /* Store the original string in the array, so we can loop the rest
     * of the algorithm. */
    tokens.push_back(str);

    // Store the split index in a 'size_t' (unsigned integer) type.
    size_t splitAt;
    // Store the size of what we're splicing out.
    size_t splitLen = splitBy.size();
    // Create a string for temporarily storing the fragment we're processing.
    std::string frag;
    // Loop infinitely - break is internal.
    while(true)
    {
        /* Store the last string in the vector, which is the only logical
         * candidate for processing. */
        frag = tokens.back();
        /* The index where the split is. */
        splitAt = frag.find(splitBy);
        // If we didn't find a new split point...
        if(splitAt == std::string::npos)
        {
            // Break the loop and (implicitly) return.
            break;
        }
        /* Put everything from the left side of the split where the string
         * being processed used to be. */
        tokens.back() = frag.substr(0, splitAt);
        /* Push everything from the right side of the split to the next empty
         * index in the vector. */
        tokens.push_back(frag.substr(splitAt+splitLen, frag.size()-(splitAt+splitLen)));
    }
}

void ReaderHelper::ReadInfrastructure(const int& index0)
{
    cout << __FUNCTION__ << ": Assign InfName " <<endl;
    const char* InfName = new char;
    string StrInfName = AddInput + "Inf-"+to_string(index0) + ".txt"; //giving the address and the name of the text file
    InfName = StrInfName.c_str();
    cout << "Inf-" << index0 << ".txt" << "\t";


    Infrastructure = fopen(InfName, "r");						//open the text file
    if (Infrastructure == NULL) {								//not found
        cout << "Cannot open Infrastructure!!" << endl;
        exit(23);
    }
    else
    {
        cout << "Infrastructure file opened successfully :D" << endl;
    }

    fscanf(Infrastructure, "%s", trash);							//each of the parameters are being read
    fscanf(Infrastructure, "%d", &LocalTripsNo);

    fscanf(Infrastructure, "%s", trash);
    fscanf(Infrastructure, "%lf", &m_SlackTime);

    fscanf(Infrastructure, "%s", trash);
    fscanf(Infrastructure, "%lf", &m_StopB);

    fscanf(Infrastructure, "%s", trash);
    fscanf(Infrastructure, "%d", &m_SAX);
    fscanf(Infrastructure, "%s", trash);
    fscanf(Infrastructure, "%d", &m_SAY);
    fscanf(Infrastructure, "%s", trash);

    fscanf(Infrastructure, "%lf", &m_Interval);
    fscanf(Infrastructure, "%s", trash);
    fscanf(Infrastructure, "%d", &m_PD);
    fscanf(Infrastructure, "%s", trash);
    fscanf(Infrastructure, "%d", &m_PND);
    fscanf(Infrastructure, "%s", trash);
    fscanf(Infrastructure, "%d", &m_NPD);
    fscanf(Infrastructure, "%s", trash);
    fscanf(Infrastructure, "%d", &m_NPND);

    fscanf(Infrastructure, "%s", trash);
    fscanf(Infrastructure, "%lf", &m_TripsGenerationTime);

    bool isTripsNotEqual = false;
    int saveTripsNo = TripsNo;
    int saveLocalTripsNo = LocalTripsNo;
    if(TripsNo > LocalTripsNo) // //if the number read trips that we read is less than the trips that we read from ins, tell us that the number of trips are not eqial in ins and inf files
    {
        isTripsNotEqual = true;
    }
    else
    {
        LocalTripsNo = TripsNo;
    }


    for(int j=0 ; j <  LocalTripsNo ; j++)					//a loop over all of the trips to read information regarding each trip
    {
        TripsInfo& ti = TripsInfoList[j];
        ti.id = j;
        fscanf(Infrastructure, "%s", trash);
        fscanf(Infrastructure, "%d", &ti.isNew);

        fscanf(Infrastructure, "%s", trash);
        fscanf(Infrastructure, "%d", &ti.VehicleIndex);

        fscanf(Infrastructure, "%s", trash);
        fscanf(Infrastructure, "%d", &ti.VehicleCapacity);
//        ti.VehicleCapacity += 200;
        m_VehicleCap = ti.VehicleCapacity;

        fscanf(Infrastructure, "%s", trash);
        fscanf(Infrastructure, "%lf", &ti.StartTime);
        ti.RealStartTime = ti.StartTime;

        fscanf(Infrastructure, "%s", trash);
        fscanf(Infrastructure, "%d", &ti.SuspendTime);

        fscanf(Infrastructure, "%s", trash);
        int emptyReq;
//        fscanf(Infrastructure, "%d", &ti.RequestNo);
        fscanf(Infrastructure, "%d", &emptyReq);
        //cout << RequestNo << "\t";

        fscanf(Infrastructure, "%s", trash);
        fscanf(Infrastructure, "%d", &ti.NodesNo);
        //cout << NodesNo << "\t";

        fscanf(Infrastructure, "%s", trash);
        fscanf(Infrastructure, "%d", &ti.C);
        m_CPno = ti.C;

        ti.MiddleNode = (ti.C - 1) / 2;				//the middle checkpoint

        int Temp;
        double Temp2;

        /// List of Nodes
        NodeCreation(TripsInfoList[j]);	// read the checkpoints coordinates.

        ////////StopsType
        fscanf(Infrastructure, "%s", trash);
        for(int i=0; i< ti.NodesNo; i++){
            fscanf(Infrastructure, "%d", &Temp);
            ti.StopsType.push_back(Temp);
        }

        ////////StopsTeta										//we are computing the theta in the code now and not using this anymore
        fscanf(Infrastructure, "%s", trash);
        for (int i = 0; i < ti.NodesNo; i++) {
            fscanf(Infrastructure, "%lf", &Temp2);
            ti.StopsTeta.push_back(Temp2);
        }

        ////////StopsB
        fscanf(Infrastructure, "%s", trash);
        for (int i = 0; i < ti.NodesNo; i++) {
            fscanf(Infrastructure, "%lf", &Temp2);
            ti.StopsB.push_back(Temp2);
        }
        VehiclesSet.insert(ti.VehicleIndex);
    }

    VehiclesNo = VehiclesSet.size();

    fclose(Infrastructure);							//close the text file


    if(isTripsNotEqual) // inf trips number smaller than ins trips number
    {
        IncreaseTrips(saveLocalTripsNo,saveTripsNo);
    }

    if(ReservedTripsNo > TripsNo)
    {
        AddReservedTrips(TripsNo,ReservedTripsNo);
    }
}

void ReaderHelper::ReadInstance(const int &index0)
{
    int Temp;
    float Temp2;
    const char* InsName = new char;
    string StrInfName = AddInput + "Ins-"+to_string(index0) + ".txt"; //giving the address and name of text file
    InsName = StrInfName.c_str();
    cout << "Ins-" << index0 << ".txt" << endl;


    Instance = fopen(InsName, "r");					///open the text file
    if (Instance == NULL) {
        cout << "Cannot open Instance!!" << endl;
        exit(23);
    }

    for(int i=0 ; i < TripsNo ; i++)					//loop over the trips
    {
        TripsInfo& ti = TripsInfoList[i];
        RequestCreation(TripsInfoList[i]);			//read the pick-up and drop-off of origin and destination of each customer

        ////////RequestsTau
        fscanf(Instance, "%s", trash);
        for (int k = 0; k < ti.RequestNo; k++) {
            fscanf(Instance, "%f", &Temp2);
            ti.RequestsTau.push_back(/*Temp2*/0.0);
        }

        GetNonCPFromInsFile(ti);			//the non-checkpoint nodes are taken from the ins text file
        /// List of Nodes

        ti.NonCPCoord.clear();
        for(set<Coord>::const_iterator j=ti.TripPickDropSetNonCP.begin() ; j!=ti.TripPickDropSetNonCP.end() ; j++)
        {					//a loop over non-checkpoints and put them into the non-checkpoint coordination set
            ti.NonCPCoord.insert(*j);
        }

        ////////W
        fscanf(Instance, "%s", trash);
        for (int i = 0; i < 3; i++) {
            fscanf(Instance, "%f", &Temp2);
            ti.W.push_back(Temp2);
        }
        m_W0 = ti.W[0];			//we read the weights of the objective functions in the ins
        m_W1 = ti.W[1];
        m_W2 = ti.W[2];
    }

    if(ReservedTripsNo > TripsNo)		//if the reserved trips are more than the trips number, the number of trips will be equal to the number of reserved trips
        TripsNo = ReservedTripsNo;
//    if(ReservedTripsNo > TripsNo)
    {
        for(int tr = 0 ; tr < TripsNo ; tr++)
        {
            TripsInfo& ti = TripsInfoList[tr];
            set<int> LocalHyb;
            vector<int> LocalHybList;
            for(int q=ti.id ; q < TripsNo ; q++)	// creation HYBR(k)
            {
                int currentID = TripsInfoList[q].id;	//
                LocalHybList.push_back(currentID);
            }
            for(int d=0 ; d < LocalHybList.size() ; d++)
            {
                ti.m_PossibleTripsForRequests.push_back(LocalHybList[d]);	// Filling HYBR(k)
            }
        }
        for(int tr = 0 ; tr < TripsNo ; tr++)		//From here this is old till 1431
        {
            TripsInfo& ti = TripsInfoList[tr];
            ti.W.resize(3);
            ti.W[0] = m_W0;
            ti.W[1] = m_W1;
            ti.W[2] = m_W2;
        }
    }

    ////////////////
    fclose(Instance);
}

void ReaderHelper::IncreaseTrips(const int &currentTripsNo, const int &targetTripsNo)
{

    for(int i = currentTripsNo ; i < targetTripsNo ; i++)
    {
        TripsInfo& ti = TripsInfoList[i];
        const TripsInfo& ti_last = TripsInfoList[i-1];
        ti = ti_last;

        double step = (double)m_SAX / (double)m_CPno;
        vector<Coord> NodeList;
        double lastX = (--ti_last.CPCoords.end())->first;
        double firstX = (ti_last.CPCoords.begin())->first;
        for(CoordSet::iterator it = ti.CPCoords.begin() ; it != ti.CPCoords.end() ; it++)
        {
            Coord cpcurrent = *it;
            cpcurrent.first += (step + lastX - firstX);
            NodeList.push_back(cpcurrent);
        }
        ti.CPCoords.clear();
        for(Index j=0 ; j < NodeList.size() ; j++)
        {
            ti.CPCoords.insert(NodeList[j]);
        }

        ti.VehicleCapacity = m_VehicleCap;
        ti.id = ti_last.id+1;
        ti.C = ti_last.CPCoords.size();
        ti.NodesNo = ti.C;

        ti.MiddleNode = (ti.C - 1) / 2;

        ti.NonCPCoord.clear();
        ti.StopsB.clear();
        ti.StopsTeta.clear();
        ti.StopsType.clear();
        ti.RequestNo = 0;
        ti.RequestsTau.clear();
        ti.RequestsType.clear();
        ti.PickStopXList.clear();
        ti.PickStopYList.clear();
        ti.DropStopXList.clear();
        ti.DropStopYList.clear();

        ti.StartTime = ti_last.StartTime + m_Interval;
        ti.RealStartTime = ti.StartTime;
    }
}

void ReaderHelper::AddReservedTrips(const int &currentTripsNo, const int &targetTripsNo)
{
    for(int i = currentTripsNo ; i < targetTripsNo ; i++)
    {
        TripsInfo ti;
        const TripsInfo& ti_last = TripsInfoList[i-1];

        double step = (double)m_SAX / (double)(m_CPno-1);
        vector<Coord> NodeList;
        double lastX = (--ti_last.CPCoords.end())->first;
        double firstX = (ti_last.CPCoords.begin())->first;
        for(CoordSet::iterator it = ti_last.CPCoords.begin() ; it != ti_last.CPCoords.end() ; it++)
        {
            Coord cpcurrent = *it;
            cpcurrent.first = doFixedPos((step + lastX - firstX) + cpcurrent.first);
            NodeList.push_back(cpcurrent);
        }
        ti.CPCoords.clear();
        for(Index j=0 ; j < NodeList.size() ; j++)
        {
            ti.CPCoords.insert(NodeList[j]);
        }

        ti.VehicleCapacity = m_VehicleCap;
        ti.id = ti_last.id+1;
        ti.C = ti_last.CPCoords.size();
        ti.NodesNo = ti.C;
        ti.MiddleNode = (ti.C - 1) / 2;

        ti.NonCPCoord.clear();
        ti.StopsB.clear();
        ti.StopsTeta.clear();
        ti.StopsType.clear();
        ti.RequestNo = 0;
        ti.RequestsTau.clear();
        ti.RequestsType.clear();
        ti.PickStopXList.clear();
        ti.PickStopYList.clear();
        ti.DropStopXList.clear();
        ti.DropStopYList.clear();

        ti.StartTime = ti_last.StartTime + m_Interval;
        ti.RealStartTime = ti.StartTime;

        TripsInfoList.push_back(ti);
    }

}

void ReaderHelper::PreProcessInstanceFile(const int &index0)           //
{																	//
    string StrInfName = AddInput + "Ins-"+to_string(index0) + ".txt";	//  we read the name of the file and address
    fstream insFile;													//
    insFile.open(StrInfName,ios_base::in | ios_base::out);				// here I open the Input text file
    int cntTrip = 0;													//
    int cntReq = 0;
    std::map<int,int> TripReqMap;										//
    if(insFile.is_open())												//if the text file is opened and we face no problem
    {
        while(!insFile.eof())											//while you have not yet reached the end of the text file
        {
            string line;
            getline(insFile,line);										//read the current line
            std::vector<std::string> nameSplit;							//
            split(line,":",nameSplit);									// the read line is separated with ":" and each item is saved
            if(nameSplit.size() == 2 && nameSplit[0] == "Requests.PickStopX") //  we check to see if there exists two parts from the above split ":", showing that each ins file starts with the Requests.PickStopX
            {																	//
                std::vector<std::string> itemSplit;								//
                split(nameSplit[1],"\t",itemSplit);								// now split using "\t" and save it in itemSplit because in the ins file we have numbers separated using \t
                itemSplit.erase(itemSplit.begin());								//
                itemSplit.erase(--itemSplit.end());								// remove the first \t and last \t (because when I generate instances there exist \t in the beginning and ending of the line)
                TripsNo++;														//
                TripReqMap[cntTrip] = itemSplit.size();							// using itemSplit.size we know the number of requests of this trip
                LocalTotalReqNo += TripReqMap[cntTrip];
            }
            if(nameSplit.size() == 2 && nameSplit[0] == "Requests.tau")         // if the read string that was separated using ":" if it was equal to 2 and the first part was the request tau we realized that all of the information of the trip is read
            {
                cntTrip++;														// thus, we add one to the counter of number of trips
            }
        }
        for(int i =0 ; i < TripsNo ; i++)										// when we finish reading the file now I know how many trips I do have and now make the data structure of tripsinfo based on the number of trips
        {
            TripsInfo ti;
            TripsInfoList.push_back(ti);										//
        }
        for(std::map<int,int>::const_iterator it = TripReqMap.begin() ; it != TripReqMap.end() ; it++) //now that we have the trips, tell me how many requests are assigned in each trip
        {
            TripsInfoList[it->first].RequestNo = it->second;					//here we fill the map of trip index and number of requests
        }
    }
    else
    {
        std::cerr << __FUNCTION__ << ":insFile Not Found" << std::endl;	//not found!
    }
    insFile.flush();
    insFile.close();
}

void ReaderHelper::NodeCreation(TripsInfo &ti)
{
    /// Set of checkpoints
    fscanf(Infrastructure, "%s", trash);
    for(int i = 0 ; i < ti.C ; i++)			//loop over the number of checkpoints
    {
        double x;
        double y;
        fscanf(Infrastructure, "%lf", &x);	//read x and y of each checkpoint
        fscanf(Infrastructure, "%lf", &y);
        pair<double,double> coord = make_pair(x,y);	//making pair of x and y
        ti.CPCoords.insert(coord);				//a set of checkpoints are being filled
    }
}

void ReaderHelper::RequestCreation(TripsInfo& ti)
{
    fscanf(Instance, "%s", trash);
    for(int i = 0 ; i < ti.RequestNo ; i++)
    {
        double px;
        fscanf(Instance, "%lf", &px);

        ti.PickStopXList.push_back(doFixedPos(px));
    }
    fscanf(Instance, "%s", trash);
    for(int i = 0 ; i < ti.RequestNo ; i++)
    {
        double py;
        fscanf(Instance, "%lf", &py);

        ti.PickStopYList.push_back(doFixedPos(py));
    }
    fscanf(Instance, "%s", trash);
    for(int i = 0 ; i < ti.RequestNo ; i++)
    {
        double dx;
        fscanf(Instance, "%lf", &dx);

        ti.DropStopXList.push_back(doFixedPos(dx));
    }
    fscanf(Instance, "%s", trash);
    for(int i = 0 ; i < ti.RequestNo ; i++)
    {
        double dy;
        fscanf(Instance, "%lf", &dy);

        ti.DropStopYList.push_back(doFixedPos(dy));
    }
}

void ReaderHelper::GetNonCPFromInsFile(TripsInfo &ti)
{
    ti.TripPickDropSetNonCP.clear();
    for(Index i=0 ; i< ti.PickStopXList.size() ; i++)	//a loop from 0 to the number of requests
    {
        Coord cp;
        Coord cd;
        cp.first  = ti.PickStopXList[i];			//making coordination of the pick-up and drop-off nodes of customers
        cp.second = ti.PickStopYList[i];
        cd.first  = ti.DropStopXList[i];
        cd.second = ti.DropStopYList[i];


        if(ti.CPCoords.find(cp) == ti.CPCoords.end())			//if the made coordinations are not equal to the checkpoints, thus they are non-checkpoints
            ti.TripPickDropSetNonCP.insert(cp);
        if(ti.CPCoords.find(cd) == ti.CPCoords.end())
            ti.TripPickDropSetNonCP.insert(cd);
    }
}
