#ifdef QT_WIDGETS_LIB
#include "GeneratorWindow.h"
#include "ui_GeneratorWindow.h"

GeneratorWindow::GeneratorWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::GeneratorWindow)
{
    ui->setupUi(this);
    rh = new ReaderHelper();

    ui->pb_Start->setDisabled(true);

    setWindowTitle("Generator");

    on_tb_create_clicked();

    rh->AddInput = "INPUT/";
    rh->AddOutput = "OUTPUT/";
#ifdef QT_CHARTS_LIB
    chart = new QChartView();
    chartWindow = new QMainWindow(this);
    chartWindow->setCentralWidget(chart);
    chartWindow->resize(400,300);
    InitGraph();
#endif


}

GeneratorWindow::~GeneratorWindow()
{
    delete ui;
    delete rh;
}



void GeneratorWindow::PrintToFile()
{
    fstream Inf;
    Inf.open("Inf_Generation.txt",ios_base::out | ios_base::trunc);
    Inf << "TripsTotal:" << "\t" << rh->TripsNo << "\n";
    for(int i=0 ; i < ui->sb_tripsNumber->value() ; i++)
    {
        int TripID = i;
        ReaderHelper::TripsInfo& ti = rh->TripsInfoList[TripID];
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
        Inf << "NonCheckpoints:" << "\n";
        for(CoordSet::iterator it = ti.NonCPCoord.begin() ; it != ti.NonCPCoord.end() ; it++)
        {
            Inf << it->first << "\t" << it->second << "\n";
        }
        Inf << "\n";
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
    Ins.open("Ins_Generation.txt",ios_base::out | ios_base::trunc);
    for(int i=0 ; i < ui->sb_tripsNumber->value() ; i++)
    {
        int TripID = i;
        ReaderHelper::TripsInfo& ti = rh->TripsInfoList[TripID];
        Ins << "Requests.PickStopX:" << "\t" ;
        for(int i=0 ; i < ti.PickStopXList.size() ; i++)
        {
            Ins << ti.PickStopXList[i] << "\t";
        }
        Ins << "\n";
        Ins << "Requests.PickStopY:" << "\t" ;
        for(int i=0 ; i < ti.PickStopYList.size() ; i++)
        {
            Ins << ti.PickStopYList[i] << "\t";
        }
        Ins << "\n";
        Ins << "Requests.DropStopX:" << "\t" ;
        for(int i=0 ; i < ti.DropStopXList.size() ; i++)
        {
            Ins << ti.DropStopXList[i] << "\t";
        }
        Ins << "\n";
        Ins << "Requests.DropStopY:" << "\t" ;
        for(int i=0 ; i < ti.DropStopYList.size() ; i++)
        {
            Ins << ti.DropStopYList[i] << "\t";
        }
        Ins << "\n";
        Ins << "Requests.tau:" << "\t" ;
        for(int i=0 ; i < ti.RequestsTau.size() ; i++)
        {
            Ins << ti.RequestsTau[i] << "\t";
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

void GeneratorWindow::PrintTripToFile(ReaderHelper *rh, const int tno_local)
{
    fstream Inf;
    Inf.open("Inf_Generation_"+QString::number(tno_local).toStdString()+".txt",ios_base::out | ios_base::trunc);
    Inf << "TripsTotal:" << "\t" << rh->TripsNo << "\n";
    for(int i=0 ; i < rh->TripsNo ; i++)
    {
        int TripID = i;
        ReaderHelper::TripsInfo& ti = rh->TripsInfoList[TripID];
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
        Inf << "NonCheckpoints:" << "\n";
        for(CoordSet::iterator it = ti.NonCPCoord.begin() ; it != ti.NonCPCoord.end() ; it++)
        {
            Inf << it->first << "\t" << it->second << "\n";
        }
        Inf << "\n";
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
    Ins.open("Ins_Generation_"+QString::number(tno_local).toStdString()+".txt",ios_base::out | ios_base::trunc);
    for(int i=0 ; i < rh->TripsNo ; i++)
    {
        int TripID = i;
        ReaderHelper::TripsInfo& ti = rh->TripsInfoList[TripID];
        Ins << "Requests.PickStopX:" << "\t" ;
        for(int i=0 ; i < ti.PickStopXList.size() ; i++)
        {
            Ins << ti.PickStopXList[i] << "\t";
        }
        Ins << "\n";
        Ins << "Requests.PickStopY:" << "\t" ;
        for(int i=0 ; i < ti.PickStopYList.size() ; i++)
        {
            Ins << ti.PickStopYList[i] << "\t";
        }
        Ins << "\n";
        Ins << "Requests.DropStopX:" << "\t" ;
        for(int i=0 ; i < ti.DropStopXList.size() ; i++)
        {
            Ins << ti.DropStopXList[i] << "\t";
        }
        Ins << "\n";
        Ins << "Requests.DropStopY:" << "\t" ;
        for(int i=0 ; i < ti.DropStopYList.size() ; i++)
        {
            Ins << ti.DropStopYList[i] << "\t";
        }
        Ins << "\n";
        Ins << "Requests.tau:" << "\t" ;
        for(int i=0 ; i < ti.RequestsTau.size() ; i++)
        {
            Ins << ti.RequestsTau[i] << "\t";
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

void GeneratorWindow::PrepareSolver()
{
    double sumObjectValue = 0.0;
    for(int i=0 ; i < ui->sb_tripsNumber->value() ; i+=2)
    {
        int TripID = i;
        ReaderHelper::TripsInfo ti = rh->TripsInfoList[TripID];
        ReaderHelper::TripsInfo ti_next = rh->TripsInfoList[TripID+1];

        Solver *solver = new Solver();
        WriteHelper *wh = new WriteHelper();
        ReaderHelper *currentRH = new ReaderHelper();
        currentRH->TripsInfoList.push_back(ti);
        currentRH->TripsInfoList.push_back(ti_next);
        currentRH->TripsNo = currentRH->TripsInfoList.size(); // always 2
        currentRH->VehiclesNo = 1; // always one vehicle for trip

        for(int j=0 ; j < currentRH->TripsNo ; j++)
        {
            ReaderHelper::TripsInfo& ti_local = currentRH->TripsInfoList[j];
            ti_local.VehicleIndex = 0;
            ///PND
            ti_local.PND.resize(ti_local.RequestNo);
            for (int k = 0; k < ti_local.RequestNo; k++) {
                ti_local.PND[k].resize(currentRH->TripsNo);
                for (int r = 0; r < currentRH->TripsNo; r++) {
                    ti_local.PND[k][r].resize(currentRH->VehiclesNo);
                }
            }


            ///NPD
            ti_local.NPD.resize(ti_local.RequestNo);
            for (int k = 0; k < ti_local.RequestNo; k++) {
                ti_local.NPD[k].resize(currentRH->TripsNo);
                for (int r = 0; r < currentRH->TripsNo; r++) {
                    ti_local.NPD[k][r].resize(currentRH->VehiclesNo);
                }
            }

            ///LocalPND
            ti_local.LocalPND.resize(ti_local.RequestNo);
            for (int k = 0; k < ti_local.RequestNo; k++) {
                ti_local.LocalPND[k].resize(currentRH->TripsNo);
                for (int r = 0; r < currentRH->TripsNo; r++) {
                    ti_local.LocalPND[k][r].resize(currentRH->VehiclesNo);
                }
            }


            ///LocalNPD
            ti_local.LocalNPD.resize(ti_local.RequestNo);
            for (int k = 0; k < ti_local.RequestNo; k++) {
                ti_local.LocalNPD[k].resize(currentRH->TripsNo);
                for (int r = 0; r < currentRH->TripsNo; r++) {
                    ti_local.LocalNPD[k][r].resize(currentRH->VehiclesNo);
                }
            }
        }

        PrintTripToFile(currentRH,i);

        solver->PreProcess(currentRH);
        solver->Init(currentRH);

        solver->doSolve(i,currentRH);
        sumObjectValue += solver->objval;
        wh->WritingData(currentRH,solver,i);
        wh->WriteTripsTime(currentRH,i);
        solver->DeletingCplexVariables(currentRH);
        solver->Free_Problem();
        solver->ClearVariables();
        currentRH->ClearData();
//        m_AllLocal_RH.push_back(currentRH);

        delete currentRH;
        delete solver;
        delete wh;
    }

    qDebug() << "\n";
    qDebug() << "--------------------------";
    qDebug() << "Sum -> All Trips Object Value = "  << sumObjectValue;
    qDebug() << "--------------------------";
}

void GeneratorWindow::on_pb_Start_clicked()
{
    //    MastTimer *simTimer = new MastTimer();
    //    simTimer->CreateTimer();
        PrepareSolver();
}

void GeneratorWindow::ComputeNonCPNumber()
{
    noPD_Req = ui->sb_pd->value() * ui->sb_reqNo->value() / 100;
    noNPD_Req = ui->sb_npd->value() * ui->sb_reqNo->value() / 100;
    noPND_Req = ui->sb_pnd->value() * ui->sb_reqNo->value() / 100;
    noNPND_Req = ui->sb_npnd->value() * ui->sb_reqNo->value() / 100;
    int allReqNo = noPD_Req+noNPD_Req+noPND_Req+noNPND_Req;
    if(allReqNo < ui->sb_reqNo->value())
    {
        int diff = ui->sb_reqNo->value() - allReqNo;
        noNPND_Req+=diff;
    }
    NonCPNo = noNPD_Req+noPND_Req+(noNPND_Req*2);
}

bool GeneratorWindow::CheckMaxReq()
{
    double VehicleCap = ui->sb_vehicleCapacity->value();
    double NodeNo = ui->sb_CheckNo->value() + NonCPNo;
    double MaxReq = (NodeNo * VehicleCap) - VehicleCap;
    double ReqNo =  ui->sb_reqNo->value();

    if(MaxReq < ReqNo)
    {
        QMessageBox::critical(this,"RequestError","Your request number ["+QString::number(ReqNo)+"] is bigger than maximum request ["+QString::number(MaxReq)+"],Please increase vehicle capacity OR decrease request number");
        return false;
    }

    qDebug() << "MaxReq = " << MaxReq;
    qDebug() << "ReqNo = " << ReqNo;

    return true;
}

void GeneratorWindow::InitNodeLoad(const ReaderHelper::TripsInfo &ti)
{
    int TripNodeSize = ti.CPCoords.size() + ti.NonCPCoord.size();
    NodeLoad.clear();
    NodeCoordMap.clear();
    for(int i=0 ; i < TripNodeSize ; i++)
    {
        NodeLoad[i] = 0;
    }
    const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
    const CoordSet& NonCPCoord = ti.NonCPCoord; ///< Non Checkpoints Coordination
    GeneratorGraphNodes.clear();
    for(CoordSet::const_iterator cit = CPCoords.begin() ; cit != CPCoords.end() ; cit++)
    {
        GeneratorGraphNodes.insert(*cit);
        ChartGraphNodes.insert(*cit);
    }
    for(CoordSet::const_iterator cit = NonCPCoord.begin() ; cit != NonCPCoord.end() ; cit++)
    {
        GeneratorGraphNodes.insert(*cit);
        ChartGraphNodes.insert(*cit);
    }

    PD_RequestCoordinations.clear();
    PND_RequestCoordinations.clear();
    NPD_RequestCoordinations.clear();
    NPND_RequestCoordinations.clear();
}

int GeneratorWindow::GetIndexFromSet(const GeneratorWindow::CoordSet &input, const GeneratorWindow::Coord &coord)
{
    int cnt = 0;
    for(CoordSet::const_iterator cit = input.begin() ; cit != input.end() ; cit++)
    {
        const Coord& _c = *cit;
        if(_c == coord)
        {
            return cnt;
        }
        cnt++;
    }
    return -1;
}
#ifdef QT_CHARTS_LIB
void GeneratorWindow::InitGraph()
{

    series0 = new QScatterSeries();
    series0->setName("Graph");
    series0->setMarkerShape(QScatterSeries::MarkerShapeCircle);
    series0->setMarkerSize(10.0);
    seriesLine = new QLineSeries();
    seriesLine->setName("Path");
    seriesLine->setColor(Qt::green);
}

void GeneratorWindow::PrepareGraph(const ReaderHelper::TripsInfo &ti)
{
    const CoordSet& CP = ti.CPCoords;
    const CoordSet& NonCP = ti.NonCPCoord;

    for(CoordSet::const_iterator cit = GeneratorGraphNodes.begin() ; cit != GeneratorGraphNodes.end() ; cit++)
    {
        const Coord& pos = *cit;
        QScatterSeries *seriesNode = new QScatterSeries();
        seriesNode->setMarkerShape(QScatterSeries::MarkerShapeCircle);
        seriesNode->setMarkerSize(10.0);
        if(CP.find(pos) != CP.end())
        {
            seriesNode->append(pos.first,pos.second);
            seriesNode->setColor(Qt::red);
        }
        else
        {
            seriesNode->append(pos.first,pos.second);
            seriesNode->setColor(Qt::blue);
        }
        if(maxRange < pos.second)
        {
            maxRange = pos.second;
        }
        chart->chart()->addSeries(seriesNode);
//        series0->append(pos.first,pos.second);
        seriesLine->append(pos.first,pos.second);
    }
}

void GeneratorWindow::DrawGraph()
{
    chart->setRenderHint(QPainter::Antialiasing);
    chart->chart()->addSeries(seriesLine);
    chart->chart()->setTitle("Generator Graph");
    chart->chart()->createDefaultAxes();
    chart->chart()->setDropShadowEnabled();
    chart->chart()->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
    chart->chart()->axisX()->setLabelsAngle(-90);
    chart->chart()->axisX()->setGridLineVisible();
    chart->chart()->axisX()->setMinorGridLineVisible();
    maxRange++;
    int y = ui->servicearea_y->value();
    int x = ui->servicearea_x->value();
    int noTrip = ui->sb_tripsNumber->value();
    chart->chart()->axisY()->setRange(-y/2.0,y/2.0);
    chart->chart()->axisX()->setRange(0,x*noTrip);

    chartWindow->show();
}
#endif
void GeneratorWindow::on_tb_create_clicked()
{
    int tripsNo = ui->sb_tripsNumber->value();
    ui->cb_trips->clear();
    for(int i=tripsNo-1 ; i >= 0 ; i--)
    {
        ui->cb_trips->insertItem(0,"Trip_"+QString::number(i));
    }
    ui->cb_trips->setCurrentIndex(0);
    rh->TripsInfoList.resize(tripsNo);
    for(int i=0 ; i < rh->TripsInfoList.size() ; i++)
    {
        rh->TripsInfoList[i].id = i;
    }
}

void GeneratorWindow::on_tb_generate_auto_clicked()
{
    ComputeNonCPNumber();
    bool isOK = CheckMaxReq();
    if(!isOK)
        return;
    rh->TripsNo = ui->sb_tripsNumber->value();
    for(int i=0 ; i < ui->sb_tripsNumber->value() ; i++)
    {
        int TripID = i;
        ReaderHelper::TripsInfo& ti = rh->TripsInfoList[TripID];
        ti.VehicleCapacity = ui->sb_vehicleCapacity->value();
        ti.RequestNo = ui->sb_reqNo->value();
        if(TripID % 2 == 0)
        {
            ti.StartTime = TotalStartTime;
            ti.SuspendTime = -1;
        }
        else
        {
            ti.StartTime = TotalStartTime;
            TotalStartTime += ui->sb_starttimeInterval->value();
            ti.SuspendTime = ui->sb_suspend->value();
        }
        ReaderHelper::TripsInfo *ti_Pre;
        if(TripID - 1 >= 0)
            ti_Pre = &rh->TripsInfoList[TripID-1];
        else
            ti_Pre = &ti;


        GenerateCheckPoints(ti,*ti_Pre);

        GenerateNonCheckPoints(ti,*ti_Pre);

        InitNodeLoad(ti);

        GenerateStops(ti,*ti_Pre);

        ComputeTripTime(ti,*ti_Pre);

        GenerateRequests(ti,*ti_Pre);

        ti.NodesNo = ti.CPCoords.size() + ti.NonCPCoord.size();
        ti.C = ti.CPCoords.size();

        for (int i = 0; i < 3; i++) {
            ti.W.push_back(1);
        }
#ifdef QT_CHARTS_LIB
        PrepareGraph(ti);
#endif
    }

    for(int i=0 ; i < ui->sb_tripsNumber->value() ; i+=2)
    {
        int TripID = i;
        ReaderHelper::TripsInfo& ti = rh->TripsInfoList[TripID];
        ReaderHelper::TripsInfo *ti_next;
        std::vector<ReaderHelper::TripsInfo> tri_pre_list;
        if(TripID + 1 < (int)rh->TripsInfoList.size())
            ti_next = &rh->TripsInfoList[TripID+1];
        else
            ti_next = &ti;
        for(int j=0 ; j < i ; j++)
        {
            tri_pre_list.push_back(rh->TripsInfoList[j]);
        }

        GenerateVehicles(ti,*ti_next,tri_pre_list); ///@todo :: fix ti_next refrence for manual generation

        rh->VehiclesSet.insert(ti.VehicleIndex);
    }

    PrintToFile();

    rh->VehiclesNo = rh->VehiclesSet.size();

    ui->pb_Start->setEnabled(true);
    ui->tb_generate->setEnabled(false);
    ui->tb_finish_generate->setEnabled(false);
#ifdef QT_CHARTS_LIB
    if(ui->cb_drawGraph->isChecked())
    {
        DrawGraph();
    }
#endif
}

void GeneratorWindow::GenerateVehicles(ReaderHelper::TripsInfo &ti, ReaderHelper::TripsInfo &ti_next,
                                       const std::vector<ReaderHelper::TripsInfo> &tri_pre_list)
{
    if(tri_pre_list.empty()) // First Trip
    {
        ReaderHelper::VehicleContainer vc(true,0);
        rh->VehicleContainerList.push_back(vc);
        ti.VehicleIndex = 0;
        ti_next.VehicleIndex = 0;
        VehicleFreeTime[0] = ti_next.RealArrivalTime;
    }
    else
    {
        bool isSameVehicle = false;
        for(std::map<int,int>::const_iterator cit = VehicleFreeTime.begin() ; cit != VehicleFreeTime.end() ; cit++)
        {
            if(ti.RealStartTime > cit->second)
            {
                ti.VehicleIndex = cit->first;
                ti_next.VehicleIndex = cit->first;
                VehicleFreeTime[cit->first] = ti_next.RealArrivalTime;
                isSameVehicle = true;
                break;
            }
        }
        if(!isSameVehicle)
        {
            ti.VehicleIndex = rh->VehicleContainerList.size();
            ti_next.VehicleIndex = rh->VehicleContainerList.size();
            VehicleFreeTime[ti.VehicleIndex] = ti_next.RealArrivalTime;
            ReaderHelper::VehicleContainer vc(true,ti.VehicleIndex);
            rh->VehicleContainerList.push_back(vc);
        }

//        double minTime = std::numeric_limits<double>::max();
//        int detectedTripIndex = -1;
//        for(int i=0 ; i < tri_pre_list.size() ; i++)
//        {
//            const ReaderHelper::TripsInfo& pre_ti = tri_pre_list[i];
//            if(pre_ti.SuspendTime != -1) // return path
//            {
//                double ArrivalTime = pre_ti.RealArrivalTime;
//                if(ArrivalTime < minTime)
//                {
//                    minTime = ArrivalTime;
//                    detectedTripIndex = i;
//                }
//            }
//        }
//        if(detectedTripIndex != -1)
//        {
//            qDebug() << "ti.RealStartTime :: " << ti.RealStartTime << " | minTime :: " << minTime;
//            if(ti.RealStartTime > minTime)
//            {
//                ti.VehicleIndex = tri_pre_list[detectedTripIndex].VehicleIndex;
//                ti_next.VehicleIndex = tri_pre_list[detectedTripIndex].VehicleIndex;
//            }
//            else
//            {
//                ti.VehicleIndex = rh->VehicleContainerList.size();
//                ti_next.VehicleIndex = rh->VehicleContainerList.size();
//                ReaderHelper::VehicleContainer vc(true,ti.VehicleIndex);
//                rh->VehicleContainerList.push_back(vc);
//            }
//        }
    }
}


void GeneratorWindow::on_tb_generate_clicked()
{
    ComputeNonCPNumber();
    rh->TripsNo = ui->sb_tripsNumber->value();
    int TripID = ui->cb_trips->currentText().split("_")[1].toInt();
    ReaderHelper::TripsInfo& ti = rh->TripsInfoList[TripID];
    ti.VehicleCapacity = ui->sb_vehicleCapacity->value();
    ti.RequestNo = ui->sb_reqNo->value();
    bool isOK = CheckMaxReq();
    if(!isOK)
        return;

    if(TripID % 2 == 0)
    {
        ti.StartTime = TotalStartTime;
        ti.SuspendTime = -1;
    }
    else
    {
        ti.StartTime = TotalStartTime;
        TotalStartTime += ui->sb_starttimeInterval->value();
        ti.SuspendTime = ui->sb_suspend->value();
    }
    ReaderHelper::TripsInfo ti_Pre;
    if(TripID - 1 >= 0)
        ti_Pre = rh->TripsInfoList[TripID-1];
    else
        ti_Pre = ti;



    GenerateCheckPoints(ti,ti_Pre);

    GenerateNonCheckPoints(ti,ti_Pre);

    InitNodeLoad(ti);

    GenerateStops(ti,ti_Pre);

    ComputeTripTime(ti,ti_Pre);

    GenerateRequests(ti,ti_Pre);

    ti.NodesNo = ti.CPCoords.size() + ti.NonCPCoord.size();
    ti.C = ti.CPCoords.size();

    for (int i = 0; i < 3; i++) {
        ti.W.push_back(1);
    }



    rh->VehiclesSet.insert(ti.VehicleIndex);
    ui->tb_generate_auto->setEnabled(false);
    ui->pb_Start->setEnabled(true);
}

void GeneratorWindow::on_tb_finish_generate_clicked()
{
    for(int i=0 ; i < ui->sb_tripsNumber->value() ; i+=2)
    {
        int TripID = i;
        ReaderHelper::TripsInfo& ti = rh->TripsInfoList[TripID];
        ReaderHelper::TripsInfo ti_next;
        std::vector<ReaderHelper::TripsInfo> tri_pre_list;
        if(TripID + 1 < (int)rh->TripsInfoList.size())
            ti_next = rh->TripsInfoList[TripID+1];
        else
            ti_next = ti;
        for(int j=0 ; j < i ; j++)
        {
            tri_pre_list.push_back(rh->TripsInfoList[j]);
        }

        GenerateVehicles(ti,ti_next,tri_pre_list);

        rh->VehiclesSet.insert(ti.VehicleIndex);
    }

    PrintToFile();

    rh->VehiclesNo = rh->VehiclesSet.size();
}


void GeneratorWindow::GenerateRequests(ReaderHelper::TripsInfo &ti, const ReaderHelper::TripsInfo &ti_pre)
{
    GeneratePDRequests(ti,ti_pre);
    GeneratePNDRequests(ti,ti_pre);
    GenerateNPDRequests(ti,ti_pre);
    GenerateNPNDRequests(ti,ti_pre);

/*    bool isRequestReassembled = false;
    int assembling = 0;
    while(!isRequestReassembled)
    {
        bool LoadFailed = false;
        for(std::map<int,int>::iterator cit = NodeLoad.begin() ; cit != NodeLoad.end() ; cit++)
        {
            int CurrentLoad      = cit->second;
            int CurrentNodeIndex = cit->first;
            int VehicleCap = ui->sb_vehicleCapacity->value();
            if(CurrentLoad > VehicleCap)
            {
                LoadFailed = true;
                break;
            }
        }
        if(!LoadFailed)
        {
            isRequestReassembled = true;
            break;
        }
        for(std::map<int,int>::iterator cit = NodeLoad.begin() ; cit != NodeLoad.end() ; cit++)
        {
            int CurrentLoad      = cit->second;
            int CurrentNodeIndex = cit->first;
            int VehicleCap = ui->sb_vehicleCapacity->value();
            if(CurrentLoad > VehicleCap)
            {

                Coord CurrentNodeCoord = *std::next(GeneratorGraphNodes.begin(),CurrentNodeIndex);
                bool isRequestChanged = false;
                for(int i=0 ; i < PD_RequestCoordinations.size() ; i++)
                {
                    Coord PickCoordPD = PD_RequestCoordinations[i].first;
                    if(PickCoordPD == CurrentNodeCoord)
                    {
                        Coord DropCoordPD = PD_RequestCoordinations[i].second;
                        bool isChanged = false;
                        for(std::map<int,int>::iterator it = NodeLoad.begin() ; it != NodeLoad.end() ; it++)
                        {
                            if(it->second < VehicleCap) //Enough Load
                            {
                                int zeroloadIndex = it->first;
                                Coord ZeroNodeCoord = *std::next(GeneratorGraphNodes.begin(),zeroloadIndex);
                                if(ZeroNodeCoord.first < DropCoordPD.first && ti.CPCoords.find(ZeroNodeCoord) != ti.CPCoords.end())
                                {
                                    PD_RequestCoordinations[i].first = ZeroNodeCoord;
                                    isChanged = true;
                                    it->second++;
                                    break;
                                }
                            }
                        }
                        if(isChanged)
                        {
                            cit->second--;
                            isRequestChanged = true;
                            break;
                        }
                    }
                }
                if(isRequestChanged)
                {
                    break;
                }
                isRequestChanged = false;
                for(int i=0 ; i < PND_RequestCoordinations.size() ; i++)
                {
                    Coord PickCoordPD = PND_RequestCoordinations[i].first;
                    if(PickCoordPD == CurrentNodeCoord)
                    {
                        Coord DropCoordPD = PND_RequestCoordinations[i].second;
                        bool isChanged = false;
                        for(std::map<int,int>::iterator it = NodeLoad.begin() ; it != NodeLoad.end() ; it++)
                        {
                            if(it->second <VehicleCap) //Enough Load
                            {
                                int zeroloadIndex = it->first;
                                Coord ZeroNodeCoord = *std::next(GeneratorGraphNodes.begin(),zeroloadIndex);
                                if(ZeroNodeCoord.first < DropCoordPD.first && ti.CPCoords.find(ZeroNodeCoord) != ti.CPCoords.end())
                                {
                                    PND_RequestCoordinations[i].first = ZeroNodeCoord;
                                    isChanged = true;
                                    it->second++;
                                    break;
                                }
                            }
                        }
                        if(isChanged)
                        {
                            cit->second--;
                            isRequestChanged = true;
                            break;
                        }
                    }
                }
                if(isRequestChanged)
                {
                    break;
                }
                isRequestChanged = false;
                for(int i=0 ; i < NPD_RequestCoordinations.size() ; i++)
                {
                    Coord PickCoordPD = NPD_RequestCoordinations[i].first;
                    if(PickCoordPD == CurrentNodeCoord)
                    {
                        Coord DropCoordPD = NPD_RequestCoordinations[i].second;
                        bool isChanged = false;
                        for(std::map<int,int>::iterator it = NodeLoad.begin() ; it != NodeLoad.end() ; it++)
                        {
                            if(it->second <VehicleCap) //Enough Load
                            {
                                int zeroloadIndex = it->first;
                                Coord ZeroNodeCoord = *std::next(GeneratorGraphNodes.begin(),zeroloadIndex);
                                if(ZeroNodeCoord.first < DropCoordPD.first && ti.NonCPCoord.find(ZeroNodeCoord) != ti.NonCPCoord.end())
                                {
                                    NPD_RequestCoordinations[i].first = ZeroNodeCoord;
                                    isChanged = true;
                                    it->second++;
                                    break;
                                }
                            }
                        }
                        if(isChanged)
                        {
                            cit->second--;
                            isRequestChanged = true;
                            break;
                        }
                    }
                }
                if(isRequestChanged)
                {
                    break;
                }
                isRequestChanged = false;
                for(int i=0 ; i < NPND_RequestCoordinations.size() ; i++)
                {
                    Coord PickCoordPD = NPND_RequestCoordinations[i].first;
                    if(PickCoordPD == CurrentNodeCoord)
                    {
                        Coord DropCoordPD = NPND_RequestCoordinations[i].second;
                        bool isChanged = false;
                        for(std::map<int,int>::iterator it = NodeLoad.begin() ; it != NodeLoad.end() ; it++)
                        {
                            if(it->second <VehicleCap) //Enough Load
                            {
                                int zeroloadIndex = it->first;
                                Coord ZeroNodeCoord = *std::next(GeneratorGraphNodes.begin(),zeroloadIndex);
                                if(ZeroNodeCoord.first < DropCoordPD.first && ti.NonCPCoord.find(ZeroNodeCoord) != ti.NonCPCoord.end())
                                {
                                    NPND_RequestCoordinations[i].first = ZeroNodeCoord;
                                    isChanged = true;
                                    it->second++;
                                    break;
                                }
                            }
                        }
                        if(isChanged)
                        {
                            cit->second--;
                            isRequestChanged = true;
                            break;
                        }
                    }
                }
                if(isRequestChanged)
                {
                    break;
                }
            }
        }
        assembling++;
        if(assembling > NodeLoad.size())
            break;
    }


    if(!isRequestReassembled)
    {
        QMessageBox::critical(this,"RequestError","For trip ["+QString::number(ti.id)+"] : Vehicle load is bigger than vehicle capacity ["+QString::number(ui->sb_vehicleCapacity->value())+
                              "],solver error occured");
    }*/

    for(int i=0 ; i < PD_RequestCoordinations.size() ; i++)
    {
        Coord start = PD_RequestCoordinations[i].first;
        Coord finish = PD_RequestCoordinations[i].second;
        ti.PickStopXList.push_back(start.first);
        ti.PickStopYList.push_back(start.second);
        ti.DropStopXList.push_back(finish.first);
        ti.DropStopYList.push_back(finish.second);
    }
    for(int i=0 ; i < PND_RequestCoordinations.size() ; i++)
    {
        Coord start = PND_RequestCoordinations[i].first;
        Coord finish = PND_RequestCoordinations[i].second;
        ti.PickStopXList.push_back(start.first);
        ti.PickStopYList.push_back(start.second);
        ti.DropStopXList.push_back(finish.first);
        ti.DropStopYList.push_back(finish.second);
    }
    for(int i=0 ; i < NPD_RequestCoordinations.size() ; i++)
    {
        Coord start = NPD_RequestCoordinations[i].first;
        Coord finish = NPD_RequestCoordinations[i].second;
        ti.PickStopXList.push_back(start.first);
        ti.PickStopYList.push_back(start.second);
        ti.DropStopXList.push_back(finish.first);
        ti.DropStopYList.push_back(finish.second);
    }
    for(int i=0 ; i < NPND_RequestCoordinations.size() ; i++)
    {
        Coord start = NPND_RequestCoordinations[i].first;
        Coord finish = NPND_RequestCoordinations[i].second;
        ti.PickStopXList.push_back(start.first);
        ti.PickStopYList.push_back(start.second);
        ti.DropStopXList.push_back(finish.first);
        ti.DropStopYList.push_back(finish.second);
    }

    qDebug() << "Request Generation Finished ... " ;
}

void GeneratorWindow::GeneratePDRequests(ReaderHelper::TripsInfo &ti, const ReaderHelper::TripsInfo &ti_pre)
{
    ///<! Generate PD Requests

    const int range_from  = 0;
    const int range_to    = ti.CPCoords.size();
//    std::vector<int> values(range_to - range_from);
//    std::cerr << __FUNCTION__ << " :Range From : " << range_from << " To : " << range_to << std::endl;

//    std::generate(values.begin(), values.end(), [value = range_from]() mutable { return value++; });
//    std::shuffle(values.begin(), values.end(), std::mt19937{std::random_device{}()});
//    values.resize(noPD_Req*2);
//    std::sort(values.begin(),values.end());

    std::vector<int> values = ShuffleGeneration(range_from,range_to,noPD_Req*2);

    CoordSet GraphNodes;
    const CoordSet& CPCoords = ti.CPCoords; ///< Checkpoints Coordination
    int step = values.size() / 2;
    for(int i=0 ; i < step ; i++)
    {
        Coord start = *std::next(CPCoords.begin(),values[i]);
        Coord finish = *std::next(CPCoords.begin(),values[i+step]);
        if(start == finish)
        {
            qDebug() << "WTH :(" ;
        }
//        ti.PickStopXList.push_back(start.first);
//        ti.PickStopYList.push_back(start.second);
//        ti.DropStopXList.push_back(finish.first);
//        ti.DropStopYList.push_back(finish.second);
        int StartNode = GetIndexFromSet(GeneratorGraphNodes,start);
        int FinishNode = GetIndexFromSet(GeneratorGraphNodes,finish);
        NodeLoad[StartNode]++;
        NodeLoad[FinishNode]--;
        NodeCoordMap[StartNode] = start;
        NodeCoordMap[FinishNode] = finish;
        ComputeTau(ti,ti_pre,start.first,start.second);
        PD_RequestCoordinations.push_back(make_pair(start,finish));
    }
}

void GeneratorWindow::GenerateNPDRequests(ReaderHelper::TripsInfo &ti, const ReaderHelper::TripsInfo &ti_pre)
{

    std::vector<int> nonCPGen = ShuffleGeneration(0,ti.NonCPCoord.size(),noNPD_Req);
    std::vector<int> CPGen = ShuffleGeneration(0,ti.CPCoords.size(),noNPD_Req);

    bool isGenerationOk = false;
    int cntError = 0;
    while(!isGenerationOk)
    {
        std::set<int> visitedIndices;
        for(int i = nonCPGen.size() -1 ; i >= 0 ; i--)
        {
//            int noncpCandidate = nonCPGen[i];
            Coord noncpCoord = *std::next(ti.NonCPCoord.begin(),nonCPGen[i]);
            bool isFounded = false;
            for(int j=CPGen.size()-1 ; j >= 0 ; j--)
            {
                Coord cpCoord = *std::next(ti.CPCoords.begin(),CPGen[j]);
                if(cpCoord.first > noncpCoord.first && visitedIndices.find(CPGen[j]) == visitedIndices.end())
                {
                    visitedIndices.insert(CPGen[j]);
                    isFounded = true;
                    break;
                }
            }
            if(!isFounded)
            {
                CPGen = ShuffleGeneration(0,ti.CPCoords.size(),noNPD_Req);
                nonCPGen = ShuffleGeneration(0,ti.NonCPCoord.size(),noNPD_Req);
                break;
            }
        }
        if(visitedIndices.size() == CPGen.size())
            isGenerationOk = true;
        cntError++;
        if(cntError > 1000)
            break;
    }
    if(cntError > 1000)
    {
        QMessageBox::critical(this,"RequestError","NPD Cannot Generate. Please Increase Checkpoint number");
        qApp->exit(0);
    }
    for(size_t i = 0 ; i < nonCPGen.size() ; i++)
    {
        Coord start = *std::next(ti.NonCPCoord.begin(),nonCPGen[i]);
        Coord finish = *std::next(ti.CPCoords.begin(),CPGen[i]);
//        ti.PickStopXList.push_back(start.first);
//        ti.PickStopYList.push_back(start.second);
//        ti.DropStopXList.push_back(finish.first);
//        ti.DropStopYList.push_back(finish.second);
        int StartNode = GetIndexFromSet(GeneratorGraphNodes,start);
        int FinishNode = GetIndexFromSet(GeneratorGraphNodes,finish);
        NodeLoad[StartNode]++;
        NodeLoad[FinishNode]--;
        NodeCoordMap[StartNode] = start;
        NodeCoordMap[FinishNode] = finish;
        ComputeTau(ti,ti_pre,start.first,start.second);
        NPD_RequestCoordinations.push_back(make_pair(start,finish));
    }
}

void GeneratorWindow::GeneratePNDRequests(ReaderHelper::TripsInfo &ti, const ReaderHelper::TripsInfo &ti_pre)
{
    std::vector<int> CPGen = ShuffleGeneration(0,ti.CPCoords.size(),noPND_Req);
    std::vector<int> nonCPGen = ShuffleGeneration(0,ti.NonCPCoord.size(),noPND_Req);


    bool isGenerationOk = false;
    int cntError = 0;
    while(!isGenerationOk)
    {
        std::set<int> visitedIndices;
        for(int i = CPGen.size() -1 ; i>= 0 ; i--)
        {
//            int cpCandidate = CPGen[i];
            Coord cpCoord = *std::next(ti.CPCoords.begin(),CPGen[i]);
            bool isFounded = false;
            for(int j=nonCPGen.size()-1 ; j>= 0 ; j--)
            {
                Coord noncpCoord = *std::next(ti.NonCPCoord.begin(),nonCPGen[j]);
                if(noncpCoord.first > cpCoord.first && visitedIndices.find(nonCPGen[j]) == visitedIndices.end())
                {
                    visitedIndices.insert(nonCPGen[j]);
                    isFounded = true;
                    break;
                }
            }
            if(!isFounded)
            {
                CPGen = ShuffleGeneration(0,ti.CPCoords.size(),noPND_Req);
                nonCPGen = ShuffleGeneration(0,ti.NonCPCoord.size(),noPND_Req);
                break;
            }
        }
        if(visitedIndices.size() == nonCPGen.size())
            isGenerationOk = true;
        cntError++;
        if(cntError > 1000)
            break;
    }
    if(cntError > 1000)
    {
        QMessageBox::critical(this,"RequestError","PND Cannot Generate. Please Increase Checkpoint number");
        qApp->exit(0);
    }
    for(size_t i = 0 ; i < CPGen.size() ; i++)
    {
        Coord start = *std::next(ti.CPCoords.begin(),CPGen[i]);
        Coord finish = *std::next(ti.NonCPCoord.begin(),nonCPGen[i]);
//        ti.PickStopXList.push_back(start.first);
//        ti.PickStopYList.push_back(start.second);
//        ti.DropStopXList.push_back(finish.first);
//        ti.DropStopYList.push_back(finish.second);
        int StartNode = GetIndexFromSet(GeneratorGraphNodes,start);
        int FinishNode = GetIndexFromSet(GeneratorGraphNodes,finish);
        NodeLoad[StartNode]++;
        NodeLoad[FinishNode]--;
        NodeCoordMap[StartNode] = start;
        NodeCoordMap[FinishNode] = finish;
        ComputeTau(ti,ti_pre,start.first,start.second);
        PND_RequestCoordinations.push_back(make_pair(start,finish));
    }
}

void GeneratorWindow::GenerateNPNDRequests(ReaderHelper::TripsInfo &ti, const ReaderHelper::TripsInfo &ti_pre)
{
    std::vector<int> values = ShuffleGeneration(0,ti.NonCPCoord.size(),noNPND_Req*2);
    CoordSet GraphNodes;
    const CoordSet& nonCPCoords = ti.NonCPCoord; ///< Non-Checkpoints Coordination
    int step = values.size() / 2;
    for(int i=0 ; i < step ; i++)
    {
        Coord start = *std::next(nonCPCoords.begin(),values[i]);
        Coord finish = *std::next(nonCPCoords.begin(),values[i+step]);
//        ti.PickStopXList.push_back(start.first);
//        ti.PickStopYList.push_back(start.second);
//        ti.DropStopXList.push_back(finish.first);
//        ti.DropStopYList.push_back(finish.second);
        int StartNode = GetIndexFromSet(GeneratorGraphNodes,start);
        int FinishNode = GetIndexFromSet(GeneratorGraphNodes,finish);
        NodeLoad[StartNode]++;
        NodeLoad[FinishNode]--;
        NodeCoordMap[StartNode] = start;
        NodeCoordMap[FinishNode] = finish;
        ComputeTau(ti,ti_pre,start.first,start.second);
        NPND_RequestCoordinations.push_back(make_pair(start,finish));
    }
}

void GeneratorWindow::GenerateStops(ReaderHelper::TripsInfo &ti, ReaderHelper::TripsInfo &ti_pre)
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
            ti.StopsB.push_back(ui->sb_stopb->value());
    }
    ComputeTeta(ti,ti_pre,NodesList);
}

void GeneratorWindow::ComputeTeta(ReaderHelper::TripsInfo &ti, ReaderHelper::TripsInfo &ti_pre, const std::vector<GeneratorWindow::Coord> &NodesList)
{
    int lastNodeIndex = 0;
    for(size_t i=0 ; i < NodesList.size()  ; i++)
    {
        const Coord& node = NodesList[i];
        if(i== 0 && ti.id == 0) // First Node Of First Trip :: 0 value for teta
        {
            ti.StopsTeta.push_back(0);
        }
        else
        {
            if(i == 0 && ti.SuspendTime != -1) // first node of retrun trip
            {
                ti.StopsTeta.push_back(ti_pre.StopsTeta[ti_pre.StopsTeta.size()-1]);
            }
            else
            {
                if(i == 0)
                {
                    ti.StopsTeta.push_back(0);
                    continue;
                }
                if(ti.CPCoords.find(node) != ti.CPCoords.end())
                {
                    double Distance = 0;
                    int sumB = 0;
                    for(size_t j=lastNodeIndex ; j <= i-1 ; j++)
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
                    double computedTeta = (computedTime+ti.StopsTeta[lastNodeIndex]+sumB);
                    double intPart = 0;
                    std::modf(computedTeta,&intPart);
                    intPart++; // Round to upper bound
                    int CT = intPart;
                    ti.StopsTeta.push_back(CT);
                    lastNodeIndex = i;
                }
                else
                {
                    ti.StopsTeta.push_back(0);
                }
            }
        }

    }
    if(ti.SuspendTime != -1)
    {
        if(ti.StopsType[0] != 2)
            ti.StopsTeta[0] += ti.SuspendTime;
        ti_pre.StopsTeta[ti_pre.StopsTeta.size()-1] += ti.SuspendTime;
    }

    if(ui->rb_maxOptimal->isChecked())
    {
        double maxTeta = std::numeric_limits<double>::min();

        for(int j=0 ; j < ti.StopsTeta.size() ; j++)
        {
            if(ti.StopsTeta[j] > maxTeta)
            {
                maxTeta = ti.StopsTeta[j];
            }
        }
        maxTeta += ui->sb_interval->value();
        qDebug() << "maxTeta :: " << maxTeta;
        double initMaxTeta = maxTeta;
        for(int j=0 ; j < ti.StopsTeta.size() ; j++)
        {
            if(j == 0 && ti.SuspendTime != -1)
            {
                ti.StopsTeta[j] = ti_pre.StopsTeta[ti_pre.StopsTeta.size()-1];
            }
            else
            {
                if(j == 0)
                    continue;
                if(ti.StopsType[j] != 2)
                {
                    ti.StopsTeta[j] = maxTeta;
                    maxTeta += initMaxTeta;
                }
            }
        }
    }
}

void GeneratorWindow::ComputeTripTime(ReaderHelper::TripsInfo& ti, const ReaderHelper::TripsInfo &ti_pre)
{
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
    if(ti.SuspendTime == -1)
    {
        ti.RealStartTime = ti.StartTime;
        ti.RealArrivalTime = ti.RealStartTime + ti.StopsTeta[ti.StopsTeta.size()-1] - ui->sb_suspend->value();
    }
    else
    {
        ti.RealStartTime = ti_pre.StopsTeta[ti_pre.StopsTeta.size()-1] + ti_pre.StartTime;
        ti.RealArrivalTime = ti.RealStartTime + ti.StopsTeta[ti.StopsTeta.size()-1];
    }
}

void GeneratorWindow::ComputeTau(ReaderHelper::TripsInfo &ti, const ReaderHelper::TripsInfo& ti_pre, const double &start_x, const double &start_y)
{
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
    bool isTauAssigned = false;
    for(size_t i=0 ; i < NodesList.size() -1 ; i++)
    {
        const Coord& node = NodesList[i];
        double x1 = node.first;
        double y1 = node.second;
        if(x1 == start_x && y1 == start_y && !isTauAssigned && i == 0)
        {
            ti.RequestsTau.push_back(0);
            isTauAssigned = true;
        }
        const Coord& nextNode = NodesList[i+1];
        double x2 = nextNode.first;
        double y2 = nextNode.second;
        double d = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        Distance += d;
        if(x2 == start_x && y2 == start_y && !isTauAssigned)
        {
            int computedTime = Distance / M_Velocity;
            double sumTeta = 0;

            double range_from  = ti.RealStartTime;
            double range_to    = ti.RealStartTime+computedTime + sumTeta;
//            std::vector<double> values(range_to - range_from);
//            std::cerr << __FUNCTION__ << ":Range From : " << range_from << " To : " << range_to << std::endl;

//            std::generate(values.begin(), values.end(), [value = range_from]() mutable { return value++; });
//            std::shuffle(values.begin(), values.end(), std::mt19937{std::random_device{}()});

//            values.resize(1);
            std::vector<double> values = ShuffleGenerationUnSortedReal(range_from,range_to,1);
//            ti.RequestsTau.push_back(values[0]);
            ti.RequestsTau.push_back(0);
            isTauAssigned = true;
        }

    }
}

void GeneratorWindow::GenerateCheckPoints(ReaderHelper::TripsInfo &ti,const ReaderHelper::TripsInfo& ti_pre)
{
    double sa_x = ui->servicearea_x->value();
    double noCP = ui->sb_CheckNo->value();
    double step = sa_x / noCP;

    double x = 0;
    double init_X = 0;
    if(ti_pre.CPCoords.empty())
        x = 0;
    else
        x = (*--ti_pre.CPCoords.end()).first;
    init_X = x;

    int y = 0;

    QVector<Coord> localCoord;
    localCoord.resize(noCP);
    pair<double,double> firstPos = make_pair(x,y);
    pair<double,double> LastPos = make_pair(sa_x+init_X,y);
    localCoord[0] = firstPos;
    localCoord[noCP-1] = LastPos;

    noCP -=2;
    step = (double)sa_x / (double)noCP;
    x+=step;
    ti.CPCoords.insert(firstPos);
    for(int i=1 ; i < localCoord.size() -1 ; i++)
    {
        pair<double,double> pos = make_pair(x,y);
        ti.CPCoords.insert(pos);
        x+=step;
    }
    ti.CPCoords.insert(LastPos);

//    int cnt = 0;
//    while(ti.CPCoords.size() != noCP)
//    {
//        if(/*ti.id == 0 && */ti.CPCoords.size() == noCP -1)
//        {
//            pair<double,double> pos = make_pair(sa_x+init_X,y);
//            ti.CPCoords.insert(pos);
//        }
//        else
//        {
//            pair<double,double> pos = make_pair(x,y);
//            ti.CPCoords.insert(pos);
//            x+=step;
//        }
//    }

}

void GeneratorWindow::GenerateNonCheckPoints(ReaderHelper::TripsInfo &ti,const ReaderHelper::TripsInfo& ti_pre)
{
    int sa_x = (*--ti.CPCoords.end()).first;
    int sa_y = ui->servicearea_y->value();
    bool isFirstTrip = false;
    int rangeFrom = 0;
    if(ti_pre.NonCPCoord.empty())
    {
        isFirstTrip = true;
        rangeFrom = 0;
    }
    if(!isFirstTrip)
    {
        rangeFrom = (*ti.CPCoords.begin()).first;
        rangeFrom++;
    }

    std::vector<int> X_NonCP;
    bool is_XnonCP_OK = false;
    while(!is_XnonCP_OK)
    {
        X_NonCP = ShuffleGenerationUnSorted(rangeFrom,sa_x,NonCPNo);
        bool isOK = true;
        for(int i=0 ; i < X_NonCP.size() ; i++)
        {
            Coord searchCoord = make_pair(X_NonCP[i],0);
            if(ti.CPCoords.find(searchCoord) != ti.CPCoords.end())
            {
                isOK = false;
                break;
            }
        }
        if(isOK)
            is_XnonCP_OK = true;
    }


    std::sort(X_NonCP.begin(),X_NonCP.end());

    std::vector<int> Y_NonCP = ShuffleGenerationUnSorted(-(sa_y/2.0),sa_y/2.0,NonCPNo);

//    Y_NonCP.resize(NonCPNo);
//    std::sort(Y_NonCP.begin(),Y_NonCP.end());



    for(int i=0 ; i < Y_NonCP.size() ; i++)
    {
        pair<double,double> pos = make_pair(X_NonCP[i],Y_NonCP[i]);
        ti.NonCPCoord.insert(pos);
    }
}

std::vector<int> GeneratorWindow::ShuffleGeneration(int from, int to, int no)
{
    const int range_from  = from;
    const int range_to    = to;
    std::vector<int> values(range_to - range_from);

    std::generate(values.begin(), values.end(), [value = range_from]() mutable { return value++; });
    std::shuffle(values.begin(), values.end(), std::mt19937{std::random_device{}()});
    if(values.size() > no)
        values.resize(no);
    std::sort(values.begin(),values.end());
    return values;


//    std::vector<int> res;
//    if(from > to)
//    {
//        qDebug() << __FUNCTION__ << "ERROR :: From:" << from << "is bigger than to:" << to;
//        qDebug() << "Generation Is Empty";
//        return res;
//    }
//    std::random_device rd;  // Will be used to obtain a seed for the random number engine
//    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
//    std::uniform_int_distribution<> dis(from, to);
//    for (int n = 0; n < no; ++n) {
//        // Use dis to transform the random unsigned int generated by gen into a
//        // double in [from, to). Each call to dis(gen) generates a new random double
//        res.push_back(dis(gen));
//    }
//    std::sort(res.begin(),res.end());
////    std::cout << '\n';
//    return res;


}

std::vector<int> GeneratorWindow::ShuffleGenerationUnSorted(int from, int to, int no)
{

    const int range_from  = from;
    const int range_to    = to;
    std::vector<int> values(range_to - range_from);

    std::generate(values.begin(), values.end(), [value = range_from]() mutable { return value++; });
    std::shuffle(values.begin(), values.end(), std::mt19937{std::random_device{}()});
    if(values.size() > no)
        values.resize(no);
//    std::sort(values.begin(),values.end());
    return values;


//    std::vector<int> res;
//    if(from > to)
//    {
//        qDebug() << __FUNCTION__ << "ERROR :: From:" << from << "is bigger than to:" << to;
//        qDebug() << "Generation Is Empty";
//        return res;
//    }
//    std::random_device rd;  // Will be used to obtain a seed for the random number engine
//    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
//    std::uniform_int_distribution<> dis(from, to);
//    for (int n = 0; n < no; ++n) {
//        // Use dis to transform the random unsigned int generated by gen into a
//        // double in [from, to). Each call to dis(gen) generates a new random double
//        res.push_back(dis(gen));
//    }
//    return res;
}

std::vector<double> GeneratorWindow::ShuffleGenerationReal(double from, double to, int no)
{

    const double range_from  = from;
    const double range_to    = to;
    std::vector<double> values(range_to - range_from);

    std::generate(values.begin(), values.end(), [value = range_from]() mutable { return value++; });
    std::shuffle(values.begin(), values.end(), std::mt19937{std::random_device{}()});
    if(values.size() > no)
        values.resize(no);
    std::sort(values.begin(),values.end());
    return values;

//    std::vector<double> res;
//    if(from > to)
//    {
//        qDebug() << __FUNCTION__ << "ERROR :: From:" << from << "is bigger than to:" << to;
//        qDebug() << "Generation Is Empty";
//        return res;
//    }
//    std::random_device rd;  // Will be used to obtain a seed for the random number engine
//    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
//    std::uniform_real_distribution<> dis(from, to);
//    for (int n = 0; n < no; ++n) {
//        // Use dis to transform the random unsigned int generated by gen into a
//        // double in [from, to). Each call to dis(gen) generates a new random double
//        res.push_back(dis(gen));
//    }
//    std::sort(res.begin(),res.end());
//    return res;
}

std::vector<double> GeneratorWindow::ShuffleGenerationUnSortedReal(double from, double to, int no)
{
    const double range_from  = from;
    const double range_to    = to;
    std::vector<double> values(range_to - range_from);

    std::generate(values.begin(), values.end(), [value = range_from]() mutable { return value++; });
    std::shuffle(values.begin(), values.end(), std::mt19937{std::random_device{}()});
    if(values.size() > no)
        values.resize(no);
//    std::sort(values.begin(),values.end());
    return values;

//    std::vector<double> res;
//    if(from > to)
//    {
//        qDebug() << __FUNCTION__ << "ERROR :: From:" << from << "is bigger than to:" << to;
//        qDebug() << "Generation Is Empty";
//        return res;
//    }
//    std::random_device rd;  // Will be used to obtain a seed for the random number engine
//    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
//    std::uniform_real_distribution<> dis(from, to);
//    for (int n = 0; n < no; ++n) {
//        // Use dis to transform the random unsigned int generated by gen into a
//        // double in [from, to). Each call to dis(gen) generates a new random double
//        res.push_back(dis(gen));
//    }
//    return res;
}

#endif
