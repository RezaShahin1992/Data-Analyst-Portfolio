#ifndef GENERATORWINDOW_H
#define GENERATORWINDOW_H

#ifdef QT_WIDGETS_LIB

#include <QMainWindow>
#include "ReaderHelper.h"
#include "Solver.h"
#include "WriteHelper.h"
#include "MastTimer.h"
#include <QFile>
#include <fstream>
#include <QMessageBox>
#ifdef QT_WIDGETS_LIB
#ifdef QT_CHARTS_LIB
#include <QScatterSeries>
#include <QLineSeries>
#include <QValueAxis>
#include <QChartView>
#include <QPoint>
#endif
#include <cmath>
#include <math.h>
#include <QtMath>

#define M_Velocity 30.0 // Unit: KM/H => Kilometer per hour
#ifdef QT_CHARTS_LIB
using namespace QtCharts;
#endif
namespace Ui {
class GeneratorWindow;
}

class GeneratorWindow : public QMainWindow
{
    Q_OBJECT

public:
    typedef set<pair<double,double>> CoordSet;
    typedef pair<double,double> Coord;

    explicit GeneratorWindow(QWidget *parent = nullptr);
    ~GeneratorWindow();

    void PrepareSolver();

    void ComputeNonCPNumber();

    bool CheckMaxReq();

    void InitNodeLoad(const ReaderHelper::TripsInfo& ti);

    int GetIndexFromSet(const CoordSet& input,const Coord& coord);
#ifdef QT_CHARTS_LIB
    void InitGraph();
    void PrepareGraph(const ReaderHelper::TripsInfo &ti);
    void DrawGraph();
#endif

    void GenerateRequests(ReaderHelper::TripsInfo& ti,const ReaderHelper::TripsInfo& ti_pre);
    void GeneratePDRequests(ReaderHelper::TripsInfo& ti,const ReaderHelper::TripsInfo& ti_pre);
    void GenerateNPDRequests(ReaderHelper::TripsInfo& ti,const ReaderHelper::TripsInfo& ti_pre);
    void GeneratePNDRequests(ReaderHelper::TripsInfo& ti,const ReaderHelper::TripsInfo& ti_pre);
    void GenerateNPNDRequests(ReaderHelper::TripsInfo& ti,const ReaderHelper::TripsInfo& ti_pre);

    void GenerateStops(ReaderHelper::TripsInfo& ti, ReaderHelper::TripsInfo &ti_pre);
    void ComputeTeta(ReaderHelper::TripsInfo& ti, ReaderHelper::TripsInfo &ti_pre,const std::vector<Coord>& NodesList);

    void ComputeTripTime(ReaderHelper::TripsInfo& ti, const ReaderHelper::TripsInfo &ti_pre);

    void GenerateVehicles(ReaderHelper::TripsInfo& ti, ReaderHelper::TripsInfo &ti_next, const std::vector<ReaderHelper::TripsInfo>& tri_pre_list);

    void ComputeTau(ReaderHelper::TripsInfo& ti,const ReaderHelper::TripsInfo& ti_pre,const double& start_x,const double& start_y);

    void GenerateCheckPoints(ReaderHelper::TripsInfo& ti,const ReaderHelper::TripsInfo& ti_pre);

    void GenerateNonCheckPoints(ReaderHelper::TripsInfo& ti,const ReaderHelper::TripsInfo& ti_pre);

    std::vector<int> ShuffleGeneration(int from,int to,int no);
    std::vector<int> ShuffleGenerationUnSorted(int from,int to,int no);
    std::vector<double> ShuffleGenerationReal(double from,double to,int no);
    std::vector<double> ShuffleGenerationUnSortedReal(double from,double to,int no);

    void PrintToFile();

    void PrintTripToFile(ReaderHelper *rh, const int tno_local);

private slots:
    void on_tb_create_clicked();

    void on_tb_generate_clicked();

    void on_tb_generate_auto_clicked();

    void on_tb_finish_generate_clicked();

    void on_pb_Start_clicked();

private:
    Ui::GeneratorWindow *ui;

    ReaderHelper* rh;

    int m_NonCPNo = 0;
    int noPD_Req = 0;
    int noNPD_Req = 0;
    int noPND_Req = 0;
    int noNPND_Req = 0;

    int NonCPNo = 0;

    int TotalStartTime = 0;

    std::map<int,int> NodeLoad;
    std::map<int,Coord> NodeCoordMap;
    std::map<int,int> VehicleFreeTime;


    CoordSet GeneratorGraphNodes;
    CoordSet ChartGraphNodes;

    QVector<std::pair<Coord,Coord>> PD_RequestCoordinations;
    QVector<std::pair<Coord,Coord>> PND_RequestCoordinations;
    QVector<std::pair<Coord,Coord>> NPD_RequestCoordinations;
    QVector<std::pair<Coord,Coord>> NPND_RequestCoordinations;

    QVector<ReaderHelper*> m_AllLocal_RH;
#ifdef QT_CHARTS_LIB
    QChartView *chart = nullptr;
    QMainWindow* chartWindow = nullptr;
    QScatterSeries *series0 = nullptr;
    QLineSeries *seriesLine = nullptr;
    double maxRange = std::numeric_limits<double>::min();
#endif

};

#endif // GENERATORWINDOW_H
#endif
#endif
