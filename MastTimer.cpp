#ifdef QT_WIDGETS_LIB
#include "MastTimer.h"
#include <QCoreApplication>


MastTimer::MastTimer(QObject *parent) : QObject(parent)
{

}

void MastTimer::AssignPointer(Solver *s, ReaderHelper *_rh, WriteHelper *_wh)
{
    solver = s;
    rh = _rh;
    wh = _wh;
}

void MastTimer::CreateTimer()
{
    timer = new QTimer();
    connect(timer,&QTimer::timeout,this,&MastTimer::timeout);
}

void MastTimer::StartTimer(int interval)
{
    if(timer != nullptr)
    {
        timer->start(interval);
        elapsed.start();
    }
}

void MastTimer::FinishSimulation()
{
    isFinished = true;
    qApp->quit();
}

bool MastTimer::IsSimulationFinished() const
{
    return isFinished;
}

void MastTimer::timeout()
{
    int elapsedSec = elapsed.elapsed() / 1000;

    const int& TripsNo = rh->TripsNo;
    if((int)visitedTrips.size() == TripsNo)
    {
        timer->stop();
        FinishSimulation();
        return;
    }
    for(int tno = 0 ; tno < TripsNo ; tno++)
    {
        const ReaderHelper::TripsInfo &ti = rh->TripsInfoList[tno];

        if(visitedTrips.find(tno) == visitedTrips.end() && ti.RealStartTime == elapsedSec)
        {
            cout << "Trip["<<tno<<"] : Started In Time : " << ti.RealStartTime << endl;
            solver->DefineConstraints(rh,tno);
//            solver->DefineConstraints(rh,tno+1);
            visitedTrips.insert(tno);
//            visitedTrips.insert(tno+1);
            break;
        }
    }
}
#endif
