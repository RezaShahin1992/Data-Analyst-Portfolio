#ifndef MASTTIMER_H
#define MASTTIMER_H

#ifdef QT_WIDGETS_LIB
#include <QObject>
#include <QElapsedTimer>
#include <QTimer>
#include <QDebug>
#include "Solver.h"
#include "WriteHelper.h"
#include "ReaderHelper.h"

class MastTimer : public QObject
{
    Q_OBJECT
public:
    explicit MastTimer(QObject *parent = nullptr);
    void AssignPointer(Solver* s, ReaderHelper* _rh,WriteHelper* _wh);
    void CreateTimer();
    void StartTimer(int interval = 0);
    void FinishSimulation();
    bool IsSimulationFinished() const;

private:
    QTimer* timer = nullptr;
    Solver* solver;
    ReaderHelper* rh;
    WriteHelper* wh;
    QElapsedTimer elapsed;

    bool isFinished = false;
    set<int> visitedTrips;

signals:

public slots:
    void timeout();

};

#endif // MASTTIMER_H
#endif
