#ifndef WRITEHELPER_H
#define WRITEHELPER_H

#include "ReaderHelper.h"
#include "Solver.h"
#include "TripsInfo.h"
#include "SolverBase.h"

class WriteHelper
{
public:
    WriteHelper();
    ~WriteHelper(){}

    void WritingData(const ReaderHelper* rh ,const SolverBase* sb,const Solver* solver,const int& ins,const string& outputname);

private:

    FILE* Out;
};

#endif // WRITEHELPER_H
