#include <iostream>
#include "Solver.h"
#include "WriteHelper.h"
#include "ReaderHelper.h"
#include "GeneratorTerminal.h"
#include "SolverBase.h"
#include "Heuristic.h"


using namespace std;

void SecondRun(int i);

int main(int argc, char *argv[])
{

    string strArguman = "";							//to initialize strArguman
    for(int i=1 ; i < argc ; i++)
    {
        strArguman += string(argv[i]);  			//argv[i] (command line argument) is converted to a string and
    }												// appended to the strArguman using +=

    if(argc > 1)
    {
        if(strArguman == "-g")						//to generate data
        {
            int cntinput = 0;
//            cout << "Enter Generator Input Number :: " << endl;		//ask us how many generator input do we have
//            cin >> cntinput;
            for(int i=0 ; i < GENERATOR_NO ; i++)						//a loop to read all of them
            {
                GeneratorTerminal gt;								//creating an object "gt" of class GeneratorTerminal
                if(gt.ReadInputFile(i))								//go to the function ReadInputFile in GeneratorTerminal and read the input of generator
                {
                    gt.StartGenerationOnTimeHorizons();				//go and generate data ...
                }
            }
        }
    }
    else
    {
        int cnt = 0;
        for(int i=0 ; i < MAIN_STEP ; i++)
        {
            SolverBase *base = new SolverBase();
            ReaderHelper *rh = new ReaderHelper();
            Solver *solver = new Solver();
            WriteHelper *wh = new WriteHelper();
            Heuristic *heu = new Heuristic(solver,base);


            bool isExist = rh->Read_Address(i);
            if(!isExist)
            {
                delete rh;
                delete solver;
                delete wh;

                continue;
            }
            cnt = i;

            double LocalSyncTime = 0.0;
            double t1 = clock();
            rh->ReadData(cnt,base);
            base->execute();
            solver->cntInputFile = cnt;
            solver->PreProcess(base);

            /**
              Heuristic Computation
              **/
            heu->doHeuristicSolve();

            heu->CleanUnusedArcsMatrix();

            heu->LoadTrips();

            solver->Init(base);

            heu->ApplyHeuToCplex();

            double t2 = clock();
            LocalSyncTime = (double)(t2-t1) / (double)CLOCKS_PER_SEC;

            t1 = clock();
            solver->doSolve(cnt,base);
            t2 = clock();
            double SolvingTime = (double)(t2-t1) / (double)CLOCKS_PER_SEC;

            fstream ConsoleOutput;
            ConsoleOutput.open("OUTPUT/Console"+to_string(cnt)+".txt",ios_base::out | ios_base::app);
            ConsoleOutput << "\n";
            ConsoleOutput << "------------RESULT--------------"<< std::endl;
            ConsoleOutput << "Sum -> All Trips Object Value = "  << solver->objval<< std::endl;
            ConsoleOutput << "Sum -> All Trips Initialization Time = "  << LocalSyncTime << std::endl;
            ConsoleOutput << "Sum -> All Trips Computation Time = "  << SolvingTime << std::endl;
            ConsoleOutput << "TotalTravelTime :: " << solver->TotalTravelTime << endl;
            ConsoleOutput << "TotalRequestTime :: " << solver->TotalRequestTime << endl;
            ConsoleOutput << "TotalRequestWaitingTime :: " << solver->TotalRequestWaitingTime << endl;
            ConsoleOutput << "TotalRejection :: " << solver->SumRej << endl;
            ConsoleOutput << "Time Horizons :: -> nbTrips = "  << base->TripsNo<< std::endl;
            ConsoleOutput << "Time -> nbRequest = "  << rh->LocalTotalReqNo<< std::endl;
            ConsoleOutput << "--------------------------"<< std::endl;
            ConsoleOutput.close();

            wh->WritingData(rh,base,solver,cnt,"OutputHeu");
//            solver->Free_Problem();

            delete solver;
            delete wh;
            delete base;
//            delete heu;

            SecondRun(i);
        }

    }

//0.03
//#ifdef QT_WIDGETS_LIB
//    return a.exec();
//#else
    return 0;
//#endif
}

void SecondRun(int i)
{
    int cnt = 0;
//    for(int i=0 ; i < MAIN_STEP ; i++)
    {
        SolverBase *base = new SolverBase();
        ReaderHelper *rh = new ReaderHelper();
        Solver *solver = new Solver();
        WriteHelper *wh = new WriteHelper();

        bool isExist = rh->Read_Address(i);
        if(!isExist)
        {
            delete rh;
            delete solver;
            delete wh;

            return;
        }
        cnt = i;

        double LocalSyncTime = 0.0;
        double t1 = clock();
        rh->ReadData(cnt,base);
        rh->ReadOutputFile(cnt);
        rh->SendOutputDataToBase(base);
        base->execute();
        solver->cntInputFile = cnt;
        solver->IsSecondRun = true;
        solver->PreProcess(base);
        solver->Init(base);
        double t2 = clock();
        LocalSyncTime = (double)(t2-t1) / (double)CLOCKS_PER_SEC;

        solver->AddDataToWarmStart(cnt,base);

        t1 = clock();
        solver->doSolve(cnt,base);
        t2 = clock();
        double SolvingTime = (double)(t2-t1) / (double)CLOCKS_PER_SEC;

        fstream ConsoleOutput;
        ConsoleOutput.open("OUTPUT/Console"+to_string(cnt)+".txt",ios_base::out | ios_base::app);
        ConsoleOutput << "\n";
        ConsoleOutput << "------------RESULT--------------"<< std::endl;
        ConsoleOutput << "Sum -> All Trips Object Value = "  << solver->objval<< std::endl;
        ConsoleOutput << "Sum -> All Trips Initialization Time = "  << LocalSyncTime << std::endl;
        ConsoleOutput << "Sum -> All Trips Computation Time = "  << SolvingTime << std::endl;
        ConsoleOutput << "TotalTravelTime :: " << solver->TotalTravelTime << endl;
        ConsoleOutput << "TotalRequestTime :: " << solver->TotalRequestTime << endl;
        ConsoleOutput << "TotalRequestWaitingTime :: " << solver->TotalRequestWaitingTime << endl;
        ConsoleOutput << "TotalRejection :: " << solver->SumRej << endl;
        ConsoleOutput << "Time Horizons :: -> nbTrips = "  << base->TripsNo<< std::endl;
        ConsoleOutput << "Time -> nbRequest = "  << rh->LocalTotalReqNo<< std::endl;
        ConsoleOutput << "--------------------------"<< std::endl;
        ConsoleOutput.close();

        wh->WritingData(rh,base,solver,cnt,"Output");
        solver->Free_Problem();

        delete solver;
        delete wh;
        delete base;

    }
}
