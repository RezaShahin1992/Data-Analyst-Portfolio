#include "WriteHelper.h"

WriteHelper::WriteHelper()
{

}

void WriteHelper::WritingData(const ReaderHelper *rh, const SolverBase *sb, const Solver *solver, const int &ins, const string &outputname)
{
    const char* OutName = new char;

    string StrInfName = rh->AddOutput+ outputname+"-"+to_string(ins) + ".txt";
    OutName = StrInfName.c_str();

    Out = fopen(OutName, "w");
    if (Out == NULL)
    {
        cout << "cannt create output!!" << endl;
        exit(12);
    }

    fprintf(Out, "Obj:	%f\n", solver->objval);

    ////////////////	 Print X-variables	 //////////////
    fprintf(Out, "%s", "X-variables:\n");
    fprintf(Out, "%s", "\n");

    int NodesNo = 0;
    int VehiclesNo = sb->VehiclesNo;
    int RequestNo = 0;
    int TripsNo = sb->TripsNo;
    vector <int> ReqType;
    for(int i=0 ; i < TripsNo ; i++)
    {
        const TripsInfo &ti = sb->TripsInfoList[i];
        RequestNo += ti.RequestNo;
        for(int j=0 ; j < ti.RequestsTau.size() ; j++)
            ReqType.push_back(ti.RequestsType[j]);
    }

    NodesNo = sb->nbPoints;

    try{
    for (int i = 0; i < NodesNo; i++)
    {
        for (int j = 0; j < NodesNo; j++)
        {
            for (int v = 0; v < VehiclesNo; v++)
            {
                if(solver->ArcsMatrix[i][j] == 0)
                    continue;
                IloNum x = solver->cplexPtr->getValue((*solver->Xval)[i][j][v]);
                if (x > 0.5)
                {
                    x = std::round(x);
                    fprintf(Out, "X[%d][%d][%d]= %f\n", i, j,v, x);
                }
            }
        }
    }
    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":X Parameters :: " << e << "\n";
    }

    try{
        for (int i = 0; i < NodesNo; i++)
        {
            fprintf(Out, "t[%d]=%f\n",i, solver->cplexPtr->getValue((*solver->Tval)[i]));
        }
    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":T Parameters :: " << e << "\n";
    }

    ////////////
    try{
        for (int i = 0; i < NodesNo; i++)
        {
            IloNumVar var = (*solver->Tbarval)[i];
            if(solver->m_Sources.find(i) == solver->m_Sources.end() && solver->m_Destinations.find(i) == solver->m_Destinations.end())
                fprintf(Out, "tbar[%d]=%f\n",i, solver->cplexPtr->getValue(var));
        }
    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":TBar Parameters :: " << e << "\n";
    }

    ////////////
    try{
        for (int k = 0; k < RequestNo; k++)
        {
            fprintf(Out, "p[%d]=%f\n",k, solver->cplexPtr->getValue((*solver->Pval)[k]));
        }
    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":P Parameters :: " << e << "\n";
    }
    ////////////
    try{
        for (int k = 0; k < RequestNo; k++)
        {
            fprintf(Out, "d[%d]=%f\n", k, solver->cplexPtr->getValue((*solver->Dval)[k]));
        }
    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":D Parameters :: " << e << "\n";
    }
    //////////
    try{
        for (int k = 0; k < RequestNo; k++)
        {
            if(ReqType[k] == PND)
            {
                for (int r = 0; r < TripsNo; r++)
                {
                    int v = sb->TripsInfoList[r].VehicleIndex;
                    if (solver->S_PND[k][r][v] != -1)
                    {
                        IloNum z = solver->cplexPtr->getValue((*solver->Zval)[k][r][v]);
                        if(z > 0.5)
                        {
                            z = std::round(z);
                            fprintf(Out, "Z[%d][%d][%d]= %f\n", k, r, v, z);
                        }

                    }
                }
            }
            if(ReqType[k] == NPD)
            {
                for (int r = 0; r < TripsNo; r++)
                {
                    int v = sb->TripsInfoList[r].VehicleIndex;

                    if (solver->S_NPD[k][r][v] != -1)
                    {
                        IloNum z = solver->cplexPtr->getValue((*solver->Zval)[k][r][v]);
                        if(z > 0.5)
                        {
                            z = std::round(z);
                            fprintf(Out, "Z[%d][%d][%d]= %f\n", k, r, v, z);
                        }
                    }
                }
            }
        }
    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":Z Parameters :: " << e << "\n";
    }
    ////////
    try{
        for(int i=0 ; i < NodesNo ; i++)
        {
            int load = solver->cplexPtr->getValue((*solver->qval)[i]);
            fprintf(Out, "q[%d]= %d\n", i, load);
        }
    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":q Parameters :: " << e << "\n";
    }
    try{
        for(int i=0 ; i < NodesNo ; i++)
        {
            for(int v=0 ; v < VehiclesNo ; v++)
            {
                int cap = solver->cplexPtr->getValue((*solver->Qval)[i][v]);
                if(cap >= 0 && cap < M - 1)
                    fprintf(Out, "Q[%d][%d]= %d\n", i,v, cap);
            }
        }
    }catch(const IloException &e)
    {
        std::cerr << __FUNCTION__ << ":Q Parameters :: " << e << "\n";
    }

    ////////////

    fclose(Out);
    cout << "objval :: " << solver->objval << endl;
}
