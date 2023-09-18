#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include "Generator.h"
#include "Initialization.h"
#include "Writer.h"
#include <vector>
#include "Macro.h"

ILOSTLBEGIN

using namespace std;

int main() {
//    Generator *gene = new Generator();
    IloEnv env;
    IloModel model(env);

    Initialization Ob1;
    Generator Gen;
    Writer wr;

    // Generate the demand  - maximum production of traditional supplier x1, and non-traditional supplier x2
    vector<vector<int>> D = Gen.DemandGenerator(500, 501, numTimePeriods, numScenarios, &Ob1, &wr);
//    vector<vector<int>> X1max = Gen.X1ProductionGenerator(50, 100, numTimePeriods, numScenarios, &Ob1, &wr);
//    vector<vector<int>> X2max = Gen.X2ProductionGenerator(50, 100, numTimePeriods, numScenarios, &Ob1, &wr);

    // Pass the generated demand into DefineConstraints

    Ob1.Init();
 //   Ob1.PrepareCplexParams(1);
    Ob1.Solve();
    wr.WritingData(&Ob1,"OutPut");

    Ob1.FreeProblem();


    return 0;
}
