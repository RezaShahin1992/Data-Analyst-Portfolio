#include "Generator.h"
#include <iostream>
#include <random>
#include <vector>
using namespace std;

random_device rd;
mt19937 gen(rd());
Generator::Generator() {


}

Generator::~Generator() {

}
vector<vector<int>> Generator::DemandGenerator(int from, int to , int time, int scenario, Initialization *Ob1, Writer *wr)
{
	uniform_int_distribution<> dis(from, to);
    //time = 5;
    //scenario = 3;


    vector<vector<int>> D(scenario, vector<int>(time));


    for (int s = 0; s < scenario; s++)
    {
          for (int t = 0; t < time; t++)
          {
              D[s][t] = dis(gen);
              cout << "Scenario " << s << ", Period " << t << ", Demand: " << D[s][t] << std::endl;
          }
      }

    // Set D in the Initialization object
    Ob1->D = D;
    wr->D=D;
    return D;

}
vector<vector<int>> Generator::X1ProductionGenerator(int from, int to , int time, int scenario, Initialization *Ob1, Writer *wr)
{
	uniform_int_distribution<> dis(from, to);
    //time = 5;
    //scenario = 3;


    vector<vector<int>> X1max(scenario, vector<int>(time));


    for (int s = 0; s < scenario; s++)
    {
          for (int t = 0; t < time; t++)
          {
        	  X1max[s][t] = dis(gen);
              cout << "Scenario " << s << ", Period " << t << ", X1max: " << X1max[s][t] << std::endl;
          }
      }

    // Set X1MAX in the Initialization object
  //  Ob1->X1max = X1max;
    wr->X1max=X1max;
    return X1max;

}
vector<vector<int>> Generator::X2ProductionGenerator(int from, int to , int time, int scenario, Initialization *Ob1, Writer *wr)
{
	uniform_int_distribution<> dis(from, to);
    //time = 5;
    //scenario = 3;


    vector<vector<int>> X2max(scenario, vector<int>(time));


    for (int s = 0; s < scenario; s++)
    {
          for (int t = 0; t < time; t++)
          {
        	  X2max[s][t] = dis(gen);
              cout << "Scenario " << s << ", Period " << t << ", X2max: " << X2max[s][t] << std::endl;
          }
      }

    // Set X2MAX in the Initialization object
 //   Ob1->X2max = X2max;
    wr->X2max=X2max;
    return X2max;

}
