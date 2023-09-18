#include "GeneratorTerminal.h"

GeneratorTerminal::GeneratorTerminal():
    GeneratorBase()
{

}

bool GeneratorTerminal::ReadInputFile(int gn)    //it reads the input of the generator
{
    Input_CntInput = gn;						//assigns the number "gn" to Input_CntInput
    fstream inputfile;							//declare inputfile of type fstream, to handle input file operations
    inputfile.open("INPUT/GeneratorInput"+to_string(gn)+".txt",ios_base::in); //open the GeneratorInput(with the value gn)
    if(inputfile.is_open())									//check if this text file is opened successfully
    {
        std::string line;									//declare variable line, this will store each line read from the input file
        while (std::getline(inputfile, line))		//a loop that starts reading each line of the input from inputfile and assigns it to line
        {
            vector<string> row_values;							//

            vector<std::string> output = GetData(line, '\t'); //call GetData function passing the line and \t as arguments
            if(output.size() == 2)				// It retrieves the data from the line by splitting it using the tab character as a delimiter. The data is stored in the output vector.
            {
                string title = output[0];
                string Strvalue = output[1];
                int value = stoi(Strvalue);
                if(title == "Input_TripsNumber:")					//starts reading the generator input information
                {
                    SetTripsNumber(value);
                }
                else if(title == "Input_CheckpointsNumber:")
                {
                    SetCheckpointsNumber(value);
                }
                else if(title == "Input_StartTimeInterval:")
                {
                    double dValue = stod(Strvalue);
                    SetStartTimeInterval(dValue);
                }
                else if(title == "Input_VehiclesCapacity:")
                {
                    SetVehiclesCapacity(value);
                }
                else if(title == "Input_StopB:")
                {
                    double dValue = stod(Strvalue);
                    SetStopB(dValue);
                }
                else if(title == "Input_RequestNumber:")
                {
                    SetRequestNumber(value);
                }
                else if(title == "Input_SuspendTime:")
                {
                    SetSuspendTime(value);
                }
                else if(title == "Input_TetaMethod:")
                {
                    SetTetaMethod(static_cast<Teta_Method>(value));
                }
                else if(title == "Input_TetaIncrementalValue:")
                {
                    SetTetaIncrementalValue(value);
                }
                else if(title == "Input_ServiceAreaX:")
                {
                    SetServiceAreaX(value);
                }
                else if(title == "Input_ServiceAreaY:")
                {
                    SetServiceAreaY(value);
                }
                else if(title == "Input_PD:")
                {
                    SetRequestPercentagePD(value);
                }
                else if(title == "Input_PND:")
                {
                    SetRequestPercentagePND(value);
                }
                else if(title == "Input_NPD:")
                {
                    SetRequestPercentageNPD(value);
                }
                else if(title == "Input_NPND:")
                {
                    SetRequestPercentageNPND(value);
                }
                else if(title == "Input_W0:")
                {
                    SetW0(value);
                }
                else if(title == "Input_W1:")
                {
                    SetW1(value);
                }
                else if(title == "Input_W2:")
                {
                    SetW2(value);
                }
                else if(title == "Normal:")
                {
                    SetNormal(value);
                }
                else if(title == "SlackTime:")
                {
                    double dValue = stod(Strvalue);
                    SetSlackTime(dValue);
                }
                else if(title == "TimeHorizon:")
                {
                    double dValue = stod(Strvalue);
                    SetTimeHorizons(dValue);
                }
            }
            else
            {
                std::cerr << "Output not equal with 2" << std::endl;
                return false;
            }
        }
        inputfile.close();
    }
    else
    {
        std::cerr << "INPUT/GeneratorInput.txt File Not Found ... :(" << std::endl;
        return false;
    }
    return true;
}

vector<std::string> GeneratorTerminal::GetData(const string &input, char delim) //
{
    std::stringstream ss;
    ss.str(input);								//assigns input string
    std::string item;							//
    vector<string> elem;						//
    while(std::getline(ss, item, delim))		//
    {
        elem.push_back(item);
    }
    return elem;
}
