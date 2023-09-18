#ifndef GENERATORTERMINAL_H
#define GENERATORTERMINAL_H
#include "GeneratorBase.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <sstream>

class GeneratorTerminal : public GeneratorBase
{
public:
    GeneratorTerminal();
    ~GeneratorTerminal(){}

    bool ReadInputFile(int gn = 0);

    vector<std::string> GetData(const std::string& input,char delim);
};

#endif // GENERATORTERMINAL_H
