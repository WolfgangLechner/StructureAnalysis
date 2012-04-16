#ifndef _PARAMETER_H
#define _PARAMETER_H

#include <iostream>
#include <exception>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

const int NUMBER_OF_PARAMETERS = 200;
const string CREATEFCCSTRING = "CREATEFCC";
const string CREATEBCCSTRING = "CREATEBCC";
const string CREATEHCPSTRING = "CREATEHCP";
const double PI = 3.14159265;
const int OrderParameterHistoSize = 1000;
const double OrderParameterHistoStart = -1.0;
const double  OrderParameterHistoStop = 1.0;
const int Q4SIZE = 200;
const int Q6SIZE = 600;
const double Q4FACTOR = 1000.0;
const double Q6FACTOR = 1000.0;
 
const int SLIQ = 0;
const int SFCC = 1;
const int SBCC = 2;
const int SHCP = 3;
const int SUND = -1;
const int SOBER = 4;


class ExceptionBadLineInParameterFile{

};

struct s_rawParameter
{
  string name;
  string value;
};

class CParameter {
  private:
    s_rawParameter rawParameter[NUMBER_OF_PARAMETERS];
  public:
    CParameter();
    virtual ~CParameter();
    void readParameter(int);
    void checkParameter();    
    //Number of Particles
    int nop;
    int baseunit;
    int minfrenkel;
    //Number of Cells
    int noc;
    //Boxsize x and y
    double boxx, boxy, boxz;
    //XYZFile
    string xyzFile;
    double neighbordistance;
 };


#endif
