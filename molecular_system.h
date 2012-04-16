#ifndef _MOLECULAR_SYSTEM_H
#define _MOLECULAR_SYSTEM_H

#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "molecule.h"
#include "parameter.h"

using namespace std;

//Value used to indicate end of cell lists
const int nilvalue = 33333333;
//Pi
const double pi = 3.141592653589793;

//Each cell has 26 neighbors (n**3)-1
struct SNeighborCells{
  int neighbors[26];
};


struct SConstants{
  //Constants, used for the calculation of the order parameters
  //see J. Chem. Phys. 129, 114707 (2008)
  int factorials[20];
  double WignerSymbol6[16];
  double WignerSymbol4[9];
  double WignerPermutations6[16];
  double WignerPermutations4[9];
  int m1[16],m2[16],m3[16];
  //q4,q6,w4,w6 values of the perfect fcc,bcc,sc,hcp,liq;
  double kfcc[4],kbcc[4],khcp[4],kliq[4],ksc[4];
};


class NumberOfParticlesNotDefinedException
{
};

class CMolecularSystem {
  private:
    //these two lists are needed for the neighbor-lists
    int *cellList;
    int *cellHead;
    //saves the 26 neighbors of each cell
    SNeighborCells *neighborcells;
    //converst int to string
    void IntToString(int , string& );
    //checks wheter cell-lists are filled
    int is_listcreated;
    //Create cell lists in the first place
    void createLists();
    //Get cell of particular particle
    int get_cellByIndex(int);
    // Fills the neighborcells[] array
    void fillNeighborCells();
    //n=1,2,... -> (nx,ny,nz)
    void convertIndex(int,int&,int&,int&);
    //(nx,ny,nz) -> n=1,2,...
    int convertIndexToCell(int,int,int);
    //x,y,z -> r, theta,phi
    void convert_SphericalCoordinates(double , double , double , double &, double &, double &);
    //Sperical harmonics
    //see Numerical Recepies for C, page 246, Chapter 6.8
    double PLM(int, int, double);
    void YLM(int , int , double , double , double &, double &);
    //q_lm see J. Chem. Phys. 129, 114707 (2008)
    void QLM(int ,int ,int ,int ,double &, double & );
    

    //Fill the constants arrays
    void fillConstants();
    //Get all neighboring particles within neighbordistance
    void get_AllNeighbors();
    //get q4,q6 from realq4[],realq6[]
    void calculate_qFromComplexVector();

    //get AQ4, AQ6
    void calculate_averagedqFromComplexVector();

    //get realQ4,realQ6
    void calculate_complexQLM();
    
    //get arealQ4,arealQ6
    void calculate_averageComplexQLM();
    //get global Q4, global Q6
    void calculate_globalqFromComplexVector();
    //Get "Frenkel-number", this is the scalar product of q_lm of
    //neighboring Particles
    double get_NumberFromBond(int,int);
    //Get "Frenkel-number" of each particle; This is the number
    //of solid bonds
    void calculate_frenkelNumbers();
    //Create the cluster recursively
    void harvestCluster(int ,int,int );
    int buildClusters(int);
    //the cluster criterium, "Is this particle a solid particle?"
    int clusterCriterium(int,int);

   //get distance between two particles with nearest image
    void get_distance(int,int,double&,double&,double&);
    //get distance between two particles with nearest image
    void get_distancePosition(int,double, double,double,double&,double&,double&);
    //get absolute distance between two particles
    double get_absDistance(int,int);
    //Applies periodic boundary conditions on particle i
    void makeperiodic(int);

    //*********************************************************
    //Miscellaneous
    //*********************************************************
    //Histogramm of local bond order parameters
    int q4Histo[OrderParameterHistoSize];
    int q6Histo[OrderParameterHistoSize];
    int aq4Histo[OrderParameterHistoSize];
    int aq6Histo[OrderParameterHistoSize];

    int q46Histo3d[OrderParameterHistoSize][OrderParameterHistoSize];
    int aq46Histo3d[OrderParameterHistoSize][OrderParameterHistoSize];
    //and its normalization
    int normq4,normq6,normaq4,normaq6;

    double QBCC[Q4SIZE][Q6SIZE];
    double QFCC[Q4SIZE][Q6SIZE];
    double QHCP[Q4SIZE][Q6SIZE];
    double QLIQ[Q4SIZE][Q6SIZE];
    
    //General purpose function to get histo bin for a given value in
    //a given interval
    int get_HistoBox(double,double,double,int);
    //inverse of get_HistoBox
    double get_xFromBox(int,double,double,int);

  public:
    //the main object where all properties of all particles are saved
    CMolecularSystem();
    virtual ~CMolecularSystem();
    //Properties of one single molecule
    CMolecule* molecules;
    CParameter* parameter;
    //all constants
    SConstants constants;
    
    //Init the system
    void initializeMolecules(int);
    void initializeMolecules();
    //and delete the System afterwards
    void deleteMolecules();

    //System can be initialized from a perfect Fcc, Bcc or Hcp
    //structure or reads in a xyz-File
    void InitializeSystem();

    void readXYZFile();
    void createFCC();
    void createBCC();
    void createHCP();

    void readParticleFile();
    void setNumberOfCells();
    //**********************************************************
    //Preparation
    //**********************************************************
    //update the Cell-Lists
    void fillLists();
    //get a random number
    double randnumber();

    //**********************************************************
    //Structure analysis
    //**********************************************************
    //calculate bond order parameters
    void calculate_bondOrderParameter();
    
    //get global Q4, global Q6
    double globalQ6,globalQ4;
    //find the greatest cluster
    int get_greatestCluster(int, int&);
    //update the q4Histo,...
    void update_bondOrderHistogramms();

    //**********************************************************
    //Output
    //**********************************************************

    //Results, structure of each particle
    void outputSimpleResults(int);


    //simplest version of xyz file
    void outputSimpleVMDXYZFile(int);

    //output q4histo.*.txt
    void outputBondOrderHistogramms(int);
    //output the OrderParameter histograms
    void writeQ4Q6Histograms();
    void clear_Q4Q6Histograms();
    void sample_Q4Q6Histograms();
    void norm_Q4Q6Histograms();


};

#endif
