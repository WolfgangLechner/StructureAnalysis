#ifndef _MOLECULE_H     
#define _MOLECULE_H

#include <iostream>
using namespace std;

const int MAXNUMBEROFNEIGHBORS = 100;

class CMolecule {
  private:
    //The position of the particle
  public:
    CMolecule();
    virtual ~CMolecule();
    double posx,posy,posz;
    double vx,vy,vz;
    double fx,fy,fz;
    double rfx,rfy,rfz;
    double potential;
    
    double sr[3];
    double sv[3];
    
    void set_position(double,double,double);
    double get_posx();
    double get_posy();
    double get_posz();
    int neighbors[MAXNUMBEROFNEIGHBORS];
    double neighbordist[MAXNUMBEROFNEIGHBORS];
    int n_neighbors;
    int cell;
    double realQ4[9],imgQ4[9];
    double realQ6[13],imgQ6[13];
    double arealQ4[9],aimgQ4[9];
    double arealQ6[13],aimgQ6[13];
    
    double Q6,Q4,W4,W6,AQ6,AQ4,AW4,AW6;
    double frenkelnumber;
    //belongs to which cluster
    int belongsto;
    int issolid;
    int structure;
};

#endif
