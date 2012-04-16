
#include "molecular_system.h"
#include "molecule.h"


CMolecularSystem::CMolecularSystem()
{
  this->fillConstants();
  this->parameter = new CParameter;
  this->parameter->nop = -1;
  this->is_listcreated = -1;
  for (int i = 0;i<OrderParameterHistoSize;i++)
  {
    this->q4Histo[i] = 0;
    this->q6Histo[i] = 0;
    this->aq4Histo[i] = 0;
    this->aq6Histo[i] = 0;
    this->normq4 = 0;
    this->normq6 = 0;
    this->normaq4 = 0;
    this->normaq6 = 0;
  }
  for (int i = 0;i<OrderParameterHistoSize;i++)
  {
     for (int j = 0;j<OrderParameterHistoSize;j++)
     {
       this->q46Histo3d[i][j] = 0;
       this->aq46Histo3d[i][j] = 0;

     }
  }
  clear_Q4Q6Histograms();
}

CMolecularSystem::~CMolecularSystem()
{
  delete parameter;
  delete [] molecules;
  delete [] cellHead;
  delete [] cellList;
}

//Initialize Simulation Box with number of Particles given
//the nop will be set to the given number
void CMolecularSystem::initializeMolecules(int numberOfParticles)
{
  this->molecules   = new CMolecule[numberOfParticles];
  this->parameter->nop = numberOfParticles;
}

//Assuming nop is allready set, the simulation box is initialized
void CMolecularSystem::initializeMolecules()
{
  try {
    if (this->parameter->nop <0) {throw NumberOfParticlesNotDefinedException();}
    else {
      this->molecules   = new CMolecule[this->parameter->nop];

    }
  } catch(NumberOfParticlesNotDefinedException& ) {}

}

//Frees the memory
void CMolecularSystem::deleteMolecules()
{
  delete [] molecules;
}


// Initialize the System depending on the choice of xyzFile.
void CMolecularSystem::InitializeSystem()
{
  if (this->parameter->xyzFile.compare("notset") == 0) {
    cerr << "xyzFilename must be set in parameter.txt\nif you want to create a fcc use CREATEFCC as filename\n";}
  else
  {
    if (this->parameter->xyzFile.compare(CREATEFCCSTRING) == 0)
    {
      this->createFCC();
    } else if (this->parameter->xyzFile.compare(CREATEBCCSTRING) == 0)
    {
      this->createBCC();
    } else if (this->parameter->xyzFile.compare(CREATEHCPSTRING) == 0)
    {
      this->createHCP();
    } else {
      this->readParticleFile();
    }


  }
}


//-------------------------------------------------------------------------------------------------------
// INITIALIZATION FROM FILE OR STRUCUTURE
//-------------------------------------------------------------------------------------------------------


//Reads configuration from file
void CMolecularSystem::readParticleFile()
{
  double posx,posy,posz;
  int nop;
  string line;
  ifstream confFile;
  confFile.open(this->parameter->xyzFile.c_str(),ifstream::in);
  if (confFile.is_open())
  {
    //the first line contains the number of particles
    getline(confFile,line);
    nop = atoi(line.c_str());
    this->parameter->nop = nop;
    this->initializeMolecules(nop);

    //Followed by boxx,boxy,boxz
    getline(confFile,line);
    this->parameter->boxx = atof(line.c_str());
    getline(confFile,line);
    this->parameter->boxy = atof(line.c_str());
    getline(confFile,line);
    this->parameter->boxz = atof(line.c_str());

    //so lets read the particles positions
    for (int ti = 0;ti<nop;ti++)
    {
      getline(confFile,line);
      string::size_type pos  = line.find(' ');
      if (pos == string::npos) break;
      posx = strtod(line.substr(0, pos ).c_str(),NULL);

      string::size_type lpos = pos + 1;
      pos  = line.find(' ', lpos);
      if (pos == string::npos) break;
      posy = strtod(line.substr(lpos, pos - lpos).c_str(),NULL);
      lpos = pos+1;

      posz = strtod(line.substr(lpos, line.length() - lpos).c_str(),NULL);

      this->molecules[ti].posx = posx;
      this->molecules[ti].posy = posy;
      this->molecules[ti].posz = posz;
    }
  }
  else
  {
    cerr << "Fatal Error : cannot open the file " <<  this->parameter->xyzFile << "\n";
  }
  setNumberOfCells();
}


//Calculates the largest possible number of cells from the neighbordistance and the boxlength
//if not set by user in the parameterfile
void CMolecularSystem::setNumberOfCells()
{
  if (this->parameter->noc < 0.0)
    {
      this->parameter->noc = int(floor((this->parameter->boxx / this->parameter->neighbordistance)))-2;
      if (this->parameter->noc < 3)
      {
        this->parameter->noc = 3;
      }
    }
}


//Creates a set of molecules in FCC structure. It assums that baseunit is set in parameter.txt
//and overrides (and overwrites) the number of particles (if set before)
void CMolecularSystem::createFCC()
{
  int t;
  int bu;
  bu = this->parameter->baseunit;
  if (bu != -1)
  {
    if (parameter->nop != -1) {cerr << "how comes that numberofparticles is already set in createfcc?\n";}
    this->parameter->nop = 4*bu*bu*bu;
    this->initializeMolecules(4*bu*bu*bu);
    for (int jz=0; jz < 2*bu ; jz++)
    {
      for (int jy=0; jy < 2*bu ; jy++)
      {
        for (int jx = 0; jx < bu; jx++)
        {
          t = 2*bu*bu*jz + bu*jy + jx;
          this->molecules[t].posx = double(jx)*this->parameter->boxx/double(bu) + this->parameter->boxx/double(2*bu) * ((jz+jy) % 2);
          this->molecules[t].posy = double(jy)*this->parameter->boxy/double(2*bu);
          this->molecules[t].posz = double(jz)*this->parameter->boxz/double(2*bu);
        }
      }
    }
  } else { cerr << "baseunit is not set, which is needed for createFCC\n";}
  setNumberOfCells();
}

//Creates a set of molecules in BCC structure. It assums that baseunit is set in parameter.txt
//and overrides (and overwrites) the number of particles (if set before)
void CMolecularSystem::createBCC()
{
  int t;
  int bu;
  bu = this->parameter->baseunit;
  if (bu != -1)
  {
    if (parameter->nop != -1) {cerr << "how comes that numberofparticles is already set in createbcc?\n";}
    this->parameter->nop = 2*bu*bu*bu;
    this->initializeMolecules(2*bu*bu*bu);
    for (int jz=0; jz < 2*bu ; jz++)
    {
      for (int jy=0; jy < bu ; jy++)
      {
        for (int jx = 0; jx < bu; jx++)
        {
          t = bu*bu*jz + bu*jy + jx;
          this->molecules[t].posx = double(jx)*this->parameter->boxx/double(bu) + this->parameter->boxx/double(2*bu) * (jz % 2);
          this->molecules[t].posy = double(jy)*this->parameter->boxy/double(bu) + this->parameter->boxy/double(2*bu) * (jz % 2);
          this->molecules[t].posz = double(jz)*this->parameter->boxz/double(2*bu);
        }
      }
    }
  } else { cerr << "baseunit is not set, which is needed for createBCC\n";}
  setNumberOfCells();
}



//Creates a set of molecules in HCP structure. It assums that baseunit is set in parameter.txt
//and overrides (and overwrites) the number of particles (if set before)
void CMolecularSystem::createHCP()
{
  double a,c,h;
  int bu,t;
  bu = this->parameter->baseunit;

  a = this->parameter->boxx/(double(bu));
  c = a*sqrt(3.0)/2.0;
  h = a*sqrt(2.0)/sqrt(3.0);
  if (bu != -1)
  {
    if (parameter->nop != -1) {cerr << "how comes that numberofparticles is already set in createHCP?\n";}
    this->parameter->nop = bu*bu*bu;
    this->initializeMolecules(bu*bu*bu);


    int jz = 0;
    for (int tz=0; tz < bu/2 ; tz++)
    {
      for (int jy=0; jy < bu ; jy++)
      {
        for (int jx=0; jx < bu ; jx++)
        {
          t  = bu*bu*jz + bu*jy + jx;
          this->molecules[t].posx = double(jx)*a + a/2.0*(jy % 2);
          this->molecules[t].posy = double(jy)*c;
          this->molecules[t].posz = double(jz)*h;
        }
      }
      jz += 1;

      for (int jy=0; jy < bu ; jy++)
      {
        for (int jx=0; jx < bu ; jx++)
        {
          t  = bu*bu*jz + bu*jy + jx;
          this->molecules[t].posx = double(jx)*a + a/2.0*((jy+1) % 2);
          this->molecules[t].posy = double(jy)*c + c/3.0;;
          this->molecules[t].posz = double(jz)*h;
        }
      }
      jz += 1;
    }
  }
  setNumberOfCells();
}



//Wigner-Symbols are hard coded
void CMolecularSystem::fillConstants()
{
  //Wigner ThreeJSymbol for l = 6
  constants.WignerSymbol6[0]  =  0.0511827; //6 -6  0
  constants.WignerSymbol6[1]  = -0.0957541; //6 -5 -1
  constants.WignerSymbol6[2]  =  0.129115;  //6 -4 -2
  constants.WignerSymbol6[3]  = -0.141438;  //6 -3 -3
  constants.WignerSymbol6[4]  =  0.127957;  //5 -5  0
  constants.WignerSymbol6[5]  = -0.106079;  //5 -4 -1
  constants.WignerSymbol6[6]  =  0.0408297; //5 -3 -2
  constants.WignerSymbol6[7]  =  0.0186119; //4 -4  0
  constants.WignerSymbol6[8]  =  0.0688184; //4 -3 -1
  constants.WignerSymbol6[9]  = -0.104459;  //4 -2 -2
  constants.WignerSymbol6[10] = -0.100039;  //3 -3 -0
  constants.WignerSymbol6[11] =  0.0452321; //3 -2 -1
  constants.WignerSymbol6[12] =  0.0511827; //2 -2  0
  constants.WignerSymbol6[13] = -0.0953576; //2 -1 -1
  constants.WignerSymbol6[14] =  0.0465298; //1 -1  0
  constants.WignerSymbol6[15] =  -0.0930595; //0  0  0

  //number of permutations for l=6
  constants.WignerPermutations6[0]  = 6.0; //6  -6   0
  constants.WignerPermutations6[1]  = 6.0; //6  -5  -1
  constants.WignerPermutations6[2]  = 6.0; //6  -4  -2
  constants.WignerPermutations6[3]  = 3.0; //6  -3  -3
  constants.WignerPermutations6[4]  = 6.0; //5  -5   0
  constants.WignerPermutations6[5]  = 6.0; //5  -4  -1
  constants.WignerPermutations6[6]  = 6.0; //5  -3  -2
  constants.WignerPermutations6[7]  = 6.0; //4  -4   0
  constants.WignerPermutations6[8]  = 6.0; //4  -3  -1
  constants.WignerPermutations6[9]  = 3.0; //4  -2  -2
  constants.WignerPermutations6[10] = 6.0; //3  -3   0
  constants.WignerPermutations6[11] = 6.0; //3  -2  -1
  constants.WignerPermutations6[12] = 6.0; //2  -2   0
  constants.WignerPermutations6[13] = 3.0; //2  -1  -1
  constants.WignerPermutations6[14] = 6.0; //1  -1   0
  constants.WignerPermutations6[15] = 1.0; //0   0   0


  //l = 4
  constants.WignerSymbol4[0]  =  0.104298; //4 -4  0
  constants.WignerSymbol4[1]  = -0.164909; //4 -3 -1
  constants.WignerSymbol4[2]  =  0.186989; //4 -2 -2
  constants.WignerSymbol4[3] =   0.156447;  //3 -3 -0
  constants.WignerSymbol4[4] =  -0.0623298; //3 -2 -1
  constants.WignerSymbol4[5] =  -0.0819482; //2 -2  0
  constants.WignerSymbol4[6] =   0.141351;  //2 -1 -1
  constants.WignerSymbol4[7] =  -0.0670485; //1 -1  0
  constants.WignerSymbol4[8] =   0.134097;  //0  0  0

  constants.WignerPermutations4[0]  = 6.0; //4  -4   0
  constants.WignerPermutations4[1]  = 6.0; //4  -3  -1
  constants.WignerPermutations4[2]  = 3.0; //4  -2  -2
  constants.WignerPermutations4[3] = 6.0; //3  -3   0
  constants.WignerPermutations4[4] = 6.0; //3  -2  -1
  constants.WignerPermutations4[5] = 6.0; //2  -2   0
  constants.WignerPermutations4[6] = 3.0; //2  -1  -1
  constants.WignerPermutations4[7] = 6.0; //1  -1   0
  constants.WignerPermutations4[8] = 1.0; //0   0   0


  constants.m1[0]  = -6; constants.m2[0]  = 6; constants.m3[0] = 0;
  constants.m1[1]  = -6; constants.m2[1]  = 5; constants.m3[1] = 1;
  constants.m1[2]  = -6; constants.m2[2]  = 4; constants.m3[2] = 2;
  constants.m1[3]  = -6; constants.m2[3]  = 3; constants.m3[3] = 3;

  constants.m1[4]  = -5; constants.m2[4]  = 5; constants.m3[4] = 0;
  constants.m1[5]  = -5; constants.m2[5]  = 4; constants.m3[5] = 1;
  constants.m1[6]  = -5; constants.m2[6]  = 3; constants.m3[6] = 2;

  constants.m1[7]  = -4; constants.m2[7]  = 4; constants.m3[7] = 0;
  constants.m1[8]  = -4; constants.m2[8]  = 3; constants.m3[8] = 1;
  constants.m1[9]  = -4; constants.m2[9]  = 2; constants.m3[9] = 2;

  constants.m1[10] = -3; constants.m2[10] = 3; constants.m3[10] = 0;
  constants.m1[11] = -3; constants.m2[11] = 2; constants.m3[11] = 1;

  constants.m1[12] = -2; constants.m2[12] = 2; constants.m3[12] = 0;
  constants.m1[13] = -2; constants.m2[13] = 1; constants.m3[13] = 1;

  constants.m1[14] = -1; constants.m2[14] = 1; constants.m3[14] = 0;
  constants.m1[15] = 0;  constants.m2[15] = 0; constants.m3[15] = 0;

  constants.factorials[0]  = 1;
  constants.factorials[1]  = 1;
  constants.factorials[2]  = 1*2;
  constants.factorials[3]  = 1*2*3;
  constants.factorials[4]  = 1*2*3*4;
  constants.factorials[5]  = 1*2*3*4*5;
  constants.factorials[6]  = 1*2*3*4*5*6;
  constants.factorials[7]  = 1*2*3*4*5*6*7;
  constants.factorials[8]  = 1*2*3*4*5*6*7*8;
  constants.factorials[9]  = 1*2*3*4*5*6*7*8*9;
  constants.factorials[10] = 1*2*3*4*5*6*7*8*9*10;
  constants.factorials[11] = 1*2*3*4*5*6*7*8*9*10*11;
  constants.factorials[12] = 1*2*3*4*5*6*7*8*9*10*11*12;

  //0 = Q4
  //1 = Q6;
  //2 = W4;
  //3 = W6;
  constants.kfcc[0] = 0.19094;
  constants.kfcc[1] = 0.57452;
  constants.kfcc[2] = -0.15932;
  constants.kfcc[3] = -0.01316;

  constants.khcp[0] = 0.09722;
  constants.khcp[1] = 0.48476;
  constants.khcp[2] = 0.13410;
  constants.khcp[3] = -0.01244;

  constants.ksc[0] = 0.76376;
  constants.ksc[1] = 0.35355;
  constants.ksc[2] = 0.15932;
  constants.ksc[3] = 0.01316;

  constants.kbcc[0] = 0.0363696;
  constants.kbcc[1] = 0.5107;
  constants.kbcc[2] = 0.15932;
  constants.kbcc[3] = 0.01316;

  constants.kliq[0] = 0.0;
  constants.kliq[1] = 0.0;
  constants.kliq[2] = 0.0;
  constants.kliq[3] = 0.0;


}


//-------------------------------------------------------------------------------------------------------
// HELPERS
//-------------------------------------------------------------------------------------------------------

//gets a random number between 0 .. 1
double CMolecularSystem::randnumber() {
return double(rand())/RAND_MAX;}

//converts a int into a string
void CMolecularSystem::IntToString(int i, string& res)
{
    ostringstream tempel;
    tempel << i;
    res = tempel.str();
}

//-------------------------------------------------------------------------------------------------------
//PREPARATION
//-------------------------------------------------------------------------------------------------------


//converts the cycling id into x, y and z id's
void CMolecularSystem::convertIndex(int index,int &idx, int &idy, int &idz)
{
  int noc;
  noc = this->parameter->noc;
  idx = index%noc;
  idy = (index/noc)%noc;
  idz = index/(noc*noc);
}

//returns the cell id of x-id, y-id and z-id
//this is used for to get the 26 neighboring cells
int CMolecularSystem::convertIndexToCell(int x, int y, int z)
{
  int noc;
  noc = this->parameter->noc;
  return  (x  + noc)%noc+ (y  + noc)%noc* noc+ (z  + noc)%noc * noc*noc;
}

//Gets all neighboring cells of each cell
void CMolecularSystem::fillNeighborCells()
{
  int x,y,z;
  int noc;
  noc = this->parameter->noc;
  for (int ci = 0;ci<noc*noc*noc;ci++)
  {
    this->convertIndex(ci,x, y, z);
    neighborcells[ci].neighbors[0] = convertIndexToCell(x-1,y-1,z-1);
    neighborcells[ci].neighbors[1] = convertIndexToCell(x  ,y-1,z-1);
    neighborcells[ci].neighbors[2] = convertIndexToCell(x+1,y-1,z-1);

    neighborcells[ci].neighbors[3] = convertIndexToCell(x-1,y  ,z-1);
    neighborcells[ci].neighbors[4] = convertIndexToCell(x  ,y  ,z-1);
    neighborcells[ci].neighbors[5] = convertIndexToCell(x+1,y  ,z-1);

    neighborcells[ci].neighbors[6] = convertIndexToCell(x-1,y+1,z-1);
    neighborcells[ci].neighbors[7] = convertIndexToCell(x  ,y+1,z-1);
    neighborcells[ci].neighbors[8] = convertIndexToCell(x+1,y+1,z-1);

    neighborcells[ci].neighbors[9]  = convertIndexToCell(x-1,y-1,z);
    neighborcells[ci].neighbors[10] = convertIndexToCell(x  ,y-1,z);
    neighborcells[ci].neighbors[11] = convertIndexToCell(x+1,y-1,z);

    neighborcells[ci].neighbors[12] = convertIndexToCell(x-1,y  ,z);
    neighborcells[ci].neighbors[13] = convertIndexToCell(x+1,y  ,z);
    neighborcells[ci].neighbors[14] = convertIndexToCell(x-1,y+1,z);
    neighborcells[ci].neighbors[15] = convertIndexToCell(x  ,y+1,z);
    neighborcells[ci].neighbors[16] = convertIndexToCell(x+1,y+1,z);

    neighborcells[ci].neighbors[17] = convertIndexToCell(x-1,y-1,z+1);
    neighborcells[ci].neighbors[18] = convertIndexToCell(x  ,y-1,z+1);
    neighborcells[ci].neighbors[19] = convertIndexToCell(x+1,y-1,z+1);

    neighborcells[ci].neighbors[20] = convertIndexToCell(x-1,y  ,z+1);
    neighborcells[ci].neighbors[21] = convertIndexToCell(x  ,y  ,z+1);
    neighborcells[ci].neighbors[22] = convertIndexToCell(x+1,y  ,z+1);

    neighborcells[ci].neighbors[23] = convertIndexToCell(x-1,y+1,z+1);
    neighborcells[ci].neighbors[24] = convertIndexToCell(x  ,y+1,z+1);
    neighborcells[ci].neighbors[25] = convertIndexToCell(x+1,y+1,z+1);
  }
}

//mallocs the lists
void CMolecularSystem::createLists()
{
  int noc;
  if (this->is_listcreated == -1)
  {
    noc = this->parameter->noc;
    cellList = new int[this->parameter->nop];
    cellHead = new int[noc*noc*noc];
    neighborcells = new SNeighborCells[noc*noc*noc];
    fillNeighborCells();
    this->is_listcreated = 1;

  }
}


//Calculates into which cell the particle i belongs
int CMolecularSystem::get_cellByIndex(int i)
{
  int cell,noc;
  double boxx,boxy,boxz;
  double x,y,z;
  boxx = this->parameter->boxx;
  boxy = this->parameter->boxy;
  boxz = this->parameter->boxz;
  noc = this->parameter->noc;

  x = this->molecules[i].posx;
  y = this->molecules[i].posy;
  z = this->molecules[i].posz;
  //cout << x << " " << y << " " << z << "\n";
  cell =  int (floor((x / boxx) * noc) +
               floor((y / boxy) * noc) * noc +
               floor((z / boxz) * noc) * noc * noc);
  return cell;
}

//Adds the particles to the cell-Lists
void CMolecularSystem::fillLists()
{
  int cell,allCells,noc;
  createLists();

  noc = this->parameter->noc;
  allCells = noc*noc*noc;
  for(int ti=0;ti<this->parameter->nop;ti++)
  {
    cellList[ti] = nilvalue;
  }
  for (int ci=0;ci<allCells;ci++)
  {
    cellHead[ci] = nilvalue;
  }
  for(int ti=0;ti<this->parameter->nop;ti++)
  {
    cell = get_cellByIndex(ti);
    this->molecules[ti].cell = cell;
    cellList[ti] = cellHead[cell];
    cellHead[cell] = ti;
  }
  
}

//Get all nearest neighbors within this->neighbordistance
//mainly used for the bond orders
void CMolecularSystem::get_AllNeighbors()
{
  double nd,d;
  nd = this->parameter->neighbordistance;
  int c,nop;
  nop = this->parameter->nop;

  for (int ti = 0;ti<nop;ti++)
  {
    for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++)
    {
      this->molecules[ti].neighbors[tn] = nilvalue;
      this->molecules[ti].neighbordist[tn] = -1.0;
    }
  }


  for (int ti = 0;ti<nop;ti++)
  {
    c = 0;
	  for (int tj = 0;tj<nop;tj++)
	  {

          if (tj != ti) { 
            d = get_absDistance(ti,tj); 
            if (d < nd) {
              this->molecules[ti].neighbors[c] = tj; 
              this->molecules[ti].neighbordist[c]=d; 
              c += 1;   
            }
         }
	  }
    this->molecules[ti].n_neighbors = c;
    cout << ti << " " << c << "\n";
  }

}

////Get all nearest neighbors within this->neighbordistance
////mainly used for the bond orders
//void CMolecularSystem::get_AllNeighbors()
//{
//  double nd,d;
//  nd = this->parameter->neighbordistance;
//  int cell, particle, maincell,c,nop;
//  nop = this->parameter->nop;
//  for (int ti = 0;ti<nop;ti++)
//  {
//    for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++)
//    {
//      this->molecules[ti].neighbors[tn] = nilvalue;
//      this->molecules[ti].neighbordist[tn] = -1.0;
//    }
//  }
//
//  for (int ti = 0;ti<nop;ti++)
//  {
//    c = 0;
//    maincell = this->molecules[ti].cell;
//    cell = this->molecules[ti].cell;
//    particle = cellHead[cell];
//    //first the main cell
//    if (particle != ti) { d = get_absDistance(ti,particle); if (d < nd) {this->molecules[ti].neighbors[c] = particle; this->molecules[ti].neighbordist[c]=d; c += 1;   }}
//    do {
//      particle = cellList[particle];
//      if (particle != ti) { d = get_absDistance(ti,particle); if (d < nd) {this->molecules[ti].neighbors[c] = particle; this->molecules[ti].neighbordist[c]=d; c += 1;   }}
//    } while (cellList[particle] != nilvalue);
//    //and all other 26 neighborcells
//    for (int ci = 0;ci<26;ci++)
//    {
//      cell = this->neighborcells[maincell].neighbors[ci];
//      particle = cellHead[cell];
//      if (particle != ti) { d = get_absDistance(ti,particle); if (d < nd) {this->molecules[ti].neighbors[c] = particle; this->molecules[ti].neighbordist[c]=d; c += 1;   }}
//      do {
//        particle = cellList[particle];
//        if (particle != ti) { d = get_absDistance(ti,particle); if (d < nd) {this->molecules[ti].neighbors[c] = particle; this->molecules[ti].neighbordist[c]=d; c += 1;   }}
//      } while (cellList[particle] != nilvalue);
//    }
//    this->molecules[ti].n_neighbors = c;
//   }
//
//
//}

//-------------------------------------------------------------------------------------------------------
//Structure analysis
//-------------------------------------------------------------------------------------------------------

//PLM from the Numerical Recipes C
//Chapter 6.8 page 246
double CMolecularSystem::PLM(int l, int m, double x)
{
  double fact,pll,pmm,pmmp1,somx2;
  int i,ll;
  pll = 0.0;
  if (m < 0 || m > l || fabs(x) > 1.0)
  cerr << "impossible combination of l and m" << "\n";
  pmm=1.0;
  if (m > 0)
  {
    somx2=sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for (i=1;i<=m;i++)
    {
      pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else
  {
    pmmp1=x*(2*m+1)*pmm;
    if (l == (m+1))
      return pmmp1;
    else
    {
      for (ll=m+2;ll<=l;ll++)
        {
          pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
          pmm=pmmp1;
          pmmp1=pll;
        }
        return pll;
     }
  }
}


//Converts (x,y,z) to (r,phi,theta)
void CMolecularSystem::convert_SphericalCoordinates(double x, double y, double z, double &r, double &phi, double &theta)
{
  //r = sqrt(x*x + y*y + z*z);
  //if (y >= 0.0) { phi = acos(x/sqrt(x*x + y*y)); } else {phi = 2.0*pi - acos(x/sqrt(x*x + y*y)); }
  //if (fabs(x)<0.00000001 && fabs(y)<0.000000001) phi = 0.0;
  //theta = pi/2 - atan(z/sqrt(x*x + y*y));
  r = sqrt(x*x+y*y+z*z);
  theta = acos(z/r);
  phi = atan2(y,x);
}

void CMolecularSystem::YLM(int l, int m, double theta, double phi, double &realYLM, double &imgYLM)
{
  double factor;
  double m_PLM;
  m_PLM = PLM(l,m,cos(theta));
  factor = ((2.0*double(l) + 1.0)*constants.factorials[l-m]) / (4.0*pi*constants.factorials[l+m]);
  realYLM = sqrt(factor) * m_PLM * cos(double(m)*phi);
  imgYLM  = sqrt(factor) * m_PLM * sin(double(m)*phi);
}

void CMolecularSystem::QLM(int l,int m,int ti,int tj,double &realYLM, double &imgYLM )
{
  double diffx, diffy, diffz;
  double r,phi,theta;
  realYLM = 0.0;
  imgYLM = 0.0;
  this->get_distance(ti,tj,diffx,diffy,diffz);
  this->convert_SphericalCoordinates(diffx, diffy, diffz, r, phi, theta);
  if (m < 0) {
    YLM(l, abs(m), theta, phi, realYLM, imgYLM);
    realYLM = pow(-1.0,m)*realYLM;
    imgYLM = pow(-1.0,m)*imgYLM;
  }
    else
  {
    YLM(l, m, theta, phi, realYLM, imgYLM);
  }
}


//Calculates the complex components of the q4m and q6m vector
void CMolecularSystem::calculate_complexQLM()
{
  //nn = number of neighbors
  int nn,nop;
  double realti,imgti;
  double realYLM,imgYLM;
  nop = this->parameter->nop;
  for (int ti= 0;ti<nop;ti++)
  {
    nn = this->molecules[ti].n_neighbors;
    for (int mi = -6;mi < 7;mi++)
    {
      realti = 0.0;
      imgti = 0.0;
      for (int ci = 0;ci<nn;ci++)
      {
        QLM(6,mi,ti,this->molecules[ti].neighbors[ci],realYLM, imgYLM);
	realti += realYLM;
	imgti += imgYLM;
      }
      realti = realti/(double(nn));
      imgti = imgti/(double(nn));
      this->molecules[ti].realQ6[mi+6] = realti;
      this->molecules[ti].imgQ6[mi+6] = imgti;
    }
    for (int mi = -4;mi < 5;mi++)
    {
      realti = 0.0;
      imgti = 0.0;
      for (int ci = 0;ci<nn;ci++)
      {
        QLM(4,mi,ti,this->molecules[ti].neighbors[ci],realYLM, imgYLM);
        realti += realYLM;
        imgti += imgYLM;
      }
      realti = realti/(double(nn));
      imgti = imgti/(double(nn));
      this->molecules[ti].realQ4[mi+4] = realti;
      this->molecules[ti].imgQ4[mi+4] = imgti;
    }
  }
}


//Calculates local Q's from the q-vector
void CMolecularSystem::calculate_averageComplexQLM()
{
  int nop;
  nop = this->parameter->nop;
  double realNQ6[13],imgNQ6[13];
  double realNQ4[9],imgNQ4[9];

  for (int ti= 0;ti<nop;ti++)
  {
    for (int mi=0;mi<13;mi++)
    {
      realNQ6[mi] = this->molecules[ti].realQ6[mi];
      imgNQ6[mi]  = this->molecules[ti].imgQ6[mi];
      for (int ci = 0;ci<this->molecules[ti].n_neighbors;ci++)
      {
        realNQ6[mi] += this->molecules[this->molecules[ti].neighbors[ci]].realQ6[mi];
        imgNQ6[mi]  += this->molecules[this->molecules[ti].neighbors[ci]].imgQ6[mi];
      }
      realNQ6[mi] = realNQ6[mi]/double(this->molecules[ti].n_neighbors+1);
      imgNQ6[mi]  = imgNQ6[mi]/double(this->molecules[ti].n_neighbors+1);

      this->molecules[ti].arealQ6[mi] = realNQ6[mi];
      this->molecules[ti].aimgQ6[mi] = imgNQ6[mi];
    }

    for (int mi=0;mi<9;mi++)
    {
      realNQ4[mi] = this->molecules[ti].realQ4[mi];
      imgNQ4[mi]  = this->molecules[ti].imgQ4[mi];
      for (int ci = 0;ci<this->molecules[ti].n_neighbors;ci++)
      {
        realNQ4[mi] += this->molecules[this->molecules[ti].neighbors[ci]].realQ4[mi];
        imgNQ4[mi]  += this->molecules[this->molecules[ti].neighbors[ci]].imgQ4[mi];
      }
      realNQ4[mi] = realNQ4[mi]/double(this->molecules[ti].n_neighbors+1);
      imgNQ4[mi]  = imgNQ4[mi]/double(this->molecules[ti].n_neighbors+1);

      this->molecules[ti].arealQ4[mi] = realNQ4[mi];
      this->molecules[ti].aimgQ4[mi] = imgNQ4[mi];
    }

  }

}


//Calculates local Q's from the q-vector
void CMolecularSystem::calculate_qFromComplexVector()
{
  int nop;
  nop = this->parameter->nop;
  double sumQ6,sumQ4;
  for (int ti= 0;ti<nop;ti++)
  {
    sumQ6 = 0.0;
    sumQ4 = 0.0;
    for (int mi = 0;mi < 13;mi++)
    {
      sumQ6 += this->molecules[ti].realQ6[mi]*this->molecules[ti].realQ6[mi] + this->molecules[ti].imgQ6[mi]*this->molecules[ti].imgQ6[mi];
      //if (ti == 0) {cout << ti << " " << mi << " " << this->molecules[ti].realQ6[mi]*this->molecules[ti].realQ6[mi] + this->molecules[ti].imgQ6[mi]*this->molecules[ti].imgQ6[mi] << "\n";}
    }
    sumQ6 = pow(((4.0*pi/13.0) * sumQ6),0.5);
    this->molecules[ti].Q6 = sumQ6;

    for (int mi = 0;mi < 9;mi++)
    {
	sumQ4 += this->molecules[ti].realQ4[mi]*this->molecules[ti].realQ4[mi] + this->molecules[ti].imgQ4[mi]*this->molecules[ti].imgQ4[mi];
    }
    sumQ4 = pow(((4.0*pi/9.0) * sumQ4),0.5);
    this->molecules[ti].Q4 = sumQ4;
  }
}


//Calculates the Averaged Q6 and Averaged Q4 from the complex vector
void CMolecularSystem::calculate_averagedqFromComplexVector()
{
  int nop;
  nop = this->parameter->nop;
  double sumAQ4,sumAQ6;

  for (int ti= 0;ti<nop;ti++)
  {
    sumAQ4 = 0.0;
    for (int hi=0;hi<9;hi++)
    {
      sumAQ4 += this->molecules[ti].arealQ4[hi]*this->molecules[ti].arealQ4[hi]+ this->molecules[ti].aimgQ4[hi]*this->molecules[ti].aimgQ4[hi];
    }
    sumAQ4  = sumAQ4 * 4.0*pi/9.0;
    sumAQ4  = pow(sumAQ4 ,0.5);

    sumAQ6 = 0.0;
    for (int hi=0;hi<13;hi++)
    {
      sumAQ6 += this->molecules[ti].arealQ6[hi]*this->molecules[ti].arealQ6[hi]+ this->molecules[ti].aimgQ6[hi]*this->molecules[ti].aimgQ6[hi];
    }
    sumAQ6  = sumAQ6 * 4.0*pi/13.0;
    sumAQ6  = pow(sumAQ6 ,0.5);

    this->molecules[ti].AQ4 = sumAQ4;
    this->molecules[ti].AQ6 = sumAQ6;
  }

}


//Calculates the Averaged Q6 and Averaged Q4 from the complex vector
void CMolecularSystem::calculate_globalqFromComplexVector()
{
  int nop;
  nop = this->parameter->nop;
  double realNQ6[13],imgNQ6[13];
  double realNQ4[9],imgNQ4[9];
  double sumGQ4,sumGQ6;
  sumGQ4 = 0.0;
  sumGQ6 = 0.0;
  for (int mi=0;mi<13;mi++)
  {
    realNQ6[mi] = 0.0;
    imgNQ6[mi] = 0.0;
  }
  for (int mi=0;mi<9;mi++)
  {
    realNQ4[mi] = 0.0;
    imgNQ4[mi] = 0.0;
  }

  for (int ti= 0;ti<nop;ti++)
  {
    for (int mi=0;mi<13;mi++)
    {
      realNQ6[mi] += this->molecules[ti].realQ6[mi]/double(nop);
      imgNQ6[mi]  += this->molecules[ti].imgQ6[mi]/double(nop);
    }

    for (int mi=0;mi<9;mi++)
    {
      realNQ4[mi] += this->molecules[ti].realQ4[mi]/double(nop);
      imgNQ4[mi]  += this->molecules[ti].imgQ4[mi]/double(nop);
    }
  }


  sumGQ6 = 0.0;
  for (int hi=0;hi<13;hi++)
  {
    sumGQ6 += realNQ6[hi]*realNQ6[hi] + imgNQ6[hi]*imgNQ6[hi];
  }
  sumGQ6  = sumGQ6 * 4.0*pi/13.0;
  sumGQ6  = pow(sumGQ6 ,0.5);

  sumGQ4 = 0.0;
  for (int hi=0;hi<9;hi++)
  {
    sumGQ4 += realNQ4[hi]*realNQ4[hi] + imgNQ4[hi]*imgNQ4[hi];
  }
  sumGQ4  = sumGQ4 * 4.0*pi/9.0;
  sumGQ4  = pow(sumGQ4 ,0.5);

  this->globalQ6 = sumGQ6;
  this->globalQ4 = sumGQ4;


}

//claculates the average of the scalar product of q6 x q6* over all neighbors
double CMolecularSystem::get_NumberFromBond(int ti,int tj)
{
  double sumSquareti,sumSquaretj;
  double realdotproduct,imgdotproduct;
  double connection;
  sumSquareti = 0.0;
  sumSquaretj = 0.0;
  realdotproduct = 0.0;
  imgdotproduct = 0.0;

  for (int mi = 0;mi < 13;mi++)
  {

    sumSquareti += this->molecules[ti].realQ6[mi]*this->molecules[ti].realQ6[mi] + this->molecules[ti].imgQ6[mi] *this->molecules[ti].imgQ6[mi];
    sumSquaretj += this->molecules[tj].realQ6[mi]*this->molecules[tj].realQ6[mi] + this->molecules[tj].imgQ6[mi] *this->molecules[tj].imgQ6[mi];
    realdotproduct += this->molecules[ti].realQ6[mi]*this->molecules[tj].realQ6[mi];
    imgdotproduct  += this->molecules[ti].imgQ6[mi] *this->molecules[tj].imgQ6[mi];
  }
  connection = (realdotproduct+imgdotproduct)/(sqrt(sumSquaretj)*sqrt(sumSquareti));
  return connection;
}


//if the scalar product of q6 x q6* is larger than 0.5 we define a bond as solid
void CMolecularSystem::calculate_frenkelNumbers()
{
  int frenkelcons;
  double scalar;
  for (int ti= 0;ti<this->parameter->nop;ti++)
  {
    frenkelcons = 0;
    for (int c = 0;c<this->molecules[ti].n_neighbors;c++)
    {
      scalar = this->get_NumberFromBond(ti,this->molecules[ti].neighbors[c]);
      if (scalar > 0.5) frenkelcons += 1;
    }
    this->molecules[ti].frenkelnumber = frenkelcons;
  }
}

void CMolecularSystem::writeQ4Q6Histograms()
{
  FILE *ofBCC,*ofFCC,*ofHCP,*ofLIQ;
  
  char IntStr[80];
  sprintf( IntStr, "./output/QBCC.dat") ;
  ofBCC = fopen(IntStr,"w");
  sprintf( IntStr, "./output/QFCC.dat") ;
  ofFCC = fopen(IntStr,"w");
  sprintf( IntStr, "./output/QHCP.dat") ;
  ofHCP = fopen(IntStr,"w");
  sprintf( IntStr, "./output/QLIQ.dat") ;
  ofLIQ = fopen(IntStr,"w");
  
  for (int j = 0;j<Q6SIZE;j++)
  {
    for (int i = 0;i<Q4SIZE;i++)
    {
      fprintf(ofBCC,"%d %d %.18lf\n",i ,j ,QBCC[i][j]);
      fprintf(ofFCC,"%d %d %.18lf\n",i ,j ,QFCC[i][j]);
      fprintf(ofHCP,"%d %d %.18lf\n",i ,j ,QHCP[i][j]);
      fprintf(ofLIQ,"%d %d %.18lf\n",i ,j ,QLIQ[i][j]);
    }
    fprintf(ofBCC,"\n");
    fprintf(ofFCC,"\n");
    fprintf(ofHCP,"\n");
    fprintf(ofLIQ,"\n");
  }
  fclose(ofBCC);
  fclose(ofFCC);
  fclose(ofHCP);
  fclose(ofLIQ);
 
}


void CMolecularSystem::clear_Q4Q6Histograms()
{
  for (int ti = 0;ti<Q4SIZE;ti++)
  {
    for (int tj = 0;tj<Q6SIZE;tj++)
    {
      this->QBCC[ti][tj] = 0.0;
      this->QFCC[ti][tj] = 0.0;
      this->QHCP[ti][tj] = 0.0;
      this->QLIQ[ti][tj] = 0.0;
    }
  }
}

void CMolecularSystem::norm_Q4Q6Histograms()
{
  double normBCC, normHCP, normFCC, normLIQ;
  normBCC = 0.0;
  normHCP = 0.0;
  normFCC = 0.0;
  normLIQ = 0.0;
  for (int ti = 0;ti<Q4SIZE;ti++)
  {
    for (int tj = 0;tj<Q6SIZE;tj++)
    {
      normBCC += QBCC[ti][tj];
      normFCC += QFCC[ti][tj];
      normHCP += QHCP[ti][tj];
      normLIQ += QLIQ[ti][tj];
    }
  }
  for (int ti = 0;ti<Q4SIZE;ti++)
  {
    for (int tj = 0;tj<Q6SIZE;tj++)
    {
      QBCC[ti][tj] /= normBCC;
      QFCC[ti][tj] /= normFCC;
      QHCP[ti][tj] /= normHCP;
      QLIQ[ti][tj] /= normLIQ;
    }
  }
}

void CMolecularSystem::sample_Q4Q6Histograms()
{
  int q4box,q6box;
  for (int ti= 0;ti<this->parameter->nop;ti++)
  {
    q4box = int(this->molecules[ti].AQ4 * Q4FACTOR);
    q6box = int(this->molecules[ti].AQ6 * Q4FACTOR);
    this->QBCC[q4box][q6box] += 1.0;
    this->QFCC[q4box][q6box] += 1.0;
    this->QHCP[q4box][q6box] += 1.0;
    this->QLIQ[q4box][q6box] += 1.0;
  }
}

//calculate histobox of value, for given start, stop, N
int CMolecularSystem::get_HistoBox(double value,double a, double b, int n)
{
  int box;
  box = 0;
  if ((value < a) || (value > b))
  {
   cerr << "value out of Histogram bound, check OrderParameterHistoStart and OrderParameterHistoStop\n";
  } else {
    box = int((value - a)/(b - a)*n);
  }
  return box;
}

//get x from box
double CMolecularSystem::get_xFromBox(int box,double a, double b, int n)
{
  double x;
  x = a + box*(b-a)/double(n);
  return x;
}


//Update alle bond order histogramms
void CMolecularSystem::update_bondOrderHistogramms()
{

  double a,b;
  int n;
  int nop;
  nop = parameter->nop;
  a = OrderParameterHistoStart;
  b = OrderParameterHistoStop;
  n = OrderParameterHistoSize;

  for (int ti= 0;ti<this->parameter->nop;ti++)
  {
    q4Histo[get_HistoBox(this->molecules[ti].Q4,a,b,n)]   += 1;
    q6Histo[get_HistoBox(this->molecules[ti].Q6,a,b,n)]   += 1;
    aq4Histo[get_HistoBox(this->molecules[ti].AQ4,a,b,n)] += 1;
    aq6Histo[get_HistoBox(this->molecules[ti].AQ6,a,b,n)] += 1;
    q46Histo3d[get_HistoBox(this->molecules[ti].Q4,a,b,n)][get_HistoBox(this->molecules[ti].Q6,a,b,n)] += 1;
    aq46Histo3d[get_HistoBox(this->molecules[ti].AQ4,a,b,n)][get_HistoBox(this->molecules[ti].AQ6,a,b,n)] += 1;
  }
  this->normq4 += nop;
  this->normq6 += nop;
  this->normaq4 += nop;
  this->normaq6 += nop;
}

//This is the main procedure that calculates the Bondorderparameters
//this one is the only public procedure to use
//first all we get all neighbors, then we caluclate the q-vector
//from this q-vector the q4,q6,aq4,aq6,w4,w6 and frenkelproduct values
//are taken
void CMolecularSystem::calculate_bondOrderParameter()
{
  //Find all particles within a radius of neighbordistance
  this->get_AllNeighbors();
  //Calculate q_lm vector
  this->calculate_complexQLM();
  //Calculate averaged version of q_lm as in J. Chem. Phys. 129, 114707 (2008)
  this->calculate_averageComplexQLM();
  //From these vectors calculate the q_4 and q_4
  this->calculate_qFromComplexVector();
  //and also aq_4 and aq_6
  this->calculate_averagedqFromComplexVector();
  //and Q_6 and Q_4
  this->calculate_globalqFromComplexVector();
  //and the number of bonds to find the largest cluster
  this->calculate_frenkelNumbers();
}


//This function can be used to change the criterium which decides whether a particle is
//solid or fluid
int CMolecularSystem::clusterCriterium(int ti,int criterium)
{
  int value;
  value = 0;
  if (criterium == 0)
  {
    if (this->molecules[ti].frenkelnumber > parameter->minfrenkel)
    { value = 1; } else {value = 0;}
  }
  return value;
}

//builds the cluster with an recursive loop over all neighbors
void CMolecularSystem::harvestCluster(int ti,int numberofCluster,int criterium)
{
  int c;
  int neighbor;
  c = 0;
  do
  {
    neighbor =  this->molecules[ti].neighbors[c];
    if ((this->molecules[neighbor].belongsto == nilvalue) && (clusterCriterium(neighbor,criterium)))
    {
      this->molecules[neighbor].belongsto = numberofCluster;
      harvestCluster(neighbor,numberofCluster,criterium);
    }
    c += 1;
    } while (this->molecules[ti].neighbors[c] != nilvalue);

}


//builds cluster, all particles that belong to the same cluster have the same .belongsto number
int CMolecularSystem::buildClusters(int criterium)
{
  int numberofCluster;
  int nop;
  nop = parameter->nop;
  numberofCluster = 0;
  for (int ti= 0;ti<nop;ti++)
  {
    this->molecules[ti].belongsto = nilvalue;
  }
  for (int ti= 0;ti<nop;ti++)
  {
    if ((this->molecules[ti].belongsto == nilvalue) && (clusterCriterium(ti,criterium)))
    {
      numberofCluster += 1;
      this->molecules[ti].belongsto = numberofCluster;
      harvestCluster(ti,numberofCluster,criterium);
    }
  }
  return numberofCluster;
}

//Searches for the greatest cluster, previously builds all .belongsto numbers
//gives back the size of the cluter with return, and the greatestbelongsto as int&
int CMolecularSystem::get_greatestCluster(int criterium,int &greatestbelongsto)
{
  int nop;
  nop = parameter->nop;
  int clusterhisto[nop];
  int greatestcluster ;
  greatestcluster = 0;

  this->buildClusters(criterium);

  for (int ti= 0;ti<nop;ti++)
  {
    clusterhisto[ti] = 0;
  }
  //Builds an histo of clustersizes
  for (int ti= 0;ti<nop;ti++)
  {
    if (this->molecules[ti].belongsto != nilvalue) clusterhisto[this->molecules[ti].belongsto] += 1;
  }
  //From the clustersize histo we can find the largest one.
  for (int ti= 0;ti<nop;ti++)
  {
    if (clusterhisto[ti] > greatestcluster)
    {
      greatestcluster = clusterhisto[ti];
      greatestbelongsto = ti;
    }
  }
  //return the size of the greatest cluster. To investigate the cluster use the variable
  //greatestbelongs to.
  return greatestcluster;
}



//-------------------------------------------------------------------------------------------------------
//DISTANCE CACULATIONS
//-------------------------------------------------------------------------------------------------------

void CMolecularSystem::makeperiodic(int ti)
{ 
  double boxx,boxy,boxz;
  boxx = this->parameter->boxx;
  boxy = this->parameter->boxy;
  boxz = this->parameter->boxz;
  
  if (this->molecules[ti].posx >=  boxx) {this->molecules[ti].posx -=boxx;}
  if (this->molecules[ti].posx <   0.0)  {this->molecules[ti].posx +=boxx;}
  if (this->molecules[ti].posy >=  boxy) {this->molecules[ti].posy -=boxy;}
  if (this->molecules[ti].posy <   0.0)  {this->molecules[ti].posy +=boxy;}
  if (this->molecules[ti].posz >=  boxz) {this->molecules[ti].posz -=boxz;}
  if (this->molecules[ti].posz <   0.0)  {this->molecules[ti].posz +=boxz;}

}


//Distance with nearest image convention
void CMolecularSystem::get_distance(int ti ,int tj ,double &diffx ,double &diffy,double &diffz)
{
  diffx = this->molecules[tj].posx - this->molecules[ti].posx;
  diffy = this->molecules[tj].posy - this->molecules[ti].posy;
  diffz = this->molecules[tj].posz - this->molecules[ti].posz;

  //nearest image
  if (diffx >  this->parameter->boxx/2.0) {diffx = diffx - this->parameter->boxx;};
  if (diffx < -this->parameter->boxx/2.0) {diffx = diffx + this->parameter->boxx;};
  if (diffy >  this->parameter->boxy/2.0) {diffy = diffy - this->parameter->boxy;};
  if (diffy < -this->parameter->boxy/2.0) {diffy = diffy + this->parameter->boxy;};
  if (diffz >  this->parameter->boxz/2.0) {diffz = diffz - this->parameter->boxz;};
  if (diffz < -this->parameter->boxz/2.0) {diffz = diffz + this->parameter->boxz;};
}


void CMolecularSystem::get_distancePosition(int ti ,double posx, double posy,double posz ,double &diffx,double &diffy,double &diffz)
{
  diffx = posx - this->molecules[ti].posx;
  diffy = posy - this->molecules[ti].posy;
  diffz = posz - this->molecules[ti].posz;

  //nearest image
  if (diffx >  this->parameter->boxx/2.0) {diffx = diffx - this->parameter->boxx;};
  if (diffx < -this->parameter->boxx/2.0) {diffx = diffx + this->parameter->boxx;};
  if (diffy >  this->parameter->boxy/2.0) {diffy = diffy - this->parameter->boxy;};
  if (diffy < -this->parameter->boxy/2.0) {diffy = diffy + this->parameter->boxy;};
  if (diffz >  this->parameter->boxz/2.0) {diffz = diffz - this->parameter->boxz;};
  if (diffz < -this->parameter->boxz/2.0) {diffz = diffz + this->parameter->boxz;};
}
    

double CMolecularSystem::get_absDistance(int ti ,int tj)
{
  double abs,diffx,diffy,diffz;
  diffx = this->molecules[tj].posx - this->molecules[ti].posx;
  diffy = this->molecules[tj].posy - this->molecules[ti].posy;
  diffz = this->molecules[tj].posz - this->molecules[ti].posz;
  //nearest image
  if (diffx >  this->parameter->boxx/2.0) {diffx = diffx - this->parameter->boxx;};
  if (diffx < -this->parameter->boxx/2.0) {diffx = diffx + this->parameter->boxx;};
  if (diffy >  this->parameter->boxy/2.0) {diffy = diffy - this->parameter->boxy;};
  if (diffy < -this->parameter->boxy/2.0) {diffy = diffy + this->parameter->boxy;};
  if (diffz >  this->parameter->boxz/2.0) {diffz = diffz - this->parameter->boxz;};
  if (diffz < -this->parameter->boxz/2.0) {diffz = diffz + this->parameter->boxz;};
  abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
  return abs;
}


//-------------------------------------------------------------------------------------------------------
//OUTPUTSECTION
//-------------------------------------------------------------------------------------------------------
//Outputs the whole system as a xyz file for VMD
//with a number for the filename

void CMolecularSystem::outputSimpleVMDXYZFile(int number)
{
  char IntStr[80];
  ofstream of;
  sprintf( IntStr, "./output/simple.%d.xyz", number);
  of.open (IntStr, ofstream::out | ofstream::app);
  if (of.is_open())
  {
    of << this->parameter->nop << "\n\n"; 
    for (int ti=0;ti<this->parameter->nop;ti++)
    {
      of << "h " << this->molecules[ti].posx  << " " << this->molecules[ti].posy  << " " <<this->molecules[ti].posz  << "\n";
    }
  }
  of.close();
}


//Output the averaged and local bond order parameters for each particle
//The columns are:
//index,x,y,z,aq4,aq6,q4,q6,numberofNeighbors
void CMolecularSystem::outputSimpleResults(int number)
{
  char IntStr[80];
  ofstream of;
  sprintf( IntStr, "./output/result.%d.dat", number);
  of.open (IntStr, ofstream::out | ofstream::trunc);
  if (of.is_open())
  {
    for (int ti=0;ti<this->parameter->nop;ti++)
    {
      of << ti << " " << this->molecules[ti].posx  << " " << this->molecules[ti].posy  << " " <<this->molecules[ti].posz
               << " " << this->molecules[ti].AQ4  << " " << this->molecules[ti].AQ6
               << " " << this->molecules[ti].Q4  << " " << this->molecules[ti].Q6
               << " " << this->molecules[ti].n_neighbors << "\n";
    }
  }
}


//Output the histograms of the bond order parameters. Make sure that
//the directory output/orderparameter exists
void CMolecularSystem::outputBondOrderHistogramms(int number)
{
  char IntStr[80];

  double a,b;
  int n;
  a = OrderParameterHistoStart;
  b = OrderParameterHistoStop;
  n = OrderParameterHistoSize;

  ofstream of;
  sprintf( IntStr, "./output/orderparameter/q4Histo.%d.txt", number) ;
  of.open (IntStr, ofstream::out | ofstream::trunc);

  if (of.is_open())
  {
    for (int i = 0;i<n;i++)
    {
      of << get_xFromBox(i,a,b,n) << " "  << this->q4Histo[i]/(double(this->normq4)*(b-a)) << "\n";
    }
  } else {cerr << "unable to open file for output, please mkdir output/orderparameter\n";}
  of.close();

  sprintf( IntStr, "./output/orderparameter/q6Histo.%d.txt", number) ;
  of.open (IntStr, ofstream::out | ofstream::trunc);

  if (of.is_open())
  {
    for (int i = 0;i<n;i++)
    {
      of << get_xFromBox(i,a,b,n) << " "  << this->q6Histo[i]/(double(this->normq6)*(b-a)) << "\n";
    }
  } else {cerr << "unable to open file for output, please mkdir output/orderparameter\n";}
  of.close();


  sprintf( IntStr, "./output/orderparameter/aq4Histo.%d.txt", number) ;
  of.open (IntStr, ofstream::out | ofstream::trunc);

  if (of.is_open())
  {
    for (int i = 0;i<n;i++)
    {
      of << get_xFromBox(i,a,b,n) << " "  << this->aq4Histo[i]/(double(this->normaq4)*(b-a)) << "\n";
    }
  } else {cerr << "unable to open file for output, please mkdir output/orderparameter\n";}
  of.close();

  sprintf( IntStr, "./output/orderparameter/aq6Histo.%d.txt", number) ;
  of.open (IntStr, ofstream::out | ofstream::trunc);

  if (of.is_open())
  {
    for (int i = 0;i<n;i++)
    {
      of << get_xFromBox(i,a,b,n) << " "  << this->aq6Histo[i]/(double(this->normaq6)*(b-a)) << "\n";
    }
  } else {cerr << "unable to open file for output, please mkdir output/orderparameter\n";}
  of.close();

  sprintf( IntStr, "./output/orderparameter/q46Histo3d.%d.txt", number) ;
  of.open (IntStr, ofstream::out | ofstream::trunc);

  if (of.is_open())
  {
    for (int i = 0;i<n;i++)
    {
      for (int j = 0;j<n;j++)
      {
        of << get_xFromBox(i,a,b,n) << " "  << get_xFromBox(j,a,b,n) << " "  <<  this->q46Histo3d[i][j]/(double(this->normaq6)*pow((b-a),2.0)) << "\n";
      }
      of << "\n";
    }

  } else {cerr << "unable to open file for output, please mkdir output/orderparameter\n";}
  of.close();


  sprintf( IntStr, "./output/orderparameter/aq46Histo3d.%d.txt", number) ;
  of.open (IntStr, ofstream::out | ofstream::trunc);

  if (of.is_open())
  {
    for (int i = 0;i<n;i++)
    {
      for (int j = 0;j<n;j++)
      {
        of << get_xFromBox(i,a,b,n) << " "  << get_xFromBox(j,a,b,n) << " "  <<  this->aq46Histo3d[i][j]/(double(this->normaq6)*pow((b-a),2.0)) << "\n";
      }
      of << "\n";
    }

  } else {cerr << "unable to open file for output, please mkdir output/orderparameter\n";}
  of.close();

}
















