///////////////////////////////////////////////////////////////////////////////////////////
//
//
// Structure Analysis Code v 1.0
// by Wolfgang Lechner, Amsterdam 2010
//
// wolfgang.lechner@gmail.com
//
//
//
// This source code is licensed under the Academic Free License
// Copyright by Wolfgang Lechner, 2010
// In addition to that I would like you to cite
// Wolfgang Lechner and Christoph Dellago, J. Chem. Phys. 129, 114707 (2008)
// and send me a short email about your project when using this code.
//
// The package shows how to use the averaged version of the local bond order parameters
// in order to determine the local structure of points or particles distributed in
// three dimensional space described in the Paper mentioned above.
// This is useful for a wide range of applications including particle simulations to
// investigate nucleation, determine the crystal structure in metals, and the study of
// glasses.
//
// This package was originally part of a particle simulation to study the nucleation
// of colloidal particles. If you are interested in a collaboration on this field feel
// free to email me about your interests.
// Also, if you have any questions or found a bug do not hesitate writing me!
//
///////////////////////////////////////////////////////////////////////////////////////////

#include "molecular_system.h"
#include "molecule.h"
#include "parameter.h"
#include <ctime>
#include <sstream>
#include <string>
#include <cstdlib>

int main(int argc, char *argv[])
{
  // The program accepts one integer as a paramter.
  // eg. executing "./main 1" means that the parameterfile
  // parameter.1.txt in the input directory is used
  // as parameter file.
  // I have included 4 different parameter files for demonstration
  // parameter.1.txt reads a configuration file. With the
  // parameters 2,3,4 perfect crystals of FCC,BCC and HCP are
  // created, respectively. If no parameter is given, parameter.txt
  // is read.
  int t_parameter = 0;
  if (argc > 1)
  {
    t_parameter = atoi(argv[1]);
  }
  // The main obejct is created. It hold all the functions and data
  // used in the analysis.
  CMolecularSystem *m_MolSys = new CMolecularSystem;
  // The parameterfile is read
  m_MolSys->parameter->readParameter(t_parameter);
  // System is initalized, memory allocated, ...
  m_MolSys->InitializeSystem();
  // Neighborlists are used, this reduces the cost of searching for neighbors
  // to O(N) instead of O(N**2).
  //m_MolSys->fillLists();
  m_MolSys->outputSimpleVMDXYZFile(t_parameter);
  
  // The structure analysis is done here, get neighbors, calculate structure
  m_MolSys->calculate_bondOrderParameter();
  // Output the results for each particle.
  m_MolSys->outputSimpleResults(t_parameter);

  // If you are interested in the histogram of the structure parameters just
  // uncomment the next 4 lines.
  //m_MolSys->update_bondOrderHistogramms();
  //m_MolSys->sample_Q4Q6Histograms();
  //m_MolSys->outputBondOrderHistogramms(t_parameter);
  //m_MolSys->outputSimpleVMDXYZFile(t_parameter);

  //Free the memory.
  m_MolSys->deleteMolecules();

  return 0;
}
