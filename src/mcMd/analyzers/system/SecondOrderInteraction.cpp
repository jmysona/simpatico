/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SecondOrderInteraction.h"
#include <mcMd/simulation/System.h>
#include <mcMd/species/Species.h>

namespace McMd
{

  using namespace Util;

  ///Constructor
  SecondOrderInteraction::SecondOrderInteraction(System& system)
   : SystemAnalyzer<System>(system),
     isInitialized_(false)
  {  setClassName("SecondOrderInteraction");}

  void SecondOrderInteraction::readParameters(std::istream& in)
  {
    readInterval(in);
    readOutputFileName(in);
    read<int>(in, "speciesId", speciesId_);
    if (speciesId_ < 0) {
      UTIL_THROW("Negative speciesId");
    }
    isInitialized_ =true; 
  }

  /*
  * Load state from an archive.
  */
  void SecondOrderInteraction::loadParameters(Serializable::IArchive& ar)
  {
    loadInterval(ar);
    loadOutputFileName(ar);
    loadParameter<int>(ar, "speciesId", speciesId_);
    // Validate
    if (speciesId_ < 0) {
      UTIL_THROW("Negative speciesId");
    }
      isInitialized_ = true;
  }
 
  /*
  *  Save state to archive
  */
  void SecondOrderInteraction::save(Serializable::OArchive& ar)
  { ar & *this; }
  
  /*
  *  Allocate matrix
  */
  void SecondOrderInteraction::setup()
  {
    if (!isInitialized_) {
      UTIL_THROW("Object is not initialized");
    }
    int speciesCapacity = system().simulation().species(speciesId_).capacity();
    Iab_.allocate(speciesCapacity, speciesCapacity);
  }

  /*
  * Calculate I_ab and other important quantities
  */
  void SecondOrderInteraction::sample(long iStep)
  { 
    DMatrix<double> CurrentIab;
    DArray<double> Ia;
    double I;
    double Isquared;
    int sampleCount = 0;
    int nMolecule = system().nMolecule(speciesId_);
    for (int i = 0; i < nMolecule; ++i) {
      for (int j = 0; j < nAtoms; ++j) {
          CurrentIab[i][j]=CurrentIab[i][j]+energy;
          Iab[i][j]=Iab[i][j]+energy;
      }
    
    
    }


  }

  void SecondOrderInteraction::output()
  {
  
  }

}


