/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McSpeciesCount.h"                  
#include <util/misc/FileMaster.h>
#include <mcMd/simulation/Simulation.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   McSpeciesCount::McSpeciesCount(McSystem& system) :
      SystemAnalyzer<McSystem>(system)
   {  setClassName("McSpeciesCount"); }

   /*
   * Read file name and open output file.
   */
   void McSpeciesCount::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
   }
 
   /*
   * Load state from an archive, and open output file.
   */
   void McSpeciesCount::loadParameters(Serializable::IArchive& ar)
   {  
      Analyzer::loadParameters(ar);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
   }

   /*
   * Save state to an archive.
   */
   void McSpeciesCount::save(Serializable::OArchive& ar)
   {  ar & *this; }

   /*
   * Evaluate energy and output to outputFile_.
   */
   void McSpeciesCount::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {
         int speciesCount;
         int nSpecies;
         nSpecies = system().simulation().nSpecies();
         for (int i =0; i < nSpecies; i++){
         speciesCount = system().nMolecule(i);
         outputFile_ << speciesCount << "     ";
         }
         outputFile_ << "\n";
      }
   }
 
   /* 
   * Summary
   */
   void McSpeciesCount::output() 
   {
      // Close *.dat file
      outputFile_.close();

      // Open and write summary file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_ << std::endl;
      outputFile_ << std::endl;
      outputFile_.close();
   }
   
}
