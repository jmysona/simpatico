/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "InterfacialLoading.h"

#include <util/global.h>

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /*
   * Constructor.
   */
   InterfacialLoading::InterfacialLoading(System& system) 
    : SystemAnalyzer<System>(system),
      outputFile_(),
      accumulator_(),
      isInitialized_(false)
   {  setClassName("InterfacialLoading"); }

   /*
   * Read parameters from file, and allocate data array.
   */
   void InterfacialLoading::readParameters(std::istream& in) 
   {
      Analyzer::readParameters(in);
      read<int>(in, "species1Id", species1Id_);
      read<int>(in, "atomType1Id", atomType1Id_);
      read<int>(in, "species2Id", species2Id_);
      read<int>(in, "atomType2Id", atomType2Id_);
      read<double>(in, "cutoff", cutoff_);
      read<int>(in, "histMax", histMax_);
      isInitialized_ = true;
   }

   /*
   * Load state from an archive.
   */
   void InterfacialLoading::loadParameters(Serializable::IArchive& ar)
   {
      Analyzer::loadParameters(ar);

      ar & accumulator_;

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void InterfacialLoading::save(Serializable::OArchive& ar)
   { ar & *this; }

   /*
   * Add particle pairs to RDF histogram.
   */
   void InterfacialLoading::setup() 
   { 
      if (!isInitialized_) {
         UTIL_THROW("Object not initialized");
      }
      accumulator_.setParam(0, histMax_); 
      accumulator_.clear(); 
   }

   /// Add particle pairs to RDF histogram.
   void InterfacialLoading::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {


      } 

   }


   /// Output results to file after simulation is completed.
   void InterfacialLoading::output() 
   {  

      // Echo parameters to a log file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_); 
      outputFile_.close();

      // Output statistical analysis to separate data file
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();

   }

}
