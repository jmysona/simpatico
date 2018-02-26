/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "InterfacialLoading.h"
#include <simp/boundary/Boundary.h>

#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <simp/species/Species.h>

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
      int atomCapacity = system().simulation().atomCapacity();
      cellList_.setAtomCapacity(atomCapacity);
   }

   /// Add particle pairs to RDF histogram.
   void InterfacialLoading::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
        //Build the cell list of the second species
        cellList_.setup(system().boundary(), cutoff_);
        System::MoleculeIterator molIter;
        Molecule::AtomIterator atomIter;
        CellList::NeighborArray neighborArray;
        Atom* otherAtom;
        double cutoffSq = cutoff_*cutoff_;
        double rsq;
        int moleculesOnInterface = 0;
        bool isMoleculeOnInterface;


        system().begin(species2Id_,molIter);
        for( ; molIter.notEnd(); ++molIter) {
          for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
            if (atomIter->typeId() == atomType2Id_) {
              system().boundary().shift(atomIter->position());
              cellList_.addAtom(*atomIter);
            }
          }
        }
        /// Cell list built, now cycle through the other species to see which are in contact
        system().begin(species1Id_,molIter);
        for( ; molIter.notEnd(); ++molIter) {
          isMoleculeOnInterface = false;
          for (molIter->begin(atomIter); atomIter.notEnd(); ++atomIter) {
            if (atomIter->typeId() == atomType1Id_) {
              cellList_.getNeighbors(atomIter->position(), neighborArray);
              for (int i = 0; i< neighborArray.size(); i++) {
                otherAtom = neighborArray[i];
                rsq = system().boundary().distanceSq(atomIter->position(),otherAtom->position());
                if (rsq < cutoffSq) {
                  isMoleculeOnInterface = true;
                  moleculesOnInterface += 1;
                  break;
                }  
              }
              if (isMoleculeOnInterface) {
                break;
              }
            }      
          }
        }
        accumulator_.sample(moleculesOnInterface);
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
