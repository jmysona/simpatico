#ifndef SPAN_DDMD_CONFIG_READER_CPP
#define SPAN_DDMD_CONFIG_READER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DdMdConfigReader.h"

#include <spAn/chemistry/Atom.h>
#include <spAn/chemistry/Group.h>
#include <spAn/chemistry/Species.h>
//#include <spAn/chemistry/MaskPolicy.h>
#include <spAn/storage/Configuration.h>

#include <util/space/Vector.h>
//#include <util/format/Int.h>
//#include <util/format/Dbl.h>

namespace SpAn
{

   using namespace Util;

   /*
   * Constructor.
   */
   DdMdConfigReader::DdMdConfigReader(Configuration& configuration, bool hasMolecules)
    : ConfigReader(configuration),
      hasMolecules_(hasMolecules)
   {  setClassName("DdMdConfigReader"); }

   /*
   * Read a configuration file.
   */
   void DdMdConfigReader::readConfig(std::ifstream& file)
   {
      // Precondition
      if (!file.is_open()) {  
            UTIL_THROW("Error: File is not open"); 
      }

      // Read and broadcast boundary
      file >> Label("BOUNDARY");
      file >> configuration().boundary();

      // Read and distribute atoms

      // Read atoms
      Atom* atomPtr;
      int atomCapacity = configuration().atoms().capacity(); // Maximum allowed id + 1
      int nAtom;          
      file >> Label("ATOMS");
      file >> Label("nAtom") >> nAtom;
      for (int i = 0; i < nAtom; ++i) {

         // Get pointer to new atom 
         atomPtr = configuration().atoms().newPtr();
 
         file >> atomPtr->id;
         if (atomPtr->id < 0) {
            std::cout << "atom id =" << atomPtr->id << std::endl;
            UTIL_THROW("Negative atom id");
         }
         if (atomPtr->id >= atomCapacity) {
            std::cout << "atom id      =" << atomPtr->id << std::endl;
            std::cout << "atomCapacity =" << atomCapacity << std::endl;
            UTIL_THROW("Invalid atom id");
         }
         file >> atomPtr->typeId;
         if (hasMolecules_) {
            file >> atomPtr->speciesId;
            if (atomPtr->speciesId < 0) {
               std::cout << "species Id  =" << atomPtr->speciesId << std::endl;
               UTIL_THROW("Negative species id");
            }
            file >> atomPtr->moleculeId; 
            if (atomPtr->moleculeId < 0) {
               std::cout << "molecule Id =" << atomPtr->moleculeId << std::endl;
               UTIL_THROW("Negative molecule id");
            }
            file >> atomPtr->atomId;
            if (atomPtr->atomId < 0) {
               std::cout << "atom id     =" << atomPtr->atomId << std::endl;
               UTIL_THROW("Negative atom id in molecule");
            }
         }
         file >> atomPtr->position;
         file >> atomPtr->velocity;

         // Finalize addition of new atom
         configuration().atoms().add();
      }

      // Read Covalent Groups
      #ifdef INTER_BOND
      if (configuration().bonds().capacity()) {
         readGroups(file, "BONDS", "nBond", configuration().bonds());
         //if (maskPolicy == MaskBonded) {
         //   setAtomMasks();
         //}
      }
      #endif

      #ifdef INTER_ANGLE
      if (configuration().angles().capacity()) {
         readGroups(file, "ANGLES", "nAngle", configuration().angles());
      }
      #endif

      #ifdef INTER_DIHEDRAL
      if (configuration().dihedrals().capacity()) {
         readGroups(file, "DIHEDRALS", "nDihedral", configuration().dihedrals());
      }
      #endif

      // Optionally add atoms to species
      if (hasMolecules_) {
         int nSpecies = configuration().nSpecies();
         if (nSpecies > 0) {
            for (int i = 0; i < nSpecies; ++i) {
               configuration().species(i).clear();
            }
            int speciesId;
            AtomStorage::Iterator iter;
            configuration().atoms().begin(iter);
            for ( ; iter.notEnd(); ++iter) {
               speciesId = iter->speciesId;
               if (speciesId < 0) {
                  UTIL_THROW("Negative speciesId");
               }
               if (speciesId >= nSpecies) {
                  UTIL_THROW("speciesId >= nSpecies");
               }
               configuration().species(speciesId).addAtom(*iter);
            }
            for (int i = 0; i < nSpecies; ++i) {
               configuration().species(i).isValid();
            }
         }
      }
   }
 
}
#endif
