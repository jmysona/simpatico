#ifndef SPAN_HOOMD_CONFIG_WRITER_CPP
#define SPAN_HOOMD_CONFIG_WRITER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HoomdConfigWriter.h"

#include <spAn/chemistry/Atom.h>
#include <spAn/chemistry/Group.h>
#include <spAn/chemistry/Species.h>
//#include <spAn/chemistry/MaskPolicy.h>
#include <spAn/storage/Configuration.h>

#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace SpAn
{

   using namespace Util;

   /*
   * Constructor.
   */
   HoomdConfigWriter::HoomdConfigWriter(Configuration& configuration)
    : ConfigReader(configuration)
   {  setClassName("HoomdConfigWriter"); }

   /* 
   * Write the configuration file.
   */
   void HoomdConfigWriter::writeConfig(std::ofstream& file)
   {
      // Precondition
      if (!file.is_open()) {  
         UTIL_THROW("Error: File is not open"); 
      }

      // Write Boundary dimensions
      file << "BOUNDARY" << std::endl << std::endl;
      file << configuration().boundary() << std::endl;
      file << std::endl;

      // Write atoms
      file << "ATOMS" << std::endl;
      int nAtom = configuration().atoms().size();
      file << "nAtom" << Int(nAtom, 10) << std::endl;
      Vector r;
      AtomStorage::Iterator iter;
      configuration().atoms().begin(iter);
      for (; iter.notEnd(); ++iter) {
         file << Int(iter->id, 10) 
              << Int(iter->typeId, 6);
         r = iter->position;
         if (hasMolecules_) {
            file << Int(iter->speciesId, 6) 
                 << Int(iter->moleculeId, 10)
                 << Int(iter->atomId, 6);
         }
         file << "\n" << r 
              << "\n" << iter->velocity << "\n";
      }

      // Write the groups
      #ifdef INTER_BOND
      if (configuration().bonds().capacity()) {
         writeGroups(file, "BONDS", "nBond", configuration().bonds());
      }
      #endif

      #ifdef INTER_ANGLE
      if (configuration().angles().capacity()) {
         writeGroups(file, "ANGLES", "nAngle", configuration().angles());
      }
      #endif

      #ifdef INTER_DIHEDRAL
      if (configuration().dihedrals().capacity()) {
         writeGroups(file, "DIHEDRALS", "nDihedral", configuration().dihedrals());
      }
      #endif

   }
 
}
#endif
