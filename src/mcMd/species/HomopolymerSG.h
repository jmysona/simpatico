#ifndef MCMD_HOMOPOLYMER_SG_H
#define MCMD_HOMOPOLYMER_SG_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SpeciesMutator.h"         // base class
#include <simp/species/Linear.h>    // base class
#include <util/containers/Pair.h>   // member

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   * A Homopolymer with a mutable type, for semigrand ensemble.
   *
   * A HomopolymerSG molecule is a Linear chain in which all atoms
   * are of the same type, as are all bonds, but in which the type
   * index for all atoms can be toggled between two values.
   *
   * \ingroup Simp_Species_Module
   */
   class HomopolymerSG : public Simp::Linear, public McMd::SpeciesMutator
   {

   public:

      /**
      * Constructor.
      */
      HomopolymerSG();

      /**
      * Destructor.
      */
      virtual ~HomopolymerSG();

      /**
      * Read parameters and initialize.
      * Unneeded
      */
      //virtual void readParameters(std::istream& in);

      /**
      * Set the type of all atoms in the molecule.
      *
      * \param molecule reference to molecule of interest.
      * \param stateId  atom type id.
      */
      virtual void setMoleculeState(Molecule& molecule, int stateId);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
   
   protected:

      /**
      * Read nAtom, a pair of atom type ids and weightRatio.
      *
      * \param in input stream
      */
      virtual void readSpeciesParam(std::istream &in);

      /**
      * Load species structure from an Archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadSpeciesParam(Serializable::IArchive &ar);

<<<<<<< HEAD

=======
>>>>>>> 0a31561d50ea070ac74ca6f560ad950c29a7029e
      /**
      * Return the same type for any particle in any chain.
      *
      * \param index atom index, in range 0,...,nAtom_ - 1
      * \return atom type index
      */
      virtual int calculateAtomTypeId(int index) const;

      /**
      * Return same bond type for any bond in any chain.
      *
      * \param index bond index, in range 0,...,nAtom_ - 1
      * \return bond type index
      */
      virtual int calculateBondTypeId(int index) const;

      #ifdef SIMP_ANGLE
      /**
      * Return same angle type for any angle in any chain.
      *
      * \param index  angle index, in range 0, ..., nAngle_ - 1
      * \return       angle type index
      */
      virtual int calculateAngleTypeId(int index) const;
      #endif

      #ifdef SIMP_DIHEDRAL
      /**
      * Return same dihedral type for any dihedral in any chain.
      *
      * \param index  dihedral index, in range 0, ..., nDihedral_ - 1
      * \return       dihedral type index
      */
      virtual int calculateDihedralTypeId(int index) const;
      #endif

   private:

      // Type index for all bonds
      int bondType_;

      #ifdef SIMP_ANGLE
      // Type index for all angles (if any)
      int angleType_;
      #endif

      #ifdef SIMP_DIHEDRAL
      // Type index for all dihedrals (if any)
      int dihedralType_;
      #endif

      // A Pair of atom type indexes.
      // Atoms in a molecule with stateId = i have type typeIds_[i]
      Pair<int> typeIds_;

      /// Ratio of statistical weights for typeId[0]/typeId[1]
      double weightRatio_;
   
   };

}
#endif
