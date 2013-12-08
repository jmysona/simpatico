#ifndef DDMD_ATOM_MAP_H
#define DDMD_ATOM_MAP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/chemistry/Atom.h>     // member template argument
#include <util/containers/DArray.h>  // member template
#include <ddMd/chemistry/Group.h>    // used in member function template
#include <util/boundary/Boundary.h>  // used in member function template
#include <util/global.h>

#ifdef UTIL_CXX11
#include <unordered_map>
#else
#include <map>
#endif

namespace Util {
   template <typename T> class ArraySet;
}

namespace DdMd
{

   using namespace Util;

   /**
   * Associative container for finding atoms identified by integer id.
   *
   * \ingroup DdMd_Storage_Module
   */
   class AtomMap 
   {

   public:

      /**
      * Constructor.
      */
      AtomMap();

      /**
      * Destructor.
      */
      ~AtomMap();

      /**
      * Set parameters, allocate memory and initialize.
      *
      * Call this or (read|load)Parameters to initialize, but not both.
      *
      * \param totalAtomCapacity max number of atoms on all processors.
      */
      void allocate(int totalAtomCapacity);

      /**
      * Add local atom.
      * 
      * \param ptr Pointer to new Atom.
      */
      void addLocal(Atom* ptr); 

      /**
      * Remove a specific Atom.
      *
      * \throw Exception if atom is not present.
      *
      * \param ptr Pointer to Atom to be removed.
      */
      void removeLocal(Atom* ptr); 

      /**
      * Add ghost atom.
      * 
      * \param ptr Pointer to new Atom.
      */
      void addGhost(Atom* ptr);

      /**
      * Remove a ghost Atom.
      *
      * This method throws an exception if no atom with this
      * id is present, but not if it does not match this pointer.
      *
      * \param ptr Pointer to Atom to be removed.
      */
      void removeGhost(Atom* ptr); 

      /**
      * Clear all ghosts from this map.  
      *
      * \param Set of all ghosts on this processor.
      */
      void clearGhosts(const ArraySet<Atom>& ghostSet);

      //@}
      /// \name Accessors 
      //@{

      /**
      * Return the number of local atoms.
      */ 
      int nLocal() const;

      /**
      * Return the number of ghosts atoms, including all images.
      */ 
      int nGhost() const;

      /**
      * Return number of ghost atom ids that are distinct from any 
      * local atom id.
      */ 
      int nGhostDistinct() const;

      /**
      * Return a pointer to an Atom with specified id.
      *
      * This method returns a pointer to an Atom with the specified id if
      * one is present, or returns a null pointer otherwise. If a local 
      * atom with the specfied id is present. If no such local atom exists,
      * but one more ghosts atoms exist, it returns a pointer to one of the
      * ghost atoms, chosen arbitrarily. 
      *
      * \param atomId integer index of desired atom
      */
      Atom* find(int atomId) const;  

      /**
      * Find image of an atom nearest a specified position.
      *
      * This function searches for an image of atom number id for which the
      * atom position is the image of itself nearest to the position parameter.
      * On return, imagePtr points to the desired atom image. If the map 
      * contains a local atom with this id whose position is not the nearest
      * image, the function returns a pointer to that local atom, and returns
      * a null pointer otherwise. Atomic coordinates must be in scaled form.
      *
      * Throws an Exception if there is no atom with specified id, or if there 
      * is no such atom with the nearest image position.
      *
      * \param atomId global identifier of required atom
      * \param position position that atom should be near
      * \param boundary Boundary, to implement nearest-image convention
      * \param imagPtr  on return, pointer to nearest image (output)
      * \return pointer to local atom that is not nearest image, or null.
      */
      Atom* findNearestImage(int atomId, const Vector& position, 
                             const Boundary& boundary, Atom*& imagePtr) const;

      /**
      * Set handles to local atoms in a Group<N> object.
      *
      * On entry, group is a Group<N> object for which the atom ids
      * for all N atoms in the Group have been set to valid values, 
      * in the range 0 <= atomId < totalAtomCapacity, but in which
      * some or all pointers have not been set. The AtomMap may not
      * contain any ghosts. 
      *
      * On exit, pointers are set correctly for all local atoms that 
      * exist in this AtomMap, or set to null for absent atoms. All 
      * old pointer values are overwritten. 
      *
      * \param group Group<N> object with known atom ids. 
      * \return number of atoms found on this processor.
      */ 
      template <int N> 
      int findGroupLocalAtoms(Group<N>& group) const;

      /**
      * Check validity of this AtomMap.
      *
      * Returns true if all is ok, or throws an Exception.
      */
      bool isValid() const;

      //@}

   private:

      #ifdef UTIL_CXX11
      typedef std::unordered_multimap<int, Atom*> GhostMap;
      #else
      typedef std::multimap<int, Atom*> GhostMap;
      #endif

      // Array of pointers to atoms, indexed by Id.
      // Elements corresponding to absent atoms hold null pointers.
      DArray<Atom*> atomPtrs_;

      // Map for extra ghost images
      GhostMap ghostMap_;

      /// Number of local atoms in this map.
      int nLocal_;

      /// Number of ghost atom with distinct ids in this map.
      int nGhostDistinct_;

      // Maximum number of atoms on all processors, maximum id + 1
      int totalAtomCapacity_;

      // Has this map been initialized (i.e., allocated)?
      bool isInitialized_;

      /*
      *  Design / invariants:
      *
      *  - If a local atom with id i is present, atomPtrs_[i]
      *    contains a pointer to that atom.
      *
      *  - If one or more ghosts with an atom id i are present,
      *    but there is no local atom with that id, atomPtrs_[i]
      *    contains a pointer to one such ghost.
      *
      *  - ghostMap_ contains pointers to all ghosts except those
      *    in atomPtrs_, stored using atom indices as keys. 
      *
      * One image of each physical atom, identified by id, is thus 
      * stored * in atomPtrs_, while ghostMap_ holds any "extra"
      * ghost images of atoms. If this processor does not contain
      * multiple images of any particle, ghostMap_ will be empty.
      */

   };

   // Inline method definitions

   /*
   * Return pointer to an Atom with specified id.
   */
   inline Atom* AtomMap::find(int atomId) const
   {  return atomPtrs_[atomId]; }

   /*
   * Return the number of local atoms.
   */ 
   inline int AtomMap::nLocal() const
   {  return nLocal_; }

   /*
   * Return the number of ghosts with distinct ids.
   */ 
   inline int AtomMap::nGhostDistinct() const
   { return nGhostDistinct_; }

   /*
   * Return the total number of ghosts, including images.
   */ 
   inline int AtomMap::nGhost() const
   { return nGhostDistinct_ + ghostMap_.size(); }

   // Template method definition

   /*
   * Set pointers to all atoms in a Group<N> object.
   */
   template <int N>
   int AtomMap::findGroupLocalAtoms(Group<N>& group) const
   {
      Atom* ptr;
      int nAtom = 0;
      for (int i = 0; i < N; ++i) {
         ptr = atomPtrs_[group.atomId(i)];
         if (ptr) {
            assert(!ptr->isGhost());
            assert(ptr->id() == group.atomId(i));
            group.setAtomPtr(i, ptr);
            ++nAtom;
         } else {
            group.clearAtomPtr(i);
         }
      }
      return nAtom;
   }

   #if 0
   /*
   * Set pointers to atoms in a Group<N> object.
   */
   template <int N>
   void AtomMap::findGroupGhostAtoms(Group<N>& group, const Boundary& boundary) 
   const
   {
      Atom* ptr;
      int nAtom = 0;
      for (int i = 0; i < N; ++i) {
         ptr = group.atomPtr(i);
         if (ptr) {
            assert(!ptr->isGhost());
            ++nAtom;
         } else {
            int atomId = group.atomId(i);
            ptr = atomPtrs_[atomId];
            if (ptr) {
               assert(ptr->isGhost());
               assert(ptr->id() == atomId);
               group.setAtomPtr(i, ptr);
               ++nAtom;
            }
         }
      }
      if (nAtom != N) {
         UTIL_THROW("Incomplete group at end of findGroupGhostAtoms");
      }
   }
   #endif

}
#endif