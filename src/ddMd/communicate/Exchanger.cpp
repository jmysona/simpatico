#ifndef EXCHANGER_CPP
#define EXCHANGER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Exchanger.h"
#include "Domain.h"
#include "Buffer.h"
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/storage/BondStorage.h>
#include <ddMd/storage/GroupIterator.h>
#include <util/mpi/MpiLogger.h>
#include <util/global.h>

#include <algorithm>
#include <string>

namespace DdMd
{
   using namespace Util;

   /*
   * Constructor.
   */
   Exchanger::Exchanger()
    : pairCutoff_(-1.0),
      ghostCapacity_(0)
   {}

   /*
   * Destructor.
   */
   Exchanger::~Exchanger()
   {}

   /*
   * Set pointers to associated objects.
   */
   void Exchanger::associate(const Domain& domain, 
                             const Boundary& boundary, 
                             AtomStorage& atomStorage, 
                             BondStorage& bondStorage, 
                             Buffer& buffer)
   {
      domainPtr_  = &domain;
      boundaryPtr_  = &boundary;
      atomStoragePtr_  = &atomStorage;
      bondStoragePtr_  = &bondStorage;
      bufferPtr_  = &buffer;
   }

   /*
   * Allocate memory.
   */
   void Exchanger::allocate()
   {
      // Preconditions
      if (!bufferPtr_->isInitialized()) {
         UTIL_THROW("Buffer must be allocated before Exchanger");
      }

      // Allocate all send and recv arrays
      ghostCapacity_ = std::max(bufferPtr_->atomCapacity(), bufferPtr_->ghostCapacity());
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < 2; ++j) {
            sendArray_(i, j).allocate(ghostCapacity_);
            recvArray_(i, j).allocate(ghostCapacity_);
         }
      }

   }

   /*
   * Set slab width used for ghosts.
   */
   void Exchanger::setPairCutoff(double pairCutoff)
   {  pairCutoff_ = pairCutoff; }

   #ifdef UTIL_MPI

   /**
   * Exchange ownership of local atoms.
   *
   * This method should be called before rebuilding the neighbor list on
   * each processor, to exchange ownership of local atoms.
   */
   void Exchanger::exchangeAtoms()
   {
      Vector lengths = boundaryPtr_->lengths();
      double bound, inner, coordinate;
      AtomIterator atomIter;
      GhostIterator ghostIter;
      GroupIterator<2> bondIter;
      Atom* atomPtr;
      int i, j, jc, k, source, dest, nSend;
      int myRank = domainPtr_->gridRank();
      int shift;
      int nLocal;
      bool choose;

      #if 0
      // Calculate atom communication plans for all local atoms
      atomStoragePtr_->begin(atomIter);
      for ( ; !atomIter.atEnd(); ++atomIter) {

         atomIter->plan().clearFlags();

         // Cartesian directions
         for (i = 0; i < Dimension; ++i) {
  
            coordinate = atomIter->position()[i];
 
            // Transmission direction
            for (j = 0; j < 2; ++j) {
   
               // j = 0 sends to lower coordinate i
               // j = 1 sends to higher coordinate i
   
               // Boundary (j=0 -> minimum, j=1 -> maximum)
               bound = domainPtr_->domainBound(i, j);

               // Decide upon plan  direction i, j
               if (j == 0) { // Communicate with lower index
                  jc = 1; // direction for reverse communication
                  choose = (coordinate < bound);
                  if (!choose) { // retain atom
                     inner = bound + pairCutoff_;
                     if (coordinate < inner) {
                        atomIter->plan().setGhost(i, j);
                     }
                  } else { // Send to lower index
                     atomIter->plan().setExchange(i, j);
                     inner = bound - pairCutoff_;
                     if (coordinate > inner) {
                        atomIter->plan().setGhost(i, jc);
                     }
                  }
               } else { // j == 1, communicate with upper index
                  jc = 0; // direction for reverse communication
                  choose = (coordinate > bound);
                  if (!choose) { // retain atom
                     inner = bound - pairCutoff_;
                     if (coordinate > inner) {
                        atomIter->plan().setGhost(i, j);
                     }
                  } else { // send to upper index
                     atomIter->plan().setExchange(i, j);
                     inner = bound + pairCutoff_;
                     if (coordinate < inner) {
                        atomIter->plan().setGhost(i, jc);
                     }
                  }
               }
   
            } // end loop j
         } // loop i

      } // end atom loop, end compute plan
      #endif

      // Clear pointers to ghosts and check local pointers in all bonds.
      int atomId;
      bondStoragePtr_->begin(bondIter);
      for ( ; !bondIter.atEnd(); ++bondIter) {
         for (k = 0; k < 2; ++k) {
            atomPtr = bondIter->atomPtr(k);
            if (atomPtr != 0) {
               if (atomPtr->isGhost()) {
                  bondIter->clearAtomPtr(k);
               } else {
                  atomId = atomPtr->id();
                  if (atomPtr != atomStoragePtr_->find(atomId)) {
                     UTIL_THROW("Error in atom pointer in bond");
                  }
               }
            }
         }
         bondIter->setPostMark(false);
      }

      atomStoragePtr_->clearGhosts();

      // Cartesian directions
      for (i = 0; i < Dimension; ++i) {

         // Transmission direction
         for (j = 0; j < 2; ++j) {

            // j = 0 sends to lower  coordinate i
            // j = 1 sends to higher coordinate i

            //Processor to receive from
            source = domainPtr_->sourceRank(i, j);

            //Processor to send to
            dest = domainPtr_->destRank(i, j);

            // Boundary (j=0 -> minimum, j=1 -> maximum)
            bound = domainPtr_->domainBound(i, j);

            // Shift due to periodic boundary conditions
            shift = domainPtr_->shift(i, j);

            sendArray_(i, j).clear();

            if (domainPtr_->grid().dimension(i) > 1) {
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::ATOM);
            }

            // Choose atoms for sending, pack and mark for removal.
            atomStoragePtr_->begin(atomIter);
            for ( ; !atomIter.atEnd(); ++atomIter) {

               atomIter->setPostMark(false);

               if (j == 0) {
                  choose = (atomIter->position()[i] < bound);
               } else {
                  choose = (atomIter->position()[i] > bound);
               }
               //assert(choose == atomIter->plan().testExchange(i, j));

               if (choose) {

                  if (domainPtr_->grid().dimension(i) > 1) {

                     sendArray_(i, j).append(*atomIter);
                     bufferPtr_->packAtom(*atomIter);
                     atomIter->setPostMark(true);

                  } else {

                     // Shift position if required by periodic b.c.
                     if (shift) {
                        atomIter->position()[i] += shift * lengths[i];
                     }
                     assert(atomIter->position()[i] > domainPtr_->domainBound(i, 0));
                     assert(atomIter->position()[i] < domainPtr_->domainBound(i, 1));

                  }
               }

            } // end atom loop

            // Send and receive only if dimension(i) > 1
            if (domainPtr_->grid().dimension(i) > 1) {

               // End atom send block
               bufferPtr_->endSendBlock();

               // Pack all bonds that contain one or more postmarked atoms.
               bufferPtr_->beginSendBlock(Buffer::GROUP, 2);
               bondStoragePtr_->begin(bondIter);
               for ( ; !bondIter.atEnd(); ++bondIter) {
                  bondIter->setPostMark(false);
                  for (k = 0; k < 2; ++k) {
                     atomPtr = bondIter->atomPtr(k);
                     if (atomPtr != 0) {
                        if (atomPtr->postMark()) {
                           bondIter->setPostMark(true);
                           break;
                        }
                     }
                  }
                  if (bondIter->postMark()) {
                     bufferPtr_->packGroup<2>(*bondIter);
                     bondIter->setPostMark(true);
                  }
               }
               bufferPtr_->endSendBlock();

               // Removal cannot be done within the above loops over atoms and 
               // groups because removal invalidates the iterator.

               // Remove chosen atoms from atomStorage
               nSend = sendArray_(i, j).size();
               for (k = 0; k < nSend; ++k) {
                  atomStoragePtr_->removeAtom(&sendArray_(i, j)[k]);
               }

               // Send to processor dest and receive from processor source
               bufferPtr_->sendRecv(domainPtr_->communicator(), source, dest);

               // Unpack atoms into atomStorage
               bufferPtr_->beginRecvBlock();
               while (bufferPtr_->recvSize() > 0) {
                  atomPtr = atomStoragePtr_->newAtomPtr();
                  bufferPtr_->unpackAtom(*atomPtr);
                  if (shift) {
                     atomPtr->position()[i] += shift * lengths[i];
                  }
                  assert(atomPtr->position()[i] > domainPtr_->domainBound(i, 0));
                  assert(atomPtr->position()[i] < domainPtr_->domainBound(i, 1));
                  atomStoragePtr_->addNewAtom();
               }
               assert(bufferPtr_->recvSize() == 0);

               // Unpack bonds into bondStorage.
               Group<2>* bondPtr;
               bufferPtr_->beginRecvBlock();
               while (bufferPtr_->recvSize() > 0) {
                  bondPtr = bondStoragePtr_->newPtr();
                  bufferPtr_->unpackGroup<2>(*bondPtr);
                  if (bondStoragePtr_->find(bondPtr->id())) {
                     bondStoragePtr_->returnPtr();
                  } else {
                     bondStoragePtr_->add();
                  }
               }
               assert(bufferPtr_->recvSize() == 0);

            }
         }
      }

      // Find all atoms for all bonds
      bondStoragePtr_->begin(bondIter); 
      for ( ; bondIter.notEnd(); ++bondIter) {
         for (k = 0; k < 2; ++k) {
            atomId = bondIter->atomId(k);
            atomPtr = atomStoragePtr_->find(atomId);
            if (atomPtr) {
               bondIter->setAtomPtr(k, atomPtr);
            } else {
               bondIter->clearAtomPtr(k);
            }
         }
      }
 
      #if 0
      // At this point:
      // All atoms on correct processor
      // No ghost atoms exist
      // All pointer to local atoms in Groups are set.
      // All pointers to ghost atoms in Groups are null.
      */
      #endif

   }

   /*
   * Exchange ghost atoms.
   *
   * Call after exchangeAtoms and before rebuilding the neighbor list on
   * time steps that require reneighboring.
   */
   void Exchanger::exchangeGhosts()
   {

      // Preconditions
      assert(bufferPtr_->isInitialized());
      assert(domainPtr_->isInitialized());
      assert(domainPtr_->hasBoundary());

      Vector        lengths = boundaryPtr_->lengths();
      double        bound, inner, coord;
      AtomIterator  localIter;
      GhostIterator ghostIter;
      Atom* atomPtr;
      int i, j, source, dest, shift;
      int myRank = domainPtr_->gridRank();
      bool choose;

      // Check that all ghosts are cleared upon entry.
      if (atomStoragePtr_->nGhost() > 0) {
         UTIL_THROW("atomStoragePtr_->nGhost() > 0");
      }

      // Cartesian directions
      for (i = 0; i < Dimension; ++i) {

         // Transmit directions
         for (j = 0; j < 2; ++j) {

            // j = 0: Send ghosts near minimum boundary to lower coordinate
            // j = 1: Send ghosts near maximum boundary to higher coordinate

            if (j == 0) {
               bound = domainPtr_->domainBound(i, 0);
               inner = bound + pairCutoff_;
            } else {
               bound = domainPtr_->domainBound(i, 1);
               inner = bound - pairCutoff_;
            }

            // Shift on receiving processor for periodic boundary conditions
            shift = domainPtr_->shift(i, j);

            if (domainPtr_->grid().dimension(i) > 1)  {
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::GHOST);
            }

            sendArray_(i, j).clear();
            recvArray_(i, j).clear();

            // Loop over local atoms on this processor
            atomStoragePtr_->begin(localIter);
            for ( ; !localIter.atEnd(); ++localIter) {

               // Decide whether to send as ghost.
               coord = localIter->position()[i];
               if (j == 0) {
                  assert(coord > bound);
                  choose = (coord < inner);
               } else {
                  assert(coord < bound);
                  choose = (coord > inner);
               }

               #if 0
               assert(choose == localIter->plan().testGhost(i, j));
               if (choose != localIter->plan().testGhost(i, j)) {
                  std::cout << "Proc " << myRank << "  "
                            << "atom " << localIter->id() << "  "
                            << "Dir  " << i << "  " << j << "  "
                            << localIter->position() << "  ";
                  for (int a= 0; a < Dimension; ++a) {
                     for (int b = 0; b < 2;  ++b) {
                        std::cout << localIter->plan().testExchange(a, b);
                     }
                  }
                  std::cout << "  ";
                  for (int a= 0; a < Dimension; ++a) {
                     for (int b = 0; b < 2;  ++b) {
                        std::cout << localIter->plan().testGhost(a, b);
                     }
                  }
                  std::cout << "  " << choose << "  " << localIter->plan().testGhost(i, j);
                  std::cout << std::endl;
                  
                  UTIL_THROW("Assert failed in exchange ghosts");
               }
               #endif

               if (choose) {

                  sendArray_(i, j).append(*localIter);

                  if (domainPtr_->grid().dimension(i) > 1)  {

                     // Pack atom for sending 
                     bufferPtr_->packGhost(*localIter);

                  } else {  // if grid dimension == 1

                     // Make a ghost copy of the local atom on this processor
                     atomPtr = atomStoragePtr_->newGhostPtr();
                     recvArray_(i, j).append(*atomPtr);
                     atomPtr->position() = localIter->position();
                     atomPtr->setTypeId(localIter->typeId());
                     atomPtr->setId(localIter->id());
                     if (shift) {
                        atomPtr->position()[i] += shift * lengths[i];
                     }
                     atomStoragePtr_->addNewGhost();

                     #ifdef UTIL_DEBUG
                     // Validate shifted positions
                     if (j == 0) {
                        assert(atomPtr->position()[i] > domainPtr_->domainBound(i, 1));
                     } else {
                        assert(atomPtr->position()[i] < domainPtr_->domainBound(i, 0));
                     }
                     #endif

                  }

               } 

            }

            // Loop over ghosts on this processor, for resending.
            atomStoragePtr_->begin(ghostIter);
            for ( ; !ghostIter.atEnd(); ++ghostIter) {

               // Decide whether to resend this ghost
               coord = ghostIter->position()[i];
               if (j == 0) {
                  choose = (coord > bound) && (coord < inner);
               } else {
                  choose = (coord < bound) && (coord > inner);
               }

               if (choose) {

                  sendArray_(i, j).append(*ghostIter);

                  if (domainPtr_->grid().dimension(i) > 1)  {
                     
                     // Pack ghost for resending
                     bufferPtr_->packGhost(*ghostIter);

                  } else {  // if grid dimension == 1

                     // Make another ghost copy on the same processor
                     atomPtr = atomStoragePtr_->newGhostPtr();
                     recvArray_(i, j).append(*atomPtr);
                     atomPtr->position() = ghostIter->position();
                     atomPtr->setTypeId(ghostIter->typeId());
                     atomPtr->setId(ghostIter->id());
                     if (shift) {
                        atomPtr->position()[i] += shift * lengths[i];
                     }
                     atomStoragePtr_->addNewGhost();

                     #ifdef UTIL_DEBUG
                     // Validate shifted position
                     if (j == 0) {
                        assert(atomPtr->position()[i] > domainPtr_->domainBound(i, 1));
                     } else {
                        assert(atomPtr->position()[i] < domainPtr_->domainBound(i, 0));
                     }
                     #endif

                  } 

               }

            }

            // Send and receive buffers
            if (domainPtr_->grid().dimension(i) > 1) {

               bufferPtr_->endSendBlock();

               source = domainPtr_->sourceRank(i, j);
               dest   = domainPtr_->destRank(i, j);
               bufferPtr_->sendRecv(domainPtr_->communicator(), source, dest);

               // Unpack ghosts and add to recvArray
               bufferPtr_->beginRecvBlock();
               while (bufferPtr_->recvSize() > 0) {

                  atomPtr = atomStoragePtr_->newGhostPtr();
                  bufferPtr_->unpackGhost(*atomPtr);
                  if (shift) {
                     atomPtr->position()[i] += shift * lengths[i];
                  }
                  recvArray_(i, j).append(*atomPtr);
                  atomStoragePtr_->addNewGhost();

                  #ifdef UTIL_DEBUG
                  // Validate ghost coordinates on receiving processor.
                  if (j == 0) {
                     assert(atomPtr->position()[i] > domainPtr_->domainBound(i, 1));
                  } else {
                     assert(atomPtr->position()[i] < domainPtr_->domainBound(i, 0));
                  }
                  #endif

               }

            }

         } // transmit direction j = 0, 1

      } // Cartesian index i = 0, ..., Dimension - 1

   }

   /*
   * Update ghost atom coordinates.
   *
   * Call on time steps for which no reneighboring is required. 
   */
   void Exchanger::updateGhosts()
   {
      Vector lengths = boundaryPtr_->lengths();
      Atom*  atomPtr;
      int    i, j, k, source, dest, size, shift;

      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < 2; ++j) {

            // Shift on receiving processor for periodic boundary conditions
            shift = domainPtr_->shift(i, j);

            if (domainPtr_->grid().dimension(i) > 1) {

               // Pack ghost positions for sending
               bufferPtr_->clearSendBuffer();
               bufferPtr_->beginSendBlock(Buffer::GHOST);
               size = sendArray_(i, j).size();
               for (k = 0; k < size; ++k) {
                  bufferPtr_->packGhost(sendArray_(i, j)[k]);
               }
               bufferPtr_->endSendBlock();
  
               // Send and receive buffers
               source = domainPtr_->sourceRank(i, j);
               dest   = domainPtr_->destRank(i, j);
               bufferPtr_->sendRecv(domainPtr_->communicator(), source, dest);
   
               // Unpack ghost positions
               bufferPtr_->beginRecvBlock();
               size = recvArray_(i, j).size();
               for (k = 0; k < size; ++k) {
                  atomPtr = &recvArray_(i, j)[k];
                  bufferPtr_->unpackGhost(*atomPtr);
                  if (shift) {
                     atomPtr->position()[i] += shift * lengths[i];
                  }
               }

            } else {

               // If grid().dimension(i) == 1, then copy positions of atoms
               // listed in sendArray to those listed in the recvArray.

               size = sendArray_(i, j).size();
               assert(size == recvArray_(i, j).size());
               for (k = 0; k < size; ++k) {
                  atomPtr = &recvArray_(i, j)[k];
                  atomPtr->position() = sendArray_(i, j)[k].position();
                  if (shift) {
                     atomPtr->position()[i] += shift * lengths[i];
                  }
               }

            }  

         } // transmit direction j = 0, 1

      } // Cartesian direction i = 0, ..., Dimension - 1

   }

   #endif

}
#endif
