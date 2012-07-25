#ifndef DDMD_GROUP_COLLECTOR_H
#define DDMD_GROUP_COLLECTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>       // base class
#include <ddMd/storage/GroupIterator.h>      // member
#include <ddMd/chemistry/Group.h>            // member
#include <util/containers/DArray.h>          // member

namespace DdMd
{

   class Domain;
   class Buffer;
   template <int N> class GroupStorage;

   using namespace Util;

   /**
   * Class for collecting Groups from processors to master processor.
   *
   * An GroupCollector collects Group objects from all processors to the
   * master processor in order to output a configuration file. 
   *
   * Usage:
   * \code
   * 
   *    GroupStorage<N> storage;
   *    Domain domain;
   *    Buffer buffer;
   *    GroupCollector collector;
   *
   *    // Initialization
   *    collector.associate(domain, storage, buffer);
   *
   *    // Communication
   *    if (domain.gridRank() == 0) {  // if master processor
   *       collector.allocate(100);
   *       collector.setup();
   *       Atom* atomPtr = collector.nextPtr();
   *       while (atomPtr) {
   *          // Write *atomPtr to file; 
   *          atomPtr = collector.nextPtr();
   *       }
   *    } else { // if not master
   *       collector.send();
   *    }
   *
   * \endcode
   *
   * \ingroup DdMd_Communicate_Module
   */
   template <int N>
   class GroupCollector : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      GroupCollector();

      /**
      * Destructor.
      */
      ~GroupCollector();

      /**
      * Initialize pointers to associated objects.
      *
      * Call on all processors, only once.
      */
      void associate(Domain& domain, GroupStorage<N>& storage, Buffer& buffer);

      /**
      * Allocate cache on master processor.
      *
      * Call only on the master processor, only once.
      */
      void allocate(int cacheSize);

      /**
      * Setup master processor for receiving.
      *
      * Call only on the master processor, just before each receive loop.
      */
      void setup();

      /**
      * Return a pointer to the next available atom, or null. 
      *
      * Call only on the master processor, within loop over groups.
      *
      * \return address of next atom, or null if no more are available.
      */
      Group<N>* nextPtr();
     
      /**
      * Send all groups to the master.
      *
      * Call on all processors except the master.
      */
      void send();

   private:

      /// Temporary array of groups, allocated only on master.
      DArray< Group<N> > recvArray_;

      /// Iterator for groups in a GroupStorage<N> (master and slaves).
      GroupIterator<N> iterator_;

      /// Pointer to associated Domain object (on master).
      Domain* domainPtr_;

      /// Pointer to associated Domain object (on master).
      GroupStorage<N>* storagePtr_;

      /// Pointer to associated Buffer object (on master).
      Buffer* bufferPtr_;

      /// Rank of processor from which groups are being received (on master).
      int source_;

      /// Number of items in receive buffer (on master).
      int recvBufferSize_;

      /// Number of items in recvArray_ (on master).
      int recvArraySize_;

      /// Index of current item in recvArray_ (on master).
      int recvArrayId_;

      /// Have all groups been processed from current source? (all).
      bool isComplete_;

   };

}
#endif