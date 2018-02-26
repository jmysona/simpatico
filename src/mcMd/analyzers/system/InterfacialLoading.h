#ifndef MCMD_INTERFACIAL_LOADING_H
#define MCMD_INTERFACIAL_LOADING_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>      // base class template
#include <mcMd/simulation/System.h>                 // base class template parameter
#include <util/accumulators/IntDistribution.h>
#include <mcMd/neighbor/CellList.h>

#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /**
   * Calculates interfacial loading based off of contact between two atom types
   */
   class InterfacialLoading : public SystemAnalyzer<System>
   {
   
   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      InterfacialLoading(System &system);

      /**
      * Read parameters from file.  
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);
  
      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /** 
      * Setup before a simulation (clear accumulator).
      */
      virtual void setup();

      /** 
      * Check the interfacial loading at the iStep and add to histogram
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);

      /** 
      * Output results to output file.
      */
      virtual void output();

   private:

      // Output file stream
      std::ofstream outputFile_;

      // InterfacialLoading statistical accumulator
      IntDistribution  accumulator_;

      /// Is this initialized (Has readParam or loadParam been called?)
      bool    isInitialized_;
      
      int species1Id_;
      int atomType1Id_;
      int species2Id_;
      int atomType2Id_;
      double cutoff_;
      int histMax_;

      CellList cellList_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void InterfacialLoading::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & accumulator_;
   }

}
#endif
