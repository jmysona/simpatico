#ifndef MCMD_SECOND_ORDER_INTERACTION_H
#define MCMD_SECOND_ORDER_INTERACTION_H


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h> // base class template

#include <util/containers/DArray.h> 
#include <util/containers/DMatrix.h> 

namespace McMd
{

  using namespace Util;

  /** 
  *  
  *  Second order interaction calculates the interaction parameters for perturbation theory of 
  *  a polymer blend to second order
  *
  **/
  class SecondOrderInteraction : public SystemAnalyzer<System>
  {
  
  public:
   
   /**
   * Constructor
   *
   *
   *
   */
   SecondOrderInteraction(System &system);
   
   /**
   * Read parameters from file.
   *
   *
   *
   **/
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
   * Clear accumulator.
   */
   virtual void setup();

   /** 
   * Calculat I_a_b pair matrix and related quantities
   *
   * \param iStep step counter
   */
   void sample(long iStep);

   /** 
   * Output results to output file.
   */
   virtual void output();

   private:

   //Output file stream
   std::ofstream outputFile_;
  
   //Species Id
   int speciesId_;

   //I_a_b matrix
   DMatrix<double> Iab_;
   
   // has read param been called
   bool isInitialized_;
  };
 /**
   * Serialize to/from an archive.
   */
   template <class Archive>
   void SecondOrderInteraction::serialize(Archive& ar, const unsigned int version)
   {
      Analyzer::serialize(ar, version);
      ar & speciesId_;

   }
}
#endif
   
