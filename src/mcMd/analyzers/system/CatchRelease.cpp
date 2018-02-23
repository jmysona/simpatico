#ifndef MCMD_CATCH_RELEASE_CPP
#define MCMD_CATCH_RELEASE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "CatchRelease.h"
#include <mcMd/simulation/System.h>
#include <mcMd/simulation/Simulation.h>
#include <simp/species/Species.h>
#include <mcMd/chemistry/Molecule.h>
#include <mcMd/chemistry/Atom.h>
#include <mcMd/species/LinearSG.h>
#include <util/misc/FileMaster.h>        
#include <util/archives/Serializable_includes.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/ioUtil.h>
#include <sstream>

namespace McMd
{

   using namespace Util;

   /// Constructor.
   CatchRelease::CatchRelease(System& system) 
    : SystemAnalyzer<System>(system),
      identifier_(system),
      hist_(),
      outputFile_(),
      speciesId_(),
      atomTypeId_(),
      cutoff_(),
      histMin_(),
      histMax_(),
      beadNumber_(),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("MicelleFlux"); }

   /// Read parameters from file, and allocate arrays.
   void CatchRelease::readParameters(std::istream& in) 
   {  
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }
      speciesPtr_ = dynamic_cast<LinearSG*>(&(system().simulation().species(speciesId_))); 
      isMutable_ = (speciesPtr_->isMutable());
      if (isMutable_) {
         read<int>(in, "speciesSubtype", speciesSubtype_);
      }
      read<int>(in, "atomTypeId", atomTypeId_);
      if (atomTypeId_ < 0) {
         UTIL_THROW("Negative atomTypeId");
      }
      if (atomTypeId_ >= system().simulation().nAtomType()) {
         UTIL_THROW("nTypeId >= nAtomType");
      }

      read<double>(in, "cutoff", cutoff_);

      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }

      read<double>(in, "radius", radius_);

      read<int>(in, "beadNumber", beadNumber_);
      read<int>(in, "tauFlip", tauFlip_);
      read<int>(in, "diblockType", diblockType_);
      read<int>(in, "homopolymerType", homopolymerType_);
      // Initialize ClusterIdentifier
      identifier_.initialize(speciesId_, atomTypeId_, cutoff_);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
      isInitialized_ = true;
      int speciesCapacity=system().simulation().species(speciesId_).capacity();
      micelleFlux_.allocate(speciesCapacity);
      priorMicelleFlux_.allocate(speciesCapacity);
      InMicelle_.allocate(speciesCapacity);
      
      for (int i = 0; i < speciesCapacity; i++) {
        priorMicelleFlux_[i]=-1;
      }
   }

   /*
   * Load state from an archive.
   */
   void CatchRelease::loadParameters(Serializable::IArchive& ar)
   {
      // Load interval and outputFileName
      Analyzer::loadParameters(ar);

      loadParameter<int>(ar,"speciesId", speciesId_);
      if (speciesId_ < 0) {
         UTIL_THROW("Negative speciesId");
      }
      if (speciesId_ >= system().simulation().nSpecies()) {
         UTIL_THROW("speciesId > nSpecies");
      }

      Species* speciesPtr;
      speciesPtr = &(system().simulation().species(speciesId_)); 
      isMutable_ = (speciesPtr->isMutable());
      if (isMutable_) {
         loadParameter<int>(ar, "speciesSubtype", speciesSubtype_);
      }

      loadParameter<int>(ar, "atomTypeId", atomTypeId_);
      if (atomTypeId_ < 0) {
         UTIL_THROW("Negative atomTypeId");
      }

      loadParameter<double>(ar, "cutoff", cutoff_);
      if (cutoff_ < 0) {
         UTIL_THROW("Negative cutoff");
      }

      loadParameter<double>(ar, "radius", radius_);
      
      loadParameter<int>(ar, "beadNumber", beadNumber_);
      identifier_.initialize(speciesId_, atomTypeId_, cutoff_);
      ar >> nSample_;

      isInitialized_ = true;
   }

   /*
   * Save state to archive.
   */
   void CatchRelease::save(Serializable::OArchive& ar)
   {
      Analyzer::save(ar);
      ar & speciesId_;
      ar & atomTypeId_;
      ar & cutoff_;
      ar & nSample_;
   }

   /*
   * Clear accumulators.
   */
   void CatchRelease::setup() 
   {  
      if (!isInitialized_) UTIL_THROW("Object is not initialized");
      nSample_ = 0;
   }

   /* 
   * Evaluate end-to-end vectors of all chains, add to ensemble.
   */
   void CatchRelease::sample(long iStep) 
   { 
     if (isAtInterval(iStep)) {
       identifier_.identifyClusters();
       ++nSample_;
       bool needsFlipping = false;
        if((iStep % tauFlip_)==0) {
          needsFlipping = true;
        }
       int moleculeToFlip;
       int bigClusterId;
       Species* speciesPtr;
       speciesPtr = &(system().simulation().species(speciesId_)); 
       int nMolecules = speciesPtr->capacity();
       ClusterLink* LinkPtr;      
       // Write the time step to the data file
       // Cycle through all the surfactant molecules; if they are in a micelle set status = 1, otherwise set status = 0;
       int numberMoleculesInMicelle=0;
       for (int i = 0; i < nMolecules; i++) {
         bool inCluster = false;
         // Get the molecule link
         LinkPtr = &(identifier_.link(i));
         // Use the link to get the cluster id
         int clusterId = LinkPtr->clusterId();
         // Determine the cluster size
         int clusterSize = identifier_.cluster(clusterId).size();
         if (clusterSize > 10) {
           InMicelle_[i]=1;
           numberMoleculesInMicelle += 1;
           bigClusterId = clusterId;
         } else {
           InMicelle_[i]=0;
         }
       }
       if (needsFlipping) {
         moleculeToFlip=system().simulation().random().uniformInt(0,numberMoleculesInMicelle);
         LinkPtr = identifier_.cluster(bigClusterId).head();
         for (int i = 0; i < moleculeToFlip; i++) {
           LinkPtr->next();
         }
         Molecule& molecule = LinkPtr->molecule();
         speciesPtr_->mutator().setMoleculeState(molecule, homopolymerType_);
       }
       Vector r;
       // Cluster COM
       micelleCOM_ = identifier_.cluster(bigClusterId).clusterCOM(atomTypeId_,system().boundary());
       //
       Vector lengths = system().boundary().lengths();
       double distance;
       for (int i = 0; i < nMolecules; i++) {
         distance = 0;
         r = system().molecule(speciesId_, i).atom(beadNumber_).position();
         Molecule& molecule = system().molecule(speciesId_, i);
         for (int j = 0; j < Dimension; j++) { 
           if (std::abs(micelleCOM_[j]-r[j]) > lengths[j]/2) {
             if ((r[j]-micelleCOM_[j]) > 0)
               distance = distance + pow(micelleCOM_[j]-r[j]+lengths[j],2);
             else 
               distance = distance + pow(micelleCOM_[j]-r[j]-lengths[j],2);
           } else {
             distance = distance + pow(micelleCOM_[j]-r[j],2);       
           } 
         } 
         distance = sqrt(distance);
         if (distance > radius_) {
           micelleFlux_[i] = 0;
           if (speciesPtr_->mutator().moleculeStateId(molecule) == homopolymerType_) {
             speciesPtr_->mutator().setMoleculeState(molecule, diblockType_);
           }
         }
         else if (InMicelle_[i] == 1) {
           micelleFlux_[i] = 1;
         } else {
           micelleFlux_[i] = priorMicelleFlux_[i];
         }
       }
       for (int i = 0; i < nMolecules; i++) {
         outputFile_ << micelleFlux_[i] << "  ";
       if (i == nMolecules -1)
         outputFile_ << "\n";
       }
       priorMicelleFlux_=micelleFlux_;
     } 
   }

   /*
   * Output results to file after simulation is completed.
   */
   void CatchRelease::output() 
   {  outputFile_.close();
      // Write parameter file
      fileMaster().openOutputFile(outputFileName(".prm"), outputFile_);
      writeParam(outputFile_);
      outputFile_.close();

      // Write histogram output
      // hist_.output(outputFile_);
   }

}
#endif 
