

#include "PointInsertionMove.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/simulation/Simulation.h>
#include <util/boundary/Boundary.h>
#include <mcMd/mcSimulation/mc_potentials.h>
#include <mcMd/species/Species.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor
   */
   PointInsertionMove::PointInsertionMove(McSystem& system)
    : SystemMove(system),
      insertionPotential_(0.0),
      speciesId_(-1)
   { setClassName("PointInsertionMove");}

   /*
   *  Read the parametrs
   */
   void PointInsertionMove::readParameters(std::istream& in)
   { 
      readProbability(in);
      read<int>(in, "speciesId", speciesId_);
      read<double>(in, "chemicalPotential", insertionPotential_);
   }
   /*
   *  Read the parametrs from archive
   */
   void PointInsertionMove::loadParameters(Serializable::IArchive &ar)
   {
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "speciesId", speciesId_);
      loadParameter<double>(ar, "chemicalPotential", insertionPotential_);
   }

   /*
   *  Save state to archive
   */
   void PointInsertionMove::save(Serializable::OArchive &ar)
   {
     McMove::save(ar);
     ar << speciesId_;
     ar << insertionPotential_;
   }
   
   /*
   *   The move itself
   */

   bool PointInsertionMove::move()
   {  
      incrementNAttempt();
      int atomMaxCap;
      int moveType_;
      currentNatom_=system().nMolecule(speciesId_);
      atomMaxCap=system().simulation().species(speciesId_).capacity();
      moveType_ = simulation().random().uniformInt(0,2);
      if (currentNatom_==0) {
         moveType_ = 0;
      }
      if (currentNatom_==atomMaxCap){
         moveType_ = 1;
      }
      if (moveType_ ==0)
      {
        return insertion();    
      }
      if (moveType_ ==1)
      {
        return deletion();
      } 
   }

   /*
   *
   */
   bool PointInsertionMove::deletion() 
   {
      Molecule* molPtr;
      molPtr = &(system().randomMolecule(speciesId_));
      Atom* theAtom;
      double V;
      theAtom = &molPtr->atom(0);
      double energy;
      #ifndef INTER_NOPAIR
      energy=system().atomPotentialEnergy(*theAtom);
      #else
      energy=0.0;
      #endif
      V = boundary().volume();
      bool accept = random().metropolis(boltzmann(-energy+insertionPotential_)*currentNatom_/V);
      if (accept) {
         incrementNAccept();
        // system().pairPotential().deleteAtom(*theAtom);
         system().removeMolecule(*molPtr);
         simulation().returnMolecule(*molPtr);
         return accept;
      } else {
         return accept;
    }
         
   }

   bool PointInsertionMove::insertion()
   {  Vector r;
      double V;
      Molecule* molPtr;
      molPtr = &(simulation().getMolecule(speciesId_));
      Atom* theAtom;
      theAtom = &molPtr->atom(0);
      system().addMolecule(*molPtr);
      boundary().randomPosition(random(), r);
      V = boundary().volume();
      boundary().shift(r);
      for (int j = 0; j < Dimension; ++j) {
         theAtom->position()[j] = r[j];
      }
      double energy;
      #ifndef INTER_NOPAIR
      energy=system().atomPotentialEnergy(*theAtom);
      #else
      energy=0.0;
      #endif
      bool accept = random().metropolis(boltzmann(energy-insertionPotential_)*V/(currentNatom_+1));
      if (accept) {
         incrementNAccept(); 
      } else {
          //system().pairPotential().deleteAtom(*theAtom);
          system().removeMolecule(*molPtr);
          simulation().returnMolecule(*molPtr);
      }     
      return accept;
   }
}

   /*
   *
   */





