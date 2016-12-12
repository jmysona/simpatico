

#include "PolymerPointExchangeMove.h"
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
   PolymerPointExchangeMove::PolymerPointExchangeMove(McSystem& system)
    : SystemMove(system),
      polymerPotential_(0.0),
      pointSpeciesId_(-1),
      polymerSpeciesId_(-1)
   { setClassName("PolymerPointExchangeMove");}

   /*
   *  Read the parametrs
   */
   void PolymerPointExchangeMove::readParameters(std::istream& in)
   { 
      readProbability(in);
      read<int>(in, "pointSpeciesId", pointSpeciesId_);
      read<int>(in, "polymerSpeciesId", polymerSpeciesId_);
      read<double>(in, "polymerChemicalPotential", polymerPotential_);
   }
   /*
   *  Read the parametrs from archive
   */
   void PolymerPointExchangeMove::loadParameters(Serializable::IArchive &ar)
   {
      McMove::loadParameters(ar);
      loadParameter<int>(ar, "pointSpeciesId", pointSpeciesId_);
      loadParameter<int>(ar, "polymerSpeciesId", polymerSpeciesId_);
      loadParameter<double>(ar, "polymerChemicalPotential", polymerPotential_);
   }

   /*
   *  Save state to archive
   */
   void PolymerPointExchangeMove::save(Serializable::OArchive &ar)
   {
     McMove::save(ar);
     ar << pointSpeciesId_;
     ar << polymerSpeciesId_;
     ar << polymerPotential_;
   }
   
   /*
   *   The move itself
   */

   bool PolymerPointExchangeMove::move()
   {  
     /* incrementNAttempt();
      int pointMaxCap;
      int polymerMaxCap;
      currentNatom_=system().nMolecule(pointSpeciesId_);
      pointMaxCap=system().simulation().species(pointSpeciesId_).capacity();
      polymerNmolecule_=system().nMolecule(polymerSpeciesId_);
      polymerMaxCap=system().simulation().species(polymerSpeciesId_).capacity();
*/

      int moveType_=0;
     // atomMaxCap=system().simulation().species(speciesId_).capacity();
     // moveType_ = simulation().random().uniformInt(0,2);
   /*   if (currentNatom_==0) {
    *     moveType_ = 0;
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
     */
       return deleteChain();
   }

   /*
   *
   */
   bool PolymerPointExchangeMove::deleteChain() 
   {
 
      Vector r;
      Molecule* molPtr;
      // Choose a chain to delete
      molPtr = &(system().randomMolecule(polymerSpeciesId_));
      Atom* theAtom;
      double V;
      for (int i = 0; i < nAtomsInMolecule; ++i) {
        theAtom = &molPtr->atom(i);
          for (int j = 0; j < Dimension; ++j) {
            r[j]=theAtom->position()[j];
          }
        double energy;
        #ifndef INTER_NOPAIR
        energy=system().atomPotentialEnergy(*theAtom);
        #else
        energy=0.0;
        #endif
      }
      V = boundary().volume();
 /*     bool accept = random().metropolis(boltzmann(-totalEnergy+insertionPotential_)*currentNatom_/V);
      if (accept) {
         incrementNAccept();
        // system().pairPotential().deleteAtom(*theAtom);
         system().removeMolecule(*molPtr);
         simulation().returnMolecule(*molPtr);
         return accept;
      } else {
         return accept;
    }
   */      
   }

   bool PolymerPointExchangeMove::insertion()
   { 
   /*   Vector r;
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
    */
   }
}

   /*
   *
   */





