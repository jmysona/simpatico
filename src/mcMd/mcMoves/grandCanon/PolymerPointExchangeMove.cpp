

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
      incrementNAttempt();
     /* int pointMaxCap;
      int polymerMaxCap;
      currentNatom_=system().nMolecule(pointSpeciesId_);
      pointMaxCap=system().simulation().species(pointSpeciesId_).capacity();
      polymerNmolecule_=system().nMolecule(polymerSpeciesId_);
      polymerMaxCap=system().simulation().species(polymerSpeciesId_).capacity();
*/    
      //Species* polymerPtr = &system.simulation().species(polymerSpeciesId_);
      //Species* solventPtr = &system.simulation().species(pointSpeciesId_);
      N_P = system().nMolecule(polymerSpeciesId_);
      N_S = system().nMolecule(pointSpeciesId_);
      N_0 = 8;
      double top = 1;
      double bottom = 1;
      for (int i = 0; i<N_0; i++) {
        top = top*(N_P*N_0-i);
        bottom = bottom*(N_S+N_0-i);
      }
      double prefactor=top/bottom;     
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
      
      Atom* atomPtr;
      float prefactor;
      //double oldEnergy = system().pairPotential().moleculeEnergy(*molPtr);
      double oldEnergy=0;
      Vector beadPositions[N_0];
      int    solventBeadIds[N_0];
      for (int i=0; i<N_0; i++) {
          atomPtr = &molPtr->atom(i);
          oldEnergy += system().bondPotential().atomEnergy(*atomPtr);
          //While getting energy, get the positions too
          beadPositions[i]=atomPtr->position();
          oldEnergy += system().pairPotential().atomEnergy(*atomPtr);
          
      }
            for (int atomId = 0; atomId < molPtr->nAtom(); ++atomId) {
                system().pairPotential().deleteAtom(molPtr->atom(atomId));
            }
      system().removeMolecule(*molPtr);
      simulation().returnMolecule(*molPtr);
      // Delete
      // Add the points and calculate energy 
      double newEnergy;
        Molecule* pointPtr;
      Atom* theAtom;
      newEnergy =0;
      for (int i=0; i<N_0; i++) {
        pointPtr = &(simulation().getMolecule(pointSpeciesId_)); 
        system().addMolecule(*pointPtr);
        theAtom =&pointPtr->atom(0);
        theAtom->position()=beadPositions[i];
        solventBeadIds[i]=system().moleculeId(*pointPtr);
        system().pairPotential().addAtom(*theAtom);
      }
      for (int i=0; i<N_0; i++) {
         pointPtr = &(system().molecule(pointSpeciesId_, solventBeadIds[i]));
         theAtom =&pointPtr->atom(0);
         newEnergy += system().pairPotential().atomEnergy(*theAtom);
      }
      //Determine if accepted or not
      bool accept = random().metropolis(prefactor*boltzmann(oldEnergy-newEnergy-polymerPotential_));
      
      if (accept) {
      } else {
        Molecule* polyPtr;        
        polyPtr = &(simulation().getMolecule(polymerSpeciesId_)); 
        system().addMolecule(*polyPtr);
        Atom* thePolymerAtom;
         
        for (int i = 0; i < N_0; i++) {
            thePolymerAtom = &polyPtr->atom(i);
            thePolymerAtom->position() = beadPositions[i];
        }
      }
      return accept;
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





