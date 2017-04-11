#ifndef MD_COULOMB_POTENTIAL_H
#define MD_COULOMB_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>                   // base class
#include <mcMd/potentials/coulomb/EwaldRSpaceAccumulator.h>  // member
#include <util/misc/Setable.h>                           // member template
#include <util/space/Tensor.h>                           // template parameter

#include <complex>

namespace McMd
{

   using namespace Util;

   /**
   * Coulomb potential for an Md simulation.
   *
   * This class computes the long-range k-space part of the
   * Coulomb forces and energies in a molecular dynamics (MD)
   * simulation, and provides accessors for both r-space and 
   * k-space contributions to the Coulomb energy and stress.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class MdCoulombPotential : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      MdCoulombPotential();

      /**
      * Destructor (does nothing).
      */
      virtual ~MdCoulombPotential();

      /// \name System energy and stress.
      //@{

      /**
      * Generate wavevectors for the current boundary.
      */
      virtual void makeWaves() = 0;

      /**
      * Current number of wavevectors.
      */
      virtual int nWave() const = 0;

      /**
      * Add k-space Coulomb forces for all atoms.
      */
      virtual void addForces() = 0;

      /**
      * Compute kspace part of Coulomb stress.
      *
      * \param stress (output) pressure
      */
      virtual void computeStress() = 0;
   
      //@}
      /// \name Energy
      //@{

      /**
      * Unset k-space energy.
      */
      virtual void unsetEnergy();

      /**
      * Calculate the long range kspace part of Coulomb energy.
      */
      virtual void computeEnergy() = 0;

      /**
      * Get long-range k-space part of Coulomb energy.
      */
      double kSpaceEnergy();

      /**
      * Return short-range r-space part of Coulomb energy.
      */
      double rSpaceEnergy();

      /**
      * Get total Coulomb energy.
      */
      double energy();

      //@}
      /// \name Stress
      //@{

      /**
      * Unset k-space stress.
      */
      virtual void unsetStress();

      /**
      * Get long-range k-space part of Coulomb stress.
      */
      Tensor kSpaceStress();

      /**
      * Return short-range r-space part of Coulomb stress.
      */
      Tensor rSpaceStress();

      /**
      * Get total Coulomb stress.
      */
      Tensor stress();

      //@}

   protected:

      /// R-space energy and stress contributions
      EwaldRSpaceAccumulator rSpaceAccumulator_;

      /// K-space part of Coulomb energy
      Setable<double> kSpaceEnergy_;

      /// K-space part of Coulomb stress.
      Setable<Tensor> kSpaceStress_;

      /// Have parameters been set?
      bool isInitialized_;

   };

} 
#endif