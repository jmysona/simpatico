/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "StressCalculator.h"
#include <util/space/Vector.h>

namespace McMd 
{

   using namespace Util;   

   /*
   * Constructor (protected).
   */
   StressCalculator::StressCalculator(bool createsStress)
    : createsStress_(createsStress)
   {  stress_.unset(); }

   /*
   * Return true iff subclass is a potential energy that creates stress.
   */
   bool StressCalculator::createsStress() const
   {  return createsStress_; }

   /*
   * Mark the stress as unknown.
   */
   void StressCalculator::unsetStress()
   {  stress_.unset(); }

   /*
   * Get the nonbonded stress tensor.
   */
   void StressCalculator::computeStress(Tensor& stress)
   {
      UTIL_CHECK(createsStress_);

      // If necessary, compute stress tensor
      if (!stress_.isSet()) {
         computeStress();
      }

      // Get full stress tensor
      stress = stress_.value();
   }

   /*
   * Get the nonbonded x, y, z pressures
   */
   void StressCalculator::computeStress(Vector& pressures)
   {
      UTIL_CHECK(createsStress_);

      // If necessary, compute stress tensor
      if (!stress_.isSet()) {
         computeStress();
      }

      // Get diagonal components of stress tensor (pressures)
      for (int i=0; i < Dimension; ++i) {
         pressures[i] = stress_.value()(i, i);
      }
   }

   /*
   * Get the isotropic pressure = Tr(stress)/3
   */
   void StressCalculator::computeStress(double& pressure)
   {
      UTIL_CHECK(createsStress_);

      // If necessary, compute stress tensor
      if (!stress_.isSet()) {
         computeStress();
      }

      // Get pressure = average of diagonal components.
      pressure = 0.0;
      for (int i=0; i < Dimension; ++i) {
         pressure += stress_.value()(i, i);
      }
      pressure = pressure/double(Dimension);
   }

}
