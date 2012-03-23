#ifndef MCMD_MPI_CPP
#define MCMD_MPI_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McMd_mpi.h"
#ifdef UTIL_MPI

namespace Util
{
   /*
   * Initialize MpiTraits< SpeciesGroup<4> >
   */
   MPI::Datatype MpiTraits< McMd::SpeciesGroup<4> >::type = MPI::BYTE;
   bool MpiTraits< McMd::SpeciesGroup<4> >::hasType = false;

   /*
   * Initialize MpiTraits< SpeciesGroup<3> >
   */
   MPI::Datatype MpiTraits< McMd::SpeciesGroup<3> >::type = MPI::BYTE;
   bool MpiTraits< McMd::SpeciesGroup<3> >::hasType = false;

   /*
   * Initialize MpiTraits< SpeciesGroup<2> >
   */
   MPI::Datatype MpiTraits< McMd::SpeciesGroup<2> >::type = MPI::BYTE;
   bool MpiTraits< McMd::SpeciesGroup<2> >::hasType = false;

   /*
   * Initialize MpiTraits< Pair<int> >
   */
   MPI::Datatype MpiTraits< Pair<int> >::type = MPI::BYTE;
   bool MpiTraits< Pair<int> >::hasType = false;

}

#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <mcMd/chemistry/SpeciesGroup.tpp>
#include <mcMd/diagnostics/util/PairSelector.h>

namespace McMd
{

   void commitMpiTypes()
   {
      Util::Vector::commitMpiType();
      Util::IntVector::commitMpiType();
      Util::Pair<int>::commitMpiType();
      SpeciesGroup<2>::commitMpiType();
      SpeciesGroup<3>::commitMpiType();
      SpeciesGroup<4>::commitMpiType();
      PairSelector::commitMpiType();
   }

}
#endif // ifdef  UTIL_MPI
#endif // ifndef MCMD_MPI_CPP
