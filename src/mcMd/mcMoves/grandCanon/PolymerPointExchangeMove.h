#ifndef MCMD_POLYMER_POINT_EXCHANGE_MOVE_H
#define MCMD_POLYMER_POINT_EXCHANGE_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h> //base class

namespace McMd
{

  using namespace Util; //???

  class McSystem;

  /**
  *  Insertion and deletion of point particle
  */

  class PolymerPointExchangeMove : public SystemMove
  {
  
  public:

     /**
     * Constructor
     */
     PolymerPointExchangeMove(McSystem& system);

     /**
     *   read the parameters
     */ 

     virtual void readParameters(std::istream& in);
     

     /**
     *   Load parameters from archive
     */

     virtual void loadParameters(Serializable::IArchive &ar);

     /**
     *  Save state to an archive
     */

     virtual void save(Serializable::OArchive &ar);

     /**
     *  Generate attempt and accept or reject the move
     */

     virtual bool move();

  private:

     /**
     *    Insertion move
     */
     
     bool insertion();

     /**
     *  Deletion move
     */

     bool deleteChain();

     int pointSpeciesId_;
     int polymerSpeciesId_;
     double polymerPotential_;
     double currentNatom_;
     int moveType_;     
     int N_0;
     int N_P;
     int N_S;
   };

}
#endif
