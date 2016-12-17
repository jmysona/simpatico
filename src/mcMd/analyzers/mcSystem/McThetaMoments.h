#ifndef MCMD_MC_THETA_MOMENTS_H
#define MCMD_MC_THETA_MOMENTS_H



#include <mcMd/analyzers/SystemAnalyzer.h> 

namespace mcMd
{

  using namespace Util;

  //Calculate moments of the quantity theta in the semigrand ensemble
  //for second order perturbation purposes

  class McThetaMoments : public SystemAnalyzer<McSystem>
  {
  
  public:

    /*
     *  Constructor
     */

     McThetaMoments(McSystem& system);
    
    /**
    *  Read parameters
    */ 
     
    virtual void readParameters(std::istream& in);

    /**
    *   Clear accumulators and other necessary setup
    */

    virtual void setup();

    /**
    *   Compute the I matrix and calculate the relevant theta quantities
    */
    virtual void sample(long iStep);

    /**
    *   Output Results at the end
    */

    virtual void output();

    /**
    * Save state to binary file archive
    *
    * \param ar binary saving (output) archive.
    */
    virtual void save(Serializable::OArchive& ar);
 
    /**
    * Load state from a binary file archive.
    *
    * \param ar binary loading (input) archive.
    */
    virtual void loadParameters(Serializable::IArchive& ar);
 
    /**
    * Serialize to/from an archive. 
    *
    * \param ar      saving or loading archive
    * \param version archive version id
    */
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version);
 
  private:
    
     // output file stream
    
     std::ofstream outputFile_;

     
 }
#endif
