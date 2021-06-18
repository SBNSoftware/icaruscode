////////////////////////////////////////////////////////////////////////                  
// \file SpaceChargeICARUS.h                                                           
//                                                                                        
// \brief header of class for storing/accessing space charge distortions for ICARUS
//                                                                                        
// \author rlazur@colostate.edu
//                                                                                        
////////////////////////////////////////////////////////////////////////

#ifndef SPACECHARGE_SPACECHARGEICARUS_H
#define SPACECHARGE_SPACECHARGEICARUS_H

// LArSoft libraries
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larcore/Geometry/Geometry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Principal/Event.h"

// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"
// c++
#include <string>
#include <map>
#include <vector>
//ROOT
#include <TGraph.h>
#include <TF1.h>
#include <TH3.h>
#include <TFile.h>

namespace spacecharge
{
    class SpaceChargeICARUS : public SpaceCharge
    {

    public:

      explicit SpaceChargeICARUS(fhicl::ParameterSet const& pset);
      SpaceChargeICARUS(SpaceChargeICARUS const&) = delete;
      virtual ~SpaceChargeICARUS() = default;

      bool Configure(fhicl::ParameterSet const& pset);
      bool Update(uint64_t ts = 0);

      // sim = forward; Cal = backward
      bool EnableSimSpatialSCE() const override; 
      bool EnableSimEfieldSCE() const override; 
      bool EnableCalSpatialSCE() const override;
      bool EnableCalEfieldSCE() const override;

      bool EnableCorrSCE() const override {return (EnableCalSpatialSCE()||EnableCalEfieldSCE()) ;}

      //Pos and Efield offsets are used primarily in larsim
      //used to calculate e-lifetime, recombination, and dedx
      geo::Vector_t GetPosOffsets(geo::Point_t const& point) const override;
      geo::Vector_t GetEfieldOffsets(geo::Point_t const& point) const override;
      //Cal offsets are for doing calibration and analysis (backwards map)
      //require TPCid to disambiguate hits that cross cathode
      geo::Vector_t GetCalPosOffsets(geo::Point_t const& point, int const& TPCid) const override;
      geo::Vector_t GetCalPosOffsets(geo::Point_t const& point, geo::TPCID const& TPCid) const;
      geo::Vector_t GetCalEfieldOffsets(geo::Point_t const& point, int const& TPCid = 1) const override { return {0.,0.,0.}; }

    private:
    protected:

      /////////////////////////////
      // DECLARE GLOBAL OBJECTS
      ////////////////////////////
      std::vector<TH3F*> SCEhistograms = std::vector<TH3F*>(9);

      //////////////////////////////
      // DECLARE FHICL PARAMETERS
      /////////////////////////////
      bool fEnableSimSpatialSCE;
      bool fEnableSimEfieldSCE;
      bool fEnableCalSpatialSCE;
      bool fEnableCalEfieldSCE;
      bool fEnableCorrSCE;
      std::string fRepresentationType;
      std::string fInputFilename;

      ////////////////////////////////
      // DECLARE SUPPLEMENTAL FUNCTIONS
      ////////////////////////////////
      void fixCoords(double* xx, double* yy, double* zz) const;
      void fixCoords(double* xx, double* yy, double* zz, int* TPCid) const;
    }; // class SpaceChargeICARUS
} //namespace spacecharge
#endif // SPACECHARGE_SPACECHARGEICARUS_H
