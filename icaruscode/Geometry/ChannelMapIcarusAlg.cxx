////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapIcarusAlg.cxx
/// \brief Interface to algorithm class for the standar, simplest detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "ChannelMapIcarusAlg.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"

#include "messagefacility/MessageLogger/MessageLogger.h" 

namespace geo{

//----------------------------------------------------------------------------
ChannelMapIcarusAlg::ChannelMapIcarusAlg(fhicl::ParameterSet const& p)
  : fSorter(geo::GeoObjectSorterStandard(p))
{
}

//----------------------------------------------------------------------------
void ChannelMapIcarusAlg::Initialize( GeometryData_t& geodata )
{
    // start over:
    Uninitialize();
  
    std::vector<geo::CryostatGeo*>& cgeo = geodata.cryostats;
    std::vector<geo::AuxDetGeo*>  & adgeo = geodata.auxDets;
  
    fNcryostat = cgeo.size();
  
    mf::LogInfo("ChannelMapIcarusAlg") << "Initializing Standard ChannelMap...";

    fSorter.SortCryostats(cgeo);
    fSorter.SortAuxDets(adgeo);
    for(size_t c = 0; c < cgeo.size(); ++c)
        cgeo[c]->SortSubVolumes(fSorter);
  
    fNTPC.resize(fNcryostat);
    fWireCounts.resize(fNcryostat);
    fNPlanes.resize(fNcryostat);
    fFirstWireProj.resize(fNcryostat);
    fOrthVectorsY.resize(fNcryostat);
    fOrthVectorsZ.resize(fNcryostat);
    fPlaneBaselines.resize(fNcryostat);
    fWiresPerPlane.resize(fNcryostat);
        fFirstChannelInNextPlane.resize(fNcryostat);
    fFirstChannelInThisPlane.resize(fNcryostat);
    fViews.clear();
    fPlaneIDs.clear();
    fTopChannel = 0;

    int RunningTotal = 0;

    for(unsigned int cs = 0; cs != fNcryostat; ++cs)
    {
        fNTPC[cs] = cgeo[cs]->NTPC();
    
        // Size up all the vectors
        fWireCounts[cs]             .resize(fNTPC[cs]);
        fFirstWireProj[cs]           .resize(fNTPC[cs]);
        fOrthVectorsY[cs]            .resize(fNTPC[cs]);
        fOrthVectorsZ[cs]            .resize(fNTPC[cs]);
        fPlaneBaselines[cs]          .resize(fNTPC[cs]);
        fWiresPerPlane[cs]           .resize(fNTPC[cs]);
        fNPlanes[cs]                     .resize(fNTPC[cs]);
        fFirstChannelInThisPlane[cs].resize(fNTPC[cs]);
        fFirstChannelInNextPlane[cs].resize(fNTPC[cs]);

        for(unsigned int TPCCount = 0; TPCCount != fNTPC[cs]; ++TPCCount)
        {
            unsigned int PlanesThisTPC = cgeo[cs]->TPC(TPCCount).Nplanes();
            fWireCounts[cs][TPCCount]   .resize(PlanesThisTPC);
            fFirstWireProj[cs][TPCCount].resize(PlanesThisTPC);
            fOrthVectorsY[cs][TPCCount] .resize(PlanesThisTPC);
            fOrthVectorsZ[cs][TPCCount] .resize(PlanesThisTPC);
            fNPlanes[cs][TPCCount]=PlanesThisTPC;
            for(unsigned int PlaneCount = 0; PlaneCount != PlanesThisTPC; ++PlaneCount)
            {
                fViews.emplace(cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).View());
                fPlaneIDs.emplace(PlaneID(cs, TPCCount, PlaneCount));
                double ThisWirePitch = cgeo[cs]->TPC(TPCCount).WirePitch(0, 1, PlaneCount);
                fWireCounts[cs][TPCCount][PlaneCount] = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Nwires();
        
                double  WireCentre1[3] = {0.,0.,0.};
                double  WireCentre2[3] = {0.,0.,0.};
        
                const geo::WireGeo& firstWire = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Wire(0);
                const double sth = firstWire.SinThetaZ(), cth = firstWire.CosThetaZ();
        
                firstWire.GetCenter(WireCentre1,0);
                cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Wire(1).GetCenter(WireCentre2,0);
        
                // figure out if we need to flip the orthogonal vector
                // (should point from wire n -> n+1)
                double OrthY = cth, OrthZ = -sth;
                if(((WireCentre2[1] - WireCentre1[1])*OrthY + (WireCentre2[2] - WireCentre1[2])*OrthZ) < 0)
                {
                    OrthZ *= -1;
                    OrthY *= -1;
                }
        
                // Overall we are trying to build an expression that looks like
                //  int NearestWireNumber = round((worldPos.OrthVector - FirstWire.OrthVector)/WirePitch);
                // That runs as fast as humanly possible.
                // We predivide everything by the wire pitch so we don't do this in the loop.
                //
                // Putting this together into the useful constants we will use later per plane and tpc:
                fOrthVectorsY[cs][TPCCount][PlaneCount] = OrthY / ThisWirePitch;
                fOrthVectorsZ[cs][TPCCount][PlaneCount] = OrthZ / ThisWirePitch;
        
                fFirstWireProj[cs][TPCCount][PlaneCount]  = WireCentre1[1]*OrthY + WireCentre1[2]*OrthZ;
                fFirstWireProj[cs][TPCCount][PlaneCount] /= ThisWirePitch;
        
                // now to count up wires in each plane and get first channel in each plane
                int WiresThisPlane = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Nwires();
                fWiresPerPlane[cs] .at(TPCCount).push_back(WiresThisPlane);
                fPlaneBaselines[cs].at(TPCCount).push_back(RunningTotal);
        
                RunningTotal += WiresThisPlane;

                fFirstChannelInThisPlane[cs].at(TPCCount).push_back(fTopChannel);
                fTopChannel += WiresThisPlane;
                fFirstChannelInNextPlane[cs].at(TPCCount).push_back(fTopChannel);
            }// end loop over planes
        }// end loop over TPCs
    }// end loop over cryostats

    // calculate the total number of channels in the detector
    fNchannels = fTopChannel;

    LOG_DEBUG("ChannelMapStandard") << "# of channels is " << fNchannels;

    return;
}
 
//----------------------------------------------------------------------------
void ChannelMapIcarusAlg::Uninitialize()
{
}

//----------------------------------------------------------------------------
std::vector<geo::WireID> ChannelMapIcarusAlg::ChannelToWire(raw::ChannelID_t channel)  const
{
    std::vector< geo::WireID > AllSegments;
    unsigned int cstat = 0;
    unsigned int tpc   = 0;
    unsigned int plane = 0;
    unsigned int wire  = 0;
  
    // first check if this channel ID is legal
    if(channel > fTopChannel)
        throw cet::exception("Geometry") << "ILLEGAL CHANNEL ID for channel " << channel << "\n";
  
    // then go find which plane, tpc and cryostat it is in from the information we stored earlier
    bool foundWid(false);
    for(unsigned int csloop = 0; csloop != fNcryostat; ++csloop){
        for(unsigned int tpcloop = 0; tpcloop != fNTPC[csloop]; ++tpcloop){
            for(unsigned int planeloop = 0; planeloop != fFirstChannelInNextPlane[csloop][tpcloop].size(); ++planeloop){
                if(channel < fFirstChannelInNextPlane[csloop][tpcloop][planeloop]){
                    cstat = csloop;
                    tpc   = tpcloop;
                    plane = planeloop;
                    wire  = channel - fFirstChannelInThisPlane[cstat][tpcloop][planeloop];
                    foundWid = true;
                    break;
                }
                if(foundWid) break;
            }// end plane loop
            if(foundWid) break;
        }// end tpc loop
        if(foundWid) break;
    }// end cryostat loop

    geo::WireID CodeWire(cstat, tpc, plane, wire);
  
    AllSegments.push_back(CodeWire);

    return AllSegments;
}

//----------------------------------------------------------------------------
raw::ChannelID_t ChannelMapIcarusAlg::Nchannels() const
{
    return fNchannels;
}

//----------------------------------------------------------------------------
unsigned int ChannelMapIcarusAlg::Nchannels
  (readout::ROPID const& ropid) const
{
    if (!HasROP(ropid)) return 0;
    // The number of channels matches the number of wires. Life is easy.
    return WireCount(FirstWirePlaneInROP(ropid));
} // ChannelMapIcarusAlg::Nchannels(ROPID)


//----------------------------------------------------------------------------
double ChannelMapIcarusAlg::WireCoordinate
  (double YPos, double ZPos, geo::PlaneID const& planeID) const
{
    // Returns the wire number corresponding to a (Y,Z) position in PlaneNo
    // with float precision.
    // B. Baller August 2014
    return    YPos*AccessElement(fOrthVectorsY, planeID)
            + ZPos*AccessElement(fOrthVectorsZ, planeID)
            - AccessElement(fFirstWireProj, planeID);
}


//----------------------------------------------------------------------------
WireID ChannelMapIcarusAlg::NearestWireID
  (const TVector3& worldPos, geo::PlaneID const& planeID) const
{
    // This part is the actual calculation of the nearest wire number, where we assume
    //  uniform wire pitch and angle within a wireplane
    // add 0.5 to have the correct rounding
    int NearestWireNumber = int(0.5 + WireCoordinate(worldPos.Y(), worldPos.Z(), planeID));
  
    // If we are outside of the wireplane range, throw an exception
    // (this response maintains consistency with the previous
    // implementation based on geometry lookup)
    if(NearestWireNumber < 0 || (unsigned int) NearestWireNumber >= WireCount(planeID))
    {
        int wireNumber = NearestWireNumber; // save for the output
    
        if(NearestWireNumber < 0 ) NearestWireNumber = 0;
        else                       NearestWireNumber = WireCount(planeID) - 1;
    
        throw InvalidWireIDError("Geometry", wireNumber, NearestWireNumber)
        << "Can't Find Nearest Wire for position ("
        << worldPos[0] << "," << worldPos[1] << "," << worldPos[2] << ")"
        << " in plane " << std::string(planeID) << " approx wire number # "
        << wireNumber << " (capped from " << NearestWireNumber << ")\n";
    }

    return geo::WireID(planeID, (geo::WireID::WireID_t) NearestWireNumber);
}

//----------------------------------------------------------------------------
// This method returns the channel number, assuming the numbering scheme
// is heirachical - that is, channel numbers run in order, for example:
//                                             (Ben J Oct 2011)                   
//                    Wire1     | 0
//           Plane1 { Wire2     | 1
//    TPC1 {          Wire3     | 2
//           Plane2 { Wire1     | 3   increasing channel number
//                    Wire2     | 4     (with no gaps)
//    TPC2 { Plane1 { Wire1     | 5
//           Plane2 { Wire1     | 6
//                    Wire2     v 7
//
raw::ChannelID_t ChannelMapIcarusAlg::PlaneWireToChannel
  (geo::WireID const& wireID) const
{
    unsigned int const* pBaseLine = GetElementPtr(fPlaneBaselines, wireID);
    // This is the actual lookup part - first make sure coordinates are legal
    if (pBaseLine) {
        // if the channel has legal coordinates, its ID is given by the wire
        // number above the number of wires in lower planes, tpcs and cryostats
        return *pBaseLine + wireID.Wire;
    }
    else{
        // if the coordinates were bad, throw an exception
        throw cet::exception("ChannelMapIcarusAlg")
            << "NO CHANNEL FOUND for " << std::string(wireID);
    }
  
    // made it here, that shouldn't happen, return raw::InvalidChannelID
    mf::LogWarning("ChannelMapIcarusAlg") << "should not be at the point in the function, returning "
                                          << "invalid channel";
    return raw::InvalidChannelID;
}


//----------------------------------------------------------------------------
SigType_t ChannelMapIcarusAlg::SignalType(raw::ChannelID_t const channel) const
{
    unsigned int nChanPerCryo = fNchannels/fNcryostat;
    unsigned int cryostat = channel / nChanPerCryo;  
    unsigned int chan_in_cryo = channel % nChanPerCryo ;
    
    unsigned int nChanPerTPC = nChanPerCryo/fNTPC[0];
    // casting wil trunc towards 0 -- faster than floor
    unsigned int tpc = chan_in_cryo / nChanPerTPC;  

    //need number of planes to know Collection 
    unsigned int PlanesThisTPC = fNPlanes[0][tpc];
    
    
    
    SigType_t sigt = geo::kMysteryType;
    if(      (channel >= fFirstChannelInThisPlane[cryostat][tpc][0]) &&
             (channel <  fFirstChannelInNextPlane[cryostat][tpc][PlanesThisTPC-2])    ){ sigt = geo::kInduction; }
    else if( (channel >= fFirstChannelInThisPlane[cryostat][tpc][PlanesThisTPC-1]) &&
             (channel <  fFirstChannelInNextPlane[cryostat][tpc][PlanesThisTPC-1])    ){ sigt = geo::kCollection; }
    else
        mf::LogWarning("BadChannelSignalType") << "Channel " << channel
                                               << " not given signal type." << std::endl;
    
    return sigt;
}


//----------------------------------------------------------------------------
View_t ChannelMapIcarusAlg::View(raw::ChannelID_t const channel) const
{
  unsigned int nChanPerCryo = fNchannels/fNcryostat;
  unsigned int cryostat = channel / nChanPerCryo;  
  unsigned int chan_in_cryo = channel % nChanPerCryo ;
  
  unsigned int nChanPerTPC = nChanPerCryo/fNTPC[0];
  // casting wil trunc towards 0 -- faster than floor
  unsigned int tpc = chan_in_cryo / nChanPerTPC;  
  
  //need number of planes to know Collection 
  unsigned int PlanesThisTPC = fNPlanes[0][tpc];

  View_t view = geo::kUnknown;
  
  
  //THIS IS A TEMPORARY THING!
  //If 4 planes, the first two are U
  //If 3 planes, the first one is U.
  
  if(PlanesThisTPC==4){
    
    //first "two planes" are with horizontal wires: give them same view
    if(      (channel >= fFirstChannelInThisPlane[cryostat][tpc][0]) &&
	     (channel <  fFirstChannelInNextPlane[cryostat][tpc][0])    ){ view = geo::kU; }
    else if( (channel >= fFirstChannelInThisPlane[cryostat][tpc][1]) &&
	     (channel <  fFirstChannelInNextPlane[cryostat][tpc][1])    ){ view = geo::kU; }
    
    else if( (channel >= fFirstChannelInThisPlane[cryostat][tpc][2]) &&
	     (channel <  fFirstChannelInNextPlane[cryostat][tpc][2])    ){ view = geo::kV; }
    else if( (channel >= fFirstChannelInThisPlane[cryostat][tpc][3]) &&
	     (channel <  fFirstChannelInNextPlane[cryostat][tpc][3])    ){ view = geo::kW; }
    else
      mf::LogWarning("BadChannelSignalType") << "Channel " << channel
					     << " not given view type.";
    
      return view;
  }
  else if(PlanesThisTPC==3){
    if(      (channel >= fFirstChannelInThisPlane[cryostat][tpc][0]) &&
	     (channel <  fFirstChannelInNextPlane[cryostat][tpc][0])    ){ view = geo::kU; }
    else if( (channel >= fFirstChannelInThisPlane[cryostat][tpc][1]) &&
	     (channel <  fFirstChannelInNextPlane[cryostat][tpc][1])    ){ view = geo::kV; }
    else if( (channel >= fFirstChannelInThisPlane[cryostat][tpc][2]) &&
	     (channel <  fFirstChannelInNextPlane[cryostat][tpc][2])    ){ view = geo::kW; }
    else
      mf::LogWarning("BadChannelSignalType") << "Channel " << channel
					     << " not given view type.";
    
    return view;
  }
  else{
    mf::LogWarning("BadChannelSignalType") << "nPLanes is weird! " << PlanesThisTPC
					   << ".";
    
    return view;
  }
}
  
  //----------------------------------------------------------------------------
std::set<View_t> const& ChannelMapIcarusAlg::Views() const
{
    return fViews;
}

//----------------------------------------------------------------------------
std::set<PlaneID> const& ChannelMapIcarusAlg::PlaneIDs() const
{
    return fPlaneIDs;
}


//----------------------------------------------------------------------------
unsigned int ChannelMapIcarusAlg::NTPCsets
  (readout::CryostatID const& cryoid) const
{
    // return the same number as the number of TPCs
    return (cryoid.isValid && cryoid.Cryostat < fNTPC.size())?fNTPC[cryoid.Cryostat]: 0;
} // ChannelMapIcarusAlg::NTPCsets()


//----------------------------------------------------------------------------
unsigned int ChannelMapIcarusAlg::MaxTPCsets() const
{
    return MaxTPCs();
} // ChannelMapIcarusAlg::MaxTPCsets()


//----------------------------------------------------------------------------
bool ChannelMapIcarusAlg::HasTPCset(readout::TPCsetID const& tpcsetid) const
{
    return tpcsetid.TPCset < NTPCsets(tpcsetid);
} // ChannelMapIcarusAlg::HasTPCset()


//----------------------------------------------------------------------------
readout::TPCsetID ChannelMapIcarusAlg::TPCtoTPCset
  (geo::TPCID const& tpcid) const
{
    return ConvertTPCtoTPCset(tpcid);
} // ChannelMapIcarusAlg::TPCtoTPCset()


//----------------------------------------------------------------------------
std::vector<geo::TPCID> ChannelMapIcarusAlg::TPCsetToTPCs
  (readout::TPCsetID const& tpcsetid) const
{
    std::vector<geo::TPCID> IDs;
    if (tpcsetid.isValid) IDs.emplace_back(ConvertTPCsetToTPC(tpcsetid));
    return IDs;
} // ChannelMapIcarusAlg::TPCsetToTPCs()


//----------------------------------------------------------------------------
geo::TPCID ChannelMapIcarusAlg::FirstTPCinTPCset
  (readout::TPCsetID const& tpcsetid) const
{
    return ConvertTPCsetToTPC(tpcsetid);
} // ChannelMapIcarusAlg::FirstTPCinTPCset()


//----------------------------------------------------------------------------
unsigned int ChannelMapIcarusAlg::MaxTPCs() const
{
    unsigned int max = 0;
    for (unsigned int nTPCs: fNTPC) if (nTPCs > max) max = nTPCs;
    return max;
} // ChannelMapIcarusAlg::MaxTPCs()


//----------------------------------------------------------------------------
unsigned int ChannelMapIcarusAlg::NROPs
    (readout::TPCsetID const& tpcsetid) const
{
    if (!HasTPCset(tpcsetid)) return 0;
    return AccessElement(fNPlanes, FirstTPCinTPCset(tpcsetid));
} // ChannelMapIcarusAlg::NROPs()


//----------------------------------------------------------------------------
unsigned int ChannelMapIcarusAlg::MaxROPs() const {
    unsigned int max = 0;
    for (auto const& cryo_tpc: fNPlanes)
        for (unsigned int nPlanes: cryo_tpc)
            if (nPlanes > max) max = nPlanes;
    return max;
} // ChannelMapIcarusAlg::MaxROPs()


//----------------------------------------------------------------------------
bool ChannelMapIcarusAlg::HasROP(readout::ROPID const& ropid) const {
    return ropid.ROP < NROPs(ropid);
} // ChannelMapIcarusAlg::HasROP()


//----------------------------------------------------------------------------
readout::ROPID ChannelMapIcarusAlg::WirePlaneToROP
  (geo::PlaneID const& planeid) const
{
    return ConvertWirePlaneToROP(planeid);
} // ChannelMapIcarusAlg::WirePlaneToROP()


//----------------------------------------------------------------------------
std::vector<geo::PlaneID> ChannelMapIcarusAlg::ROPtoWirePlanes
  (readout::ROPID const& ropid) const
{
    std::vector<geo::PlaneID> IDs;
    if (ropid.isValid) IDs.emplace_back(FirstWirePlaneInROP(ropid));
    return IDs;
} // ChannelMapIcarusAlg::ROPtoWirePlanes()


//----------------------------------------------------------------------------
std::vector<geo::TPCID> ChannelMapIcarusAlg::ROPtoTPCs
  (readout::ROPID const& ropid) const
{
  std::vector<geo::TPCID> IDs;
  // we take the TPC set of the ROP and convert it straight into a TPC ID
  if (ropid.isValid) IDs.emplace_back(ConvertTPCsetToTPC(ropid.asTPCsetID()));
  return IDs;
} // ChannelMapIcarusAlg::ROPtoTPCs()


//----------------------------------------------------------------------------
readout::ROPID ChannelMapIcarusAlg::ChannelToROP
  (raw::ChannelID_t channel) const
{
    // which wires does the channel cover?
    std::vector<geo::WireID> wires = ChannelToWire(channel);
  
    // - none:
    if (wires.empty()) return {}; // default-constructed ID, invalid
  
    // - one: maps its plane ID into a ROP ID
    return WirePlaneToROP(wires[0]);
} // ChannelMapIcarusAlg::ROPtoTPCs()


//----------------------------------------------------------------------------
raw::ChannelID_t ChannelMapIcarusAlg::FirstChannelInROP
  (readout::ROPID const& ropid) const
{
    if (!ropid.isValid) return raw::InvalidChannelID;
    return (raw::ChannelID_t)
    AccessElement(fPlaneBaselines, ConvertROPtoWirePlane(ropid));
} // ChannelMapIcarusAlg::FirstChannelInROP()


//----------------------------------------------------------------------------
geo::PlaneID ChannelMapIcarusAlg::FirstWirePlaneInROP
  (readout::ROPID const& ropid) const
{
    return ConvertROPtoWirePlane(ropid);
} // ChannelMapIcarusAlg::FirstWirePlaneInROP()


//----------------------------------------------------------------------------
readout::TPCsetID ChannelMapIcarusAlg::ConvertTPCtoTPCset
  (geo::TPCID const& tpcid)
{
    if (!tpcid.isValid) return {}; // invalid ID, default-constructed
    return {
        (readout::CryostatID::CryostatID_t) tpcid.Cryostat,
        (readout::TPCsetID::TPCsetID_t) tpcid.TPC
    };
} // ChannelMapIcarusAlg::ConvertTPCtoTPCset()


//----------------------------------------------------------------------------
geo::TPCID ChannelMapIcarusAlg::ConvertTPCsetToTPC
  (readout::TPCsetID const& tpcsetid)
{
    if (!tpcsetid.isValid) return {};
    return {
        (geo::CryostatID::CryostatID_t) tpcsetid.Cryostat,
        (geo::TPCID::TPCID_t) tpcsetid.TPCset
    };
} // ChannelMapIcarusAlg::ConvertTPCsetToTPC()


//----------------------------------------------------------------------------
readout::ROPID ChannelMapIcarusAlg::ConvertWirePlaneToROP
  (geo::PlaneID const& planeid)
{
    if (!planeid.isValid) return {}; // invalid ID, default-constructed
    return {
        (readout::CryostatID::CryostatID_t) planeid.Cryostat,
        (readout::TPCsetID::TPCsetID_t) planeid.TPC,
        (readout::ROPID::ROPID_t) planeid.Plane
    };
  
} // ChannelMapIcarusAlg::ConvertWirePlaneToROP()


//----------------------------------------------------------------------------
geo::PlaneID ChannelMapIcarusAlg::ConvertROPtoWirePlane
  (readout::ROPID const& ropid)
{
    if (!ropid.isValid) return {};
    return {
        (geo::CryostatID::CryostatID_t) ropid.Cryostat,
        (geo::TPCID::TPCID_t) ropid.TPCset,
        (geo::PlaneID::PlaneID_t) ropid.ROP
    };
} // ChannelMapIcarusAlg::ConvertROPtoWirePlane()


//----------------------------------------------------------------------------
  
} // namespace
