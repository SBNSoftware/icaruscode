////////////////////////////////////////////////////////////////////////
//
// ROIFinder class - An ROI finding module for complete deconvolved waveforms
//
// usher@slac.stanford.edu
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <utility> // std::pair<>
#include <memory> // std::unique_ptr<>
#include <iomanip>
#include <fstream>
#include <random>

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h" 
#include "canvas/Utilities/Exception.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "icaruscode/TPC/SignalProcessing/RecoWire/ROITools/IROILocator.h"

#include "icarus_signal_processing/WaveformTools.h"
#include "icarus_signal_processing/Denoising.h"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_arena.h"
#include "tbb/spin_mutex.h"
#include "tbb/concurrent_hash_map.h"


namespace {

  /// Helper: lazily returns the expanded content of a set of `recob::Wires`.
struct PlaneWireData 
{
    std::size_t size() const { return wires.size(); }
    void resize(std::size_t nWires)
    { wires.clear(); wires.resize(nWires, nullptr); }
    void addWire(std::size_t iWire, recob::Wire const& wire)
    { wires.at(iWire) = &wire; }
    icarus_signal_processing::ArrayFloat operator() () const
    {
        icarus_signal_processing::ArrayFloat data;
        data.resize(wires.size());
        for (auto [ iWire, wire ]: util::enumerate(wires))
        {
          if (wire) data[iWire] = wire->Signal();
          else      data[iWire] = std::vector<float>(4096,0.);
        }
      return data;
}
private:
    std::vector<recob::Wire const*> wires;
}; // PlaneWireData
  
} // local namespace


///creation of calibrated signals on wires
namespace caldata {

    tbb::spin_mutex roifinderSpinMutex;

class ROIFinder : public art::EDProducer
{
public:
// create calibrated signals on wires. this class runs 
// an fft to remove the electronics shaping.     
    explicit ROIFinder(fhicl::ParameterSet const& pset);
    virtual ~ROIFinder();
    void     produce(art::Event& evt); 
    void     beginJob(); 
    void     endJob();                 
    void     reconfigure(fhicl::ParameterSet const& p);
private:
    using PlaneIDToDataPair    = std::pair<std::vector<raw::ChannelID_t>,PlaneWireData>;
    using PlaneIDToDataPairMap = std::map<geo::PlaneID,PlaneIDToDataPair>;
    using PlaneIDVec           = std::vector<geo::PlaneID>;

    // Define a class to handle processing for individual threads
    class multiThreadDeconvolutionProcessing 
    {
    public:
        multiThreadDeconvolutionProcessing(ROIFinder const&            parent,
                                           art::Event&                 event,
                                           const PlaneIDVec&           planeIDVec,
                                           const PlaneIDToDataPairMap& planeIDToDataPairMap, 
                                           std::vector<recob::Wire>&   wireColVec,
                                           std::vector<recob::Wire>&   morphedVec)
            : fROIFinder(parent),
              fEvent(event),
              fPlaneIDVec(planeIDVec),
              fPlaneIDToDataPairMap(planeIDToDataPairMap),
              fWireColVec(wireColVec),
              fMorphedVec(morphedVec)
        {}
        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t idx = range.begin(); idx < range.end(); idx++)
                fROIFinder.processPlane(idx, fEvent, fPlaneIDVec, fPlaneIDToDataPairMap, fWireColVec, fMorphedVec);
        }
    private:
        const ROIFinder&            fROIFinder;
        art::Event&                 fEvent;
        const PlaneIDVec&           fPlaneIDVec;
        const PlaneIDToDataPairMap& fPlaneIDToDataPairMap;
        std::vector<recob::Wire>&   fWireColVec;
        std::vector<recob::Wire>&   fMorphedVec;
    };

    // Function to do the work
    void  processPlane(size_t,
                       art::Event&,
                       const PlaneIDVec&,
                       const PlaneIDToDataPairMap&, 
                       std::vector<recob::Wire>&,
                       std::vector<recob::Wire>&) const;

    // This is for the baseline...
    float getMedian(const icarus_signal_processing::VectorFloat, const unsigned int) const;

    std::vector<art::InputTag>                                 fWireModuleLabelVec;         ///< vector of modules that made digits
    bool                                                       fOutputMorphed;              ///< Output the morphed waveforms
    bool                                                       fDiagnosticOutput;           ///< secret diagnostics flag
    size_t                                                     fEventCount;                 ///< count of event processed
    
    std::map<size_t,std::unique_ptr<icarus_tool::IROILocator>> fROIToolMap;

    const geo::GeometryCore*                                   fGeometry = lar::providerFrom<geo::Geometry>();
    
}; // class ROIFinder

DEFINE_ART_MODULE(ROIFinder)

//-------------------------------------------------
ROIFinder::ROIFinder(fhicl::ParameterSet const& pset) : EDProducer{pset}
{
    this->reconfigure(pset);

    for(const auto& wireLabel : fWireModuleLabelVec)
    {
        produces< std::vector<recob::Wire>>(wireLabel.instance());

        if (fOutputMorphed) produces<std::vector<recob::Wire>>(wireLabel.instance() + "M");
    }
}

//-------------------------------------------------
ROIFinder::~ROIFinder()
{
}

//////////////////////////////////////////////////////
void ROIFinder::reconfigure(fhicl::ParameterSet const& pset)
{
    // Recover the parameters
    fWireModuleLabelVec  = pset.get<std::vector<art::InputTag>>("WireModuleLabelVec",  std::vector<art::InputTag>()={"decon1droi"});
    fOutputMorphed       = pset.get< bool                     >("OutputMorphed",                                              true);
    fDiagnosticOutput    = pset.get< bool                     >("DaignosticOutput",                                          false);
    
    // Recover the list of ROI finding tools
    const fhicl::ParameterSet& roiFinderTools = pset.get<fhicl::ParameterSet>("ROIFinderToolVec");

    // Make a mapping between plane id and a plane's ROI finder
    // This allows different ROI finders per plane but, more important, different parameters
    for(const std::string& roiFinderTool : roiFinderTools.get_pset_names())
    {
        const fhicl::ParameterSet& roiFinderToolParamSet = roiFinderTools.get<fhicl::ParameterSet>(roiFinderTool);
        size_t                     planeIdx              = roiFinderToolParamSet.get<size_t>("Plane");
        
        fROIToolMap[planeIdx] = art::make_tool<icarus_tool::IROILocator> (roiFinderToolParamSet);
    }
    
    return;
}

//-------------------------------------------------
void ROIFinder::beginJob()
{
    fEventCount = 0;
} // beginJob

//////////////////////////////////////////////////////
void ROIFinder::endJob()
{
}

//////////////////////////////////////////////////////
void ROIFinder::produce(art::Event& evt)
{
    // We need to loop through the list of Wire data we have been given
    for(const auto& wireLabel : fWireModuleLabelVec)
    {
        // make a collection of Wires
        std::unique_ptr<std::vector<recob::Wire>> wireCol(new std::vector<recob::Wire>);
    
        std::unique_ptr<std::vector<recob::Wire>> morphedCol(new std::vector<recob::Wire>);
    
        std::cout << "ROIFinder, looking for decon1droi data at " << wireLabel << std::endl;

        mf::LogInfo("ROIFinder") << "ROIFinder, looking for decon1droi data at " << wireLabel << std::endl;
    
        // Read in the collection of full length deconvolved waveforms
        // Note we assume this list is sorted in increasing channel number!
        art::Handle< std::vector<recob::Wire>> wireVecHandle;
        
        evt.getByLabel(wireLabel, wireVecHandle);
    
        std::cout << "Recovered Wire data, size: " << wireVecHandle->size() << std::endl;

        mf::LogInfo("ROIFinder") << "--> Recovered wire data, size: " << wireVecHandle->size() << std::endl;
    
        if (!wireVecHandle->size())
        {
            evt.put(std::move(wireCol), wireLabel.instance());
            fEventCount++;
            
            return;
        }
    
        // The first step is to break up into groups by logical TPC/plane in order to do the parallel loop
        PlaneIDToDataPairMap planeIDToDataPairMap;
        PlaneIDVec           planeIDVec;
        size_t               numChannels(0);
    
        for(const auto& wire : *wireVecHandle)
        {
            raw::ChannelID_t channel = wire.Channel();
            
            std::vector<geo::WireID> wireIDVec = fGeometry->ChannelToWire(channel);

            if (wireIDVec.empty()) continue;
    
            for(const auto& wireID : wireIDVec)
            {
                const geo::PlaneID& planeID = wireID.planeID();
    
                PlaneIDToDataPairMap::iterator mapItr = planeIDToDataPairMap.find(planeID);
    
                // Make sure the array is initialized
                if (mapItr == planeIDToDataPairMap.end())
                {
                    unsigned int nWires = fGeometry->Nwires(planeID);
    
                    std::pair<PlaneIDToDataPairMap::iterator,bool> mapInsert = planeIDToDataPairMap.insert({planeID,PlaneIDToDataPair()});
    
                    if (!mapInsert.second) std::cout << "Failed to insert, is this possible?" << std::endl;
    
                    mapItr = mapInsert.first;
    
                    mapItr->second.first.resize(nWires);
                    mapItr->second.second.resize(nWires);
    
                    planeIDVec.emplace_back(planeID);
                }
    
                // Add waveform to the 2D array
                mapItr->second.first[wireID.Wire] = channel;
                mapItr->second.second.addWire(wireID.Wire, wire);
                numChannels++;
            }
        }

        // We might need this... it allows a temporary wire object to prevent crashes when some data is missing
        std::vector<recob::Wire> tempWireVec(numChannels/4);

        // Check integrity of map
        for(auto& mapInfo : planeIDToDataPairMap)
        {
            const std::vector<raw::ChannelID_t>&        channelVec = mapInfo.second.first;
            const icarus_signal_processing::ArrayFloat& dataArray  = mapInfo.second.second();

            for(size_t idx = 0; idx < channelVec.size(); idx++)
            {
                if (dataArray[idx].size() < 100) 
                {
                    mf::LogInfo("ROIFinder") << "  **> Found truncated wire, size: " << dataArray[idx].size() << ", channel: " << channelVec[idx] << std::endl;

                    std::vector<float>               zeroVec(4096,0.);
                    recob::Wire::RegionsOfInterest_t ROIVec;

                    ROIVec.add_range(0, std::move(zeroVec));

                    std::vector<geo::WireID> wireIDVec = fGeometry->ChannelToWire(channelVec[idx]);

                    // Given channel a large number so we know to not save
                    mapInfo.second.first[idx] = 100000 + idx;

                    tempWireVec.emplace_back(recob::WireCreator(std::move(ROIVec),idx,fGeometry->View(wireIDVec[0].planeID())).move());

                    mapInfo.second.second.addWire(idx,tempWireVec.back());
               }
            }
        }
   
        // Reserve the room for the output
        wireCol->reserve(wireVecHandle->size());
    
        // ... Launch multiple threads with TBB to do the deconvolution and find ROIs in parallel
        multiThreadDeconvolutionProcessing deconvolutionProcessing(*this, evt, planeIDVec, planeIDToDataPairMap, *wireCol, *morphedCol);
    
        tbb::parallel_for(tbb::blocked_range<size_t>(0, planeIDVec.size()), deconvolutionProcessing);
        
        // Time to stroe everything
        if(wireCol->size() == 0) mf::LogWarning("ROIFinder") << "No wires made for this event.";

        evt.put(std::move(wireCol), wireLabel.instance());

        if (fOutputMorphed) evt.put(std::move(morphedCol), wireLabel.instance()+"M");
    }

  fEventCount++;

  return;
} // produce

void  ROIFinder::processPlane(size_t                      idx,
                              art::Event&                 event,
                              const PlaneIDVec&           planeIDVec,
                              const PlaneIDToDataPairMap& planeIDToDataPairMap, 
                              std::vector<recob::Wire>&   wireColVec,
                              std::vector<recob::Wire>&   morphedVec) const
{
    // Recover the planeID for this thread
    const geo::PlaneID& planeID = planeIDVec[idx];

    // And the data array 
    PlaneIDToDataPairMap::const_iterator  mapItr = planeIDToDataPairMap.find(planeID);

    if (mapItr == planeIDToDataPairMap.end())
    {
        std::cout << "We know this cannot happen" << std::endl;
        return;
    }

    const PlaneIDToDataPair& planeIDToDataPair = mapItr->second;

    const icarus_signal_processing::ArrayFloat& dataArray  = planeIDToDataPair.second();
    const std::vector<raw::ChannelID_t>&        channelVec = planeIDToDataPair.first;

    // Keep track of our selected values
    icarus_signal_processing::ArrayFloat outputArray(dataArray.size(),icarus_signal_processing::VectorFloat(dataArray[0].size(),0.));
    icarus_signal_processing::ArrayBool  selectedVals(dataArray.size(),icarus_signal_processing::VectorBool(dataArray[0].size(),false));

    fROIToolMap.at(planeID.Plane)->FindROIs(event, dataArray, mapItr->first, outputArray, selectedVals);

//    icarus_signal_processing::ArrayFloat morphedWaveforms(dataArray.size());

    // Copy the "morphed" array
    if (fOutputMorphed)
    {
        for(size_t waveIdx = 0; waveIdx < outputArray.size(); waveIdx++)
        {
            // skip if a bad channbel
            if (channelVec[idx] >= 100000) continue;

            // First get a lock to make sure we don't conflict
            tbb::spin_mutex::scoped_lock lock(roifinderSpinMutex);

            recob::Wire::RegionsOfInterest_t ROIVec;

            ROIVec.add_range(0, std::move(outputArray[waveIdx]));

            raw::ChannelID_t channel = planeIDToDataPair.first[waveIdx];
            geo::View_t      view    = fGeometry->View(channel);

            morphedVec.push_back(recob::WireCreator(std::move(ROIVec),channel,view).move());
        }
    }

    // Ok, now go through the refined selected values array and find ROIs
    // Define the ROI and its container
    using CandidateROI    = std::pair<size_t, size_t>;
    using CandidateROIVec = std::vector<CandidateROI>;

    size_t leadTrail(1);

    for(size_t waveIdx = 0; waveIdx < selectedVals.size(); waveIdx++)
    {
        // Skip if a bad channel
        if (channelVec[waveIdx] >= 100000) continue;

        // Set up an object... 
        CandidateROIVec candidateROIVec;

        // Search for ROIs in current waveform
        icarus_signal_processing::VectorBool& selVals = selectedVals[waveIdx];

        for(size_t idx = 0; idx < selVals.size(); idx++)
        {
            if (selVals[idx])
            {
                Size_t startTick = idx >= leadTrail ? idx - leadTrail : 0;

                while(idx < selVals.size() && selVals[idx]) idx++;

                size_t stopTick  = idx < selVals.size() - leadTrail ? idx + leadTrail : selVals.size();

                candidateROIVec.emplace_back(startTick, stopTick);
            }
        }

        // merge overlapping (or touching) ROI's
        if(candidateROIVec.size() > 1)
        {
            // temporary vector for merged ROIs
            CandidateROIVec tempRoiVec;

            // Loop through candidate roi's
            size_t startRoi = candidateROIVec.front().first;
            size_t stopRoi  = startRoi;

            for(auto& roi : candidateROIVec)
            {
                // Should we merge roi's?
                if (roi.first <= stopRoi + 50)
                { 
                    // Make sure the merge gets the right start/end times
                    startRoi = std::min(startRoi,roi.first);
                    stopRoi  = std::max(stopRoi,roi.second);
                }
                else
                {
                    tempRoiVec.emplace_back(startRoi,stopRoi);

                    startRoi = roi.first;
                    stopRoi  = roi.second;
                }
            }

            // Make sure to get the last one
            tempRoiVec.push_back(CandidateROI(startRoi,stopRoi));

            candidateROIVec = tempRoiVec;
        }

        // vector that will be moved into the Wire object
        recob::Wire::RegionsOfInterest_t ROIVec;

        const icarus_signal_processing::VectorFloat& waveform = dataArray[waveIdx];
        icarus_signal_processing::VectorFloat holder;

        // We need to copy the deconvolved (and corrected) waveform ROI's
        for(const auto& candROI : candidateROIVec)
        {
            // First up: copy out the relevent ADC bins into the ROI holder
            size_t roiLen = candROI.second - candROI.first;

            holder.resize(roiLen);

            std::copy(waveform.begin()+candROI.first, waveform.begin()+candROI.second, holder.begin());

            // Now we do the baseline determination and correct the ROI
            // For now we are going to reset to the minimum element
            // Get slope/offset from first to last ticks
            if (holder.size() < 40)
            {
                float dADC   = (holder.back() - holder.front()) / float(holder.size());
                float offset = holder.front();

                for(auto& adcVal : holder)
                {
                    adcVal -= offset;
                    offset += dADC;
                }
            }

            // add the range into ROIVec
            ROIVec.add_range(candROI.first, std::move(holder));
        }

        // Check for emptiness
        if (!ROIVec.empty())
        {
            // First get a lock to make sure we don't conflict
            tbb::spin_mutex::scoped_lock lock(roifinderSpinMutex);

            raw::ChannelID_t channel = planeIDToDataPair.first[waveIdx];
            geo::View_t      view    = fGeometry->View(channel);

            // Since we process logical TPC images we need to watch for duplicating entries 
            // We can do that by checking to see if a channel has already been added...
            std::vector<geo::WireID> wireIDVec = fGeometry->ChannelToWire(channel);

            if (wireIDVec.size() > 1)
            {
                std::vector<recob::Wire>::iterator wireItr = std::find_if(wireColVec.begin(),wireColVec.end(),[channel](const auto& wire){return wire.Channel() == channel;});

                if (wireItr != wireColVec.end())
                {
                    if (fDiagnosticOutput) std::cout << "******************* Found duplicate entry for channel " << channel << " ************************" << std::endl;
                    continue;
                }
            }

            wireColVec.push_back(recob::WireCreator(std::move(ROIVec),channel,view).move());
        }
    }

    return;
}

float ROIFinder::getMedian(icarus_signal_processing::VectorFloat vals, const unsigned int nVals) const
{
    float median(0.);

    if (nVals > 2) 
    {
        if (nVals % 2 == 0) 
        {
            const auto m1 = vals.begin() + nVals / 2 - 1;
            const auto m2 = vals.begin() + nVals / 2;
            std::nth_element(vals.begin(), m1, vals.begin() + nVals);
            const auto e1 = *m1;
            std::nth_element(vals.begin(), m2, vals.begin() + nVals);
            const auto e2 = *m2;
            median = (e1 + e2) / 2.0;
        } 
        else 
        {
            const auto m = vals.begin() + nVals / 2;
            std::nth_element(vals.begin(), m, vals.begin() + nVals);
            median = *m;
        }
    }

    return median;
}

} // end namespace caldata
