/**
 * @file   icaruscode/Filters/FilterDirts_module.cc
 * @brief  Select only neutrino events with the primary vertex outside the active volume and that deposit energy in the active volume (dirt events)
 * @authors Christian Farnese (farnese@pd.infn.it),
 *
 * @date   Jan 10, 2024
 * 
 *
 */


// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardataalg/MCDumpers/MCDumperUtils.h" // sim::TruthXXXName()
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/ROOTGeometryNavigator.h"
#include "larcorealg/Geometry/GeoNodePath.h"
#include "larcorealg/Geometry/TransformationMatrix.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::...::makeFromCoords()
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/CoreUtils/StdUtils.h" // util::begin(), util::end()
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h" // geo::Point_t
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimChannel.h"

// nutools libraries
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ standard libraries
#include <regex>
#include <algorithm> // std::sort(), std::binary_search()
#include <vector>
#include <string>
#include <atomic>
#include <utility> // std::move()
#include <cmath> // std::abs()

#include "canvas/Utilities/InputTag.h"

// -----------------------------------------------------------------------------
namespace icarus::simfilter { class FilterDirts; }


/**
 * The filter selects "qualifying" events that:
 * 
 *  * have interaction vertex outside the active volume or outside specific volumes that can be indicated using the configuration parameters
 *  * deposit an energy > 0 ( the filter can be improved introducing a threshold that can be set using a configuration parameter )
 * 
 * The event is kept if there is _at least_ one qualifying interaction.
 * In that case, the whole event is passed (including any other interaction).
 * Please note: this class has been created starting from the `FilterNeutrinoActive` filter
 *
 *
 *
 * Input
 * ======
 * 
 * In the input file, there should be the truth MC information and the Geant4 propagation: so run both gen and g4 stages before this filter.
 * 
 * 
 * Configuration
 * ==============
 * 
 *
 * Filtering
 * ----------
 * 
 * This is a filter module: its configuration must appear in the `filters`
 * table and to be included in a path in `trigger_path`, and the output module
 * must be told to include that path in `SelectEvents` configuration key; e.g.:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * physics: {
 *   filters: {
 *     eventActive: @local::icarus_FilterDirts
 *     # ...
 *   } # filters
 *   # ...
 *   
 *   appliedFilters: [ eventActive ]
 *   
 *   trigger_paths: [ appliedFilters ]
 *   
 * } # physics
 * 
 * 
 * outputs: {
 *   rootoutput: {
 *     module_type:  RootOutput
 *     fileName:    "%ifb_%tc-%p.root"
 *     SelectEvents: [ appliedFilters ]
 *   }
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 */
class icarus::simfilter::FilterDirts: public art::EDFilter {
    
public:
    
    /// Configuration of box volume geometry.
    struct BoxCoordConfig {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        
        fhicl::Atom<double> Xmin {
            Name("Xmin"),
            Comment("minimum x coordinate of the volume [cm]")
        };
        fhicl::Atom<double> Ymin {
            Name("Ymin"),
            Comment("minimum y coordinate of the volume [cm]")
        };
        fhicl::Atom<double> Zmin {
            Name("Zmin"),
            Comment("minimum z coordinate of the volume [cm]")
        };
        fhicl::Atom<double> Xmax {
            Name("Xmax"),
            Comment("maximum x coordinate of the volume [cm]")
        };
        fhicl::Atom<double> Ymax {
            Name("Ymax"),
            Comment("maximum y coordinate of the volume [cm]")
        };
        fhicl::Atom<double> Zmax {
            Name("Zmax"),
            Comment("maximum z coordinate of the volume [cm]")
        };
        
    }; // struct BoxCoordConfig
    
    
    /// Configuration parameter structure.
    struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        
        fhicl::Atom<bool> inActive {
            Name("inActive"),
            Comment("rejects events with interactions in TPC active volumes"),
            false // default
        };

        fhicl::Atom<art::InputTag> TPCchannelTag {
            Name("TPCchannelTag"),
            Comment("data product with simulated TPC channel charge"),
            "larg4intime" // default
        };
        
        fhicl::OptionalSequence<fhicl::Table<BoxCoordConfig>> volumeBoxes {
            Name("volumeBoxes"),
            Comment("interactions in box volumes by world coordinates [cm]")
        };
        
        fhicl::Sequence<std::string> volumeNames {
            Name("volumeNames"),
            Comment("interactions in box volumes by name (std::regex pattern)"),
            std::vector<std::string>{} // default value
        };
        
        fhicl::Atom<std::string> logCategory {
            Name("logCategory"),
            Comment("category name for message facility message stream"),
            "FilterDirts" // default
        };
        
    }; // Config
    
    using Parameters = art::EDFilter::Table<Config>;
    
    
    /// Constructor: reads configuration and extracts information from geometry.
    explicit FilterDirts(Parameters const& config);
    
    /// Framework hook: applies the filter.
    virtual bool filter(art::Event& event) override;
    
    /// Framework hook: prints the summary of the passed events.
    virtual void endJob() override;
    
private:
    
    // --- BEGIN -- Configuration parameters -----------------------------------
    
    /// Volumes for qualifying interactions.
    std::vector<geo::BoxBoundedGeo> fVolumes;
    art::InputTag const fTPCchannelTag; ///< `sim::SimChannel` data product tag.   
    std::string const fLogCategory; ///< Category name for the output stream.
    
    // --- END -- Configuration parameters -------------------------------------
    
    
    // --- BEGIN -- Counters ---------------------------------------------------
    
    std::atomic<unsigned int> fNObserved { 0U }; ///< Number of observed events.
    std::atomic<unsigned int> fNPassed { 0U }; ///< Number of passed events.
    
    // --- END -- Counters -----------------------------------------------------
    
    /// Adds all active volumes of detector into the qualifying volume list.
    void addActiveVolumes();
    
    /// Adds the specified volumes into the qualifying volume list.
    void addVolumeBoxes
    (fhicl::OptionalSequence<fhicl::Table<BoxCoordConfig>> const& boxConfig);
    
    /// Adds all `volName` from geometry into the qualifying volume list.
    unsigned addVolumeByName(std::string const& volumeName);
    
    
    /// Returns whether the interaction described in `truth` qualifies.
    bool qualifying(simb::MCTruth const& truth) const;
    
    /// Returns whether the location is among the accepted ones.
    bool qualifyingLocation(geo::Point_t const& location) const;
    
    
    /// Returns a sorted copy of the specified collection.
    template <typename Coll>
    static Coll sorted(Coll const& coll);
    
    
}; // icarus::simfilter::FilterDirts



// -----------------------------------------------------------------------------
// ---  Implementation
// -----------------------------------------------------------------------------
icarus::simfilter::FilterDirts::FilterDirts
(Parameters const& config)
: art::EDFilter(config)
, fTPCchannelTag(config().TPCchannelTag())
, fLogCategory(config().logCategory())
{
    
    { // local scope
        mf::LogInfo log(fLogCategory);
        
       
        log << "\nConfiguration of qualifying volumes:";
        
    }
    
    //
    // load volumes from the different sources
    //
    if (config().inActive()) addActiveVolumes();
    
    addVolumeBoxes(config().volumeBoxes);
    
    for (std::string const& volName: config().volumeNames())
        addVolumeByName(volName);
    
    //
    // check that we are doing at least something
    //
    if (fVolumes.empty()) {
        
        throw art::Exception(art::errors::Configuration)
        << "No filtering action specified (volume).\n"
        ;
        
    } // if no filter
    
} // icarus::simfilter::FilterDirts::FilterDirts()


// -----------------------------------------------------------------------------
bool icarus::simfilter::FilterDirts::filter(art::Event& event) {
    
    ++fNObserved;
    
    /*
     * Consider all truth information available in the event.
     * Any record of any truth data product will be enough to pass the event.
     */
    //std::vector<art::Handle<std::vector<simb::MCTruth>>> allTruth;
    //event.getManyByType(allTruth);
    
    std::vector<sim::SimChannel> const& charge   = *(event.getValidHandle<std::vector<sim::SimChannel>>(fTPCchannelTag));
    
    float total_quenched_energy=0;

	for (sim::SimChannel const& chargechannel: charge) //loop on SimChannel
    	{	
        for (sim::TDCIDE const& tdcide: chargechannel.TDCIDEMap()) //loop on TDC
        {
            for (sim::IDE const& ide: tdcide.second) //loop on IDE
            {
                total_quenched_energy += ide.energy;
            }    //loop on IDE
        }     //loop on TDC
    	}//loop on SimChannel


    mf::LogDebug(fLogCategory) << "Total energy: " << total_quenched_energy;
   
    if(total_quenched_energy==0)return false;
    
    bool outside_volume = false;
    auto allTruth = event.getMany<std::vector<simb::MCTruth>>();
    
    if (allTruth.empty()) { // is this real data?
        throw art::Exception(art::errors::ProductNotFound)
        << event.id() << " has no truth information!\n";
    } // if no truth
    
    mf::LogDebug(fLogCategory)
    << "Event " << event.id() << " (#" << fNObserved << ") has "
    << allTruth.size() << " truth data products.";
    
    for (auto const& handle: allTruth) {
        
        art::InputTag const& tag [[maybe_unused]] = handle.provenance()->inputTag();
        
        std::vector<simb::MCTruth> const& truths = *handle;
        if (truths.empty()) {
            mf::LogTrace(fLogCategory)
            << "No truth records from " << tag.encode() << ": skipped.";
            continue;
        } // if no truth
        
        for (auto const& [ iTruth, truth ]: util::enumerate(truths)) {
            
            mf::LogTrace(fLogCategory)
            << "Processing record [" << (iTruth + 1U) << "/" << truths.size()
            << "] from " << tag.encode();
            
            if (!qualifying(truth))
		{
		outside_volume = true;
		break;
 		}           
        } // for truth record
        
    } // for truth data product
    
    //mf::LogTrace(fLogCategory) << "Event " << event.id() << " (#" << fNObserved
    //<< ")  does not pass the filter (" << fNPassed << "/" << fNObserved
    //<< " passed so far).";
    return outside_volume;
    
} // icarus::simfilter::FilterDirts::filter()


// -----------------------------------------------------------------------------
void icarus::simfilter::FilterDirts::endJob() {
    
    mf::LogInfo log(fLogCategory);
    log
    << "FilterDirts: passed " << fNPassed << " / " << fNObserved
    << " events";
    if (fNObserved > 0U)
        log << " (" << (float(fNPassed) * 100. / fNObserved) << "%)";
    
} // icarus::simfilter::FilterDirts::endJob()


// -----------------------------------------------------------------------------
void icarus::simfilter::FilterDirts::addActiveVolumes() {
    
    geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());
    
    for (geo::TPCGeo const& TPC: geom.Iterate<geo::TPCGeo>()) {
        
        geo::BoxBoundedGeo const& box = TPC.ActiveBoundingBox();
        
        mf::LogVerbatim(fLogCategory)
        << "[volume #" << fVolumes.size() << "] active volume from " << TPC.ID()
        << ": [ " << box.Min() << " -- " << box.Max() << " ]";
        
        fVolumes.push_back(box);
        
    } // for all TPCs
    
} // icarus::simfilter::FilterDirts::addActiveVolumes()


// -----------------------------------------------------------------------------
void icarus::simfilter::FilterDirts::addVolumeBoxes
(fhicl::OptionalSequence<fhicl::Table<BoxCoordConfig>> const& boxConfig)
{
    std::vector<BoxCoordConfig> boxParams;
    
    if (!boxConfig(boxParams)) return;
    
    for (auto const& [ iBox, boxParam ]: util::enumerate(boxParams)) {
        
        geo::BoxBoundedGeo box {
            boxParam.Xmin(), boxParam.Xmax(),
            boxParam.Ymin(), boxParam.Ymax(),
            boxParam.Zmin(), boxParam.Zmax()
        };
        
        mf::LogVerbatim(fLogCategory)
        << "[volume #" << fVolumes.size() << "] box coordinates #" << iBox
        << ": [ " << box.Min() << " -- " << box.Max() << " ]";
        
        fVolumes.push_back(std::move(box));
        
    } // for boxes
    
} // icarus::simfilter::FilterDirts::addVolumeBoxes()


// -----------------------------------------------------------------------------
unsigned int icarus::simfilter::FilterDirts::addVolumeByName
(std::string const& volumePattern)
{
    
    geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());
    
    //
    // find the path of all volumes matching the given pattern
    //
    std::regex const namePattern { volumePattern };
    std::vector<geo::GeoNodePath> volumePaths;
    auto const matchMe = [&pattern=namePattern](std::string const& s)
    { std::smatch match; return (std::regex_match(s, match, pattern)); };
    auto const findVolume = [&volumePaths, &patternMatcher=matchMe](auto& path)
    {
        if (patternMatcher(path.current().GetVolume()->GetName()))
            volumePaths.push_back(path);
        return true;
    };
    
    geo::ROOTGeometryNavigator navigator { *(geom.ROOTGeoManager()) };
    
    navigator.apply(findVolume);
    
    if (volumePaths.empty()) {
        throw art::Exception(art::errors::Configuration)
        << "No volume matching '" << volumePattern
        << "' has been found in the detector '" << geom.DetectorName()
        << "'.\n";
    }
    
    //
    // convert each volume into world coordinates and add it to the list
    //
    for (auto const& [ iVolume, path ] : util::enumerate(volumePaths)) {
        
        //
        // find the coordinates of the volume in local coordinates
        //
        TGeoShape const* pShape = path.current().GetVolume()->GetShape();
        auto pBox = dynamic_cast<TGeoBBox const*>(pShape);
        if (!pBox) {
            throw cet::exception("FilterDirts") << "Volume '"
            << path.current().GetName() << "' is a " << pShape->IsA()->GetName()
            << ", not a TGeoBBox.\n";
        }
        
        geo::Point_t const origin
        = geo::vect::makeFromCoords<geo::Point_t>(pBox->GetOrigin());
        geo::Vector_t const diag = {
            std::abs(pBox->GetDX()), std::abs(pBox->GetDY()), std::abs(pBox->GetDZ())
        };
        
        //
        // convert to world coordinates
        //
        
        auto const trans
        = path.currentTransformation<geo::TransformationMatrix>();
        
        geo::Point_t min, max;
        trans.Transform(origin - diag, min);
        trans.Transform(origin + diag, max);
        
        //
        // add to the coordinates
        //
        geo::BoxBoundedGeo box { min, max };
        
        mf::LogVerbatim(fLogCategory)
        << " c* [volume #" << fVolumes.size() << "] volume box '"
        << path.current().GetVolume()->GetName()
        << "' [(" << (iVolume + 1U) << "/" << volumePaths.size()
        << "): [ " << box.Min() << " -- " << box.Max() << " ]";
        
        fVolumes.push_back(std::move(box));
        
    } // for all volume paths
    
    return volumePaths.size();
} // icarus::simfilter::FilterDirts::addVolumeByName()


// -----------------------------------------------------------------------------
bool icarus::simfilter::FilterDirts::qualifying
(simb::MCTruth const& truth) const
{
    /*
     * Apply all the needed cuts:
     * * interaction type
     * * current type
     * * location
     *
     */
    
    //
    // only neutrino record types may qualify:
    //
    if (!truth.NeutrinoSet()) {
        mf::LogTrace(fLogCategory)
        << "Interaction does not qualify because it is not tagged as neutrino.";
        return false;
    }
    
    simb::MCNeutrino const& nuInfo = truth.GetNeutrino();
    simb::MCParticle const& nu = nuInfo.Nu();
    
    //
    //
    // location
    //
    if (!fVolumes.empty() && !qualifyingLocation({ nu.Vx(), nu.Vy(), nu.Vz() }))
        return false;
    
    // success, after all
    return true;
    
} // icarus::simfilter::FilterDirts::qualifying()


// -----------------------------------------------------------------------------
bool icarus::simfilter::FilterDirts::qualifyingLocation
(geo::Point_t const& location) const
{
    
    mf::LogTrace log(fLogCategory);
    log << "Interaction location: " << location << " cm";
    
    for (auto const& [ iBox, box ]: util::enumerate(fVolumes)) {
        if (!box.ContainsPosition(location)) continue;
        
        log << " => in volume #" << iBox
        << " [ " << box.Min() << " -- " << box.Max() << " ] => :-)";
        
        return true;
    } // for
    
    log << " => :-(";
    return false;
    
} // icarus::simfilter::FilterDirts::qualifyingLocation()


// -----------------------------------------------------------------------------
template <typename Coll>
Coll icarus::simfilter::FilterDirts::sorted(Coll const& coll) {
    
    // copy, sort, return
    auto sortedColl { coll };
    std::sort(util::begin(sortedColl), util::end(sortedColl));
    return sortedColl;
    
} // icarus::simfilter::FilterDirts::sorted()



// -----------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::simfilter::FilterDirts)


// -----------------------------------------------------------------------------
