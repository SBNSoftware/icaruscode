////////////////////////////////////////////////////////////////////////
// Class:       SimInfoOverlayFilter
// Plugin Type: filter (art v2_11_03)
// File:        SimInfoOverlayFilter_module.cc
//
// Generated at Sat Nov 24 16:09:21 2018 by Wesley Ketchum using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Provenance/ProductID.h"
#include "canvas/Persistency/Provenance/canonicalProductName.h"
#include "canvas/Utilities/FriendlyName.h"
#include "art/Persistency/Common/PtrMaker.h"

#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <unordered_map>

//gallery...
#include "gallery/Event.h"
#include "gallery/Handle.h"
#include "gallery/ValidHandle.h"

//association includes
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Assns.h"

//generator products
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/BeamGateInfo.h"

//G4 products
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "lardataobj/MCBase/MCParticleLite.h"

//pot product
#include "larcoreobj/SummaryData/POTSummary.h"


namespace mix {
  class SimInfoOverlayFilter;
}


class mix::SimInfoOverlayFilter : public art::EDFilter {
public:
  explicit SimInfoOverlayFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimInfoOverlayFilter(SimInfoOverlayFilter const &) = delete;
  SimInfoOverlayFilter(SimInfoOverlayFilter &&) = delete;
  SimInfoOverlayFilter & operator = (SimInfoOverlayFilter const &) = delete;
  SimInfoOverlayFilter & operator = (SimInfoOverlayFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  bool beginRun(art::Run & r) override;
  bool beginSubRun(art::SubRun & sr) override;
  void endJob() override;
  bool endRun(art::Run & r) override;
  bool endSubRun(art::SubRun & sr) override;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const  &fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;

private:

  //Infrastructure to read the event...
  std::vector<std::string> fSimInputFileNames;
  gallery::Event gEvent;

  //input tag lists...
  //note, if these are empty, that's ok: then we won't produce them on the event

  std::vector<art::InputTag>    fMCTruthInputModuleLabels;
  std::vector<art::InputTag>    fMCFluxInputModuleLabels;
  std::vector<art::InputTag>    fGTruthInputModuleLabels;
  std::vector<art::InputTag>    fBeamGateInputModuleLabels;

  std::vector<art::InputTag>    fMCParticleInputModuleLabels;
  std::vector<art::InputTag>    fSimEnergyDepositInputModuleLabels;
  std::vector<art::InputTag>    fAuxDetSimChannelInputModuleLabels;
  std::vector<art::InputTag>    fMCParticleLiteInputModuleLabels;
  std::vector<art::InputTag>    fSimChannelInputModuleLabels;
  std::vector<art::InputTag>    fSimPhotonsInputModuleLabels;

  //now for the associations...
  std::vector<art::InputTag>    fMCTruthMCFluxAssnsInputModuleLabels;
  std::vector<art::InputTag>    fMCTruthGTruthAssnsInputModuleLabels;
  std::vector<art::InputTag>    fMCTruthMCParticleAssnsInputModuleLabels;

  // Added by C Thorpe: New fhicl parameters needed to handle separate gen/g4 stages in auxiallary file
  bool fAuxHasSeparateGenG4;
  std::vector<art::InputTag>    fMCTruthMCParticleAssnsMCTruthLookupLabel;

  template<class T>
  using CollectionMap = std::unordered_map< std::string, std::unique_ptr<T> >;

  CollectionMap< std::vector<simb::MCTruth> >     fMCTruthMap;
  CollectionMap< std::vector<simb::MCFlux> >      fMCFluxMap;
  CollectionMap< std::vector<simb::GTruth> >      fGTruthMap;
  CollectionMap< std::vector<sim::BeamGateInfo> > fBeamGateInfoMap;

  CollectionMap< art::Assns<simb::MCTruth,simb::MCFlux,void> > fMCTruthMCFluxAssnsMap;
  CollectionMap< art::Assns<simb::MCTruth,simb::GTruth,void> > fMCTruthGTruthAssnsMap;

  CollectionMap< std::vector<simb::MCParticle> >      fMCParticleMap;

  CollectionMap< std::vector<sim::SimEnergyDeposit> > fSimEnergyDepositMap;
  CollectionMap< std::vector<sim::AuxDetSimChannel> > fAuxDetSimChannelMap;
  CollectionMap< std::vector<sim::MCParticleLite> >   fMCParticleLiteMap;//Ivan
  CollectionMap< std::vector<sim::SimChannel> >       fSimChannelMap;
  CollectionMap< std::vector<sim::SimPhotons> >       fSimPhotonsMap;

  CollectionMap< art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo> > fMCTruthMCParticleAssnsMap;
  //verbose
  int fVerbosity;

  //pot accounting
  bool fFillPOTInfo;
  art::InputTag fPOTSummaryTag;
  double fPOTSum_totpot;
  double fPOTSum_totgoodpot;
  double fPOTSum_totspills;
  double fPOTSum_goodspills;

  void MakePOTMap();
  std::map<unsigned int,double> fSR_POTPerEvent;
  std::map<unsigned int,double> fSR_GoodPOTPerEvent;
  std::map<unsigned int,double> fSR_SpillsPerEvent;
  std::map<unsigned int,double> fSR_GoodSpillsPerEvent;

  void FillInputModuleLabels(fhicl::ParameterSet const &,
			     std::string,
			     std::vector<art::InputTag> &);
  void FillInputModuleLabels(fhicl::ParameterSet const &);

  template<class T>
  void DeclareProduces(std::vector<art::InputTag> const&);
  void DeclareProduces();

  template<class T>
  void InitializeCollectionMap(std::vector<art::InputTag> const&,
			       CollectionMap<T> &);
  void InitializeCollectionMaps();

  template<class T>
  void PutCollectionOntoEvent(art::Event &,
			      CollectionMap<T> &);
  void PutCollectionsOntoEvent(art::Event &);

  template<class T>
  std::map< std::pair<art::ProductID,size_t> , art::Ptr<T> >
  FillCollectionMap(std::vector<art::InputTag> const&,
		    CollectionMap< std::vector<T> > &,
		    std::string,
		    art::Event &);

  template< class T, class U >
  void FillAssnsCollectionMap(std::vector<art::InputTag> const&,
			      CollectionMap< art::Assns<T,U,void> > &,
			      std::map< std::pair<art::ProductID,size_t> , art::Ptr<T> > const&,
			      std::map< std::pair<art::ProductID,size_t> , art::Ptr<U> > const&);
  template< class T, class U, class D>
  void FillAssnsCollectionMap(std::vector<art::InputTag> const&,
			      CollectionMap< art::Assns<T,U,D> > &,
			      std::map< std::pair<art::ProductID,size_t> , art::Ptr<T> > const&,
			      std::map< std::pair<art::ProductID,size_t> , art::Ptr<U> > const&);
};



mix::SimInfoOverlayFilter::SimInfoOverlayFilter(fhicl::ParameterSet const & p)
  : EDFilter{p},
  fSimInputFileNames(p.get<std::vector<std::string>>("SimInputFileNames")),
  gEvent(fSimInputFileNames),
  fAuxHasSeparateGenG4(p.get<bool>("AuxHasSeparateGenG4",false)),
  fVerbosity(p.get<int>("Verbosity",-1))
{
  FillInputModuleLabels(p);
  DeclareProduces();

  fFillPOTInfo = p.get<bool>("FillPOTInfo",true);
  if(fFillPOTInfo){
    fPOTSummaryTag = p.get<art::InputTag>("POTSummaryTag");
    fPOTSum_totpot = 0.0;
    fPOTSum_totgoodpot = 0.0;
    fPOTSum_totspills = 0.0;
    fPOTSum_goodspills = 0.0;

    this->produces<sumdata::POTSummary,art::InSubRun>();
    MakePOTMap();
  }

}

void mix::SimInfoOverlayFilter::FillInputModuleLabels(fhicl::ParameterSet const & p,
						      std::string s,
						      std::vector<art::InputTag> & labels)
{
  labels = p.get< std::vector<art::InputTag> >(s,std::vector<art::InputTag>());
}

void mix::SimInfoOverlayFilter::FillInputModuleLabels(fhicl::ParameterSet const & p)
{
  FillInputModuleLabels(p,"MCTruthInputModuleLabels",fMCTruthInputModuleLabels);
  FillInputModuleLabels(p,"MCFluxInputModuleLabels",fMCFluxInputModuleLabels);
  FillInputModuleLabels(p,"GTruthInputModuleLabels",fGTruthInputModuleLabels);
  FillInputModuleLabels(p,"BeamGateInputModuleLabels",fBeamGateInputModuleLabels);

  FillInputModuleLabels(p,"MCParticleInputModuleLabels",fMCParticleInputModuleLabels);
  FillInputModuleLabels(p,"SimEnergyDepositInputModuleLabels",fSimEnergyDepositInputModuleLabels);
  FillInputModuleLabels(p,"AuxDetSimChannelInputModuleLabels",fAuxDetSimChannelInputModuleLabels);
  FillInputModuleLabels(p,"MCParticleLiteInputModuleLabels",fMCParticleLiteInputModuleLabels); //Ivan
  FillInputModuleLabels(p,"SimChannelInputModuleLabels",fSimChannelInputModuleLabels);
  FillInputModuleLabels(p,"SimPhotonsInputModuleLabels",fSimPhotonsInputModuleLabels);

  FillInputModuleLabels(p,"MCTruthMCFluxAssnsInputModuleLabels",fMCTruthMCFluxAssnsInputModuleLabels);
  FillInputModuleLabels(p,"MCTruthGTruthAssnsInputModuleLabels",fMCTruthGTruthAssnsInputModuleLabels);
  FillInputModuleLabels(p,"MCTruthMCParticleAssnsInputModuleLabels",fMCTruthMCParticleAssnsInputModuleLabels);


  // New fhicl parameter needed to handle separate gen/g4 passes in the aux file, crashes when trying to construct
  // assns otherwise

  if(fAuxHasSeparateGenG4){
    FillInputModuleLabels(p,"MCTruthMCParticleAssnsMCTruthLookupLabel",fMCTruthMCParticleAssnsMCTruthLookupLabel);
    if(!fMCTruthMCParticleAssnsMCTruthLookupLabel.size())
      throw cet::exception("SimInfoOverlayFilter")
        << "Expecting aux file produced with separate gen/g4 passes, MCTruthMCParticleAssnsMCTruthLookupLabel in config needs to point to gen info in aux file";
  }
}

template<class T>
void mix::SimInfoOverlayFilter::DeclareProduces(std::vector<art::InputTag> const& labels)
{
  for(auto label : labels)
    this->produces<T>(label.instance());
}

void mix::SimInfoOverlayFilter::DeclareProduces()
{
  DeclareProduces< std::vector<simb::MCTruth>     >(fMCTruthInputModuleLabels);
  DeclareProduces< std::vector<simb::MCFlux>      >(fMCFluxInputModuleLabels);
  DeclareProduces< std::vector<simb::GTruth>      >(fGTruthInputModuleLabels);
  DeclareProduces< std::vector<sim::BeamGateInfo> >(fBeamGateInputModuleLabels);

  DeclareProduces< std::vector<simb::MCParticle>      >(fMCParticleInputModuleLabels);
  DeclareProduces< std::vector<sim::SimEnergyDeposit> >(fSimEnergyDepositInputModuleLabels);
  DeclareProduces< std::vector<sim::AuxDetSimChannel> >(fAuxDetSimChannelInputModuleLabels);
  DeclareProduces< std::vector<sim::MCParticleLite>      >(fMCParticleLiteInputModuleLabels); //Ivan
  DeclareProduces< std::vector<sim::SimChannel>       >(fSimChannelInputModuleLabels);
  DeclareProduces< std::vector<sim::SimPhotons>       >(fSimPhotonsInputModuleLabels);

  DeclareProduces< art::Assns<simb::MCTruth,simb::MCFlux,void> >
    (fMCTruthMCFluxAssnsInputModuleLabels);
  DeclareProduces< art::Assns<simb::MCTruth,simb::GTruth,void> >
    (fMCTruthGTruthAssnsInputModuleLabels);
  DeclareProduces< art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo> >
    (fMCTruthMCParticleAssnsInputModuleLabels);
}

template<class T>
void mix::SimInfoOverlayFilter::InitializeCollectionMap(std::vector<art::InputTag> const& labels,
							CollectionMap<T> & colmap)
{
  for(auto const& label : labels)
    colmap[label.instance()] = std::make_unique<T>();
}

void mix::SimInfoOverlayFilter::InitializeCollectionMaps()
{
  InitializeCollectionMap(fMCTruthInputModuleLabels,fMCTruthMap);
  InitializeCollectionMap(fMCFluxInputModuleLabels,fMCFluxMap);
  InitializeCollectionMap(fGTruthInputModuleLabels,fGTruthMap);
  InitializeCollectionMap(fBeamGateInputModuleLabels,fBeamGateInfoMap);

  InitializeCollectionMap(fMCTruthMCFluxAssnsInputModuleLabels,fMCTruthMCFluxAssnsMap);
  InitializeCollectionMap(fMCTruthGTruthAssnsInputModuleLabels,fMCTruthGTruthAssnsMap);

  InitializeCollectionMap(fMCParticleInputModuleLabels,fMCParticleMap);
  InitializeCollectionMap(fSimEnergyDepositInputModuleLabels,fSimEnergyDepositMap);
  InitializeCollectionMap(fAuxDetSimChannelInputModuleLabels,fAuxDetSimChannelMap);
  InitializeCollectionMap(fMCParticleLiteInputModuleLabels,fMCParticleLiteMap); //Ivan
  InitializeCollectionMap(fSimChannelInputModuleLabels,fSimChannelMap);
  InitializeCollectionMap(fSimPhotonsInputModuleLabels,fSimPhotonsMap);

  std::cout << "Filling fMCTruthMCParticleAssnsMap" << std::endl;
  std::cout << "fMCTruthMCParticleAssnsInputModuleLabels.size()=" <<  fMCTruthMCParticleAssnsInputModuleLabels.size() <<  std::endl;
  if(fMCTruthMCParticleAssnsInputModuleLabels.size()){
    std::cout << "Using label " << fMCTruthMCParticleAssnsInputModuleLabels.at(0) << std::endl;
    auto canonical_product_name =  art::canonicalProductName(
        art::friendlyname::friendlyName("art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo>"),
        fMCTruthMCParticleAssnsInputModuleLabels.at(0).label(),
        fMCTruthMCParticleAssnsInputModuleLabels.at(0).instance(),
        fMCTruthMCParticleAssnsInputModuleLabels.at(0).process()
        );
    auto product_id = art::ProductID(canonical_product_name);
    std::cout << "Equivalent ProductID: " << product_id << std::endl;
  }
  InitializeCollectionMap(fMCTruthMCParticleAssnsInputModuleLabels,fMCTruthMCParticleAssnsMap);
  std::cout << "fMCTruthMCParticleAssnsMap.size()=" << fMCTruthMCParticleAssnsMap.size() << std::endl;

  if(fVerbosity>1)
    std::cout << "Finished initializing collection maps." << std::endl;
}

template<class T>
void mix::SimInfoOverlayFilter::PutCollectionOntoEvent(art::Event & e,
						       CollectionMap<T> & colmap)
{
  for(auto & cp : colmap)
    e.put(std::move(cp.second),cp.first);
}

void mix::SimInfoOverlayFilter::PutCollectionsOntoEvent(art::Event & e)
{
  PutCollectionOntoEvent(e,fMCTruthMap);
  PutCollectionOntoEvent(e,fMCFluxMap);
  PutCollectionOntoEvent(e,fGTruthMap);
  PutCollectionOntoEvent(e,fBeamGateInfoMap);

  PutCollectionOntoEvent(e,fMCTruthMCFluxAssnsMap);
  PutCollectionOntoEvent(e,fMCTruthGTruthAssnsMap);

  PutCollectionOntoEvent(e,fMCParticleMap);
  PutCollectionOntoEvent(e,fSimEnergyDepositMap);
  PutCollectionOntoEvent(e,fAuxDetSimChannelMap);
  PutCollectionOntoEvent(e,fMCParticleLiteMap);
  PutCollectionOntoEvent(e,fSimChannelMap);
  PutCollectionOntoEvent(e,fSimPhotonsMap);

  PutCollectionOntoEvent(e,fMCTruthMCParticleAssnsMap);

  if(fVerbosity>1)
    std::cout << "Finished putting collections onto the event." << std::endl;
}

template<class T>
std::map< std::pair<art::ProductID,size_t> , art::Ptr<T> >
mix::SimInfoOverlayFilter::FillCollectionMap(std::vector<art::InputTag> const& labels,
					     CollectionMap< std::vector<T> > & colmap,
					     std::string type_name,
					     art::Event & e)
{
  std::map< std::pair<art::ProductID,size_t> , art::Ptr<T> >
    orig_artptr_lookup;

  for(auto label : labels){
    auto & out_ptrvec = colmap[label.instance()];
    art::PtrMaker<T> makeArtPtr(e, label.instance());

    gallery::Handle< std::vector<T> > handle;
    if(!gEvent.getByLabel< std::vector<T> >(label,handle)) continue;

    //in the future, this will be the way we do it
    //auto product_id = handle.id();

    //but, that's only available in gallery v1_10_00 and up.
    //For now, we gotta do it ourselves!
    if(label.process().size()==0)
      throw cet::exception("SimInfoOverlayFilter")
	<< "ERROR! All InputTags must be FULLY specified until gallery v1_10_00 is available."
	<< " InputTag: " << label << " not fully specified!";
    auto canonical_product_name =  art::canonicalProductName(art::friendlyname::friendlyName(type_name),
							     label.label(),
							     label.instance(),
							     label.process());
    auto product_id = art::ProductID(canonical_product_name);

    if(fVerbosity>2)
      std::cout << "Reading gallery handle. Canonical name = " << canonical_product_name << std::endl;

    auto const& vec(*handle);
    for(size_t i_obj=0; i_obj<vec.size(); ++i_obj){
      out_ptrvec->emplace_back(vec[i_obj]);
      orig_artptr_lookup[std::make_pair(product_id,i_obj)] = makeArtPtr(out_ptrvec->size()-1);
    }
  }

  if(fVerbosity>2)
    std::cout << "\tArtPtr lookup map has " << orig_artptr_lookup.size() << " entries." << std::endl;

  return orig_artptr_lookup;

}

template< class T, class U >
void mix::SimInfoOverlayFilter::FillAssnsCollectionMap(std::vector<art::InputTag> const& labels,
						       CollectionMap< art::Assns<T,U,void> > & colmap,
						       std::map< std::pair<art::ProductID,size_t> , art::Ptr<T> > const&
						       orig_artptr_lookup_left,
						       std::map< std::pair<art::ProductID,size_t> , art::Ptr<U> > const&
						       orig_artptr_lookup_right)
{

  for(auto label : labels){
    auto & out_ptrAssns = colmap[label.instance()];

    gallery::Handle< art::Assns<T,U,void> > handle;
    if(!gEvent.getByLabel< art::Assns<T,U,void> >(label,handle)) continue;

    auto const& assns(*handle);
    for(size_t i_assn=0; i_assn<assns.size(); ++i_assn){
      auto id_left  = std::make_pair(assns[i_assn].first.id(),assns[i_assn].first.key());
      auto id_right = std::make_pair(assns[i_assn].second.id(),assns[i_assn].second.key());

      auto newptr_left = orig_artptr_lookup_left.at(id_left);
      auto newptr_right = orig_artptr_lookup_right.at(id_right);

      out_ptrAssns->addSingle(newptr_left,newptr_right);
    }
  }

}

template< class T, class U, class D>
void mix::SimInfoOverlayFilter::FillAssnsCollectionMap(std::vector<art::InputTag> const& labels,
						       CollectionMap< art::Assns<T,U,D> > & colmap,
						       std::map< std::pair<art::ProductID,size_t> , art::Ptr<T> > const&
						       orig_artptr_lookup_left,
						       std::map< std::pair<art::ProductID,size_t> , art::Ptr<U> > const&
						       orig_artptr_lookup_right)
{

  std::cout << "Starting FillAssnsCollectionMap" << std::endl;

  // Notes on the inputs:
  // colmap = fMCTruthMCParticleAssnsMap, CollectionMap filled with art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo>
  // orig_artptr_lookup_left, mctruth_artptr_lookup, map between product id and art::Ptr<simb::MCTruth> filled before calling this function
  // orig_artptr_lookup_right, mcpart_artptr_lookup, map between product id and art::Ptr<simb::MCParticle> filled at the start of filter

  for(auto label : labels){
    std::cout << "label=" << label << std::endl;
    auto & out_ptrAssns = colmap[label.instance()];

    auto canonical_product_name =  art::canonicalProductName(
        art::friendlyname::friendlyName("art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo>"),
        label.label(),
        label.instance(),
        label.process()
        );
    auto product_id = art::ProductID(canonical_product_name);
    std::cout << "Equivalent ProductID: " << product_id << std::endl;

    gallery::Handle< art::Assns<T,U,D> > handle;
    if(!gEvent.getByLabel< art::Assns<T,U,D> >(label,handle)) continue;

    auto const& assns(*handle);
    for(size_t i_assn=0; i_assn<assns.size(); ++i_assn){
      auto id_left  = std::make_pair(assns[i_assn].first.id(),assns[i_assn].first.key());
      auto id_right = std::make_pair(assns[i_assn].second.id(),assns[i_assn].second.key());
      auto const& data = assns.data(i_assn);

      std::cout << "Looking for Product ID " << id_left.first << " in mctruth_artptr_lookup" << std::endl;
      auto newptr_left = orig_artptr_lookup_left.at(id_left);

      std::cout << "Looking for Product ID " << id_right.first << " in mcpart_artptr_lookup" << std::endl;
      auto newptr_right = orig_artptr_lookup_right.at(id_right);

      out_ptrAssns->addSingle(newptr_left,newptr_right,data);
    }

  }

}

bool mix::SimInfoOverlayFilter::filter(art::Event & e)
{

  std::cout << std::endl << "STARTING FILTER" << std::endl;

  InitializeCollectionMaps();

  //check if we have exhausted our simulation events. If so, we return false.
  if(gEvent.atEnd()) {
    if(fVerbosity>0)
      std::cout << "We've reached the end of the simulation event stream. Returning false." << std::endl;

    PutCollectionsOntoEvent(e);
    return false;
  }

  if(fVerbosity>1){
    std::cout << "Processing input event: "
	      << "Run " << e.run() << ", "
	      << "Event " << e.event() << std::endl;
    std::cout << "Processing simulation event: "
	      << "Run " << gEvent.eventAuxiliary().run() << ", "
	      << "Event " << gEvent.eventAuxiliary().event() << std::endl;
  }

  if(fFillPOTInfo){
    /*
       auto eventsInGalleryFile = gEvent.numberOfEventsInFile();
       gallery::Handle< sumdata::POTSummary > potsum_handle;
       if(!gEvent.getByLabel<sumdata::POTSummary>(fPOTSummaryTag,potsum_handle))
       throw cet::exception("SimInfoOverlayFilter") << "No POTSummary object with tag " << fPOTSummaryTag;

       auto const& potsum(*potsum_handle);
       fPOTSum_totpot += potsum.totpot/eventsInGalleryFile;
       fPOTSum_totgoodpot += potsum.totgoodpot/eventsInGalleryFile;
       fPOTSum_totspills += double(potsum.totspills)/eventsInGalleryFile;
       fPOTSum_goodspills += double(potsum.goodspills)/eventsInGalleryFile;
       */

    fPOTSum_totpot += fSR_POTPerEvent[gEvent.eventAuxiliary().subRun()];
    fPOTSum_totgoodpot += fSR_GoodPOTPerEvent[gEvent.eventAuxiliary().subRun()];
    fPOTSum_totspills += fSR_SpillsPerEvent[gEvent.eventAuxiliary().subRun()];
    fPOTSum_goodspills += fSR_GoodSpillsPerEvent[gEvent.eventAuxiliary().subRun()];

  }

  auto mctruth_artptr_lookup = FillCollectionMap<simb::MCTruth>(fMCTruthInputModuleLabels,
								fMCTruthMap,
								"std::vector<simb::MCTruth>",
								e);
  auto mcflux_artptr_lookup = FillCollectionMap<simb::MCFlux>(fMCFluxInputModuleLabels,
							      fMCFluxMap,
							      "std::vector<simb::MCFlux>",
							      e);
  auto gtruth_artptr_lookup = FillCollectionMap<simb::GTruth>(fGTruthInputModuleLabels,
							      fGTruthMap,
							      "std::vector<simb::GTruth>",
							      e);
  FillCollectionMap<sim::BeamGateInfo>(fBeamGateInputModuleLabels,
				       fBeamGateInfoMap,
				       "std::vector<sim::BeamGateInfo>",
				       e);

  FillAssnsCollectionMap<simb::MCTruth,simb::MCFlux>(fMCTruthMCFluxAssnsInputModuleLabels,
						     fMCTruthMCFluxAssnsMap,
						     mctruth_artptr_lookup,
						     mcflux_artptr_lookup);

  FillAssnsCollectionMap<simb::MCTruth,simb::GTruth>(fMCTruthGTruthAssnsInputModuleLabels,
						     fMCTruthGTruthAssnsMap,
						     mctruth_artptr_lookup,
						     gtruth_artptr_lookup);

  std::cout << "Filling mcpart_artptr_lookup" << std::endl;
  std::cout << "fMCParticleInputModuleLabels.size()= " << fMCParticleInputModuleLabels.size() << std::endl;
  if(fMCParticleInputModuleLabels.size()){
    std::cout << "Using label " << fMCParticleInputModuleLabels.at(0) << std::endl;
    auto canonical_product_name =  art::canonicalProductName(
        art::friendlyname::friendlyName("std::vector<simb::MCParticle>"),
        fMCParticleInputModuleLabels.at(0).label(),
        fMCParticleInputModuleLabels.at(0).instance(),
        fMCParticleInputModuleLabels.at(0).process()
        );
    auto product_id = art::ProductID(canonical_product_name);
    std::cout << "Equivalent ProductID: " << product_id << std::endl;
  }

  auto mcpart_artptr_lookup = FillCollectionMap<simb::MCParticle>(fMCParticleInputModuleLabels,
							          fMCParticleMap,
							          "std::vector<simb::MCParticle>",
							          e);
  std::cout << "mcpart_artptr_lookup.size()=" << mcpart_artptr_lookup.size() << std::endl;


  FillCollectionMap<sim::SimEnergyDeposit>(fSimEnergyDepositInputModuleLabels,
				       fSimEnergyDepositMap,
				       "std::vector<sim::SimEnergyDeposit>",
				       e);

  FillCollectionMap<sim::AuxDetSimChannel>(fAuxDetSimChannelInputModuleLabels,
				       fAuxDetSimChannelMap,
				       "std::vector<sim::AuxDetSimChannel>",
				       e);
  
  FillCollectionMap<sim::MCParticleLite>(fMCParticleLiteInputModuleLabels,
				       fMCParticleLiteMap,
				       "std::vector<sim::MCParticleLite>",
				       e);

  FillCollectionMap<sim::SimChannel>(fSimChannelInputModuleLabels,
				       fSimChannelMap,
				       "std::vector<sim::SimChannel>",
				       e);

  FillCollectionMap<sim::SimPhotons>(fSimPhotonsInputModuleLabels,
				       fSimPhotonsMap,
				       "std::vector<sim::SimPhotons>",
				       e);
  
  // this could be empty if running a "downstream" SimInfoMixer instance
  if(fMCTruthInputModuleLabels.empty() && !fMCTruthMCParticleAssnsInputModuleLabels.empty()) {
    // assume "generator" module generated the mctruth already
    auto mclists = e.getMany<std::vector<simb::MCTruth>>();
    for(size_t mcl = 0; mcl < mclists.size(); ++mcl){
      art::Handle< std::vector<simb::MCTruth> > mclistHandle = mclists[mcl];
      std::string process_name = mclistHandle.provenance()->processName();
      std::string instance_name = mclistHandle.provenance()->productInstanceName();
      std::string module_name = mclistHandle.provenance()->moduleLabel();
      // make sure we are only using products added in this session
      if(process_name != moduleDescription().processName()) {
        continue;
      }
      std::cout << "Making mctruth_artptr_lookup" << std::endl;
      std::cout << "Using label " << module_name << " " << instance_name << " " << fMCTruthMCParticleAssnsInputModuleLabels.begin()->process() << std::endl;
      auto canonical_product_name =  art::canonicalProductName(
          art::friendlyname::friendlyName("std::vector<simb::MCTruth>"),
	  module_name,
	  instance_name,
	  fMCTruthMCParticleAssnsInputModuleLabels.begin()->process());
      auto product_id = art::ProductID(canonical_product_name);

      // If the auxillary file was made with separate Gen/G4 steps - the wrong productID is assigned here
      // causing a crash when tying to construct the assns, extra fhicl parameter to give the correct label
      if(fAuxHasSeparateGenG4){
        std::cout << "Config indicates separate gen/g4 passes used in the aux file" << std::endl;
        canonical_product_name =  art::canonicalProductName(
            art::friendlyname::friendlyName("std::vector<simb::MCTruth>"),
            fMCTruthMCParticleAssnsMCTruthLookupLabel.at(0).label(),
            fMCTruthMCParticleAssnsMCTruthLookupLabel.at(0).instance(),
            fMCTruthMCParticleAssnsMCTruthLookupLabel.at(0).process());
          product_id = art::ProductID(canonical_product_name);
      }
      std::cout << "ProductID: " << product_id << std::endl;
      for(size_t m = 0; m < mclistHandle->size(); ++m){
        art::Ptr<simb::MCTruth> mct(mclistHandle, m);

        mctruth_artptr_lookup[std::make_pair(product_id,m)] = mct;
      }
    }
  }

  std::cout << "List of ProductIDs in mctruth_artptr_lookup:" << std::endl;
  std::map< std::pair<art::ProductID,size_t> , art::Ptr<simb::MCTruth> >::iterator it_T;
  for(it_T = mctruth_artptr_lookup.begin();it_T != mctruth_artptr_lookup.end();it_T++)
        std::cout << it_T->first.first << std::endl;

  std::cout << "List of ProductIDs in mcpart_artptr_lookup:" << std::endl;
  std::map< std::pair<art::ProductID,size_t> , art::Ptr<simb::MCParticle> >::iterator it_U;
  for(it_U = mcpart_artptr_lookup.begin();it_U != mcpart_artptr_lookup.end();it_U++){
        std::cout << it_U->first.first << std::endl;
        break;
  }

  std::cout << "Calling FillAssnsCollectionMap" << std::endl;
  FillAssnsCollectionMap<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo>
    (fMCTruthMCParticleAssnsInputModuleLabels,
     fMCTruthMCParticleAssnsMap,
     mctruth_artptr_lookup,
     mcpart_artptr_lookup);
  std::cout << "Finished FillAssnsCollectionMap" << std::endl;

  //put onto event and loop the gallery event
  PutCollectionsOntoEvent(e);
  gEvent.next();
  return true;
}

// POT map creation updated to handle multiple subruns in auxillary event file
void mix::SimInfoOverlayFilter::MakePOTMap(){

   gEvent.first();

   // Reset everything
   fSR_POTPerEvent.clear();
   fSR_GoodPOTPerEvent.clear();
   fSR_SpillsPerEvent.clear();
   fSR_GoodSpillsPerEvent.clear();

   unsigned int current_sr = 0;
   unsigned int events_in_sr = 0;

   double current_sr_POT = 0;
   double current_sr_GoodPOT = 0;
   int current_sr_Spills = 0;
   int current_sr_GoodSpills = 0;

   // Iterate through all events in auxillary file, count how many are in each subrun
   while(!gEvent.atEnd()){

      if(gEvent.eventAuxiliary().subRun() != current_sr){
         if(current_sr != 0){
            fSR_POTPerEvent[current_sr] = current_sr_POT/events_in_sr;
            fSR_GoodPOTPerEvent[current_sr] = current_sr_GoodPOT/events_in_sr;
            fSR_SpillsPerEvent[current_sr] = (double)current_sr_Spills/events_in_sr;
            fSR_GoodSpillsPerEvent[current_sr] = (double)current_sr_GoodSpills/events_in_sr;
         }
         events_in_sr = 0;
         current_sr = gEvent.eventAuxiliary().subRun();
      }

      gallery::Handle< sumdata::POTSummary > potsum_handle;
      if(!gEvent.getByLabel<sumdata::POTSummary>(fPOTSummaryTag,potsum_handle))
         throw cet::exception("SimInfoOverlayFilterGenG4") << "No POTSummary object with tag " << fPOTSummaryTag;

      auto const& potsum(*potsum_handle);
      current_sr_POT = potsum.totpot;
      current_sr_GoodPOT = potsum.totgoodpot;
      current_sr_Spills = potsum.totspills;
      current_sr_GoodSpills = potsum.goodspills;

      events_in_sr++;
      gEvent.next();
   }

   // Add the last sr
   if(current_sr != 0){
      fSR_POTPerEvent[current_sr] = current_sr_POT/events_in_sr;
      fSR_GoodPOTPerEvent[current_sr] = current_sr_GoodPOT/events_in_sr;
      fSR_SpillsPerEvent[current_sr] = (double)current_sr_Spills/events_in_sr;
      fSR_GoodSpillsPerEvent[current_sr] = (double)current_sr_GoodSpills/events_in_sr;
   }

   // std::map<unsigned int,double>::iterator it;
   // for (it = fSR_POTPerEvent.begin(); it != fSR_POTPerEvent.end(); it++) std::cout << "SR=" << it->first << " POT/Event=" << it->second << std::endl;

   gEvent.first();
}

void mix::SimInfoOverlayFilter::beginJob()
{
}

bool mix::SimInfoOverlayFilter::beginRun(art::Run & r)
{
  return true;
}

bool mix::SimInfoOverlayFilter::beginSubRun(art::SubRun & sr)
{
  if(fFillPOTInfo){
    fPOTSum_totpot = 0.0;
    fPOTSum_totgoodpot = 0.0;
    fPOTSum_totspills = 0.0;
    fPOTSum_goodspills = 0.0;
  }
  return true;
}

void mix::SimInfoOverlayFilter::endJob()
{
}

bool mix::SimInfoOverlayFilter::endRun(art::Run & r)
{
  return true;
}

bool mix::SimInfoOverlayFilter::endSubRun(art::SubRun & sr)
{
  if(fFillPOTInfo){
    std::unique_ptr<sumdata::POTSummary> srpot_ptr(new sumdata::POTSummary());
    srpot_ptr->totpot = fPOTSum_totpot;
    srpot_ptr->totgoodpot = fPOTSum_totgoodpot;
    srpot_ptr->totspills = (int)fPOTSum_totspills;
    srpot_ptr->goodspills = (int)fPOTSum_goodspills;
    sr.put(std::move(srpot_ptr), art::subRunFragment());
  }
  return true;
}

void mix::SimInfoOverlayFilter::respondToCloseInputFile(art::FileBlock const & fb)
{
}

void mix::SimInfoOverlayFilter::respondToCloseOutputFiles(art::FileBlock const & fb)
{
}

void mix::SimInfoOverlayFilter::respondToOpenInputFile(art::FileBlock const  &fb)
{
}

void mix::SimInfoOverlayFilter::respondToOpenOutputFiles(art::FileBlock const & fb)
{
}

DEFINE_ART_MODULE(mix::SimInfoOverlayFilter)