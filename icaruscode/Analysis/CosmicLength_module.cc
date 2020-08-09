/**
 * @file   CosmicLength_module.cc
 * @brief  Access CRT data and reco products and compare to MCTruth info 
 * @author Ryan Howell (ryan.howell@rochester.edu)
 * 
 * The last revision of this code was done in August 2020.
 */

// LArSoft includes
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/Simulation/LArG4Parameters.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TROOT.h"

// C++ includes
#include <map>
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <iostream>
#include <utility>
#include <array>

using std::string;
using std::vector;
using std::map;
using std::set;
using std::pair;

//using cmath::pow;
//using cmath::sqrt;

namespace icarus {
    namespace crt {

        //-----------------------------------------------------------------------
        //-----------------------------------------------------------------------
        // class definition
        /**
         * @brief Example analyzer
         * 
         * This class extracts information from the generated and reconstructed
         * particles.
         *
         * It produces histograms for the simulated particles in the input file:
         * - PDG ID (flavor) of all particles
         * - momentum of the primary particles selected to have a specific PDG ID
         * - length of the selected particle trajectory
         * 
         * It also produces two ROOT trees.
         *
         * The first ROOT tree contains information on the simulated
         * particles, including "dEdx", a binned histogram of collected
         * charge as function of track range.
         * 
         * Configuration parameters
         * =========================
         * 
         * - *SimulationLabel* (string, default: "largeant"): tag of the input data
         *   product with the detector simulation information (typically an instance
         *   of the LArG4 module)
         *
         */
        class CosmicLength: public art::EDAnalyzer {
            public:

                struct Config {

                    // Save some typing:
                    using Name = fhicl::Name;
                    using Comment = fhicl::Comment;

                    // One Atom for each parameter
                    fhicl::Atom<art::InputTag> SimulationLabel {
                        Name("SimulationLabel"),
                            Comment("tag of the input data product with the detector simulation information")
                    };
                }; // Config

            using Parameters = art::EDAnalyzer::Table<Config> ;

            // -------------------------------------------------------------------
            // -------------------------------------------------------------------
            // Standard constructor for an ART module with configuration validation;
            // we don't need a special destructor here.

            /// Constructor: configures the module (see the Config structure above)
            explicit CosmicLength(Parameters
                const &config);

            virtual void beginJob() override;
            virtual void beginRun(const art::Run &run) override;
            virtual void analyze(const art::Event &event) override;

            private:

            void fillHistograms(art::Handle<vector<simb::MCTruth>> genHandle, art::Handle<vector<sim::SimEnergyDeposit>> depositHandle, art::Handle<vector<simb::MCParticle>> particleHandle);
            void fillTree(art::Handle<vector<simb::MCTruth>> genHandle, art::Handle<vector<simb::MCParticle>> particleHandle);

            // The parameters we'll read from the .fcl file.
            art::InputTag fSimulationProducerLabel; ///< The name of the producer that tracked simulated particles through the detector

            // The n-tuples we'll create.
            TTree *fSingleCosmicMC;

            TH1D *fTotalTimeDiff;
            TH1D *fPrimaryStartTime;
            TH1I *pdgs;
            TH1I *nGen;

            TH1D *hDistance; 

            TH2I *pdgDaughters;
            TH2I *pdgTime;

            // The comment lines with the @ symbols define groups in doxygen. 
            /// @name The variables that will go into both n-tuples.
            /// @{
            int fEvent; ///< number of the event being processed
            int fRun; ///< number of the run being processed
            int fSubRun; ///< number of the sub-run being processed
            /// @}

            /// @name The variables that will go into the Simulation n-tuple.
            /// @{

            int fNGen;
            int fNDaught; ///< number of daughters belonging to this MCParticle
            uint32_t fSimHits; ///< number of trajectory points for each MCParticle

            double fPrimaryTimeDiff;
            double firstTime;
            double lastTime;

            double fTotalDiff;
            double totalFirstTime;
            double totalLastTime;

            double fEntry[4];
            double fExit[4];
            double fDistance;

            int fPDG;
            int fContained;

            int fFullyContained;
            bool fill;

            std::vector<int>fMotherTracks;
            std::map<int,
            int>genPdgTrackMap;
            std::map<int,
            sim::SimEnergyDeposit>depositMap;

            const std::map<int,
                const char*>binToName = {
                    {
                        1,
                        "neutron"
                    },
                    {
                        2,
                        "proton"
                    },
                    {
                        3,
                        "electon"
                    },
                    {
                        4,
                        "positron"
                    },
                    {
                        5,
                        "muon"
                    },
                    {
                        6,
                        "anti-muon"
                    },
                    {
                        7,
                        "gamma"
                    }
                };
            const std::map<int, int>pdgToBin = {
                {
                    2112,
                    1
                },
                {
                    2212,
                    2
                },
                {
                    11,
                    3
                },
                {
                    -11,
                    4
                },
                {
                    13,
                    5
                },
                {
                    -13,
                    6
                },
                {
                    22,
                    7
                }
            };
            const Int_t nBins = 7;

            // Other variables that will be shared between fTotalTimeDifferent methods.
            geo::GeometryCore
            const *fGeometryService; ///< pointer to Geometry provider
            int fTriggerOffset; ///< (units of ticks) time of expected neutrino event

        }; // class CosmicLength

        //-----------------------------------------------------------------------
        //-----------------------------------------------------------------------
        // class implementation

        //-----------------------------------------------------------------------
        // Constructor
        // 
        // Note that config is a Table<Config>, and to access the Config
        // value we need to use an operator: "config()". In the same way,
        // each element in Config is an Atom<Type>, so to access the type we
        // again use the call operator, e.g. "SimulationLabel()".

        CosmicLength::CosmicLength(Parameters
            const &config): EDAnalyzer(config), fSimulationProducerLabel(config().SimulationLabel()) {
            // Get a pointer to the geometry service provider.
            fGeometryService = lar::providerFrom<geo::Geometry> ();
            // The same for detector TDC clock services.
            //fTimeService = lar::providerFrom<detinfo::DetectorClocksService>();
            // Access to detector properties.
            const detinfo::DetectorProperties *detprop = lar::providerFrom<detinfo::DetectorPropertiesService> ();
            fTriggerOffset = detprop->TriggerOffset();
        }

        //-----------------------------------------------------------------------
        void CosmicLength::beginJob() {
            std::cout << "starting analysis job" << std::endl;

            // Access ART's TFileService, which will handle creating and writing
            // histograms and n-tuples for us. 
            art::ServiceHandle<art::TFileService> tfs;

            // Define our n-tuples
            fSingleCosmicMC = tfs->make<TTree> ("cosmicMC", "cosmicLengthSim");
            nGen = tfs->make<TH1I>("nGen", "Number of Generated Particles per Cosmic Sample", 60, 0, 5000);

            fTotalTimeDiff = tfs->make<TH1D>("fTotalTimeDiff", "Primary Cosmic Interaction Time in Detector", 160, 3e5, 1e7);
            fTotalTimeDiff->GetXaxis()->SetTitle("Primary Cosmic Last Interaction Time - First Interaction Time (ns)");

            pdgs = tfs->make<TH1I>("pdgs", "Primary Cosmic PDGs", 2450, -200, 2250);
            fPrimaryStartTime = tfs->make<TH1D>("fPrimaryStartTime", "Primary Cosmic Start Time", 60, -1150000, 1150000);

            pdgDaughters = tfs->make<TH2I>("pdgDaughters", "Number of Daughters per Primary Cosmic", 50, 0, 200, 7, 0.5, 7.5);
            pdgTime = tfs->make<TH2I>("pdgTime", "Length of Detector Interactions per Primary Cosmic", 300, 0, 300, 7, 0.5, 7.5);

            hDistance = tfs->make<TH1D>("hDistance","Distance of Throughgoing Particles",50,0,1000);

            pdgDaughters->GetXaxis()->SetTitle("number of daughter particles that entered a cryostat");
            pdgTime->GetXaxis()->SetTitle("interaction length in detector (ns)");

            for (int i = 1; i<8; i++) {
                const char * hLabel = binToName.at(i);
                pdgDaughters->GetYaxis()->SetBinLabel(i, hLabel);
                pdgTime->GetYaxis()->SetBinLabel(i, hLabel);
            }

            // Define the branches of our simulation n-tuple
            fSingleCosmicMC->Branch("event", &fEvent, "event/I");
            fSingleCosmicMC->Branch("entryPos", &fEntry, "entryPos[4]/F");
            fSingleCosmicMC->Branch("exitPos", &fExit, "exitPos[4]/F");
            fSingleCosmicMC->Branch("distance", &fDistance, "distance/D");
            fSingleCosmicMC->Branch("pdg", &fPDG, "pdg/I");
            fSingleCosmicMC->Branch("fullyContained", &fFullyContained, "fullyContained/I");

        }

        void CosmicLength::beginRun(const art::Run &/*run*/ ) {
            //art::ServiceHandle<sim::LArG4Parameters> larParameters;
            //fElectronsToGeV = 1./larParameters->GeVToElectrons();
            std::cout << "beginning run" << std::endl;
        }

        //-----------------------------------------------------------------------
        void CosmicLength::fillHistograms(art::Handle<vector<simb::MCTruth>> genHandle, art::Handle<vector<sim::SimEnergyDeposit>> depositHandle, art::Handle<vector<simb::MCParticle>> particleHandle) {

            //get TPC objects
            geo::CryostatGeo const &cryo0 = fGeometryService->Cryostat(0);
            geo::CryostatGeo const &cryo1 = fGeometryService->Cryostat(1);
            
            geo::TPCGeo const &tpc00 = cryo0.TPC(0);
            geo::TPCGeo const &tpc01 = cryo0.TPC(1);
            //geo::TPCGeo const &tpc02 = cryo0.TPC(2);
            //geo::TPCGeo const &tpc03 = cryo0.TPC(3);

            geo::TPCGeo const &tpc10 = cryo1.TPC(0);
            geo::TPCGeo const &tpc11 = cryo1.TPC(1);
            //geo::TPCGeo const &tpc12 = cryo1.TPC(2);
            //geo::TPCGeo const &tpc13 = cryo1.TPC(3);
            
            auto const &truth = (*genHandle)[0];
            fNGen = truth.NParticles();
            nGen->Fill(fNGen);

            for (int i = 0; i<fNGen; i++) {
                auto const &part = truth.GetParticle(i); //simb::MCParticle
                genPdgTrackMap.insert({part.TrackId(),part.PdgCode()});
                fPrimaryStartTime->Fill(part.T(0));
            }

            for (auto const &deposit: (*depositHandle)) {
                depositMap.insert({deposit.TrackID(),deposit});
            }

            totalFirstTime = 1e19;
            totalLastTime = -1e19;
            
            //loop over primary particles
            for (auto const &pair: genPdgTrackMap) {
                int motherID = pair.first;
                int genPdg = pair.second;

                fill = false;

                fNDaught = 0;
                firstTime = 1e19;
                lastTime = -1e19;

                for (auto const &particle: (*particleHandle)) {
                    sim::SimEnergyDeposit energyDeposit = depositMap.find(particle.TrackId())->second;

                    if (particle.Mother() == motherID && energyDeposit.Energy()>.05) {
                        fNDaught++;

                        int fSimHits = particle.NumberTrajectoryPoints();

                        bool wasInCryo = false;

                        //loop over trajectory points
                        for (int i = 0; i<fSimHits; i++) {
                            const TLorentzVector &pos = particle.Position(i); // 4-position in World coordinates

                            const double point[3] = {
                                pos.X(),
                                pos.Y(),
                                pos.Z()
                            };

                            if (tpc00.ContainsPosition(point) || tpc01.ContainsPosition(point) || tpc10.ContainsPosition(point) || tpc11.ContainsPosition(point)) {
                                fill = true;
                                wasInCryo = true;
                            }

                            //finds the earliest time a primary's particle was in the cryostat
                            if (pos.T()<firstTime && wasInCryo) {
                                firstTime = pos.T();
                            }

                            //finds the earliest time anything was in the cryostat
                            if (pos.T()<totalFirstTime && wasInCryo) {
                                totalFirstTime = pos.T();
                            }

                            //finds the last time a primary's particle was in the cryostat
                            if (pos.T()>lastTime && wasInCryo) {
                                lastTime = pos.T();
                            }

                            //finds the last time anything was in the cryostat
                            if (pos.T()>totalLastTime && wasInCryo) {
                                totalLastTime = pos.T();
                            }

                            if (!tpc00.ContainsPosition(point) && !tpc01.ContainsPosition(point) && !tpc10.ContainsPosition(point) && !tpc11.ContainsPosition(point)) {
                                wasInCryo = false;
                            }
                        } //for trajectory points
                    } //if particle is daughter of the primary cosmic
                } //for particles

                if (fill) {
                    fPrimaryTimeDiff = lastTime - firstTime;

                    pdgDaughters->Fill(fNDaught, pdgToBin.at(genPdg));
                    pdgTime->Fill(fPrimaryTimeDiff, pdgToBin.at(genPdg));
                    pdgs->Fill(genPdg);
                }

            } //for primary particles

            fTotalTimeDiff->Fill(totalLastTime - totalFirstTime);
        }

        void CosmicLength::fillTree(art::Handle<vector<simb::MCTruth>> genHandle, art::Handle<vector<simb::MCParticle>> particleHandle) {
            //get TPC objects
            geo::CryostatGeo const &cryo0 = fGeometryService->Cryostat(0);
            geo::CryostatGeo const &cryo1 = fGeometryService->Cryostat(1);
            
            geo::TPCGeo const &tpc00 = cryo0.TPC(0);
            geo::TPCGeo const &tpc01 = cryo0.TPC(1);
            //geo::TPCGeo const &tpc02 = cryo0.TPC(2);
            //geo::TPCGeo const &tpc03 = cryo0.TPC(3);

            geo::TPCGeo const &tpc10 = cryo1.TPC(0);
            geo::TPCGeo const &tpc11 = cryo1.TPC(1);
            //geo::TPCGeo const &tpc12 = cryo1.TPC(2);
            //geo::TPCGeo const &tpc13 = cryo1.TPC(3);
            
            auto const &truth = (*genHandle)[0];
            fNGen = truth.NParticles();

            for (int i = 0; i<fNGen; i++) {
                auto const &part = truth.GetParticle(i); //simb::MCParticle
                int fSimHits = part.NumberTrajectoryPoints();

                float firstTime = 1e19;
                float lastTime = -1e19;

                bool wasInCryo = false;

                fFullyContained = 1;

                for (int i = 0; i<fSimHits; i++) {
                    const TLorentzVector &pos = part.Position(i); // 4-position in World coordinates
                    const double point[3] = {pos.X(),pos.Y(),pos.Z()};

                    if (tpc00.ContainsPosition(point) || tpc01.ContainsPosition(point) || tpc10.ContainsPosition(point) || tpc11.ContainsPosition(point)) {
                        fPDG = part.PdgCode();
                        fill = true;
                        wasInCryo = true;
                    }

                    if (pos.T()<firstTime && wasInCryo) {
                        firstTime = pos.T();
                        fEntry[0] = pos.X();
                        fEntry[1] = pos.Y();
                        fEntry[2] = pos.Z();
                        fEntry[3] = pos.T();
                    }

                    //finds the last time a primary's particle was in the cryostat
                    if (pos.T()>lastTime && wasInCryo) {
                        lastTime = pos.T();
                        fExit[0] = pos.X();
                        fExit[1] = pos.Y();
                        fExit[2] = pos.Z();
                        fExit[3] = pos.T();
                    }

                    if (wasInCryo && !tpc00.ContainsPosition(point) && !tpc01.ContainsPosition(point) && !tpc10.ContainsPosition(point) && !tpc11.ContainsPosition(point)) {
                        wasInCryo = false;
                        fFullyContained = 0;
                    }
                }

                if (wasInCryo) {
                    fDistance = sqrt(pow((fExit[0]-fEntry[0]),2) + pow((fExit[1]-fEntry[1]),2) + pow((fExit[2]-fEntry[2]),2));
                    hDistance -> Fill(fDistance);
                    fSingleCosmicMC->Fill();
                }
            }

            for (auto const &part: (*particleHandle)) {
                int fSimHits = part.NumberTrajectoryPoints();

                float firstTime = 1e19;
                float lastTime = -1e19;

                bool wasInCryo = false;

                fFullyContained = 1;

                for (int i = 0; i<fSimHits; i++) {
                    const TLorentzVector &pos = part.Position(i); // 4-position in World coordinates
                    const double point[3] = {pos.X(),pos.Y(),pos.Z()};

                    if (tpc00.ContainsPosition(point) || tpc01.ContainsPosition(point) || tpc10.ContainsPosition(point) || tpc11.ContainsPosition(point)) {

                        fPDG = part.PdgCode();
                        fill = true;
                        wasInCryo = true;
                    }

                    if (pos.T()<firstTime && wasInCryo) {
                        firstTime = pos.T();
                        fEntry[0] = pos.X();
                        fEntry[1] = pos.Y();
                        fEntry[2] = pos.Z();
                        fEntry[3] = pos.T();
                    }

                    //finds the last time a primary's particle was in the cryostat
                    if (pos.T()>lastTime && wasInCryo) {
                        lastTime = pos.T();
                        fExit[0] = pos.X();
                        fExit[1] = pos.Y();
                        fExit[2] = pos.Z();
                        fExit[3] = pos.T();
                    }

                    if (wasInCryo && !tpc00.ContainsPosition(point) && !tpc01.ContainsPosition(point) && !tpc10.ContainsPosition(point) && !tpc11.ContainsPosition(point)) {
                        wasInCryo = false;
                        fFullyContained = 0;
                    }
                }

                if (wasInCryo) {
                    fDistance = sqrt(pow((fExit[0]-fEntry[0]),2) + pow((fExit[1]-fEntry[1]),2) + pow((fExit[2]-fEntry[2]),2));
                    hDistance -> Fill(fDistance);
                    fSingleCosmicMC->Fill();
                }
            }
        }

        //-----------------------------------------------------------------------
        void CosmicLength::analyze(const art::Event &event) {
            MF_LOG_DEBUG("CRT") << "beginning analyis" << '\n';

            // Start by fetching some basic event information for our n-tuple.
            fEvent = event.id().event();
            fRun = event.run();
            fSubRun = event.subRun();

            // Define "handle" to Generator level MCTruth objects
            art::Handle<vector<simb::MCTruth>> genHandle;

            // Define a "handle" to point to a vector of MCParticle objects.
            art::Handle<vector<simb::MCParticle>> particleHandle;

            // Define a "handle" to point to a vector of SimEnergyDeposit objects
            art::Handle<vector<sim::SimEnergyDeposit>> depositHandle;

            if (!event.getByLabel("generator", genHandle)) {
                std::cout << "could not get handle to gen objects!!!" << std::endl;
            }

            if (!event.getByLabel("ionization", depositHandle)) {
                std::cout << "could not get handle to ionization objects!!!" << std::endl;
            }

            // If there aren't any simb::MCParticle object art will 
            // display a "ProductNotFound" exception message and may skip
            // all processing for the rest of this event or stop the execution.
            if (!event.getByLabel(fSimulationProducerLabel, particleHandle)) {
                // If we have no MCParticles at all in an event, then we're in
                // big trouble. Throw an exception.
                throw cet::exception("CosmicLength") <<
                    " No simb::MCParticle objects in this event - " <<
                    " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
            }

            if ((*genHandle).size()>1) throw cet::exception("CosmicLength") << "gen stage MCParticle vector has more than 1 entry!" << std::endl;

            fillHistograms(genHandle, depositHandle, particleHandle);
            fillTree(genHandle, particleHandle);

            /*
            const TGeoVolume* tpc00Active = tpc00.ActiveVolume();
            const TGeoVolume* tpc01Active = tpc01.ActiveVolume();
            const TGeoVolume* tpc02Active = tpc02.ActiveVolume();
            const TGeoVolume* tpc03Active = tpc03.ActiveVolume();

            const TGeoVolume* tpc10Active = tpc10.ActiveVolume();
            const TGeoVolume* tpc11Active = tpc11.ActiveVolume();
            const TGeoVolume* tpc12Active = tpc12.ActiveVolume();
            const TGeoVolume* tpc13Active = tpc13.ActiveVolume();
            */

        } // CosmicLength::analyze()

        DEFINE_ART_MODULE(CosmicLength)

    } // namespace crt
} // namespace icarus