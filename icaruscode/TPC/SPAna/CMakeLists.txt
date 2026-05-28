cet_build_plugin(SPAna art::module
                       LIBRARIES
                            art::Framework_Core
                            art::Framework_Services_Registry
                            art_root_io::TFileService_service
                            art::Framework_Principal
                            art::Persistency_Common
                            art::Utilities canvas::canvas
                            cetlib::cetlib cetlib_except::cetlib_except
                            ROOT::X3d
                            art_root_io::tfile_support
                            art_root_io::art_root_io
                            art_root_io::dict
                            lardataobj::AnalysisBase
                            lardataobj::RecoBase
                            larcorealg::Geometry
                            larcore::Geometry_Geometry_service
                            larcorealg::GeoAlgo
                            sbnobj::Common_Reco
                            larcorealg::GeoAlgo
                            lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
                            lardataalg::DetectorInfo
                            ROOT::Minuit
                            sbnobj::Common_Calibration_dict
                            larevt::SpaceCharge
			    nug4::ParticleNavigation
			    larsim::MCCheater_BackTrackerService_service
			    larsim::MCCheater_ParticleInventoryService_service
			    larsim::Simulation
)

install_headers()
install_source()
install_fhicl()
