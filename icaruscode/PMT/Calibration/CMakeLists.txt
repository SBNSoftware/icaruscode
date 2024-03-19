add_subdirectory(CaloTools)
add_subdirectory(fcl)

cet_build_plugin(PMTLaserCalibration art::module
                    LIBRARIES
                icaruscode_PMT_Algorithms
                icaruscode_PMT_Calibration_CaloTools
                icaruscode_IcarusObj
		lardataobj::RawData
                lardataobj::Simulation
                lardataobj::RecoBase
                art::Framework_Core
                art::Framework_Principal
                art::Framework_Services_Registry
                art_root_io::TFileService_service
                art_root_io::tfile_support
                art::Persistency_Provenance
                canvas::canvas
                messagefacility::MF_MessageLogger
                fhiclcpp::fhiclcpp
                cetlib::cetlib
                cetlib_except::cetlib_except
                Boost::system
              )

cet_build_plugin(PMTBackgroundphotonsCalibration art::module
                    LIBRARIES
                icaruscode_PMT_Algorithms
                icaruscode_PMT_Calibration_CaloTools
                lardataobj::RawData
                lardataobj::Simulation
                lardataobj::RecoBase
                art::Framework_Core
                art::Framework_Principal
                art::Framework_Services_Registry
                art_root_io::TFileService_service
                art_root_io::tfile_support
                art::Persistency_Provenance
                canvas::canvas
                messagefacility::MF_MessageLogger
                fhiclcpp::fhiclcpp
                cetlib::cetlib
                cetlib_except::cetlib_except
                Boost::system
              )


art_make_exec(  NAME 
                    fitPulseDistribution 
                SOURCE 
                    fitPulseDistribution.cc
                LIBRARIES
                    icaruscode_PMT_Algorithms
                    icaruscode_PMT_Calibration_CaloTools
              )


install_headers()
install_fhicl()
install_source()

