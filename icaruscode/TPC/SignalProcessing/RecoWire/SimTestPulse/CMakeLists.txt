

art_make(	LIB_LIBRARIES
			lardataobj::RawData
			lardataobj::RecoBase
			lardata::Utilities
			fhiclcpp::fhiclcpp
			cetlib::cetlib

		MODULE_LIBRARIES
			icaruscode::TPC_SignalProcessing_RecoWire_SimTestPulse
			icaruscode::TPC_Utilities_SignalShapingICARUSService_service
			lardata::Utilities
			lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
			lardataobj::Simulation
			larcoreobj::SummaryData
			larcorealg::Geometry
			larevt::CalibrationDBI_IOVData
			larevt::CalibrationDBI_Providers
			larreco::HitFinder
			lardataobj::RecoBase
			larcorealg::Geometry
			larcore::Geometry_Geometry_service
			lardata::Utilities
			nurandom::RandomUtils_NuRandomService_service
			lardata::Utilities
			lardataobj::RawData
			larevt::CalibrationDBI_IOVData
			larevt::CalibrationDBI_Providers
			lardataobj::RecoBase
			lardata::ArtDataHelper
			art::Framework_Core
			art::Framework_Principal
			art::Framework_Services_Registry
			art_root_io::tfile_support ROOT::Core
			art_root_io::TFileService_service
			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
			canvas::canvas
			messagefacility::MF_MessageLogger
			messagefacility::headers
			fhiclcpp::fhiclcpp
			cetlib::cetlib
			ROOT::FFTW
			ROOT::Tree
)

install_fhicl()

