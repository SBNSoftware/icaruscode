cet_make(
	SUBDIRS
		"details"
	LIBRARIES
		range-v3::range-v3
		art::Framework_Principal
		sbnobj::ICARUS_PMT_Trigger_Data
		icaruscode::PMT_Algorithms
		lardataalg::DetectorInfo
		lardataalg::MCDumpers
		lardataobj::RawData
		icarusalg::Utilities
		nusimdata::SimulationBase
		messagefacility::MF_MessageLogger
		fhiclcpp::fhiclcpp
	)

install_headers(SUBDIRS "details")
install_source(SUBDIRS "details")
install_fhicl(SUBDIRS "details")
