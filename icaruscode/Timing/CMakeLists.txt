cet_enable_asserts()

set(    MODULE_LIBRARIES
        lardataobj::RecoBase
        larcore::Geometry_Geometry_service
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
	art_root_io::TFileService_service
        art_root_io::tfile_support
        lardataobj::RawData
        art::Framework_Core
        art::Framework_Principal
        art::Utilities
        canvas::canvas
        messagefacility::MF_MessageLogger
        fhiclcpp::fhiclcpp
        cetlib::cetlib
        cetlib_except::cetlib_except
        ROOT::Tree
)
set(	LIB_LIBRARIES
        art::Framework_Services_Registry
        messagefacility::MF_MessageLogger
        lardataobj::RecoBase
        lardata::Utilities
        icaruscode_IcarusObj
        sbnobj::Common_PMT_Data
	larevt::CalibrationDBI_IOVData
	larevt::CalibrationDBI_Providers
)
set(	SERVICE_LIBRARIES
	icaruscode_Timing
        icaruscode_IcarusObj
        larcore::Geometry_Geometry_service
        lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
)

file(GLOB lib_srcs *.cxx)

art_make_library( SOURCE ${lib_srcs} LIBRARIES PUBLIC ${LIB_LIBRARIES})

cet_build_plugin( OpHitTimingCorrection art::module LIBRARIES PUBLIC ${MODULE_LIBRARIES})

cet_build_plugin( PMTTimingCorrectionService art::service LIBRARIES PUBLIC ${SERVICE_LIBRARIES})

cet_build_plugin( PMTBeamSignalsExtractor art::producer LIBRARIES PUBLIC ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()
