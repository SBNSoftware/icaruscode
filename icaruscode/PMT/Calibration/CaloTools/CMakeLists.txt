cet_enable_asserts()

add_definitions(-DEIGEN_FFTW_DEFAULT)

art_make(
          LIB_LIBRARIES
                        larcorealg::Geometry
                        larcore::Geometry_Geometry_service
                        lardata::Utilities
                        larevt::Filters
                        lardataobj::RawData
                        lardata::ArtDataHelper
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        art_root_io::TFileService_service
                        art_root_io::tfile_support
                        art::Persistency_Provenance
                        art::Persistency_Common
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        Eigen3::Eigen
                        FFTW3::FFTW3
                        Boost::system
			ROOT::Core
			ROOT::Gpad
        )

install_headers()
install_fhicl()
install_source()
