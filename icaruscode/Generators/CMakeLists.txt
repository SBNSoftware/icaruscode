#find_ups_product(ifdhc)
#include_directories( $ENV{IFDHC_FQ_DIR}/inc )
#set( IFDH_AS_SERVICE_LIB_LIST
#     art_Framework_Services_Registry
#     ${Boost_SYSTEM_LIBRARY}
#     ${IFDH_SERVICE} )
#include_directories ( . )

cet_make_library(
      LIBRARY_NAME     
                  icaruscode_hepmc
      LIB_LIBRARIES
                  larcorealg::Geometry
                  larcore::Geometry_Geometry_service
                  lardataobj::RecoBase
                  lardataobj::AnalysisBase 
                  lardataobj::RawData
                  lardata::RecoObjects 
                  lardata::Utilities
                  lardataobj::Simulation
                  nusimdata::SimulationBase
                  larcoreobj::SummaryData
                  lardataobj::RawData
                  lardataobj::RecoBase
                  lardataobj::AnalysisBase
                  ifdh_art::IFDH_service
                  art::Persistency_Provenance
                  art::Persistency_Common
                  art::Utilities
                  canvas::canvas
                  cetlib::cetlib
                  cetlib_except::cetlib_except
                  ROOT::Core
                  ROOT::Core ROOT::Physics
)

cet_build_plugin(HepMCFileGen art::module 
         SOURCE HepMCFileGen_module.cc 
      LIBRARIES larcorealg::Geometry
                larcore::Geometry_Geometry_service
                lardataobj::RecoBase
                lardataobj::AnalysisBase
                lardataobj::RawData
                lardata::RecoObjects
                lardata::Utilities
                lardataobj::Simulation
                nusimdata::SimulationBase
                larcoreobj::SummaryData
                lardataobj::RawData
                lardataobj::RecoBase
                lardataobj::AnalysisBase
                ifdh_art::IFDH_service
                art::Persistency_Provenance
                art::Persistency_Common
                art::Utilities
                canvas::canvas
                cetlib::cetlib
                cetlib_except::cetlib_except
                ROOT::Core
                ROOT::Physics
)
     
install_headers()
install_source()
install_fhicl()
