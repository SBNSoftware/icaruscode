cet_build_plugin(HepMCFileGen art::module
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

cet_build_plugin(SelectCathodeCrossingGenParticles art::module
      LIBRARIES sbnalg::Utilities
                larcore::Geometry_Geometry_service
                larcore::headers
                lardataalg::MCDumpers
                larcorealg::Geometry
                larcorealg::headers
                nusimdata::SimulationBase
)

install_headers()
install_source()
install_fhicl()
