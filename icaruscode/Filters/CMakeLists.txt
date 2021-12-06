art_make( 
          MODULE_LIBRARIES
                        icarusalg_Utilities
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation 
                        nug4::ParticleNavigation
                        lardataobj_Simulation
                        lardata_Utilities
                        larevt_Filters
                        lardataalg_MCDumpers
                        lardataobj_RawData
                        lardataobj_RecoBase
                        lardata_RecoObjects
                        nusimdata::SimulationBase
                        art::Framework_Core
                        art::Framework_Principal
                        art::Framework_Services_Registry
                        ROOT::Core
                        art::Persistency_Common 
                        art::Persistency_Provenance 
                        art::Utilities 
                        messagefacility::MF_MessageLogger
                        
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib
                        cetlib_except::cetlib_except
                        ROOT::EG
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
                        ${ROOT_GENVECTOR}
                        ${ROOT_BASIC_LIB_LIST}
        )

# install_headers()
install_fhicl()
install_source()

