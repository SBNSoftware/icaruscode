
art_make(
  SERVICE_LIBRARIES
      icarusalg_Geometry
      larcorealg_Geometry
      art::Utilities
      messagefacility::MF_MessageLogger
      fhiclcpp::fhiclcpp
      cetlib::cetlib
      cetlib_except::cetlib_except
      ROOT::Geom
      ROOT::XMLIO
      ROOT::Gdml
      ${ROOT_BASIC_LIB_LIST}
  TOOL_LIBRARIES
      icarusalg_Geometry
      larcorealg_Geometry
      art::Utilities
      
      fhiclcpp::fhiclcpp
      cetlib::cetlib
      cetlib_except::cetlib_except
      ROOT::Geom
      ROOT::XMLIO
      ROOT::Gdml
      ${ROOT_BASIC_LIB_LIST}    
  )

