# Build the module
art_make( MODULE_LIBRARIES
           larcorealg_Geometry
           larreco_Calorimetry
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
           ${MF_MESSAGELOGGER}
         )

install_headers()
install_fhicl()
install_source()

