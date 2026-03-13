
set( IFDH_ART_DIR $ENV{IFDH_ART_DIR} )

add_subdirectory(Timing)
add_subdirectory(IcarusObj) 
add_subdirectory(Analysis)
add_subdirectory(CRT)
add_subdirectory(Decode)
add_subdirectory(Filters)
add_subdirectory(Generators)
add_subdirectory(Geometry)
# LegacyLArG4 is not available when building with geant 4.11
#add_subdirectory(LArG4)
add_subdirectory(PMT)
add_subdirectory(RecoUtils)
add_subdirectory(TPC)
add_subdirectory(Utilities)
add_subdirectory(Overlays)
#add_subdirectory(Supera)
