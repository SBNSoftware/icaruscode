# Install fcl files in /job subdirectory.

install_fhicl()

# Also put a copy in the source tree.

FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )

# add subdirectories
add_subdirectory(FlashMatchSimple)
add_subdirectory(ForCITests)
#add_subdirectory(archive)
add_subdirectory(Definitions)
add_subdirectory(larcv)
add_subdirectory(ntuple)
add_subdirectory(Stage0)
add_subdirectory(Stage1)
add_subdirectory(trigger)

