# Install fcl files in /job subdirectory.

install_fhicl()

# Also put a copy in the source tree.

FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )

# add subdirectories
add_subdirectory(FlashMatchSimple)
add_subdirectory(overburden)
add_subdirectory(archive)

