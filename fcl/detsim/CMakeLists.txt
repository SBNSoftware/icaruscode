# add underlying subfolders
add_subdirectory(SBNNoise)
add_subdirectory(commissioning)
add_subdirectory(archive)

# Install fcl files in /job subdirectory.
install_fhicl()

# Also put a copy in the source tree.

FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )
