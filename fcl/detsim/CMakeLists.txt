# add underlying subfolders
add_subdirectory(overburden)
add_subdirectory(SBNNov17)
add_subdirectory(SBNMar18)
add_subdirectory(SBNNoise)

# Install fcl files in /job subdirectory.
install_fhicl()

# Also put a copy in the source tree.

FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )
