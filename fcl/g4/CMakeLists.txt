# Install fcl files in /job subdirectory.

install_fhicl()

# Also put a copy in the source tree.

FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )

# Add the underlying folders

add_subdirectory(modifiers)
add_subdirectory(2d_drift_sim)

