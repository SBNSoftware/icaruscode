# Install fcl files in /job subdirectory.

install_fhicl()

# Also put a copy in the source tree.
FILE(GLOB fcl_files *.fcl)
install_source( EXTRAS ${fcl_files} )


# Also install the xml file(s) for this
# They need to be found via FW_SEARCH_PATH

# install xml files
file(GLOB xml_file *.xml)
install_fw(LIST ${xml_file})

