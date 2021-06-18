# Wire-Cell configuration files, required at runtime.
FILE(GLOB cfg_file *.jsonnet *.fcl )
message(STATUS "Configuration files ${cfg_file}")
#install(FILES ${data} DESTINATION #${product}/${version}/WireCellData/pgrapher/experiment/icarus)
install_wp(LIST ${cfg_file} SUBDIRNAME pgrapher/experiment/icarus)

#add_subdirectory(fhicl)

