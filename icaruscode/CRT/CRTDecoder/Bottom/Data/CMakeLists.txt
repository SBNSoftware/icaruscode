# project name
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
#include_directories( "$ENV{DUNEPDSPRCE_INC}" )
#cet_find_library( RCEDAMLIB NAMES protodune-dam PATHS ENV DUNEPDSPRCE_LIB NO_DEFAULT_PATH )
#project(dune-raw-data_Overlays)


art_make_library(
  LIBRARIES
    artdaq_core::artdaq-core_Data
  )
art_dictionary(DICTIONARY_LIBRARIES icaruscode::CRT_CRTDecoder_Bottom_Data)
install_headers()
install_source()
