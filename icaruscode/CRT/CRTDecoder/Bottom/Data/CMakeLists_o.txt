# project name
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
#include_directories( "$ENV{DUNEPDSPRCE_INC}" )
#cet_find_library( RCEDAMLIB NAMES protodune-dam PATHS ENV DUNEPDSPRCE_LIB NO_DEFAULT_PATH )
#project(dune-raw-data_Overlays)
art_make(BASENAME_ONLY
               LIBRARY_NAME ICARUS-raw-data_CRT
               LIB_LIBRARIES
               artdaq-core_Data
  ${ARTDAQ_DAQDATA}
  ${CETLIB_EXCEPT}
  ${Boost_SYSTEM_LIBRARY}
  ${FHICLCPP}
  ${MF_MESSAGELOGGER}
  ${MF_MESSAGEUTILITIES}
  ${RCEDAMLIB}
  ${GALLERY_LIB}
  gallery
  canvas
  ${ROOT_BASIC_LIB_LIST}
  ${ART_LIB}
  art_Framework_Principal
  pthread
  z
  )
install_headers()
install_source()
