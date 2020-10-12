# ======================================================================
#  icaruscode main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsoft/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 3.12-3.14
cmake_minimum_required (VERSION 3.12...3.14 FATAL_ERROR)

project(icaruscode LANGUAGES CXX)

project(icaruscode)

message(STATUS
  "\n-- ============================================================================="
  "\n-- === *** Start of configuration of ${PROJECT_NAME} ***"
  "\n-- ===")

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH
  $ENV{CANVAS_ROOT_IO_DIR}/Modules
  $ENV{CANVAS_DIR}/Modules
  $ENV{CETBUILDTOOLS_DIR}/Modules
  ${CMAKE_MODULE_PATH}
  )

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -I $ENV{CLANG_FQ_DIR}/include/c++/v1
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larcore )
find_ups_product( larcoreobj )
find_ups_product( lardata )
find_ups_product( lardataobj )
find_ups_product( larevt )
find_ups_product( larsim )
find_ups_product( larreco )
find_ups_product( larpandora )
find_ups_product( larana )
find_ups_product( larpandora )
find_ups_product( nutools )
find_ups_product( nug4 )
find_ups_product( nurandom )
find_ups_product( art )
find_ups_product( art_root_io)
find_ups_product( cetbuildtools )
find_ups_product( libwda )
find_ups_product( postgresql )
find_ups_product( eigen )
find_ups_product( fftw )
find_ups_boost( )
find_ups_product( ifdh_art )
find_ups_product( guideline_sl )
find_ups_product( libwda )
find_ups_product( sbnobj )
find_ups_product( sbncode )
find_ups_product( icarus_signal_processing )
find_ups_product( icarusalg )
find_ups_product( icarusutil )
find_ups_product( icarus_data )

find_ups_product( artdaq_core )
find_ups_product( sbndaq_artdaq_core )

cet_find_library( PQ NAMES pq PATHS ENV POSTGRESQL_LIBRARIES NO_DEFAULT_PATH )
cet_find_library( NURANDOM_RANDOMUTILS_NURANDOMSERVICE_SERVICE 
                  NAMES nurandom_RandomUtils_NuRandomService_service
                  PATHS ENV NURANDOM_LIB NO_DEFAULT_PATH )

# deal with fftw
cet_find_library(ICARUS_FFTW3_LIBRARY  NAMES fftw3  PATHS $ENV{FFTW_LIBRARY} )
cet_find_library(ICARUS_FFTW3F_LIBRARY NAMES fftw3f PATHS $ENV{FFTW_LIBRARY} )
cet_find_library(ICARUS_FFTW3L_LIBRARY NAMES fftw3l PATHS $ENV{FFTW_LIBRARY} )
cet_find_library(ICARUS_FFTW3Q_LIBRARY NAMES fftw3q PATHS $ENV{FFTW_LIBRARY} )
set(ICARUS_FFTW_LIBRARIES ${ICARUS_FFTW3_LIBRARY} ${ICARUS_FFTW3F_LIBRARY} ${ICARUS_FFTW3L_LIBRARY} ${ICARUS_FFTW3Q_LIBRARY})

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(icaruscode)
add_subdirectory(scripts)
add_subdirectory(doc)

# subdirectory for test code
add_subdirectory(test)

# Top level fcl files.
add_subdirectory(fcl)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)

message(STATUS
  "\n-- ==="
  "\n-- === *** End of configuration of ${PROJECT_NAME} ***"
  "\n-- ============================================================================="
  "\n-- "
  )

