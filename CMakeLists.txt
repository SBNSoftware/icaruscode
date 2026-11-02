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


cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
# cmake_policy(VERSION 3.19)

find_package(cetmodules 3.24.01 REQUIRED)
#project(icaruscode VERSION 10.25.00 LANGUAGES CXX)

set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 10.27.00)
project(icaruscode LANGUAGES CXX)

message(STATUS
  "\n-- ============================================================================="
  "\n-- === *** Start of configuration of ${PROJECT_NAME} ${${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING} ***"
  "\n-- ==="
  )

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic
              $<$<CXX_COMPILER_ID:GNU>:-Wno-maybe-uninitialized>
              $<$<CXX_COMPILER_ID:Clang>:-Wno-return-stack-address>
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package(nufinder REQUIRED ) # necessary for CMake to include sbncode

# in addition to the ones already EXPORTed by icarusalg
# (CLHEP/ROOT/GSL, CET/canvas, nutools/LArSoft/SBN data and algorithms)

find_package(Boost COMPONENTS system program_options REQUIRED)
find_package(FFTW3 REQUIRED )
find_package(FFTW3f REQUIRED )
find_package(FFTW3q REQUIRED )
find_package(FFTW3l REQUIRED )
find_package(Range-v3 REQUIRED )
find_package(Eigen3 REQUIRED )
find_package(SQLite3 REQUIRED )
find_package(hep_concurrency REQUIRED )
find_package(TorchScatter REQUIRED )

find_package(libwda REQUIRED )
find_package(ifdhc REQUIRED )
find_package(art REQUIRED )
find_package(art_root_io REQUIRED )
find_package(messagefacility REQUIRED )
find_package(ifdh_art REQUIRED )

find_package(artdaq_core REQUIRED )

find_package(nurandom REQUIRED )
find_package(nug4 REQUIRED )

find_package(larcore REQUIRED )
find_package(lardata REQUIRED )
find_package(larsim REQUIRED )
find_package(larevt REQUIRED )
find_package(larreco REQUIRED )
find_package(larrecodnn REQUIRED )
find_package(larana REQUIRED )
find_package(larpandoracontent REQUIRED )
find_package(larpandora REQUIRED )

find_package(sbnobj REQUIRED )
find_package(sbnalg REQUIRED )
find_package(sbndaq_artdaq_core REQUIRED )
find_package(sbncode REQUIRED )

find_package(WireCell)

find_package(icarusalg REQUIRED )
find_package(icarus_signal_processing REQUIRED )

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

message(STATUS
  "\n-- ==="
  "\n-- === *** End of configuration of ${PROJECT_NAME} ***"
  "\n-- ============================================================================="
  "\n-- "
  )

cet_cmake_config()
