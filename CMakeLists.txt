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


cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

find_package(cetmodules 3.19.02 REQUIRED)
project(icaruscode VERSION 09.60.01 LANGUAGES CXX)

message(STATUS
  "\n-- ============================================================================="
  "\n-- === *** Start of configuration of ${PROJECT_NAME} ***"
  "\n-- ===")

# cetbuildtools contains our cmake modules
#find_package(cetbuildtools REQUIRED)

#list(APPEND CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules $ENV{ART_DIR}/Modules )

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -I $ENV{CLANG_FQ_DIR}/include/c++/v1
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

#these are minimum required versions,not thea ctual product versions
find_package(larsoft)
find_package(libwda)
find_package(PostgreSQL)
find_package(Eigen3)
find_package(ifbeam)
find_package(ifdhc)
find_package(ifdh_art)
find_package(guideline_sl)
find_package(range-v3)
find_package(sbnobj)
find_package(sbndaq-artdaq-core)
find_package(icarus_signal_processing)
find_package(icarusalg)
find_package(icarusutil)
find_package(icarus_data)
find_package(artdaq_core)
find_package(larcoreobj)
find_package(larcorealg)
find_package(lardataobj)
find_package(lardataalg)
find_package(larcore)
find_package(lardata)
find_package(larsim)
find_package(larevt)
find_package(larreco)
find_package(larana)
find_package(larpandora)
find_package(larpandoracontent)
find_package(nusimdata)
find_package(nurandom)
find_package(nug4)
find_package(art)
find_package(canvas)
find_package(messagefacility)
find_package(FFTW3)
find_package(libwda)
find_package(ifdh_art)
find_package(ROOT)
find_package(Boost COMPONENTS system program_options REQUIRED)
find_package(protobuf)
find_package(PandoraSDK)
find_package(PandoraMonitoring)
find_package(SQLite3)
find_package(artdaq_core)
find_package(sbndaq_artdaq_core)
find_package(sbndaq-artdaq-core)
find_package(PostgreSQL)

cet_find_library(PQ NAMES pq )

#deal with fftw
find_package( PkgConfig REQUIRED)
pkg_search_module(FFTW3 REQUIRED IMPORTED_TARGET GLOBAL fftw3f fftw3 fftw3q fftw3l fftw3_threads fftw3f_threeads fftw3q_threads)
set(ICARUS_FFTW_LIBRARIES PkgConfig::FFTW3)


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
#add_subdirectory(ups)

# packaging utility
#include(UseCPack)

message(STATUS
  "\n-- ==="
  "\n-- === *** End of configuration of ${PROJECT_NAME} ***"
  "\n-- ============================================================================="
  "\n-- SQLITE3: ${SQLITE3}"
  "\n-- "
  )

cet_cmake_config()
