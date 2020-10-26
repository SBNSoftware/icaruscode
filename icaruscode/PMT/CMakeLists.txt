add_subdirectory(OpticalTools)
add_subdirectory(OpReco)
add_subdirectory(Algorithms)
add_subdirectory(LibraryMappingTools)
add_subdirectory(Trigger)
add_subdirectory(Decode)
add_subdirectory(Calo)

# Removing AVX instructions because assuming the grid nodes
# being less than 6 year old proved to be pretentious.
# AVX is not so relevant after the last commit anyway
# cet_set_compiler_flags(EXTRA_CXX_FLAGS -mavx)

simple_plugin(SimPMTIcarus "module"
  icaruscode_PMT_Algorithms
  lardataobj_RawData
  lardataobj_Simulation
  nurandom_RandomUtils_NuRandomService_service
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  art_Persistency_Provenance
  ${CANVAS}
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${CETLIB}
  cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  ${Boost_SYSTEM_LIBRARY}
  )

simple_plugin(PMTWaveformBaselines "module"
  icaruscode_PMT_Algorithms
  sbnobj_ICARUS_PMT_Data
  lardataobj_RawData
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  )


simple_plugin(AsymGaussPulseFunctionTool "tool"
  icaruscode_PMT_Algorithms
  )

simple_plugin(AsymExpPulseFunctionTool "tool"
  icaruscode_PMT_Algorithms
  )

simple_plugin(CustomPulseFunctionTool "tool"
  icaruscode_PMT_Algorithms
  ROOT::Hist
  )


art_make(
  EXCLUDE
    "SimPMTIcarus_module.cc"
    "PMTWaveformBaselines_module.cc"
    "AsymGaussPulseFunctionTool_tool.cc"
    "AsymExpPulseFunctionTool_tool.cc"
    "CustomPulseFunctionTool_tool.cc"
          LIB_LIBRARIES 
            lardataobj_RawData
            lardataobj_RecoBase
            lardata_Utilities
            ${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
			${Boost_SYSTEM_LIBRARY}
	        MODULE_LIBRARIES larcorealg_Geometry
	  		larcore_Geometry_Geometry_service
	        lardata_Utilities
			larevt_Filters
			lardataobj_RawData
			larevt_CalibrationDBI_IOVData
			larevt_CalibrationDBI_Providers
			lardataobj_RecoBase
			lardata_ArtDataHelper
			nusimdata_SimulationBase
            nurandom_RandomUtils_NuRandomService_service
			${CLHEP}
            ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	        ${FHICLCPP}
			${CETLIB}
			cetlib_except
            ${ROOT_GEOM}
            ${ROOT_XMLIO}
            ${ROOT_GDML}
			${ROOT_FFTW}
			${ROOT_BASIC_LIB_LIST}
			${Boost_SYSTEM_LIBRARY}
        )

install_headers()
install_fhicl()
install_source()
install_scripts(
  SUBDIRS scripts
  EXTRAS
    "scripts/MergePhotonLibrary.C"
    "scripts/RecompressROOTtrees.C"
  )

