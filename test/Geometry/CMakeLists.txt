
# geometry test on ICARUS geometries

cet_transitive_paths(FHICL_DIR BINARY IN_TREE)
cet_test_env_prepend(FHICL_FILE_PATH "." ${TRANSITIVE_PATHS_WITH_FHICL_DIR} ${CMAKE_CURRENT_LIST_DIR} ${icarusalg_BINARY_DIR}/job)
cet_transitive_paths(LIBRARY_DIR BINARY IN_TREE)
cet_test_env_prepend(CET_PLUGIN_PATH ${TRANSITIVE_PATHS_WITH_LIBRARY_DIR})
cet_transitive_paths(GDML_DIR BINARY IN_TREE)
cet_test_env_prepend(FW_SEARCH_PATH ${TRANSITIVE_PATHS_WITH_GDML_DIR})
# these tests run a FCL file and fail only if lar exits with a bad exit code;
# they are fast and always run.
cet_test(geometry_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_geometry_icarus.fcl
#  OPTIONAL_GROUPS Broken
)

cet_test(dump_channel_map_icarus_test HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./dump_icarus_channelmap.fcl
)

install_fhicl()
