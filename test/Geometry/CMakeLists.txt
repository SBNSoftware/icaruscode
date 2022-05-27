
# geometry test on ICARUS geometries

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
