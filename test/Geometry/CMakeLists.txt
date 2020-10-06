
# geometry test on ICARUS geometries

# these tests run a FCL file and fail only if lar exits with a bad exit code;
# they are fast and always run.
# FCL files need to be copied to the test area (DATAFILES directive) since they
# are not installed.
cet_test(geometry_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_geometry_icarus.fcl
#  OPTIONAL_GROUPS Broken
)

# FCL files need to be copied to the test area (DATAFILES directive) since they
# are not installed.
cet_test(dump_channel_map_icarus_test HANDBUILT
  DATAFILES dump_icarus_channelmap.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./dump_icarus_channelmap.fcl
)

# geometry dump for geometry _single_induction_nooverburden
cet_test(geometry_single_induction_nooverburden_geometry_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_single_induction_nooverburden_icarus_geometry.fcl
)

# geometry test for geometry _single_induction_nooverburden
cet_test(geometry_single_induction_nooverburden_geometry_test_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_geometry_single_induction_nooverburden_icarus.fcl
)

# channel mapping dump for geometry _single_induction_nooverburden
cet_test(geometry_single_induction_nooverburden_channel_mapping_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_single_induction_nooverburden_icarus_channelmap.fcl
)

# geometry dump for geometry _single_induction_overburden
cet_test(geometry_single_induction_overburden_geometry_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_single_induction_overburden_icarus_geometry.fcl
)

# geometry test for geometry _single_induction_overburden
cet_test(geometry_single_induction_overburden_geometry_test_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_geometry_single_induction_overburden_icarus.fcl
)

# channel mapping dump for geometry _single_induction_overburden
cet_test(geometry_single_induction_overburden_channel_mapping_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_single_induction_overburden_icarus_channelmap.fcl
)

# geometry dump for geometry _single_induction
cet_test(geometry_single_induction_geometry_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_single_induction_icarus_geometry.fcl
)

# geometry test for geometry _single_induction
cet_test(geometry_single_induction_geometry_test_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_geometry_single_induction_icarus.fcl
)

# channel mapping dump for geometry _single_induction
cet_test(geometry_single_induction_channel_mapping_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_single_induction_icarus_channelmap.fcl
)

# geometry dump for geometry _split_induction_nooverburden
cet_test(geometry_split_induction_nooverburden_geometry_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_split_induction_nooverburden_icarus_geometry.fcl
)

# geometry test for geometry _split_induction_nooverburden
cet_test(geometry_split_induction_nooverburden_geometry_test_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_geometry_split_induction_nooverburden_icarus.fcl
)

# channel mapping dump for geometry _split_induction_nooverburden
cet_test(geometry_split_induction_nooverburden_channel_mapping_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_split_induction_nooverburden_icarus_channelmap.fcl
)

# geometry dump for geometry _split_induction_overburden
cet_test(geometry_split_induction_overburden_geometry_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_split_induction_overburden_icarus_geometry.fcl
)

# geometry test for geometry _split_induction_overburden
cet_test(geometry_split_induction_overburden_geometry_test_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_geometry_split_induction_overburden_icarus.fcl
)

# channel mapping dump for geometry _split_induction_overburden
cet_test(geometry_split_induction_overburden_channel_mapping_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_split_induction_overburden_icarus_channelmap.fcl
)

# geometry dump for geometry _split_induction
cet_test(geometry_split_induction_geometry_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_split_induction_icarus_geometry.fcl
)

# geometry test for geometry _split_induction
cet_test(geometry_split_induction_geometry_test_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_geometry_split_induction_icarus.fcl
)

# channel mapping dump for geometry _split_induction
cet_test(geometry_split_induction_channel_mapping_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_split_induction_icarus_channelmap.fcl
)

# geometry dump for geometry _nooverburden
cet_test(geometry_nooverburden_geometry_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_nooverburden_icarus_geometry.fcl
)

# geometry test for geometry _nooverburden
cet_test(geometry_nooverburden_geometry_test_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_geometry_nooverburden_icarus.fcl
)

# channel mapping dump for geometry _nooverburden
cet_test(geometry_nooverburden_channel_mapping_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_nooverburden_icarus_channelmap.fcl
)

# geometry dump for geometry _overburden
cet_test(geometry_overburden_geometry_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_overburden_icarus_geometry.fcl
)

# geometry test for geometry _overburden
cet_test(geometry_overburden_geometry_test_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config test_geometry_overburden_icarus.fcl
)

# channel mapping dump for geometry _overburden
cet_test(geometry_overburden_channel_mapping_dump_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config dump_overburden_icarus_channelmap.fcl
)


install_headers()
install_fhicl()
install_source()
