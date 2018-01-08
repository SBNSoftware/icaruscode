
# geometry test on ICARUS geometries

# these tests run a FCL file and fail only if lar exits with a bad exit code;
# they are fast and always run.
# FCL files need to be copied to the test area (DATAFILES directive) since they
# are not installed.
cet_test(geometry_icarus HANDBUILT
  DATAFILES test_geometry_icarus.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./test_geometry_icarus.fcl
#  OPTIONAL_GROUPS Broken
)


#
# the following are tests equivalent to the module ones, but with no framework
#
# unit test (use the hard-coded channel mapping for ICARUS geometry)
cet_test(geometry_icarus_test
  SOURCES geometry_icarus_test.cxx
  DATAFILES test_geometry_icarus.fcl
  TEST_ARGS ./test_geometry_icarus.fcl
  LIBRARIES icarus_Geometry
            larcorealg_Geometry
            GeometryTestLib
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            cetlib_except
	    ${ROOT_CORE}
#  OPTIONAL_GROUPS Broken
)


# unit test (use the hard-coded configuration for ICARUS geometry)
# this uses BOOST for the test
cet_test(geometry_iterator_icarus_test
  SOURCES geometry_iterator_icarus_test.cxx
  DATAFILES ./test_geometry_iterators_icarus.fcl
  TEST_ARGS ./test_geometry_iterators_icarus.fcl
  LIBRARIES icarus_Geometry
            larcorealg_Geometry
            GeometryTestLib
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            cetlib_except
	    ${ROOT_CORE}
  USE_BOOST_UNIT
)


# unit test (use the hard-coded configuration for ICARUS geometry)
cet_test(geometry_iterator_loop_icarus_test
  SOURCES geometry_iterator_loop_icarus_test.cxx
  DATAFILES ./test_geometry_iterators_icarus.fcl
  TEST_ARGS ./test_geometry_iterators_icarus.fcl
  LIBRARIES icarus_Geometry
            larcorealg_Geometry
            GeometryTestLib
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            cetlib_except
	    ${ROOT_CORE}
)


# FCL files need to be copied to the test area (DATAFILES directive) since they
# are not installed.
cet_test(dump_channel_map_icarus_test HANDBUILT
  DATAFILES dump_icarus_channelmap.fcl
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ./dump_icarus_channelmap.fcl
)

