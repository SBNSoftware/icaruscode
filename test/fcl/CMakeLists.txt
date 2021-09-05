
# enable this test when all FHiCL configuration in ICARUS is fixed

cet_test( testFHiCLfiles_icarus.sh PREBUILT
  TEST_ARGS "--quiet" "--no-color" "--exclude=[Dd]eprecated" "--exclude=[Aa]rchived?" "--exclude-from=SkipTesting.txt"
  DATAFILES testFHiCLfiles.sh "${PROJECT_SOURCE_DIR}/fcl/SkipTesting.txt"
  )

install_scripts()
