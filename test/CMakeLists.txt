# cmake driver file for testing from CET build tools
include(CetTest)

# Always enable asserts (for tests only)
cet_enable_asserts()

# test directories
add_subdirectory(Geometry)
add_subdirectory(fcl)
add_subdirectory(PMT)
add_subdirectory(Decode)

# Continuous Integration tests
add_subdirectory(ci)

