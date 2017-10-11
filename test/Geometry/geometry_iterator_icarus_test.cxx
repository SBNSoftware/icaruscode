/**
 * @file   geometry_iterator_icarus_test.cxx
 * @brief  Unit test for geometry iterators on ICARUS detector.
 * @date   October 10, 2017
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * 
 * Usage: just run the executable.
 * Boost unit testing environment keeps the arguments secret anyway.
 */

// Boost test libraries; defining this symbol tells boost somehow to generate
// a main() function; Boost is pulled in by geometry_boost_unit_test_base.h
#define BOOST_TEST_MODULE GeometryIteratorTestICARUS

// ICARUS libraries
#include "icaruscode/Geometry/ChannelMapIcarusAlg.h"

// LArSoft libraries
#include "test/Geometry/geometry_unit_test_icarus.h"
#include "test/Geometry/GeometryIteratorTestAlg.h"
#include "larcorealg/TestUtils/boost_unit_test_base.h"

//------------------------------------------------------------------------------
//---  The test environment
//---

// we define here all the configuration that is needed;
// in the specific, the type of the channel mapping and a proper test name,
// used for output only; we use IcarusGeometryEnvironmentConfiguration
// as base class, that is already configured to use ICARUS geometry.
// We wrap it in testing::BoostCommandLineConfiguration<> so that it can learn
// the configuration file name from the command line.
struct IcarusGeometryConfiguration:
  public testing::BoostCommandLineConfiguration<
    icarus::testing::IcarusGeometryEnvironmentConfiguration
      <geo::ChannelMapIcarusAlg>
    >
{
  /// Constructor: overrides the application name; ignores command line
  IcarusGeometryConfiguration()
    { SetApplicationName("GeometryIteratorUnitTest"); }
}; // class IcarusGeometryConfiguration

/*
 * Our fixture is based on GeometryTesterFixture, configured with the object
 * above.
 * It provides:
 * - `Tester`, a configured instance of the test algorithm.
 */
class IcarusGeometryIteratorTestFixture:
  private testing::GeometryTesterEnvironment<IcarusGeometryConfiguration>
{
    public:
  geo::GeometryIteratorTestAlg Tester;
  
  /// Constructor: initialize the tester with the Geometry from base class
  IcarusGeometryIteratorTestFixture(): Tester(TesterParameters())
    { Tester.Setup(*Geometry()); }

}; // class IcarusGeometryIteratorTestFixture



//------------------------------------------------------------------------------
//---  The tests
//---

BOOST_FIXTURE_TEST_SUITE
  (GeometryIteratorsIcarus, IcarusGeometryIteratorTestFixture)
// BOOST_GLOBAL_FIXTURE(IcarusGeometryIteratorTestFixture);


BOOST_AUTO_TEST_CASE( AllTests )
{
  Tester.Run();
} // BOOST_AUTO_TEST_CASE( AllTests )

/*
BOOST_AUTO_TEST_CASE( CryostatIteratorsTest )
{
  Tester.CryostatIteratorsTest();
} // BOOST_AUTO_TEST_CASE( CryostatIteratorsTest )



BOOST_AUTO_TEST_CASE( TPCIteratorsTest )
{
  Tester.TPCIteratorsTest();
} // BOOST_AUTO_TEST_CASE( TPCIteratorsTest )



BOOST_AUTO_TEST_CASE( PlaneIteratorsTest )
{
  Tester.PlaneIteratorsTest();
} // BOOST_AUTO_TEST_CASE( PlaneIteratorsTest )



BOOST_AUTO_TEST_CASE( WireIteratorsTest )
{
  Tester.WireIteratorsTest();
} // BOOST_AUTO_TEST_CASE( WireIteratorsTest )
*/

BOOST_AUTO_TEST_SUITE_END()

