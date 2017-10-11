/**
 * @file   geometry_unit_test_icarus.h
 * @brief  Class for objects initializing ICARUS geometry
 * @date   October 10, 2017
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * 
 * Provides an environment for easy set up of ICARUS-aware tests.
 * Keep in mind that the channel mapping algorithm must be hard-coded and, if
 * using Boost unit test, the configuration file location must be hard coded too
 * (or you can use the provided configuration).
 * 
 * For an example of usage, see icaruscode/test/Geometry/geometry_test_icarus.cxx
 */

#ifndef TEST_GEOMETRY_UNIT_TEST_ICARUS_H
#define TEST_GEOMETRY_UNIT_TEST_ICARUS_H

// LArSoft libraries
#include "test/Geometry/geometry_unit_test_base.h"
#include "icaruscode/Geometry/ChannelMapIcarusAlg.h"

// C/C++ standard libraries
#include <string>


namespace icarus {
  
  /// Namespace including ICARUS-specific testing
  namespace testing {
    
    /** ************************************************************************
     * @brief Class holding the configuration for a ICARUS fixture
     * @tparam CHANNELMAP the class used for channel mapping
     * @see BasicGeometryEnvironmentConfiguration
     *
     * This class needs to be fully constructed by the default constructor
     * in order to be useful as Boost unit test fixture.
     * It is supposed to be passed as a template parameter to another class
     * that can store an instance of it and extract configuration information
     * from it.
     * 
     * This class should be used with ChannelMapIcarusAlg.
     * 
     * We reuse BasicGeometryEnvironmentConfiguration as base class and then we
     * fix its setup.
     */
    template <typename CHANNELMAP = geo::ChannelMapIcarusAlg>
    struct IcarusGeometryEnvironmentConfiguration:
      public ::testing::BasicGeometryEnvironmentConfiguration<CHANNELMAP>
    {
      // remember that BasicGeometryEnvironmentConfiguration is not polymorphic
      using base_t
        = ::testing::BasicGeometryEnvironmentConfiguration<CHANNELMAP>;
      
      /// Default constructor
      IcarusGeometryEnvironmentConfiguration() { IcarusDefaultInit(); }
      
      /// Constructor; accepts the name as parameter
      IcarusGeometryEnvironmentConfiguration(std::string name):
        IcarusGeometryEnvironmentConfiguration()
        { base_t::SetApplicationName(name); }
      
      
        private:
      void IcarusDefaultInit()
        {
          // overwrite the configuration that happened in the base class:
          base_t::SetApplicationName("IcarusGeometryTest");
          base_t::SetDefaultGeometryConfiguration(R"(
            SurfaceY:          130e2 # in cm, vertical distance to the surface
            Name:             "icarus"
            GDML:             "test_27Sep.gdml"
            ROOT:             "test_27Sep.gdml"
            SortingParameters: {}
            )");
        }
    }; // class IcarusGeometryEnvironmentConfiguration<>
    
    
  } // namespace testing
} // namespace icarus

#endif // TEST_GEOMETRY_UNIT_TEST_ICARUS_H
