/**
 * @file    icaruscode/Utilities/RandFastGauss.h
 * @brief   Fast approximate Gaussian random translator.
 * @date    February 15, 2020
 */

#ifndef ICARUSCODE_UTILITIES_RANDFASTGAUS_H
#define ICARUSCODE_UTILITIES_RANDFASTGAUS_H


// ICARUS libraries
#include "icaruscode/Utilities/FastAndPoorGauss.h"

// CLHEP
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RandomEngine.h"

// C/C++ standard libraries
#include <memory> // std::shared_ptr


// -----------------------------------------------------------------------------
namespace util { class RandFastGauss; }

/**
 * @brief Normal distribution focussing on speed.
 * 
 * The random number is generated via `util::FastAndPoorGauss<double>`.
 * 
 * @note This class is incomplete.
 */
class util::RandFastGauss: public CLHEP::HepRandom {

    public:

  /// Constructor: borrows an engine but does not manage it.
  RandFastGauss(
    CLHEP::HepRandomEngine& anEngine,
    double mean = 0.0, double stdDev = 1.0
    );
  
  /**
   * @brief Constructor: takes ownership of the engine.
   * 
   * The ownership of the specified engine is transferred to this object, which
   * will dispose of it at the end of its life.
   */
  RandFastGauss(
    CLHEP::HepRandomEngine* anEngine,
    double mean = 0.0, double stdDev = 1.0
    );

  // --- BEGIN -- Not implemented ----------------------------------------------
  static double shoot();

  static double shoot(double mean, double stdDev);

  static void shootArray
    (const int size, double* vect, double mean = 0.0, double stdDev = 1.0);

  static double shoot(CLHEP::HepRandomEngine* anEngine);

  static double shoot
    (CLHEP::HepRandomEngine* anEngine, double mean, double stdDev);

  static void shootArray(CLHEP::HepRandomEngine* anEngine,
    const int size, double* vect, double mean = 0.0, double stdDev = 1.0);

  // --- END -- Not implemented ------------------------------------------------
  
  /// Extracts a single normal value under the default distribution.
  double fire() { return fTransform(normal()); }

  /// Extracts a single normal value under the specified normal distribution.
  double fire(double mean, double stdDev) { return normal() * stdDev + mean; }
  
  // --- BEGIN -- Not implemented ----------------------------------------------
  void fireArray(const int size, double* vect);
  void fireArray(const int size, double* vect, double mean, double stdDev);
  // --- END -- Not implemented ------------------------------------------------

  /// Extracts a single normal value under the default distribution.
  virtual double operator()() override { return fire(); }
  
  /// Extracts a single normal value under the specified normal distribution.
  double operator()(double mean, double stdDev) { return fire(mean, stdDev); }

  /// Returns the name of the distribution.
  virtual std::string name() const override { return distributionName(); }
  
  /// Returns the default random generator engine.
  virtual CLHEP::HepRandomEngine& engine() override { return *localEngine; }

  /// Returns the name of the distribution.
  static std::string distributionName() { return "RandFastGauss"; }
  
  // --- BEGIN -- Not implemented ----------------------------------------------
  virtual std::ostream& put(std::ostream& os) const override
    { return os; }
  virtual std::istream& get(std::istream& is) override
    { return is; }
  
  // Methods overriding the base class static saveEngineStatus ones,
  // by adding extra data so that save in one program, then further gaussians,
  // will produce the identical sequence to restore in another program, then 
  // generating gaussian randoms there 

  static void saveEngineStatus( const char filename[] = "Config.conf" );
  // Saves to file the current status of the current engine.

  static void restoreEngineStatus( const char filename[] = "Config.conf" );
  // Restores a saved status (if any) for the current engine.

  static std::ostream& saveFullState ( std::ostream & os );
  // Saves to stream the state of the engine and cached data.

  static std::istream& restoreFullState ( std::istream & is );
  // Restores from stream the state of the engine and cached data.

  static std::ostream& saveDistState ( std::ostream & os );
  // Saves to stream the state of the cached data.

  static std::istream& restoreDistState ( std::istream & is );
  // Restores from stream the state of the cached data.
  // --- END -- Not implemented ------------------------------------------------


protected:
  
  util::GaussianTransformer<double> fTransform;

  std::shared_ptr<CLHEP::HepRandomEngine> localEngine;

  double normal() { return fToGauss(localEngine->flat()); }

private:
  
  /// Translates uniform number in [ 0, 1 ] into a Gaussian number.
  util::FastAndPoorGauss<32768U> fToGauss; // TODO make it static
  

}; // util::RandFastGauss


// -----------------------------------------------------------------------------
// ---  inline implementation
// -----------------------------------------------------------------------------
util::RandFastGauss::RandFastGauss(
  CLHEP::HepRandomEngine& anEngine,
  double mean /* = 0.0 */, double stdDev /* = 1.0 */
)
  : HepRandom()
  , fTransform(mean, stdDev)
  , localEngine(&anEngine, [](void const*){})
  {}


// -----------------------------------------------------------------------------
util::RandFastGauss::RandFastGauss(
  CLHEP::HepRandomEngine* anEngine,
  double mean /* = 0.0 */, double stdDev /* = 1.0 */
)
  : HepRandom()
  , fTransform(mean, stdDev)
  , localEngine(anEngine)
  {}


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_RANDFASTGAUS_H
