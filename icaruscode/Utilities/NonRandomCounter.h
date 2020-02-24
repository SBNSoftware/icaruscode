/**
 * @file   icaruscode/Utilities/NonRandomCounter.h
 * @brief  Non-random number engine for profiling purposes.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 28, 2020
 *
 */

#ifndef ICARUSCODE_UTILITIES_NONRANDOMCOUNTER_H
#define ICARUSCODE_UTILITIES_NONRANDOMCOUNTER_H

// CLHEP
#include "CLHEP/Random/defs.h"
#include "CLHEP/Random/RandomEngine.h"

// C/C++ standard library
#include <istream>
#include <string>
#include <limits> // std::numeric_limits<>


// -----------------------------------------------------------------------------
namespace util { class NonRandomCounter; }

/**
 * @brief Fast random engine which returns sequential numbers.
 *
 * This generator does **not** produce pseudorandom numbers.
 * It is used only as a replacement of a real engine for profiling purposes.
 * 
 * The range of values spans [ 0, 1 [.
 */
class util::NonRandomCounter: public CLHEP::HepRandomEngine {

    public:

  NonRandomCounter() = default;

  NonRandomCounter(long seed): count(seed) {}

  NonRandomCounter(std::istream& is): NonRandomCounter(readLong(is)) {}

  virtual double flat() override { return doFlat(); }

  virtual void flatArray(const int size, double* vect) override
    { double* end = vect + size; while (vect != end) *(vect++) = doFlat(); }

  virtual void setSeed(long seed, int) override
    { count = static_cast<unsigned long>(seed); }

  virtual void setSeeds(long const* seeds, int _) override
    { setSeed(*seeds, _); }

  virtual void saveStatus(const char filename[] = "NonRandomCounter.conf") const
    override;

  virtual void restoreStatus(const char filename[] = "NonRandomCounter.conf")
    override;

  virtual void showStatus() const override;

  virtual std::string name() const override { return "NonRandomCounter"; }

    private:

  unsigned long count = 0U;

  double doFlat()
    {
      return static_cast<double>(++count)
        / std::numeric_limits<unsigned long>::max();
    }


  static long readLong(std::istream& is) { unsigned long l; is >> l; return l; }

}; // util::NonRandomCounter


// -----------------------------------------------------------------------------
// ---  inline implementation
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_NONRANDOMCOUNTER_H
