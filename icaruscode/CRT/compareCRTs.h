#ifndef ICARUS_compareCRTs_h
#define ICARUS_compareCRTs_h

#include <utility>

namespace icarus{

  //----------------------------------------------------------------------------
  // Sorting CRT co-ordinates: decreasing vertical coordinate (Y),
  // next increasing beam coordinate (Z), and next increasing drift direction (X).

  template <typename T>
  bool compareCRTs(T const& t1, T const& t2)
  {
    // Sort by decreasing y, then increasing z, then increasing x
    auto const& [c1, c2] = std::pair{t1.GetCenter(), t2.GetCenter()};
    if (c1.Y() != c2.Y()) return c1.Y() > c2.Y();
    if (c1.Z() != c2.Z()) return c1.Z() < c2.Z();
    return c1.X() < c2.X();
  }

} //namespace icarus

#endif // ICARUS_compareCRTs_h
