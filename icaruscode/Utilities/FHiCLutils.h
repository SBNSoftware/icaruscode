/**
 * @file   icaruscode/Utilities/FHiCLutils.h
 * @brief  Plots to inform trigger design decisions.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   August 8, 2019
 * 
 * This library is header only.
 */

#ifndef ICARUSCODE_UTILITIES_FHICLUTILS_H
#define ICARUSCODE_UTILITIES_FHICLUTILS_H

// C/C++ standard libraries
#include <vector>


namespace util::fhicl {
  
  //--------------------------------------------------------------------------
  /*
   * Currently util::Quantity objects have serious issues with FHiCL validation,
   * so we are not reading them directly as such. This pulls in a number of
   * workarounds for:
   *  * every `Quantity` parameter: its FHiCL parameter must be declared as
   *    the base type `Quantity::value_t`
   *  * every optional parameter must be then read indirectly since a reference
   *    to the exact type of the parameter is required in such reading
   *  * for sequences, direct vector assignment
   *    (`vector<Quantity> = vector<Quantity::value_t>`) won't work, and since
   *    `Sequence::operator()` returns a temporary, this needs to be wraped too
   */
  
  /**
   * @brief Helper to return a converted sequence from FHiCL configuration.
   * @tparam SeqValueType type of the sequence value being extracted
   * 
   * This class can be used to apply implicit conversion of a `fhicl::Sequence`
   * of objects into a `std::vector` of another type, as sometimes
   * `fhicl::Sequence` does not support the desired type directly.
   * For example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct MyData {
   *   int fValue;
   *   
   *   MyData(int value): fValue(value) {}
   * };
   * 
   * struct Config {
   *   
   *   fhicl::Sequence<int> Data {
   *     fhicl::Name("Data"),
   *     fhicl::Comment("The data"),
   *     { 0, 1, 2 }
   *     };
   *   
   * };
   * 
   * std::vector<MyData> fData;
   * 
   * MyClass(Config const& config)
   *   : fData(util::FHiCLsequenceWrapper(config.Data()))
   *   {}
   * 
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename SeqValueType>
  struct SequenceWrapper;
  
  
  //--------------------------------------------------------------------------
  /**
   * @brief Returns the value of an optional parameter as `std::optional`.
   * @tparam Optional FHiCL optional class (e.g. `fhicl::OptionalAtom`)
   * @param parameter the optional FHiCL parameter
   * @return the value of the parameter
   *
   * This utility allows to carry the information whether an optional parameter
   * was specified or not.
   * It may also help with single-line initialization of data members from
   * FHiCL optional parameters, e.g.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct Config {
   *
   *   fhicl::Atom<unsigned int> A {
   *     fhicl::Name("A"),
   *     fhicl::Comment("parameter A")
   *     };
   *
   *   fhicl::OptionalAtom<unsigned int> B {
   *     fhicl::Name("B"),
   *     fhicl::Comment("parameter B (default: same as A)")
   *     };
   *
   * };
   *
   * unsigned int fA, fB;
   *
   * MyClass(Config const& config)
   *   : fA(config.A())
   *   , fB(util::fhicl::getOptionalValue(config.B).value_or(fA))
   *   {}
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * This feature will be rendered obsolete by Redmine Issue #23653.
   *
   */
  template <typename Optional>
  std::optional<typename Optional::value_type> getOptionalValue
    (Optional const& parameter);


  //--------------------------------------------------------------------------

  /**
   * @brief Returns the value of an optional parameter, or a default value.
   * @tparam T type of the value being returned
   * @tparam Optional FHiCL optional class (e.g. `fhicl::OptionalAtom`)
   * @param parameter the optional FHiCL parameter
   * @param defValue default value to be returned
   * @return the value of the parameter, or `defValue` if absent
   *
   * It is usually preferable to use a non-optional FHiCL parameter with a
   * default value, but some types are not well suited for this.
   * Also, this function does not buy much unless the parameter reading needs to
   * be done in a single statement, like in the initialization list of a
   * constructor:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * struct Config {
   *
   *   fhicl::OptionalSequenceTable<BoardCfg> Boards {
   *     fhicl::Name("MyStruct"),
   *     fhicl::Comment("data configuration elements")
   *     };
   *
   * };
   *
   * using AllBoardCfg_t = std::vector<BoardCfg>;
   *
   * AllBoardCfg_t fBoards;
   *
   *
   * MyClass(Config const& config)
   *   : fBoards
   *     (util::fhicl::getOptionalValue<BoardCfg>(config.Boards, AllBoardCfg_t{}))
   *   {}
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   */
  template <typename T, typename Optional>
  T getOptionalValue(Optional const& parameter, T defValue);


  //--------------------------------------------------------------------------


} // namespace util


//--------------------------------------------------------------------------
/*
 * Currently util::Quantity objects have serious issues with FHiCL validation,
 * so we are not reading them directly as such. This pulls in a number of
 * workarounds for:
 *  * every `Quantity` parameter: its FHiCL parameter must be declared as
 *    the base type `Quantity::value_t`
 *  * every optional parameter must be then read indirectly since a reference
 *    to the exact type of the parameter is required in such reading
 *  * for sequences, direct vector assignment
 *    (`vector<Quantity> = vector<Quantity::value_t>`) won't work, and since
 *    `Sequence::operator()` returns a temporary, this needs to be wrapped too
 */
template <typename SeqValueType>
struct util::fhicl::SequenceWrapper {
  SeqValueType seqValue;

  SequenceWrapper(SeqValueType seqValue): seqValue(seqValue) {}

  template <typename T>
  operator std::vector<T>() const
    { return { seqValue.cbegin(), seqValue.cend() }; }

}; // util::fhicl::SequenceWrapper


//--------------------------------------------------------------------------
template <typename Optional>
std::optional<typename Optional::value_type>
util::fhicl::getOptionalValue(Optional const& parameter) {

  using Value_t = typename Optional::value_type;

  if (!parameter.hasValue()) return std::nullopt;

  Value_t value;
  parameter(value);
  return { value };

} // util::fhicl::getOptionalValue(Optional& parameter)


//--------------------------------------------------------------------------
template <typename T, typename Optional>
T util::fhicl::getOptionalValue(Optional const& parameter, T defValue)
  { parameter(defValue); return defValue; }


//--------------------------------------------------------------------------


#endif // ICARUSCODE_UTILITIES_FHICLUTILS_H
