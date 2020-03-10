/**
 * @file   icaruscode/Utilities/WaveformOperations.h
 * @brief  Operations on waveform samples.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   June 27, 2019
 *
 * This is a header-only library.
 */

#ifndef ICARUSCODE_UTILITIES_WAVEFORMOPERATIONS_H
#define ICARUSCODE_UTILITIES_WAVEFORMOPERATIONS_H

/**
 * @brief Functions to manipulate waveform sample values.
 *
 * This namespace provides trivial functions to manage operations on waveform
 * samples given a certain polarity.
 *
 *
 */
namespace icarus::waveform_operations {

  namespace details {

    // -------------------------------------------------------------------------
    template <typename Sample, int Polarity, typename = void>
    struct FlipImpl;

    template <typename Sample, int Polarity>
    struct FlipImpl<Sample, Polarity, std::enable_if_t<(Polarity > 0)>> {
      static constexpr Sample transform(Sample sample) { return sample; }
    };

    template <typename Sample, int Polarity>
    struct FlipImpl<Sample, Polarity, std::enable_if_t<(Polarity < 0)>> {
      static constexpr Sample transform(Sample sample) { return -sample; }
    };

    template <int Polarity, typename Sample>
    constexpr Sample flip(Sample sample)
      { return FlipImpl<Sample, Polarity>::transform(sample); }


    // -------------------------------------------------------------------------
    /**
     * @brief Operations on a waveform with a fixed baseline.
     * @tparam Sample type of ADC counts in the waveform
     * @tparam Transform transformation to apply on the samples
     *
     * This object provides several functions to operate on waveforms
     * independently of their "polarity".
     *
     * It can be used as an object that stores the waveform baseline, in which
     * case a few member functions provide operations with respect to that
     * baseline. Static functions independent of the baseline are also offered.
     *
     * The type of operations offered here are:
     *  * shifts with respect to a baseline
     *  * distance from the baseline
     *  * relative comparisons
     *
     */
    template <typename Sample, Sample Transform(Sample)>
    struct WaveformTransformedOperations {

      using Sample_t = Sample; ///< Type of ADC samples.

      // --- BEGIN --- Constructors --------------------------------------------
      /// Constructor: sets the baseline to `0`.
      WaveformTransformedOperations() = default;

      /// Constructor: sets the `baseline`.
      constexpr WaveformTransformedOperations(Sample_t baseline)
        : fBaseline(baseline) {}

      // --- END --- Constructors ----------------------------------------------


      // --- BEGIN --- Operations relative to the baseline ---------------------
      /// @name Operations relative to the baseline
      /// @{

      /// Shift (addition) of an offset `shift` to the baseline.
      constexpr Sample_t shiftFromBaseline(Sample_t shift) const
        { return shiftBy(fBaseline, shift); }

      /// Distance of `sample` from the baseline.
      constexpr Sample_t subtractBaseline(Sample_t sample) const
        { return subtractBaseline(sample, fBaseline); }

      /// @}
      // --- END --- Operations relative to the baseline -----------------------


      // --- BEGIN --- Operations ----------------------------------------------
      /// @name Operations involving two samples.
      /// @{

      /// Difference between two samples (`to - from`).
      static constexpr Sample_t distance(Sample_t from, Sample_t to)
        { return transform(to - from); }

      /// Shift (addition) of an offset `shift` to a baseline.
      static constexpr Sample_t shiftBy(Sample_t baseline, Sample_t shift)
        { return baseline + transform(shift); }

      /// Distance of `sample` from the baseline: just mnemonic for `distance`.
      static constexpr Sample_t subtractBaseline
        (Sample_t sample, Sample_t baseline)
        { return distance(baseline, sample); }

      /// @}
      // --- END --- Operations ------------------------------------------------


      // --- BEGIN --- Comparisons ---------------------------------------------
      /// @name Comparisons
      /// @{

      static constexpr bool lessThan(Sample_t a, Sample_t b)
        { return transform(a) < transform(b); }
      static constexpr bool greaterThan(Sample_t a, Sample_t b)
        { return transform(a) > transform(b); }
      static constexpr bool noLessThan(Sample_t a, Sample_t b)
        { return transform(a) >= transform(b); }
      static constexpr bool noGreaterThan(Sample_t a, Sample_t b)
        { return transform(a) <= transform(b); }

      /// @}
      // --- END --- Comparisons -----------------------------------------------


        private:

      Sample_t fBaseline { 0 }; ///< Waveform baseline [ADC counts]

      // this is just to make it explicit that this is constexpr.
      static constexpr Sample_t transform(Sample_t sample)
        { return Transform(sample); }

    }; // WaveformTransformedOperations


    // -------------------------------------------------------------------------


  } // namespace details


  // ---------------------------------------------------------------------------
  /**
   * @brief Waveform operations of waveforms with specified polarity.
   * @tparam Sample type of ADC count in the waveform
   * @tparam Polarity whether the positive signal develops above the baseline
   *                  (positive polarity, +1) or below it (negative polarity,
   *                  -1)
   * @see details::WaveformTransformedOperations
   *
   * This is an alias of `details::WaveformTransformedOperations`.
   */
  template <typename Sample, int Polarity>
  using Operations = details::WaveformTransformedOperations
    <Sample, details::flip<Polarity, Sample>>;

  /// Waveform operations for positive polarity waveforms.
  template <typename Sample>
  using PositivePolarityOperations = Operations<Sample, +1>;

  /// Waveform operations for negative polarity waveforms.
  template <typename Sample>
  using NegativePolarityOperations = Operations<Sample, -1>;


  // ---------------------------------------------------------------------------


} // namespace icarus::waveform_operations


#endif // ICARUSCODE_UTILITIES_WAVEFORMOPERATIONS_H
