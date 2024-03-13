#ifndef artdaq_core_Data_dictionary_control_hh
#define artdaq_core_Data_dictionary_control_hh

// This header defines the CPP symbol HIDE_FROM_ROOT. This symbol
// can be used to hide code from Root dictionaries, as shown below
//
// Usage:
//
//   #if HIDE_FROM_ROOT
//      void function_that_you_want_to_hide() { ... };
//      void another_function_to_hide() {...};
//   #endif
//

#undef HIDE_FROM_ROOT
#if !defined(__GCCXML__) && defined(__GXX_EXPERIMENTAL_CXX0X__)
#define HIDE_FROM_ROOT 1
#endif

#endif /* artdaq_core_Data_dictionary_control_hh */
