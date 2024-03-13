//File: CRTTrigger.h
//Brief: A CRT::Trigger stores the activity in a single ProtoDUNE-SP CRT module when the module triggered itself to be read out.  CRT modules 
//       trigger a readout when there is activity on a channel that is above some threshold.  When a CRT module triggers a readout, all channels' 
//       activity is read out, regardless of ADC counts.  A channel with above-baseline(?) activity at the time of readout is represented as a 
//       CRT::Hit stored by the owned CRT::Trigger.  So, a CRT::Trigger can only have 0 CRT::Hits if it was default-constructed.  
//
//       I don't yet know whether it is possible to recover which hit caused the readout.  Based 
//       this preliminary CRT data format on Matt Strait's artdaq-core code for handling the Double Chooz outer veto modules as the ProtoDUNE-SP 
//       Cosmic Ray Tagger: https://github.com/straitm/artdaq_core_demo/blob/crt/artdaq-core-demo/Overlays/CRTFragment.hh
//Author: Andrew Olivier aolivier@ur.rochester.edu

//Include guards so that this file doesn't interfere with itself when included multiple times
#ifndef CRT_TRIGGER_H
#define CRT_TRIGGER_H

//c++ includes
#include <cstdint>
#include <vector>
#include <limits>
#include <cstddef>

namespace CRT
{
  class Hit
  {
    public:
      //Constructor from information in CRT::Fragment plus detector name for LArSoft
      Hit(uint8_t channel, /*const std::string& detName,*/ uint16_t adc): fChannel(channel), /*fAuxDetname(detName),*/ fADC(adc)
      {
      }

      //Default constructor to make ROOT happy.  Set ADC peak to -1 and everything else to largest possible value so that default-constructed 
      //hits are obvious.  
      Hit(): fChannel(std::numeric_limits<decltype(fChannel)>::max()), /*fAuxDetName("NULL"),*/ fADC(std::numeric_limits<decltype(fADC)>::max()) {}

      //No resources managed by a CRT::Hit, so use default destructor
      virtual ~Hit() = default; 
  
      //Public access to stored information
      inline size_t Channel() const { return fChannel; }

      //inline std::string AuxDetName() const { return fAuxDetName; }

      inline short ADC() const { return fADC; }

      //Check whether this Hit was default-constructed
      inline bool IsDefault() const { return fADC == std::numeric_limits<decltype(fADC)>::max(); }

      //Print out to some ostream-like object for debugging
      template <class STREAM>
      void dump(STREAM& stream) const
      {
        stream << "CRT::Hit Dump:\n"
               << "Channel: " << fChannel << "\n"
               //<< "AuxDetName: " << fAuxDetName << "\n"
               << "ADC: " << fADC << "\n"
               << "Was this CRT::Hit default-constructed? " << (IsDefault()?"true":"false") << "\n";
        //return stream;
      }

    private: //The last thing I read about LArSoft data products and polymorphism told me not to make polymorphic data products.  
      //Information for identifying which CRT module strip this hit was recorded on 
      //TODO: Should there be some kind of AuxDetID for storing this detector-sensitive pair in the Geometry service? Probably not worth the effort 
      //      just for this class.  
      //TODO: Modify/write AuxDetGeo channel mapping such that AuxDetGeo index IS CRT module label from documentation
      size_t fChannel; //The index of the AuxDetSensitiveGeo that represents this strip in a CRT module
      //std::string fAuxDetName; //The name of the volume from the geometry that represents the CRT module strip this hit was recorded in. 
                               //LArSoft sometimes needs this information to look up AuxDetSensitiveGeos.   

      short fADC; //Baseline-subtracted Analog to Digital Converter value at time when this hit was read out.  
                  //TODO: Is this ADC value the integral over some time period?  The result of a fit?  
  };

  //Also overload operator << for std::ostream to make debugging CRT::Hits as easy as possible.  Really just call dump() internally.
  template <class STREAM>
  STREAM& operator << (STREAM& lhs, const CRT::Hit& hit)
  {
    return hit.dump(lhs);
  }

  class Trigger
  {
    public:
      //Constructor from information that comes from Matt's artdaq::Fragments.  Takes ownership of the CRT::Hits given.  Note that this 
      //is the only time that CRT::Hits can be added to a Trigger in accordance with some data product design notes I found.   
      //TODO: Should timestamp assembly be handled here or elsewhere?
      Trigger(const unsigned short channel, /*const std::string& detName,*/ const unsigned long long timestamp, 
              std::vector<CRT::Hit>&& hits): fChannel(channel), /*fDetName(detName),*/ fTimestamp(timestamp), fHits(hits) 
      {}

      Trigger(): fChannel(std::numeric_limits<decltype(fChannel)>::max()), /*fDetName(""),*/ 
                 fTimestamp(std::numeric_limits<decltype(fTimestamp)>::max()), fHits() {} //Default constructor to satisfy ROOT.  

      //User access to stored information.  See member variables for explanation
      inline unsigned short Channel() const { return fChannel; }
      //const std::string& DetName() const { return fDetName; }
      inline unsigned long long Timestamp() const { return fTimestamp; }
      inline const std::vector<CRT::Hit>& Hits() const { return fHits; }

      //Check whether this Hit was default-constructed
      inline bool IsDefault() const { return fChannel == std::numeric_limits<decltype(fChannel)>::max(); }

      //Dumping an object is also helpful for debugging.  
      template <class STREAM>
      void dump(STREAM& stream) const
      {
        stream << "CRT::Trigger dump:\n"
               << "Module: " << fChannel << "\n"
               //<< "Detector Name: " << fDetName << "\n"
               << "Timestamp: " << fTimestamp << "\n"
               << "Hits->  ";
        for(const auto& hit: fHits) hit.dump(stream);

        //return stream;
      }

    private:
      unsigned short fChannel; //Mapping to CRT module that was Triggered.  Index into the Geometry service's array of AuxDetGeos
                               //that returns the representation of this module.  Module AuxDetGeo shown contain strips as AuxDetSensitiveGeos 
                               //(see CRT::Hit::fChannel).  
                               //TODO: Ideally, this will correspond to CRT electronics channel one day.  
      //std::string fDetName; //Name of the geometry representation for the module that Triggered itself.  Needed sometimes by the Geometry service.  

      unsigned long long fTimestamp; //Timestamp when this Trigger occurred.  First 32 bits are Linux time since epoch; last 32 bits are time in 
                                     //ns. 
      std::vector<CRT::Hit> fHits; //All activity in CRT strips within this module when it was read out 
  };

  //ostream operator overload so users can do mf::LogWarning("CRTTriggerError") << /*stuff*/ << trigger << "\n". 
  template <class STREAM>
  STREAM& operator <<(STREAM& lhs, const CRT::Trigger& trigger)
  {
    return trigger.dump(lhs);
  }
}

#endif //CRT_TRIGGER_H
