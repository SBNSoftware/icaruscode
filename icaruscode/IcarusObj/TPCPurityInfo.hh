/**
 * \class TPCPurityInfo
 *
 * \ingroup anab
 *
 * \brief TPCPurity Analysis Info
 *
 */

#ifndef TPCPurityInfo_hh_
#define TPCPurityInfo_hh_


namespace anab {

  struct TPCPurityInfo{
    
    unsigned int Run;
    unsigned int Subrun;
    unsigned int Event;

    unsigned int Cryostat;
    unsigned int TPC;

    double Attenuation;
    double FracError;

    double Wires;
    double Ticks;

    TPCPurityInfo();
    void Print() const;
  };

}


#endif
