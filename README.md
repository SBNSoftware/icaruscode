![ICARUS logo](http://icarus.lngs.infn.it/img/n3.jpg)

ICARUS simulation and reconstruction software
==============================================

--Ivan's Branch

The documentation of ICARUS is spread mainly between:

* [ICARUS wiki in GitHub](https://sbnsoftware.github.io/icaruscode_wiki/Wiki) at https://sbnsoftware.github.io/icaruscode_wiki/Wiki;
  there is also an [older Fermilab Redmine wiki](https://cdcvs.fnal.gov/redmine/projects/icaruscode/wiki)
* automatically generated (Doxygen) [code documentation](https://icarus-exp.fnal.gov/at_work/software/doc/icaruscode/latest) at https://icarus-exp.fnal.gov/at_work/software/doc/icaruscode/versionlist.html
  (requires authentication)

ICARUS software is spread across several repositories under [GitHub SBNSoftware](https://github.com/SBNSoftware). Help yourself:
* [`icaruscode`](https://github.com/SBNSoftware/icaruscode) is the main entry point (code based on _art_)
* [`icarusalg`](https://github.com/SBNSoftware/icarusalg) includes _art_-independent code
* [`icarus_signal_processing`](https://github.com/SBNSoftware/icarus_signal_processing) includes _art_-independent code
* `icarus_data` distributes large data files
* [`icarusutil`](https://github.com/SBNSoftware/icarusutil) mostly contains experiment customization of [`larbatch`](https://github.com/LArSoft/larbatch) utilities
* [`sbncode`](https://github.com/SBNSoftware/sbncode) is code shared with other SBN experiments
* [`sbnobj`](https://github.com/SBNSoftware/sbnobj) contains data object definitions shared with SBN
* [`sbndaq_artdaq_core`](https://github.com/SBNSoftware/sbndaq_artdaq_core) interfaces with DAQ and holds data objects definitions for ICARUS "raw" data

ICARUS software is in great part based on LArSoft. Here are some shortcuts to its GitHub repositories:
* [`larsoft`](https://github.com/LArSoft/larsoft) (umbrella repository)
* data product repositories:
    * [`larcoreobj`](https://github.com/LArSoft/larcoreobj), basic geometry data products
    * [`lardataobj`](https://github.com/LArSoft/larcorealg), most simulation and reconstruction data products (not the `simb` namespace)
* [`larcorealg`](https://github.com/LArSoft/larcorealg), geometry code and some bsic utilities
* [`lardataalg`](https://github.com/LArSoft/lardataalg), some basic service providers (`LArProperties`, `DetectorClocks`, `DetectorProperties`) and utilities
* [`larevt`](https://github.com/LArSoft/larevt), with database code
* [`larsim`](https://github.com/LArSoft/larsim), detector simulation (esp. GEANT4)
* [`larana`](https://github.com/LArSoft/larana), including calorimetry and optical reconstruction
* [`larreco`](https://github.com/LArSoft/larreco), most reconstruction algorithms and modules
* [`art`](https://github.com/art-framework-suite/art), framework underpinning data processing (many support libraries are also in [art-framework-suite](https://github.com/art-framework-suite))
