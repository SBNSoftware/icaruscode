![ICARUS logo](http://icarus.lngs.infn.it/img/n3.jpg)

ICARUS simulation and reconstruction software
==============================================

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
