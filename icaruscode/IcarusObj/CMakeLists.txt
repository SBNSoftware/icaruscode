
##cet_make_library(SOURCE 
##                 ChannelROI.cxx
##                 PMTWaveformTimeCorrection.cxx
##                 LIBRARIES
##                 PRIVATE
##                 messagefacility::MF_MessageLogger
##                )
cet_make (
           NO_DICTIONARY
         )

art_dictionary(DICTIONARY_LIBRARIES
                 icaruscode::IcarusObj
                 lardataobj_RecoBase
                 canvas::canvas 
              )

install_headers()
install_source()
