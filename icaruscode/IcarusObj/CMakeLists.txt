
cet_make_library(SOURCE 
                 ChannelROI.cxx
                 LIBRARIES
                 PRIVATE
                 messagefacility::MF_MessageLogger
                )

art_dictionary(DICTIONARY_LIBRARIES
                 icaruscode::IcarusObj
                 lardataobj_RecoBase
                 canvas::canvas 
              )

install_headers()
install_source()
