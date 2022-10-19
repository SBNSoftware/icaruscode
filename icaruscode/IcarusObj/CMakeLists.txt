
cet_make_library(SOURCE
                 ChannelROI.cxx
                 LIBRARIES
                 PRIVATE
                 messagefacility::MF_MessageLogger
                 larcoreobj::headers
                 lardataobj::headers
                )

art_dictionary(DICTIONARY_LIBRARIES
                 icaruscode::IcarusObj
                 lardataobj::RecoBase
                 lardataalg::headers
                 canvas::canvas
              )

install_headers()
install_source()
