cet_make_library(SOURCE 
                 ChannelROI.cxx
                 LIBRARIES
                 PRIVATE
                 messagefacility::MF_MessageLogger
                )

##art_dictionary()
build_dictionary(DICTIONARY_LIBRARIES
                 icaruscode::IcarusObj
                 canvas::canvas 
                )

install_headers()
install_source()
