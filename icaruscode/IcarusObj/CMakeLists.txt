file(GLOB lib_srcs *.cxx)
cet_make_library(SOURCE
                 ${lib_srcs}
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
                 sbnobj::ICARUS_PMT_Trigger_Data
              )

install_headers()
install_source()
