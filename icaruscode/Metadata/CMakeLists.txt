GenerateRepoVersionSource(${CMAKE_PROJECT_NAME})

message(DEBUG "Generated source file '${${CMAKE_PROJECT_NAME}_GIT_REPO_VERSION_SOURCE}.h'")

art_make_library()

cet_build_plugin(icaruscodeRepositoryVersion art::tool
  LIBRARIES
    icaruscode::Metadata
    icarusalg::Metadata
    sbncode::Metadata
    sbnobj::Common_Metadata
  )

install_headers(LIST "${CMAKE_CURRENT_BINARY_DIR}/${${CMAKE_PROJECT_NAME}_GIT_REPO_VERSION_SOURCE}.h")
install_source(LIST "${CMAKE_CURRENT_BINARY_DIR}/${${CMAKE_PROJECT_NAME}_GIT_REPO_VERSION_SOURCE}.cxx")

