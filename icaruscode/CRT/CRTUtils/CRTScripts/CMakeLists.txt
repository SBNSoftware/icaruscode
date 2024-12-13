# The txt files here are used to load Top CRT corrections
# They need to be found via FW_SEARCH_PATH

# install txt files
file(GLOB txt_file *.txt)
install_fw(LIST ${txt_file})