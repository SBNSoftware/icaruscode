cet_enable_asserts()

# Gather all source files in this Tools subdirectory.
file(GLOB tool_srcs *.cc)

# Build a library target from these sources.
# This will create a target whose name is based on the directory (or you can force it with a variable).
art_make_library(
  LIBRARIES PUBLIC ROOT::Core ROOT::Gpad ROOT::Geom
  SOURCE ${tool_srcs}
)

install_headers()
install_source()
