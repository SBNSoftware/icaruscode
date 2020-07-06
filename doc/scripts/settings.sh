#!/usr/bin/env bash
#
# Settings for the scripts in this directory.
# 
# Usage:
# 
# # -- BEGIN -- boilerplate settings and library loading -------------------------
# SCRIPTNAME="$(basename "$0")"
# SCRIPTDIR="$(dirname "$0")"
# 
# declare LibraryToLoad
# for LibraryToLoad in 'settings.sh' 'utilities.sh' ; do
#   
#   source "${SCRIPTDIR%/}/${LibraryToLoad}" || exit $?
#   
# done
# unset LibraryToLoad
# # -- END -- boilerplate settings and library loading ---------------------------
#

declare -r DefaultExperimentName='ICARUS'

# name for a log directory (scripts may set it where they want)
declare -r LogDir='logs'

# the GIT branch or commit to check out before generating the documentation.
declare -r DefaultBranch='master'

# the subdirectory under the GIT repository where documentation generation might be found
declare -r RepoDocSubdir='doc'

# where the content is published
declare -r PublishBaseDir='/web/sites/i/icarus-exp.fnal.gov/htdocs/at_work/software/doc'

# directory of Doxygen generation metadata, relative to Doxygen output directory
declare -r MetadataFileRelPath='meta/doxygen'

# name of the index file for package software versions (must match the hard-coded name from the web site)
declare -r VersionListFile='versionlist.html'

# name of the current version of the software
declare -r LatestLinkName='latest'

