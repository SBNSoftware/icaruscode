#!/usr/bin/env bash
#
# General settings for the scripts in this directory.
# Override the settings in `experiment_settings.sh`.
# 
# Usage:
# 
# # -- BEGIN -- boilerplate settings and library loading -------------------------
# SCRIPTNAME="$(basename "$0")"
# SCRIPTDIR="$(dirname "$0")"
# 
# declare LibraryToLoad
# for LibraryToLoad in 'settings.sh' 'experiment_settings.sh' 'utilities.sh' ; do
#   
#   source "${SCRIPTDIR%/}/${LibraryToLoad}" || exit $?
#   
# done
# unset LibraryToLoad
# # -- END -- boilerplate settings and library loading ---------------------------
#

declare DefaultExperimentName='ICARUS'

# name for a log directory (scripts may set it where they want)
declare LogDir='logs'

# the GIT branch or commit to check out before generating the documentation.
declare DefaultBranch='master'

# the subdirectory under the GIT repository where documentation generation might be found
declare RepoDocSubdir='doc'

# where the content is published
declare PublishBaseDir='/web/sites/i/icarus-exp.fnal.gov/htdocs/at_work/software/doc'

# directory of Doxygen generation metadata, relative to Doxygen output directory
declare MetadataFileRelPath='meta/doxygen'

# name of the index file for package software versions (must match the hard-coded name from the web site)
declare VersionListFile='versionlist.html'

# name of the current version of the software
declare LatestLinkName='latest'

# GitHub group name for fetching code
declare GitHubExperimentGroup='SBNSoftware'
