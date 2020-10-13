#!/usr/bin/env bash
#
# Custom settings for the scripts in this directory.
# General default values are in 'settings.sh';
# the custom ones for the experiment are here.
#

declare -r DefaultExperimentName='ICARUS'

# other repositories to check out for a complete documentation
declare -ra AdditionalRepoNames=( 'sbnobj' 'sbncode' 'icarusalg' 'icarus_signal_processing' )

# where the content is published
declare -r PublishBaseDir='/web/sites/i/icarus-exp.fnal.gov/htdocs/at_work/software/doc'

declare GitHubExperimentGroup='SBNSoftware'
