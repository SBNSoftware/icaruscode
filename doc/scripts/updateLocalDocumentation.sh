#!/usr/bin/env bash
#
# Generation of icaruscode Doxygen documentation:
#  * fetches the latest Doxygen tag files for LArSoft and ROOT
#    (at most once per day)
#  * sets up the latest UPS version of `icaruscode` (needed for depedencies)
#  * checks out the latest tag of `icaruscode` `master` branch
#  * runs Doxygen (output is expected in the ${OutputDir} directory)
#
# -----------------------------------------------------
# Implementation note
# -----------------------------------------------------
#
# This script currently does not support generating documentation from multiple
# repositories. To implement such support, a possible way is to have the user
# specify a list of repositories including a main one (say, the first one)
# which is where the Doxygen configuration file is. Then that configuration file
# is processed to include as INPUT all the repositories, instead of just one.
# Some additional work is needed to achieve that.
#
#
# Changes
# --------
#
# 20200426 (petrillo@slac.stanford.edu) [v1.0]
#   first version in `icaruscode`; metadata version: 1
# 20200520 (petrillo@slac.stanford.edu) [v1.1]
#   using GitHub sources
#

# -- BEGIN -- boilerplate settings and library loading -------------------------
SCRIPTNAME="$(basename "$0")"
SCRIPTDIR="$(dirname "$0")"

declare LibraryToLoad
for LibraryToLoad in 'settings.sh' 'utilities.sh' ; do
  
  source "${SCRIPTDIR%/}/${LibraryToLoad}" || exit $?
  
done
unset LibraryToLoad
# -- END -- boilerplate settings and library loading ---------------------------

# ------------------------------------------------------------------------------
SCRIPTVERSION="1.1"


declare -ri MetadataVersion='1' # integral numbers please
declare -ri GitHubExperimentGroup='SBNSoftware'


# ------------------------------------------------------------------------------
function printHelp() {
  
  cat <<EOH

Generates the documentation of ${ExperimentName}.

Usage:  ${SCRIPTNAME}  [options]

Supported options:
--experiment=EXPERIMENTNAME [autodetect]
    the name of the experiment (e.g. 'MicroBooNE')
--reponame=REPONAME [autodetect]
    the name of the repository (e.g. 'uboonecode')
--branch=BRANCH , --tag=TAG , --commit=COMMIT
    use this branch, tag or commit instead of branch '${DefaultBranch}';
    the content of the repository will be changed
--currentbranch , -c
    do not change the branch of the GIT repository; if the repository is
    downloaded, its default branch is used ('${DefaultBranch}');
    otherwise, the current status is used
--clean
    removes the old documentation if present; by default, this does not happen
--force-tags
    forces the update of Doxygen tag files even if they have been downloaded
    since less than one day
--version , -V
    prints the script version and exits (hint: it's version ${SCRIPTVERSION}
--help , -h , -?
    prints this usage message and exits

EOH

} # printHelp()


function printVersion() {
  echo "${SCRIPTNAME} version ${SCRIPTVERSION}."
} # printVersion()


# ------------------------------------------------------------------------------

function UpdateTag() {
  #
  # Downloads the specified tag file and renames it as directed, if older than Age
  #
  # Usage:  UpdateTag  TagFile URL Description
  #
  local -r Age='yesterday'
  
  local TagFile="$1"
  local URL="$2"
  shift 2
  local Description="$*"
  
  local TempFile="$(mktemp --dry-run )"
  touch --date="$Age" "$TempFile"
  
  if [[ "$TagFile" -nt "$TempFile" ]] && isFlagUnset ForceTags ; then
    echo "Keeping '${TagFile}' (run with \`--force-tags\` to override)."
  else
    echo "Fetching tag file for ${Description} ('${TagFile}' from '${URL}')"
    wget --output-document "$TagFile" "$URL"
    touch "$TagFile"
  fi
  
  rm "$TempFile"

} # UpdateTag()


function UpdateROOTtag() {

  local ROOTVersion
  ROOTVersion="$(DetectROOTversion)"
  LASTFATAL "Can't detect ROOT version!" # did set up fail?
  
  local ROOTVersionTag="${ROOTVersion%/*}"
  ROOTVersionTag="${ROOTVersionTag//.}"
  
  echo "Fetching Doxygen tag for ROOT version ${ROOTVersion}"
  UpdateTag 'ROOT.tag' "https://root.cern/doc/v${ROOTVersionTag}/ROOT.tag" "ROOT"

} # UpdateROOTtag()


function UpdateLArSoftTag() {

#   local LArSoftVersion="$(UPSversion 'larsoft')"

  # Only one documented version: hope you like it
  UpdateTag 'LArSoft.tag' 'http://nusoft.fnal.gov/larsoft/doxsvn/html/doxytags-larsoft.xml' "LArSoft"

} # UpdateLArSoftTag()


function PrepareDoxyfile() {
  local -r DoxyfileTemplate="$1"
  local -r Doxyfile="$2"
  
  ReplaceEnvironmentVariables "$DoxyfileTemplate" "$Doxyfile"
} # PrepareDoxyfile()


function RestoreGITrepository() {
  
  local Path="$1"
  local Commit="$2"
  local -i Stashed="${3:-0}"
  
  # restore the status of the repository
  if [[ -n "$Commit" ]]; then
    echo "Restoring the repository (to ${Commit})"
    git -C "$Path" checkout "$Commit"
  fi
  if [[ -n "$Commit" ]]; then
    echo "Restoring unstaged content of repository"
    isFlagSet Stashed && git -C "$Path" stash pop --quiet
  fi
  
} # RestoreGITrepository()


################################################################################
function RunDoxygen() {
  
  #
  # find a Doxygen configuration file
  #
  local -r ExperimentName="$1"
  
  local MasterDoxyfile
  MasterDoxyfile="$(FindDoxyfile "$ExperimentName")"
  LASTFATAL "Can't find a doxygen configuration file for '${ExperimentName}'"
  
  local ExperimentCodeName="$(ExperimentCodeProduct "$ExperimentName" )"
  
  # in the future this constraint can be relaxed:
  local ExperimentCodeRepoPath="$ExperimentCodeName"
  
  local DateTag="$(date '+%Y%m%d')"
  local DoxygenLog="${LogDir%/}/$(basename "${MasterDoxyfile%.doxy*}")-${DateTag}.log"
  
  local -i res=0
  local NewGITclone=0
  if [[ ! -r "${ExperimentCodeRepoPath}/.git" ]]; then
    local -a GitSources=(
      $(GitHubRemoteURL "$ExperimentCodeName" "$GitHubExperimentGroup" )
      $(RedmineGITremoteURL "$ExperimentCodeName")
      )
    
    for GitSource in "${GitSources[@]}" ; do
      echo "Attempting to clone the GIT source tree of '${ExperimentCodeName}' from '${GitSource}'"
      git clone "$GitSource" "$ExperimentCodeName" >&2
      res=$?
      [[ $res == 0 ]] && break
      echo "  (failed)"
    done
    
    [[ -r "${ExperimentCodeRepoPath}/.git" ]] || FATAL 1 "Could not download '${ExperimentCodeName}' code!"
    
    NewGITclone=1
  fi
  
  local -i GITstashed=0
  local GIToldCommit

  if [[ "$RequestedExperimentCodeBranch" != 'HEAD' ]]; then
    
    local -r ExperimentCodeBranch="${RequestedExperimentCodeBranch:-$DefaultBranch}"
    
    if isFlagUnset NewGITclone ; then
      echo "Fetching updates for the existing GIT source tree of '${ExperimentCodeRepoPath}'"
      git -C "$ExperimentCodeRepoPath" fetch
    fi
    
    # remove everything that could bother our update;
    # git stash does not tell a script if a stash was created
    local -i nStashes="$(cd "$ExperimentCodeRepoPath" && git stash list | wc -l)"
    git -C "$ExperimentCodeRepoPath" stash save --quiet -- "${SCRIPTNAME} temporary for Doxygen generation - $(date)"
    LASTFATAL "Failed to save local changes in '${ExperimentCodeRepoPath}'. Giving up."
    local -i nNewStashes="$(cd "$ExperimentCodeRepoPath" && git stash list | wc -l)"
    if [[ "$nNewStashes" -gt "$nStashes" ]]; then
      echo "Changes in GIT repository temporary stashed."
      GITstashed=1
    fi
    
    # also keep track of where we are
    GIToldCommit="$(cd "$ExperimentCodeRepoPath" && git reflog -n 1 --format='format:%H')"
    local OldScriptChecksum="$(md5sum < "$0")"
    
    echo "Updating the GIT branch..."
    git -C "$ExperimentCodeRepoPath" checkout "$ExperimentCodeBranch"
    git -C "$ExperimentCodeRepoPath" rebase
    local -i res=$?
    if [[ $res != 0 ]]; then
      RestoreGITrepository "$ExperimentCodeRepoPath" "$GIToldCommit" "$GITstashed"
      FATAL "$res" "Failed to check out branch '${ExperimentCodeBranch}' of GIT repository ${ExperimentCodeRepoPath}."
    fi
    
    local GITnewCommit="$(cd "$ExperimentCodeRepoPath" && git reflog -n 1 --format='format:%H')"
    [[ "$GITnewCommit" == "$GIToldCommit" ]] && unset GIToldCommit # no change, no need to track
    local NewScriptChecksum="$(md5sum < "$0")"
    
    if [[ "$OldScriptChecksum" != "$NewScriptChecksum" ]]; then
      FATAL 1 "The attempt to update the repository '${ExperimentCodeRepoPath}' on the fly caused a change in this script ('${SCRIPTNAME}'). Please update the repository manually and try again."
      RestoreGITrepository "$ExperimentCodeRepoPath" "$GIToldCommit" "$GITstashed"
    fi
    
  fi
  
  local -r ExperimentCodeVersion="$(cd "$ExperimentCodeRepoPath" && git describe --abbrev=0)" # use as tag
  
  # special: if no branch was requested in particular, we move to the tagged version of it
  [[ -z "$RequestedExperimentCodeBranch" ]] && git -C "$ExperimentCodeRepoPath" checkout "$ExperimentCodeVersion"
  
  local GitDescription
  GitDescription="$(cd "$ExperimentCodeRepoPath" && git describe 2> /dev/null)"
  [[ $? == 0 ]] && echo "Repository status: ${GitDescription}"
  
  mkdir -p "$LogDir"
  
  local -r Doxyfile="${TMPDIR:-/tmp}/${ExperimentCodeName}-${ExperimentCodeVersion}.doxy"
  echo "Creating Doxygen configuration from '${MasterDoxyfile}'..."
  PrepareDoxyfile "$MasterDoxyfile" "$Doxyfile"
  
  cat <<EOB

Now starting Doxygen (log file: '${DoxygenLog}')
--------------------------------------------------------------------------------

EOB
  
  local ConfiguredOutputDir MetadataFile
  ConfiguredOutputDir="$(ExtractValueFromDoxyfile 'HTML_OUTPUT' "$Doxyfile")"
  if [[ $? == 0 ]]; then
    echo "Directory for HTML output: '${ConfiguredOutputDir}'"
    MetadataFile="${ConfiguredOutputDir%/}/${MetadataFileRelPath}"
    local OldRunDate OldExperimentCodeName OldVersion
    OldRunDate="$(ExtractFromMetadata 'Date' "$MetadataFile")"
    OldExperimentCodeName="$(ExtractFromMetadata 'ExperimentCode' "$MetadataFile")"
    OldVersion="$(ExtractFromMetadata 'CodeVersion' "$MetadataFile")"
    [[ -n "$OldRunDate" ]] && echo "Existing output created or updated on ${OldRunDate} for ${OldExperimentCodeName} ${OldVersion}."
  else
    echo "Could not detect the HTML output directory from '${Doxyfile}'"
    ConfiguredOutputDir='' # probably redundant
  fi
  
  if [[ -n "$ConfiguredOutputDir" ]] && [[ -d "$ConfiguredOutputDir" ]]; then
    if [[ -n "$OldVersion" ]] && [[ -n "$OldExperimentCodeName" ]] \
      && [[ "$OldVersion" == "$ExperimentCodeVersion" ]] \
      && [[ "$OldExperimentCodeName" == "$ExperimentCodeName" ]] \
      && isFlagUnset CleanOldDocs
    then
      echo "Existing output not removed as it is from the same version (run with \`--clean\` to override)."
    else
      echo "Removing old output ('${ConfiguredOutputDir}')..."
      rm -rf "$ConfiguredOutputDir"
    fi
  fi
  
  local -r DoxygenName='doxygen'
  local doxygen="$(which "${DoxygenName[@]}")"
  
  local -a Cmd=( $doxygen "$Doxyfile" "&>" "$DoxygenLog" )
  echo "${Cmd[@]}"
  eval "${Cmd[@]}"
  
  local -i ExitCode=$?
  
  if [[ -d "$ConfiguredOutputDir" ]]; then
    mv "$Doxyfile" "$ConfiguredOutputDir" && echo "Doxygen configuration copied into '${ConfiguredOutputDir}'."
  fi
  
  # restore the status of the repository
  RestoreGITrepository "$ExperimentCodeRepoPath" "$GIToldCommit" "$GITstashed"
  
  
  if [[ -d "$ConfiguredOutputDir" ]]; then
    mkdir -p "$(dirname "$MetadataFile")"
    cat <<EOB

Writing metadata into '${MetadataFile}'.

EOB
  ### BEGIN metadata version 1:
    cat <<EOM > "$MetadataFile"
MetadataVersion = '${MetadataVersion}'
Experiment      = '${ExperimentName}'
Date            = '$(date)'
Host            = '$(hostname)'
User            = '$(whoami)'
Configuration   = '${Doxyfile}'
ExperimentCode  = '${ExperimentCodeName}'
CodeVersion     = '${ExperimentCodeVersion:-"n/a"}'
GITreference    = '${GitDescription:-"n/a"}'
WordDirectory   = '$(pwd)'
CommandLine     = '${Cmd[@]}'
DoxygenPath     = '${doxygen}'
DoxygenVersion  = '$($doxygen --version)'
LogFile         = '${DoxygenLog}'
ExitCode        = '${ExitCode}'
EOM
  fi
  ### END metadata version 1
  
  
  cat <<EOB

${ExperimentCodeVersion:+"This was for ${ExperimentCodeName} ${ExperimentCodeVersion}."}

--------------------------------------------------------------------------------
Doxygen completed its execution.

EOB
  
  return $ExitCode

} # RunDoxygen()



################################################################################
###  Start here!
################################################################################
#
# parameter parsing
#
declare -i NoMoreOptions
declare -i DoVersion=0 DoHelp=0
declare -i CleanOldDocs=0 ForceTags=0 UpdateBranch=1
declare RequestedExperimentCodeBranch
declare ExitWithCode
declare Param
for (( iParam = 1 ; iParam <= $# ; ++iParam )); do
  Param="${!iParam}"
  if isFlagUnset NoMoreOptions && [[ "${Param:0:1}" == '-' ]]; then
    case "$Param" in
      ( '--experiment='* )         ExperimentName="${Param#--*=}" ;;
      ( '--reponame='* )           ExperimentCodeName="${Param#--*=}" ;;
      ( '--branch='* | '--tag='* | '--commit='* )
                                   RequestedExperimentCodeBranch="${Param#--*=}" ;;
      ( '--currentbranch' | '-c' ) RequestedExperimentCodeBranch='HEAD' ;;
      ( '--clean' | '-C' )         CleanOldDocs=1 ;;
      ( '--force-tags' )           ForceTags=1 ;;
      ( '--version' | '-V' )       DoVersion=1 ;;
      ( '--help' | '-h' | '-?' )   DoHelp=1 ;;
      ( '-' | '--' )               NoMoreOptions=1 ;;
      ( * )
        ERROR "Unknown option: '${Param}'."
        ExitWithCode=1
    esac
  else
    PositionalArguments+=( "$Param" )
  fi
done

if [[ "${#PositionalArguments[@]}" -gt 0 ]]; then
  FATAL 1 "Found ${#PositionalArguments[@]} spurious arguments on command line: ${PositionalArguments[@]}."
fi


: ${ExperimentName:="$(FindExperiment)"}
LASTFATAL "Can't detect the experiment name."

: ${ExperimentCodeName:="$(ExperimentCodeProduct "$ExperimentName")"}
LASTFATAL "Can't put together the experiment code repository name for '${ExperimentName}'."


if isFlagSet DoVersion ; then
  printVersion
  [[ -z "$ExitWithCode" ]] && ExitWithCode=0
fi
if isFlagSet DoHelp ; then
  printHelp
  [[ -z "$ExitWithCode" ]] && ExitWithCode=0
fi

[[ -n "$ExitWithCode" ]] && exit "$ExitWithCode"


################################################################################
declare -i nErrors
cat <<EOB

================================================================================
Set up
--------------------------------------------------------------------------------

EOB

GITrepoName="$(GITrepositoryName)"
[[ $? == 0 ]] && [[ "$GITrepoName" == "$ExperimentCodeName" ]] && FATAL 1 "Please do not run this script from within a GIT repository '${ExperimentCodeName}'."

ExperimentSetup "$ExperimentName"
LASTFATAL "Set up not completed, can't proceed."

cat <<EOB

================================================================================
ROOT tag
--------------------------------------------------------------------------------

EOB
UpdateROOTtag || let ++nErrors

cat <<EOB

================================================================================
LArSoft tag
--------------------------------------------------------------------------------

EOB
UpdateLArSoftTag || let ++nErrors

cat <<EOB

================================================================================
All setup done${nErrors:+" (with ${nErrors})"}
--------------------------------------------------------------------------------

EOB

[[ -z "$nErrors" ]] || exit 1


cat <<EOB

================================================================================
Now running Doxygen
--------------------------------------------------------------------------------

EOB

RunDoxygen "$ExperimentName"
