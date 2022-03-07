#!/usr/bin/env bash
#
# Ran without arguments, will publish the current Doxygen output.
#
#
# Changes
# --------
#
# 20200520 (petrillo@slac.stanford.edu) [v1.1]
#   * added script argument parsing
#   * added `--update` option
#   * the option to remove the existing data (`CleanOutputArea`) is not enabled
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


declare -r RepoDir="./"


# ------------------------------------------------------------------------------
function printHelp() {
  
  cat <<EOH

Copies the documentation into the public web page storage area.

Usage:  ${SCRIPTNAME}  [options]

Supported options:
--srcdocdir=SOURCEDIR [autodetect]
    directory where the HTML documentation has been rendered
--experiment=EXPERIMENTNAME [autodetect]
    the name of the experiment (e.g. 'MicroBooNE')
--codeversion=VERSION [from metadata]
    the version of the code being documented
--update
    if there is already documentation for the same version, update it
    instead of refusing to proceed; the existing directory is not removed
    before publishing
--version , -V
    prints the script version and exits (hint: it's version ${SCRIPTVERSION}
--help , -h , -?
    prints this usage message and exits

EOH

} # printHelp()


function printVersion() {
  echo "${SCRIPTNAME} version ${SCRIPTVERSION}."
} # printVersion()


################################################################################
###  Start here!
################################################################################
#
# parameter parsing
#
declare -i NoMoreOptions
declare -i DoVersion=0 DoHelp=0
declare -i CleanOutputArea=0 UpdateOutputArea=0
declare ExitWithCode
declare Param
for (( iParam = 1 ; iParam <= $# ; ++iParam )); do
  Param="${!iParam}"
  if isFlagUnset NoMoreOptions && [[ "${Param:0:1}" == '-' ]]; then
    case "$Param" in
      ( '--experiment='* )         ExperimentName="${Param#--*=}" ;;
      ( '--srcdocdir='* )          SourceDir="${Param#--*=}" ;;
      ( '--codeversion='* )        Version="${Param#--*=}" ;;
      ( '--update' )               UpdateOutputArea=1 ;;
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

[[ "${#PositionalArguments[@]}" -gt 0 ]] && FATAL 1 "Found ${#PositionalArguments[@]} spurious arguments on command line: ${PositionalArguments[@]}."

: ${ExperimentName="$(FindExperiment)"}
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


# ------------------------------------------------------------------------------
#
# find the GIT repository (may contain doxygen file!)
#
declare GitRepoPath
GitRepoPath="$(FindGitRepository "$ExperimentCodeName")"
if [[ $? == 0 ]] && [[ -n "$GitRepoPath" ]]; then
  echo "GIT repository: '${GitRepoPath}'"
else
  echo "GIT repository not found."
  GitRepoPath=""
fi

#
# find the Doxygen configuration file
#
declare Doxyfile
Doxyfile="$(FindDoxyfile "$ExperimentName")"
if [[ $? == 0 ]] && [[ -n "$Doxyfile" ]]; then
  echo "Doxygen configuration: '${Doxyfile}'"
else
  echo "Doxygen configuration file not found."
  Doxyfile=""
fi

#
# finally find the directory of Doxygen output (hint: it's "html")
#
if [[ -z "$SourceDir" ]]; then
  [[ -r "$Doxyfile" ]] || FATAL 1 "The source directory must be specified."
  SourceDir="$(ExtractValueFromDoxyfile 'HTML_OUTPUT' "$Doxyfile")"
  [[ -n "$SourceDir" ]] || FATAL 1 "Can't extract the output directory from the Doxygen configuration file."
  echo "Doxygen output directory detected: '${SourceDir}'."
  MetadataFile="${SourceDir%/}/${MetadataFileRelPath}"
  if [[ -r "$MetadataFile" ]]; then
    MetadataVersion="$(ExtractFromMetadata 'MetadataVersion' "$MetadataFile")"
    echo "  => metadata (version ${MetadataVersion}) found at: '${MetadataFile}'"
  fi
fi
[[ -d "$SourceDir" ]] || FATAL 2 "Can't find the source directory '${SourceDir}'."

#
# with metadata we know everything about generation of documentation
#
if [[ "${MetadataVersion:-0}" -ge 1 ]]; then
  Version="$(ExtractFromMetadata 'CodeVersion' "$MetadataFile")"
  [[ $? == 0 ]] && [[ -n "$Version" ]] && echo "Code version from metadata: ${Version}."
fi

[[ -z "$Version" ]] && FATAL 1 "The version of the documentation must be specified."

#
# transfer the content
#
declare TransferLogFile="${LogDir:+${LogDir%/}/}transfer-${ExperimentCodeName}-${Version}.log"
declare -r DestVersionsDir="${PublishBaseDir:+${PublishBaseDir%/}/}${ExperimentCodeName}"
mkdir -p "$DestVersionsDir"
declare -r DestDir="${DestVersionsDir}/${Version}"
if [[ -d "$DestDir" ]]; then
  if isFlagSet CleanOutputArea ; then
    echo "The existing content in '${DestDir}' will be removed before proceeding."
    rm -Rf "$DestDir"
  elif isFlagSet UpdateOutputArea ; then
    echo "The existing content in '${DestDir}' will be updated."
    # no action actually needed
  else
    FATAL 1 "Output directory '${DestDir}' already exists: remove it first!"
  fi
fi
declare -a Cmd=( rsync -av "${SourceDir%/}/" "${DestDir%/}/" ">&" "$TransferLogFile")

echo "Copying the documentation content:"
echo "${Cmd[@]}"
eval "${Cmd[@]}"
LASTFATAL "Error while copying the content from '${SourceDir}' to '${DestDir}'."

#
# update the link
#
echo "Updating latest version link:"
declare -r LatestLinkPath="${DestVersionsDir}/${LatestLinkName}"
[[ -h "$LatestLinkPath" ]] && rm "$LatestLinkPath"
ln -sfv "$(basename "$DestDir")" "$LatestLinkPath"
