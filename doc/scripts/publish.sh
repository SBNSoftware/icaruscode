#!/usr/bin/env bash
#
# Ran without arguments, will publish the current Doxygen output.
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
SCRIPTVERSION="1.0"


declare -r RepoDir="./"


# ------------------------------------------------------------------------------
declare Version="$1"
declare SourceDir="$2"

ExperimentName="$(FindExperiment)"
LASTFATAL "Can't detect the experiment name."

declare -r ExperimentCodeName="$(ExperimentCodeProduct "$ExperimentName" )"

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
[[ -d "$DestDir" ]] && FATAL 1 "Output directory '${DestDir}' already exists: remove it first!"
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
