#!/usr/bin/env bash
#
# Run to update the HTML page with the available versions.
# 
# Used HTML classes:
# 
# aside         : for messages of lesser importance
# reponame      : a code repository name
# GITbranch     : GIT hash or branch tag
# version       : a code version
# currentversion: the current code version
# updateline    : (block) the "last updated" line
# warning       : a message about possible failure
# 
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

# additional settings
declare -r RepoDir="./"

# template for the version list HTML file, relative to the script directory
declare -r TemplateFileName='versionlisttemplate.html'

# the first line with this tag in the template will be replaced with the list
declare -r VersionListTag='VERSIONLISTTAG'
# the first line with this tag in the template will be replaced with update info
declare -r UpdateTag='UPDATETAG'

# ------------------------------------------------------------------------------
function FillVersionsList() {
  
  local -r DocDir="${1%/}/"
  
  [[ -d "$DocDir" ]] || FATAL 2 "Can't find documentation base directory '${DocDir}'."
  
  DocVersions=() # global
  
  local VersionCandidatePath
  for VersionCandidatePath in "${DocDir}"* ; do
    [[ -d "$VersionCandidatePath" ]] || continue
    local VersionCandidate="$(basename "$VersionCandidatePath")"
    [[ "$VersionCandidate" =~ ^v[[:digit:]]{2}_[[:digit:]]{2}_[[:digit:]]{2}(_[[:digit:]]{2})?$ ]] || continue
    DocVersions+=( "$VersionCandidate" )
    
  done
  
} # FillVersionsList()


function MakeVersionItem() {
  
  local -r DocVersion="$1"
  local -r DocPath="$2"
  local -r Indent="${3:-"    "}"

  local ExpCodeName="$ExperimentCodeName"
  local Version="$DocVersion"
  local DocDate
  local DocGenExitCode
  
  local -r MetadataFile="${DocPath%/}/${MetadataFileRelPath}"
  if [[ -r "$MetadataFile" ]]; then
    # metadata file is currently (version 1.0) created by UpdateFiles.sh
    local MetadataVersion
    MetadataVersion="$(ExtractFromMetadata 'MetadataVersion' "$MetadataFile")"
    
    ExpCodeName="$(ExtractFromMetadata 'ExperimentCode' "$MetadataFile")" # maybe overwrite
    Version="$(ExtractFromMetadata 'CodeVersion' "$MetadataFile")" # maybe overwrite
    DocDate="$(ExtractFromMetadata 'Date' "$MetadataFile")" # maybe overwrite
    DocGenExitCode="$(ExtractFromMetadata 'ExitCode' "$MetadataFile")" # maybe overwrite
    DocCodeGITreference="$(ExtractFromMetadata 'GITreference' "$MetadataFile")" # maybe overwrite
    
  fi
  
  local Msg="<span class=\"reponame\">${ExpCodeName}</span> <span class=\"version\">${Version}</span>"
  if [[ -z "$DocGenExitCode" ]]; then
    Msg+=' <span class="aside">(no additional information available)</span>'
  else
    Msg+=" (documentation from ${DocDate}"
    if [[ -n "$DocCodeGITreference" ]] && [[ "$DocCodeGITreference" != "$Version" ]]; then
      Msg+=" for GIt reference <span class=\"GITbranch\">${DocCodeGITreference}</span>'"
    fi
    Msg+=')'
    if [[ "$DocGenExitCode" != 0 ]]; then
      Msg+=" <span class=\"warning\">[it may be incomplete: generation exited with code ${DocGenExitCode}]</span>"
    fi
  fi
  
  echo "${Indent}<li><a href=\"${DocVersion}\">${Msg}</a></li>"
  
} # MakeVersionItem()


function MakeHTMLlist() {
  
  local -r DocDir="$1"
  
  [[ -d "$DocDir" ]] || FATAL 2 "Can't find documentation base directory '${DocDir}'."
  
  # check for the special link
  if [[ -h "${DocDir%/}/${LatestLinkName}" ]]; then
    local -r LatestVersion="$(readlink "${DocDir%/}/${LatestLinkName}")"
    echo "    <li><a href=\"${LatestLinkName}\" class=\"currentversion\">latest documented version (${LatestVersion})</a></li>"
  fi
  
  # look for the actual versions
  local VersionCandidatePath
  for VersionCandidatePath in "${DocDir%/}/"* ; do
    
    [[ -d "$VersionCandidatePath" ]] || continue
    local VersionCandidate="$(basename "$VersionCandidatePath")"
    [[ "$VersionCandidate" =~ ^v[[:digit:]]{2}_[[:digit:]]{2}_[[:digit:]]{2}(_[[:digit:]]{2})?$ ]] || continue
    
    MakeVersionItem "$VersionCandidate" "$VersionCandidatePath"
    
  done
  
} # MakeHTMLlist()


function UpdatingRemarks() {
  local -r Indent="${1:-"  "}"
  cat <<EOF
${Indent}
${Indent}<hr>
${Indent}<div class="updateline">Last list refresh: $(date)</div>
${Indent}
EOF
} # UpdatingRemarks()


function MakeVersionListFile() {
  # Writes a HTML file on screen
  
  local -r TemplateFilePath="$1"
  local -r DocDir="$2"
  
  #
  # write a header comment
  #
  cat <<EOC
<!--
  This HTML file was automatically generated:
  
  Date:               $(date)
  Host:               $(hostname)
  Working directory: '$(pwd)'
  Command:            ${0} ${@}
  Script version:     ${SCRIPTVERSION}
  
  -->
EOC
  
  #
  # "parse" the template
  #
  # the horror...
  local Line
  local OldIFS="$IFS"
  IFS=''
  while read Line ; do
    if   [[ "$Line" =~ ${VersionListTag} ]]; then # replace this line
      MakeHTMLlist "$DocDir"
    elif [[ "$Line" =~ ${UpdateTag} ]]; then # replace this line
      UpdatingRemarks
    else
      # trying hard to avoid any reinterpretation; `echo` *might* have sufficed
      cat <<< "$Line"
    fi
  done < "$TemplateFilePath"
  IFS="$OldIFS"
  
} # MakeVersionListFile()


# ------------------------------------------------------------------------------

#
# argument "parsing"
#
declare DestFile="$1" # default set later

#
# checks
#
ExperimentName="$(FindExperiment)"
LASTFATAL "Can't detect the experiment name."

declare TemplateFilePath="${SCRIPTDIR%/}/${TemplateFileName}"
[[ -r "$TemplateFilePath" ]] || FATAL 2 "Can't find template file '${TemplateFileName}' in '${SCRIPTDIR}'."

declare -r ExperimentCodeName="$(ExperimentCodeProduct "$ExperimentName")"

declare DocDir="${PublishBaseDir%/}/${ExperimentCodeName}"
[[ -r "$DocDir" ]] || FATAL 2 "Can't find publishing base directory '${DocDir}'."

# set a default destination if needed
declare -r DefaultDestFile="${DocDir%/}/${VersionListFile}"
: ${DestFile:="$DefaultDestFile"}

#
# write the HTML file
#
if [[ -e "$DestFile" ]]; then
  echo "Overwriting: '${DestFile}'"
else
  echo "Writing: '${DestFile}'"
fi

MakeVersionListFile "$TemplateFilePath" "$DocDir" > "$DestFile"

