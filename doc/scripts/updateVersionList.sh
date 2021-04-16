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
# Changes
# --------
# 
# 202004xx (petrillo@slac.stanford.edu) [v1.0]
#   original version
# 20210412 (petrillo@slac.stanford.edu) [v1.1]
#   added support for command line options
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

# additional settings
declare -r RepoDir="./"

# template for the version list HTML file, relative to the script directory
declare -r DefaultTemplateFileName='versionlisttemplate.html'

# the first line with this tag in the template will be replaced with the list
declare -r DefaultVersionListTag='VERSIONLISTTAG'
# the first line with this tag in the template will be replaced with update info
declare -r DefaultUpdateTag='UPDATETAG'


# ------------------------------------------------------------------------------
function printHelp() {
  
  cat <<EOH

Updates the HTML file presenting the list of available documentation of
${ExperimentName:-the currently configured experiment}.

Usage:  ${SCRIPTNAME}  [options]

Supported options:
--experiment=EXPERIMENTNAME [${ExperimentName:-autodetect}]
    the name of the experiment (e.g. 'MicroBooNE')
--outputfile=OUTPUTHTML
    sets the list file to be created (by default, the repository name under the
    configured PublishBaseDir${PublishBaseDir:+", i.e. '${PublishBaseDir}'"})
--reponame=REPONAME [${ExperimentCodeName:-autodetect}]
    the name of the repository (e.g. 'uboonecode')
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
# argument parsing
#
declare -i NoMoreOptions
declare -i DoVersion=0 DoHelp=0
declare -i CleanOutputArea=0 UpdateOutputArea=0
declare TemplateFileName="$DefaultTemplateFileName"
declare VersionListTag="$DefaultVersionListTag"
declare UpdateTag="$DefaultUpdateTag"


declare ExitWithCode
declare Param
for (( iParam = 1 ; iParam <= $# ; ++iParam )); do
  Param="${!iParam}"
  if isFlagUnset NoMoreOptions && [[ "${Param:0:1}" == '-' ]]; then
    case "$Param" in
      ( '--experiment='* )         ExperimentName="${Param#--*=}" ;;
      ( '--outputfile='* )         DestFile="${Param#--*=}" ;;
      ( '--reponame='* )           ExperimentCodeName="${Param#--*=}" ;;
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

: ${ExperimentName:="$(FindExperiment)"}
: ${ExperimentCodeName:="$(ExperimentCodeProduct "$ExperimentName")"}

if isFlagSet DoVersion ; then
  printVersion
  [[ -z "$ExitWithCode" ]] && ExitWithCode=0
fi
if isFlagSet DoHelp ; then
  printHelp
  [[ -z "$ExitWithCode" ]] && ExitWithCode=0
fi

[[ -n "$ExitWithCode" ]] && exit "$ExitWithCode"


#
# checks
#
# (some checks are delayed because they are not fatal when just printing help)
# 
[[ -n "$ExperimentName" ]] || FATAL 1 "Can't detect the experiment name."

[[ -n "$ExperimentCodeName" ]] || FATAL 1 "Can't put together the experiment code repository name for '${ExperimentName}'."

declare TemplateFilePath="${SCRIPTDIR%/}/${TemplateFileName}"
[[ -r "$TemplateFilePath" ]] || FATAL 2 "Can't find template file '${TemplateFileName}' in '${SCRIPTDIR}'."

# not controlled by option because used only to compose `DestFile`:
PublishDir="${PublishBaseDir%/}/${ExperimentCodeName}"
[[ -r "$PublishDir" ]] || FATAL 2 "Can't find publishing base directory '${PublishDir}'."

# set a default destination if needed
: ${DestFile:="${PublishDir%/}/${VersionListFile}"}

#
# write the HTML file
#
if [[ -e "$DestFile" ]]; then
  echo "Overwriting: '${DestFile}'"
else
  echo "Writing: '${DestFile}'"
fi

MakeVersionListFile "$TemplateFilePath" "$PublishDir" > "$DestFile"

