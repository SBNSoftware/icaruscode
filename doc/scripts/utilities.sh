#!/usr/bin/env bash
#
# Common utilities for the scripts in this directory.
#

# ------------------------------------------------------------------------------
function isFlagSet() {
  local -r VarName="$1"
  [[ -n "${!VarName//0}" ]]
} # isFlagSet()

function isFlagUnset() {
  local -r VarName="$1"
  [[ -z "${!VarName//0}" ]]
} # isFlagUnset()

function STDERR() { echo "$*" >&2 ; }
function ERROR() { STDERR "ERROR: $*" ; }
function FATAL() {
  local -i Code="$1"
  shift
  STDERR "FATAL (${Code}): $*"
  exit "$Code"
} # FATAL()
function LASTFATAL() {
  local -i Code="$?"
  [[ "$Code" == 0 ]] || FATAL "$Code" "$@"
} # LASTFATAL()


# ------------------------------------------------------------------------------
function isRelativePath() {
  local -r Path="$1"
  [[ "${Path:0:1}" != '/' ]]
} # isRelativePath()


function simplifyPath() {
  # Replaces '.' and '..' path components by plain text substitution.
  local Path="$(sed -E 's@/+@/@g' <<< "${1}/")" # replace double / with single ones
  
  # build a stack of directories, from the root on
  local Dir OtherPath
  local -a Stack
  local -i nParents=0 # number of parent directories at root
  while true ; do
    
    Dir="${Path%%/*}"
    OtherPath="${Path#${Dir}/}"
    case "$Dir" in
      ( '.' )
        ;;
      ( '..' )
        local -i lastElem=$((${#Stack[@]} - 1))
        if [[ $lastElem -ge 0 ]] && [[ "${Stack[$lastElem]}" != '..' ]]; then
          # if the element is empty then it's root (should be lastElem=0)
          # and we can't escalate to the parent
          [[ -n "${Stack[$lastElem]}" ]] && unset Stack[$lastElem]
        else
          let ++nParents
        fi
        ;;
      ( * )
        Stack+=( "$Dir" ) ;;
    esac
    
    [[ "${Dir}/" == "$Path" ]] && break
    Path="$OtherPath"
  done
  
  # merge the elements left
  local ReconstitutedPath
  while [[ $nParents -gt 0 ]]; do
    ReconstitutedPath+='../'
    let --nParents
  done
  for Elem in "${Stack[@]}" ; do
    ReconstitutedPath+="${Elem}/"
  done
  
  [[ "$ReconstitutedPath" != '/' ]] && ReconstitutedPath="${ReconstitutedPath%/}"
  echo "$ReconstitutedPath"
  
} # simplifyPath()


function makeAbsolutePath() {
  local Path="$1"
  local RelativeTo="$2"
  
  if isRelativePath "$Path" ; then
    [[ -z "$RelativeTo" ]] && RelativeTo="$(pwd)"
    Path="${RelativeTo%/}/${Path}"
  fi
  simplifyPath "$Path"
  
} # makeAbsolutePath()


# ------------------------------------------------------------------------------
function isExperiment() {
  #
  # Prints the name of the experiment as deduced.
  # Returns 0 on success, 1 if deduction failed, and 2 if more than one
  # experiments match.
  #
  
  # trick: look for the BlueArc mounts
  
  local CandidateData
  local -a Experiments
  for CandidateData in /*/data ; do
    local Candidate="$(dirname "${CandidateData#/}")"
    
    # is there an 'app' directory too?
    [[ -d "/${Candidate}/app" ]] || continue
    Experiments+=( "$Candidate" )
  done
  
  case "${#Experiments[@]}" in
    ( 0 ) return 1 ;;
    ( 1 ) ;;
    ( * ) return 2 ;;
  esac
  echo "${Experiments[0]^^}"
  
} # isExperiment()


function FindExperiment() {
  #
  # Prints the name of the experiment the directory or machine belongs to.
  #
  local -r ExperimentNameFile='ExperimentName'
  
  [[ -r "$ExperimentNameFile" ]] && cat "$ExperimentNameFile" && return
  
  local ExperimentName
  ExperimentName="$(isExperiment 2> /dev/null)"
  
  [[ $? != 0 ]] && ExperimentName="$(basename "$(pwd)")"
  
  echo "$ExperimentName"
  [[ -n "$ExperimentName" ]]
  
} # FindExperiment()


function FindLatestUPS() {
  local -r ProductName="$1"
  local -r Qualifiers="$2" # optional
  
  ups list -aKVERSION "$ProductName" ${Qualifiers:+-q "$Qualifiers"} | tr -d '" ' | sort -u | tail -n 1
  
} # FindLatestUPS()


function PickAnyQualifiersFor() {
  local -r ProductName="$1"
  local -r Version="$2"
  
  ups list -aKQUALIFIERS "$ProductName" "$Version" | tr -d '" ' | sort -u | tail -n 1
  
} # PickAnyQualifiersFor()


function UPSversion() {
  # Product must already be set up
  local -r Product="$1"
  
  local VarName
  VarName="${Product^^}_DIR"
  [[ -n "${!VarName}" ]]
  declare -p "$VarName" >& /dev/null || return 2
  
  VarName="${Product^^}_VERSION"
  declare -p "$VarName" >& /dev/null || return 1
  
  echo "${!VarName}"
} # UPSversion()


function NeedSetup() {
  local Experiment="$1"
  
  local ExperimentCodeName="$(ExperimentCodeProduct "$Experiment" )"
  
  local VarName="${ExperimentCodeName}_DIR"
  [[ -z "${!VarName}" ]]
  
} # NeedSetup()


function ExperimentSetup() {
  
  local ExperimentName="$1"
  
  local -r ExperimentCodeName="$(ExperimentCodeProduct "$ExperimentName")"
  
  if ! NeedSetup "$ExperimentName" ; then
    echo "${ExperimentName} software already set up (${ExperimentCodeName} $(UPSversion "$ExperimentCodeName")"
    return 0
  fi
  
  local -r experiment="${ExperimentName,,}"
  local -a SetupScriptCandidates=(
    "/cvmfs/${experiment}.opensciencegrid.org/products/${experiment}/setup_${experiment}.sh"
    "/cvmfs/${experiment}.opensciencegrid.org/products/setup_${experiment}.sh" 
  )
  
  local SetupScript
  for SetupScript in "${SetupScriptCandidates[@]}" "" ; do
    [[ -r "$SetupScript" ]] && break
  done
  
  [[ -r "$SetupScript" ]] || FATAL 1 "Can't find setup script for experiment '${Experiment}'."
  
  echo "Setting up the UPS area for ${Experiment}..."
  source "$SetupScript"
  
  local -r LatestVersion=$(FindLatestUPS "$ExperimentCodeName")
  
  local -r Qualifiers=$(PickAnyQualifiersFor "$ExperimentCodeName" "$LatestVersion")
  
  echo "Setting up: ${ExperimentCodeName} ${LatestVersion}${Qualifiers:+" (${Qualifiers})"}"
  source "$(ups setup "$ExperimentCodeName" "$LatestVersion" ${Qualifiers:+-q "$Qualifiers"})"
  
  source "$(ups setup git)" # failed? doesn't matter
  
} # ExperimentSetup()


function DetectROOTversion() {
  
  root-config --version
  
} # DetectROOTversion()


# ------------------------------------------------------------------------------
function ExperimentCodeProduct() {
  
  local Experiment="$1"
  
  # won't work for lariatsoft and dunetpc
  echo "${Experiment,,}code"
  
} # ExperimentCodeProduct()


function RedmineGITremoteURL() {
  local -r RepoName="$1"
  echo "http://cdcvs.fnal.gov/projects/${RepoName}"
} # RedmineGITremoteURL()


function FindGitRepository() {
  
  local RepoName="$1"
  
  local CandidateRepoDir
  for CandidateRepoDir in "${RepoDir:+${RepoDir%/}}/${RepoName}" ; do
    [[ -r "${CandidateRepoDir}/.git" ]] && break
  done
  
  [[ -r "$CandidateRepoDir" ]] || LASTFATAL "Could not find the local GIT repository for '${RepoName}'."
  echo "$CandidateRepoDir"
  
} # FindGitRepository()


function isInGITrepository() {
  
  local Path="$(makeAbsolutePath "${1:-.}")"
  
  while true ; do
    [[ -d "${Path%/}/.git" ]] && return 0
    [[ "$Path" == '/' ]] && break
    Path="$(dirname "$Path")"
  done
  return 1
} # isInGITrepository()


function GITrepositoryPath() {
  
  local Path="$(makeAbsolutePath "${1:-.}")"
  
  while true ; do
    if [[ -d "${Path%/}/.git" ]]; then
      echo "${Path%/}"
      return 0
    fi
    [[ "$Path" == '/' ]] && break
    Path="$(dirname "$Path")"
  done
  return 1
} # GITrepositoryPath()


function GITrepositoryName() {
  
  local Path
  Path="$(GITrepositoryPath "$1")" && basename "$Path"
  
} # GITrepositoryName()


function hasGITtag() {
  
  local -r Tag="$1"
  local -r Directory="$2"
  
  git${Directory:+ -C "$Directory"} rev-parse "$Tag" >& /dev/null

} # hasGITtag()


function checkoutGITtagIfExists() {
  
  local -r Tag="$1"
  local -r Directory="$2"
  
  hasGITtag && git${Directory:+ -C "$Directory"} checkout "$Tag"
  
  git${Directory:+ -C "$Directory"} rebase

} # checkoutGITtagIfExists()


# ------------------------------------------------------------------------------
function FindDoxyfile() {
  
  local -r ExptName="$1"
  
  local -r ExptCodeName="$(ExperimentCodeProduct "$ExptName" )"
  
  local -a CandidateDoxyfiles
  local GITrepoPath
  for GITrepoPath in "$(GITrepositoryPath)" "$(GITrepositoryPath "$0")" ; do
    [[ -n "$GITrepoPath" ]] || continue
    CandidateDoxyfiles+=(
      "${GITrepoPath}/${RepoDocSubdir}/${ExptCodeName}.doxy"
      "${GITrepoPath}/${RepoDocSubdir}/Doxyfile.${ExptCodeName}"
      "${GITrepoPath}/${RepoDocSubdir}/${ExptName}.doxy"
      "${GITrepoPath}/${RepoDocSubdir}/Doxyfile.${ExptCode}"
      "${GITrepoPath}/${RepoDocSubdir}/Doxyfile"
    )
  done
  
  CandidateDoxyfiles+=(
    "${ExptCodeName}/docs/Doxyfile"
    "Doxyfile.${ExptCodeName}"
    "Doxyfile.${ExptName}"
    "${ExptCodeName}.doxy"
    "${ExptName}.doxy"
    'Doxyfile'
    )
  
  local Doxyfile
  for Doxyfile in "${CandidateDoxyfiles[@]}" "" ; do
    [[ -r "$Doxyfile" ]] || continue
    echo "$Doxyfile"
    return
  done
  
  return 1
} # FindDoxyfile()


function Unquote() {
  
  local String="$1"
  String="${String#\"}"
  String="${String%\"}"
  String="${String#\'}"
  String="${String%\'}"
  echo $String
  
} # Unquote()


function ExtractFromMetadata() {
  # 
  # Prints the value of metadata value Key in the specified metadata file.
  # Metadata files are created by this same script.
  # 
  # An error code (1) is returned if the variable is not defined.
  # 
  # Usage:  ExtractFromMetadata Key MetadataFile
  # 
  # This works only for variables with values defined on a single line.
  # Also, it does not do any decent parsing.
  # 
  local -r Key="$1"
  local -r DataFile="$2"
  
  [[ -r "$DataFile" ]] || return 1 # no file, no key
  
  # note that we are avoiding using `eval` for security reasons
  local DefLine
  DefLine="$(grep -E "^[[:blank:]]*${Key}[[:blank:]]*=[[:blank:]]*" "$DataFile" | sed -e 's/[[:blank:]]*=[[:blank:]]*/=/' | sed -e 's/[[:blank:]]*#.*$//g')"
  local -i res=$?
  [[ $res != 0 ]] && return $res
  [[ -z "$DefLine" ]] && return 1
  
  local "$DefLine" || FATAL 1 "Could not assign '${Key}' with the line \`local \"${DefLine}\"\`."
  
  Unquote "${!Key}"
  
  return 0
} # ExtractFromMetadata()


function ExtractValueFromDoxyfile() {
  # 
  # Prints the value assigned to VarName in the specified file.
  # An error code (1) is returned if the variable is not defined.
  # 
  # Usage:  ExtractValueFromDoxyfile VarName DefinitionFile
  # 
  # This works only for variables with values defined on a single line.
  # Also, it does not do any decent parsing.
  # 
  local -r VarName="$1"
  local -r Doxyfile="$2"
  
  # since our metadata has a format compatible with Doxygen, we use the same tool
  ExtractFromMetadata "$VarName" "$Doxyfile"
  
} # ExtractValueFromDoxyfile()


# ------------------------------------------------------------------------------
function FindFirstVariableName() {
  #
  # Prints the name of the first "variable" in the input stream, introduced as
  # "${VarName}", with VarName a valid bash variable identifier
  #
  
  #
  # this means:
  # * do not print anything by default
  # * look for a line with pattern ${VarName} (VarName made of letters, underscores and, except the first character, numbers
  # * when found, execute the three commands:
  #     * match the variable name and replace the whole line with it
  #     * print the new line
  #     * quit
  # 
  sed -n -E '/\$\{[[:alpha:]_][[:alnum:]_]*\}/{s/^.*\$\{([[:alpha:]_][[:alnum:]_]*)\}.*$/\1/;p;q}'
} # FindFirstVariableName()


function EscapeSedReplacementLiteralPattern() {
  local Value="$1"
  local Escape="$2"
  Value="${Value//\\/\\\\}"
  [[ "$Escape" != '\' ]] && Value="${Value//${Escape}/\\${Escape}}"
  echo "$Value"
} # EscapeSedReplacementLiteralPattern()


function ReplaceEnvironmentVariables() {
  #
  # Using `sed`, which means a lot of potential problems...
  #
  local -r TemplateFile="$1"
  local -r DestFile="$2"
  
  if [[ "$TemplateFile" -ef "$DestFile" ]]; then
    FATAL 1 "ReplaceEnvironmentVariables does not support in-place changes ('${TemplateFile}' and '${DestFile}' are the same file)."
  fi 
  
  cp "$TemplateFile" "$DestFile"
  LASTFATAL "ReplaceEnvironmentVariables: can't create file '${DestFile}'."
  
  local VarName LastVarName
  while true ; do
    VarName="$(FindFirstVariableName < "$DestFile")"
    [[ $? == 0 ]] || break # done
    [[ -z "$VarName" ]] && break # done
    
    # avoid infinite loops (from bugs, presumably)
    [[ "$LastVarName" == "$VarName" ]] && FATAL 1 "Internal error: variable '${VarName}' still present in '${DestFile}' after substitution."
    LastVarName="$VarName"
    
    local Value
    Value="${!VarName}"
    LASTFATAL "ReplaceEnvironmentVariables: error evaluating the variable '${VarName}'."
    
    echo "Replacing: '\${${VarName}}' => '${Value}'"
    sed -i -E "s@\\$\\{${VarName}\\}@$(EscapeSedReplacementLiteralPattern "$Value" '@')@g" "$DestFile"
    LASTFATAL "ReplaceEnvironmentVariables: error replacing variable '${VarName}' with value '${Value}'."
    
  done
  return 0
} # ReplaceEnvironmentVariables()


# ------------------------------------------------------------------------------
# for test: bash utilities.sh testCommand args
if [[ "${BASH_SOURCE[0]}" == "$0" ]] && [[ $# -gt 0 ]]; then
  echo "CMD> $*"
  "$@"
else
  true
fi

# ------------------------------------------------------------------------------
