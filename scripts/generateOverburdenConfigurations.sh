#!/usr/bin/env bash

declare -r SCRIPTNAME="$(basename "$0")"
declare -r SCRIPTDIR="$(dirname "$0")"
declare -r SCRIPTVERSION="1.0"

declare -Ar GenerationScriptNames=(
  ['overburden']='makeOverburdenVersionOf.sh'
  ['nooverburden']='makeNoOverburdenVersionOf.sh'
  )

declare -r FHiCLbaseDir="${ICARUSCODE_DIR}/fcl"

declare -r JobListFileName="OverburdenConfigurations.txt"


# ------------------------------------------------------------------------------
declare -r EndLine=$'\n'
function STDERR() { echo "$*" >&2 ; }
function ERROR() { STDERR "ERROR: $*" ; }
function FATAL() {
  local Code="$1"
  shift
  STDERR "FATAL(${Code}): $*"
  exit "$Code"
} # FATAL()


# ------------------------------------------------------------------------------
function PrintHelp() {
  cat <<EOH
Generates overburden and no-overburden versions of job file specified in
'${JobListFileName}'.

Usage:  ${SCRIPTNAME}

This script must be run in a MRB development area with icarus_code checked out.

EOH

} # help()


# ------------------------------------------------------------------------------
function MakeUpName() {
  local TemplateName="${1%.fcl}"
  local Tag="$2"
  
  local Name="$TemplateName"
  
  while true ; do
  
    Name="${Name/_standard/_${Tag}}"
    [[ "$Name" != "$TemplateName" ]] && break
    
    Name="${Name/standard_/${Tag}_}"
    [[ "$Name" != "$TemplateName" ]] && break
    
    Name+="_${Tag}"
    
    break
  done
  
  echo "${Name}.fcl"
  
} # MakeUpName()


function GenerateConfigurationsFor() {
  local FHiCLname="${1%.fcl}.fcl"
  
  #
  # find the original FHiCL file (we need the full path)
  #
  local -a Matches=( $(find "$FHiCLbaseDir" -name "$FHiCLname" ) )
  case "${#Matches[@]}" in
    ( 0 )
      ERROR "No configuration file '${FHiCLname}' found in '${FHiCLbaseDir}'"
      return 2
      ;;
    ( 1 )
      FHiCLpath="${Matches[0]}"
      ;;
    ( * )
      ERROR "Multiple (${#Matches[@]}) configuration files '${FHiCLbaseDir}' found in '${FHiCLbaseDir}':${Matches/%/${EndLine}  }"
      return 2
  esac
  
  #
  # find the original FHiCL file (we need the full path)
  #
  local ScriptTag Script ConfName
  local -a Cmd
  local -i nErrors=0
  for ScriptTag in "${!GenerationScriptPaths[@]}" ; do
    Script="${GenerationScriptPaths["$ScriptTag"]}"
    
    ConfName="$(MakeUpName "$FHiCLname" "$ScriptTag")"
    
    echo "'${FHiCLname}' => '${ConfName}' ['${ScriptTag}']"
    Cmd=( "$Script" "$FHiCLpath" "$ConfName" )
    if [[ -n "${FAKE//0}" ]]; then
      echo "$SHELL" "${Cmd[@]}"
    else
      "$SHELL" "${Cmd[@]}"
    fi
    if [[ $? != 0 ]]; then
      ERROR "Failed to run ${ScriptTag} script on '${FHiCLpath}'"
      let ++nErrors
      continue
    fi
  done
  
  [[ $nErrors == 0 ]] # return value
} # GenerateConfigurationsFor()


# ------------------------------------------------------------------------------
#
# setup check
#

# find scripts
declare -A GenerationScriptPaths
for ScriptTag in "${!GenerationScriptNames[@]}" ; do
  ScriptName="${GenerationScriptNames[$ScriptTag]}"
  
  while true ; do
    
    ScriptPath="$(which "$ScriptName" 2> /dev/null)"
    [[ $? == 0 ]] && break
    
    ScriptPath="${SCRIPTDIR}/${ScriptName}"
    [[ $? == 0 ]] && break
    
    FATAL 2 "Could not locate the script '${ScriptName}'!"
  done
  GenerationScriptPaths+=( ["$ScriptTag"]="$ScriptPath" )
done

# find job list
declare JobListFilePath
for CandidateDir in "$FHiCLbaseDir" "." ; do
  JobListFilePath="${CandidateDir}/${JobListFileName}"
  [[ -r "$JobListFilePath" ]] && break
done
[[ -r "$JobListFilePath" ]] || FATAL 2 "Couldn't find job list file '${JobListFileName}'"


# extract job list
declare Jobs
while read Line ; do
  [[ -z "$Line" ]] && continue
  [[ "${Line:0:1}" == "#" ]] && continue
  FHiCLname="${Line%.fcl}.fcl"
  Jobs+=( "$FHiCLname" )
done < "$JobListFilePath"


cat <<EOH
--- Setup ----------------------------------------------------------------------
Configuration base directory: '${FHiCLbaseDir}'
Job name list:                '${JobListFileName}'
Generation scripts (${#GenerationScriptPaths[@]}):${GenerationScriptPaths[@]/#/${EndLine}  }
Job configuration (${#Jobs[@]}):${Jobs[@]/#/${EndLine}  }
--------------------------------------------------------------------------------
EOH

# we need the full path of the FHiCL file
declare -i nErrors=0
for Job in "${Jobs[@]}"; do
  GenerateConfigurationsFor "$Job" || let ++nErrors
done

[[ "$nErrors" -gt 0 ]] && FATAL 1 "Detected ${nErrors} failures."


# ------------------------------------------------------------------------------





exit
# ------------------------------------------------------------------------------




# script configuration
declare -r SubdirName='overburden'
declare -r DropInConfigurationName='use_overburden_geometry_icarus.fcl'


# ------------------------------------------------------------------------------

if [[ $# -lt 2 ]]; then
  help
  exit
fi

declare -r SourceFile="$1"
declare -r DestName="$(basename "${2%.fcl}.fcl")"

declare -r SourceDir="$(dirname "$SourceFile")"
declare -r SourceName="$(basename "${SourceFile%.fcl}.fcl")"
declare -r DestDir="${SourceDir%/}/${SubdirName}"
declare -r DestFile="${DestDir:+"${DestDir%/}/"}${DestName}"

if [[ "$SourceName" == "$DestName" ]]; then
  echo "Cowardly refusing to create a new FHiCL file with the same name as the original one ('${SourceName}')." >&2
  exit 1
fi

echo "'${SourceName}' => '${DestName}'"

[[ -d "$DestDir" ]] || mkdir "$DestDir"

cat <<EOF > "$DestFile"
#
# File:    ${DestName}
# Purpose: Overburden version of \`${SourceName}\`.
# Author:  $(printUserName)
# Date:    $(date)
# 
# This file was automatically generated by ${SCRIPTNAME} version ${SCRIPTVERSION}.
#

#include "${SourceName}"

# turn to overburden geometry:
#include "${DropInConfigurationName}"

EOF
