#!//usr/bin/env bash
#
# Run with `--help` for terse help (or see `PrintHelp()` function below).
#

SCRIPTNAME="$(basename "$0")"


# ------------------------------------------------------------------------------
function STDERR() { echo "$@" >&2 ; }
function ERROR() { STDERR "ERROR: $*" ; }
function FATAL() {
  local Code="$1"
  shift
  STDERR "FATAL(${Code}): $*"
  exit "$Code"
} # FATAL()


function isFlagSet() {
  local VarName="$1"
  [[ -n "${!VarName//0}" ]]
} # isFlagSet()


# ------------------------------------------------------------------------------
function PrintHelp() {
  cat <<EOH
Checks the exit code of each of the jobs from the specified XML files.
The XML files are listed in the file pointed by the \`XMLfileList\` argument.

Usage:  ${SCRIPTNAME}  [options]  XMLfileList

The jobs must have already been checked with \`project.py --check\` (but the
check may have reported failure because the check on no-input, no-ROOT-output
jobs is not well supported).


Options:

--help , -h , -?
    prints these usage instructions

EOH
} # PrintHelp()


### ----------------------------------------------------------------------------
function XMLoutputDir() {
  local XMLfile="$1"
  
  "$ProjectPY" --xml "$XMLfile" --outdir | grep -v '^Stage'
  
} # XMLoutputDir()


### ----------------------------------------------------------------------------
function ProcessXMLfile() {
  
  local -r ExitCodeFileName='larStage0.stat'

  local XMLfile="$1"
  
  local OutputDir="$(XMLoutputDir "$XMLfile" )"
  if [[ ! -d "$OutputDir" ]]; then
    ERROR "${XMLfile} : output directory not found ('${OutputDir}')"
    return 1
  fi
  
  local JobIDlist="${OutputDir}/jobids.list"
  if [[ ! -r "$JobIDlist" ]]; then
    ERROR "${XMLfile} : job ID list not found; have you run \`project.py --check\`? ('${JobIDlist}')"
    return 1
  fi
  
  local JobID
  local -i nJobs=0
  local -i nErrors=0
  while read JobID ; do
    let ++nJobs
    
    local JobNo="${JobID%@*}"
    local JobOutDir="${OutputDir}/${JobNo/./_}"
    if [[ ! -d "$JobOutDir" ]]; then
      ERROR "${XMLfile} ${JobID} : job output directory not found; might it be still running? ('${JobOutDir}')"
      let ++nErrors
      continue
    fi
    
    local ExitCodeFile="${JobOutDir}/${ExitCodeFileName}"
    if [[ ! -r "$ExitCodeFile" ]]; then
      ERROR "${XMLfile} ${JobID} : exit code file not found; might it be still running? ('${ExitCodeFile}')"
      let ++nErrors
      continue
    fi
    
    local JobExitCode="$(< "$ExitCodeFile")"
    if [[ $JobExitCode != "0" ]]; then
      ERROR "${XMLfile} ${JobID} : job exited with error code '${JobExitCode}'"
      let ++nErrors
    fi
    
    echo "${XMLfile} ${JobID}: exit code ${JobExitCode}"
  
  done < "$JobIDlist"
  
  if [[ $nJobs -gt 1 ]]; then
    if [[ $nErrors == 0 ]]; then
      echo "${XMLfile}: all ${nJobs} succeeded."
    else 
      ERROR "${XMLfile}: ${nErrors}/${nJobs} jobs present issues."
    fi
  fi
  
  [[ $nErrors == 0 ]]
} # ProcessXMLfile()


### ----------------------------------------------------------------------------
###
### argument parser
###

declare -i DoHelp=0
declare -i NoMoreOptions=0
declare -a Arguments
for (( iParam = 1 ; iParam <= $#; ++iParam )); do
  Param="${!iParam}"
  if isFlagSet NoMoreOptions || [[ "${Param:0:1}" != "-" ]]; then
    Arguments+=( "$Param" )
  else
    case "$Param" in
      ( '--help' | '-h' | '-?' ) DoHelp=1 ;;
      ( '--' | '-' ) NoMoreOptions=1 ;;
      ( * )
        FATAL 1 "Unsupported option #${iParam} ('${Param}'); run with \`--help\` for usage instructions."
        ;;
    esac
  fi
done

if isFlagSet DoHelp ; then
  echo
  PrintHelp
  exit 0
fi


[[ ${#Arguments[@]} == 0 ]] && FATAL 1 "The name of XML file list is mandatory; run with \`--help\` for usage instructions."

XMLfileList="${Arguments[0]}"

if [[ ${#Arguments[@]} -gt 1 ]]; then
  FATAL 1 "Arguments not supported: ${Arguments[@]}."
fi



###
### setup
###
ProjectPY="$(which 'project.py')"
[[ -x "$ProjectPY" ]] || FATAL 1 "Can't find \`project.py\` script."


###
### checks
###
[[ -r "$XMLfileList" ]] || FATAL 2 "Can't open the XML file list '${XMLfileList}'."

declare -i nXML=0
declare -i nErrors=0
while read XMLfile ; do
  let ++nXML
  ProcessXMLfile "$XMLfile" || let ++nErrors
  
done < "$XMLfileList"

if [[ $nErrors == 0 ]]; then
  echo "All ${nXML} job sets terminated without error."
else
  STDERR "${nErrors}/${nXML} job sets included jobs with errors (note: stderr stream contains only the error messages)."
  exit 1
fi

exit 0
# ------------------------------------------------------------------------------
