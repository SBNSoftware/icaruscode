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

function isFlagUnset() {
  local VarName="$1"
  [[ -z "${!VarName//0}" ]]
} # isFlagUnset()


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

--goodlist=LISTNAME [from the XML file list name]
--badlist=LISTNAME [from the XML file list name]
    create a file list of all job configuration XML files which succeeded
    and one of the jobs which included at least one failure
--outputlist [from the XML file list name]
    create a file list of all output files detected from the job which fully
    succeeded (the ones which would end in the \`--goodlist\` list)
--nogoodlist , -G
--nobadlist , -B
--nooutputlist , -O
    do not fill a list of good jobs, bad jobs, or output files (see options
    \`--goodlist\`, \`--badlist\` and \`--outputlist\` respectively
--skipgood[=LISTNAME] [same as in \`--goodlist\`]
--skipbad[=LISTNAME] [same as in \`--badlist\`]
    if a XML file is already in an existing good or bad list, the check is not
    repeated (but the output file list is still regenerated)

--maxjobs=LIMIT
    if specified, processes at most LIMIT input jobs (mostly debugging purpose)

--help , -h , -?
    prints these usage instructions

EOH
} # PrintHelp()


### ----------------------------------------------------------------------------
declare -a TemporaryFileList=( )

function CreateTemporaryFile() {
  # 
  # CreateTemporaryFile  VarName [Tag]
  # 
  # creates and registers a temporary file; the file name is stored in a
  # environment variable with name VarName; registered temporary files
  # will be removed on exit
  # 
  local VarName="$1"
  local -r Tag="$2"
  local TempFile="$(mktemp --tmpdir "${SCRIPTNAME%.sh}${Tag:+"-${Tag}"}-temp.XXXXXX")" || return $?
  
  eval "$VarName"="$TempFile"
  
  echo "${!VarName}"
  TemporaryFileList+=( "${!VarName}" )
  return 0
  
} # CreateTemporaryFile()


function CleanTemporaryFiles() {
  # removes all registered temporary files, no error if removal fails
  
  rm -f "${TemporaryFileList[@]}"
  
} # CleanTemporaryFiles()

trap CleanTemporaryFiles EXIT


### ----------------------------------------------------------------------------
function XMLoutputDir() {
  local XMLfile="$1"
  
  "$ProjectPY" --xml "$XMLfile" --outdir | grep -v '^Stage'
  
} # XMLoutputDir()


### ----------------------------------------------------------------------------
function FindOutputFiles() {
  local -r JobOutDir="$1"
  
  # we want the TFileService file; heuristic here is that that one is closed
  # after the RootOutput files (which should not be there anyway)
  ls -t "${JobOutDir}/"*.root 2> /dev/null | head -n 1
  
  return ${PIPESTATUS[0]} # error if `ls` failed
} # FindOutputFiles()


### ----------------------------------------------------------------------------
function ExtractCampaignName() {

  local -r XMLfileList="$1"
  
  local temp="$(basename "$XMLfileList")"
  temp="${temp%.list}"
  temp="${temp%.filelist}"
  temp="${temp%-xml}"
  echo "$temp"
  
} # ExtractCampaignName()


### ----------------------------------------------------------------------------
function IsJobKnownGood() {
  local XMLfile="$1"
  local GoodXMLfile="$1"
  for GoodXMLfile in "${KnownGoodJobs[@]}" ; do
    [[ "$XMLfile" == "$GoodXMLfile" ]] && return 0
  done
  return 1
} # IsJobKnownGood()


function IsJobKnownBad() {
  local XMLfile="$1"
  local BadXMLfile="$1"
  for BadXMLfile in "${KnownBadJobs[@]}" ; do
    [[ "$XMLfile" == "$BadXMLfile" ]] && return 0
  done
  return 1
} # IsJobKnownBad()


### ----------------------------------------------------------------------------
function ProcessXMLfile() {
  
  local -r ExitCodeFileName='larStage0.stat'

  local XMLfile="$1"
  
  local GoodBoy=0
  while true ; do
    if IsJobKnownBad "$XMLfile" ; then
      ERROR "${XMLfile}: already known to have failed."
      [[ -n "$FailureXMLlist" ]] && echo "$XMLfile" >> "$FailureXMLlist"
      break
    fi
    
    IsJobKnownGood "$XMLfile" && GoodBoy=1
    
    local OutputDir="$(XMLoutputDir "$XMLfile" )"
    if isFlagUnset GoodBoy && [[ ! -d "$OutputDir" ]]; then
      ERROR "${XMLfile} : output directory not found ('${OutputDir}')"
      break
    fi
    
    local JobIDlist="${OutputDir}/jobids.list"
    if isFlagUnset GoodBoy && [[ ! -r "$JobIDlist" ]]; then
      ERROR "${XMLfile} : job ID list not found; have you run \`project.py --check\`? ('${JobIDlist}')"
      break
    fi
    
    local JobID
    local -a ROOToutputFiles
    local -i nJobs=0
    local -i nErrors=0
    while read JobID ; do
      let ++nJobs
      
      local JobExitCode
      local JobNo="${JobID%@*}"
      local JobOutDir="${OutputDir}/${JobNo/./_}"
      if isFlagUnset GoodBoy ; then
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
        
        JobExitCode="$(< "$ExitCodeFile")"
      else
        JobExitCode=0
      fi
      
      if [[ $JobExitCode != "0" ]]; then
        ERROR "${XMLfile} ${JobID} : job exited with error code '${JobExitCode}'"
        let ++nErrors
      else
        
        if [[ -n "$OutputFileList" ]]; then
          ROOToutputFiles=( "$(FindOutputFiles "$JobOutDir")" )
          if [[ $? == 0 ]] && [[ "${#ROOToutputFiles[@]}" -gt 0 ]]; then
            for ROOToutputFile in "${ROOToutputFiles[@]}" ; do
              echo "$ROOToutputFile"
            done >> "$OutputFileList"
          else
            # if output file list filling is requested, we require an output file
            # to be detected, or else it's an error no matter the exit code value
            ERROR "${XMLfile} ${JobID} : no suitable output files detected in '${JobOutDir}'"
            let ++nErrors
          fi
        fi # if output file list requested
      fi # if exit code != 0 ... else
      
      echo "${XMLfile} ${JobID}: exit code ${JobExitCode}${OutputFileList:+" (${#ROOToutputFiles[@]} output files)"}"
      
    done < "$JobIDlist"
    
    if [[ $nErrors == 0 ]]; then
      GoodBoy=1
      if [[ $nJobs -gt 1 ]]; then
        echo "${XMLfile}: all ${nJobs} succeeded."
      fi
    else 
      GoodBoy=0
      if [[ $nJobs -gt 1 ]]; then
        ERROR "${XMLfile}: ${nErrors}/${nJobs} jobs present issues."
      fi
      
    fi
    
    break
  done # fake loop
  
  if isFlagSet GoodBoy ; then
    [[ -n "$SuccessXMLlist" ]] && echo "$XMLfile" >> "$SuccessXMLlist"
  else
    [[ -n "$FailureXMLlist" ]] && echo "$XMLfile" >> "$FailureXMLlist"
  fi
  
  [[ $nErrors == 0 ]]
} # ProcessXMLfile()


### ----------------------------------------------------------------------------
###
### argument parser
###

declare -i DoHelp=0 NoSuccessXMLlist=0 NoFailureXMLlist=0 NoOutputFileList=0
declare -i SkipGoodJobs=0 SkipBadJobs=0
declare -i MaxProcessedJobs=0
declare -i NoMoreOptions=0
declare -a Arguments
for (( iParam = 1 ; iParam <= $#; ++iParam )); do
  Param="${!iParam}"
  if isFlagSet NoMoreOptions || [[ "${Param:0:1}" != "-" ]]; then
    Arguments+=( "$Param" )
  else
    case "$Param" in
      ( '--goodlist='* )   FinalSuccessXMLlist="${Param=--*=}" ;;
      ( '--badlist='* )    FinalFailureXMLlist="${Param=--*=}" ;;
      ( '--outputlist='* ) FinalOutputList="${Param=--*=}" ;;
      ( '--nogoodlist' | '-G' )   NoSuccessXMLlist=1 ;;
      ( '--nobadlist' | '-B' )    NoFailureXMLlist=1 ;;
      ( '--nooutputlist' | '-O' ) NoOutputFileList=1 ;;
      ( '--skipgood='* )   OldGoodJobList="${Param#--*=}" ; SkipGoodJobs=1 ;;
      ( '--skipgood' )     SkipGoodJobs=1 ;;
      ( '--skipbad='* )    OldBadJobList="${Param#--*=}" ; SkipBadJobs=1 ;;
      ( '--skipbad' )      SkipBadJobs=1 ;;
      ( '--maxjobs='* )    MaxProcessedJobs="${Param#--*=}" ;;
      ( '--help' | '-h' | '-?' )  DoHelp=1 ;;
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

# lists are created in temporary files and then copied to the final destination
# at the end, to avoid problems if the final destination is on a file system
# which does not support appending text to a file (like dCache's)

declare -r XMLfileListDir="$(dirname "$XMLfileList")"
declare -r CampaignName="$(ExtractCampaignName "$XMLfileList")"

declare -r DefaultSuccessXMLlist="${XMLfileListDir}/${CampaignName}-goodxml.list"
declare -r DefaultFailureXMLlist="${XMLfileListDir}/${CampaignName}-badxml.list"

if isFlagSet NoSuccessXMLlist ; then
  [[ -z "$FinalSuccessXMLlist" ]] || FATAL 1 "Do you want an good job list, or not?!? (requested both \`--goodlist=\"${FinalSuccessXMLlist}\"\` and \`--nogoodlist\`)"
elif [[ -z "$SuccessXMLlist" ]]; then
  FinalSuccessXMLlist="$DefaultSuccessXMLlist"
fi
[[ -n "$FinalSuccessXMLlist" ]] && CreateTemporaryFile SuccessXMLlist 'goodxml'

if isFlagSet NoFailureXMLlist ; then
  [[ -z "$FinalFailureXMLlist" ]] || FATAL 1 "Do you want an bad job list, or not?!? (requested both \`--badlist=\"${FinalFailureXMLlist}\"\` and \`--nobadlist\`)"
elif [[ -z "$FailureXMLlist" ]]; then
  FinalFailureXMLlist="$DefaultFailureXMLlist"
fi
[[ -n "$FinalFailureXMLlist" ]] && CreateTemporaryFile FailureXMLlist 'badxml'

if isFlagSet NoOutputFileList ; then
  [[ -z "$FinalOutputList" ]] || FATAL 1 "Do you want a output file list, or not?!? (requested both \`--outputlist=\"${FinalOutputList}\"\` and \`--nooutputlist\`)"
elif [[ -z "$OutputFileList" ]]; then
  FinalOutputList="${XMLfileListDir}/${CampaignName}-files.list"
fi
[[ -n "$FinalOutputList" ]] && CreateTemporaryFile OutputFileList 'files'


declare -a KnownGoodJobs
if isFlagSet SkipGoodJobs ; then
  [[ -n "$OldGoodJobList" ]] || OldGoodJobList="${FinalSuccessXMLlist:-${DefaultSuccessXMLlist}}"
  [[ -r "$OldGoodJobList" ]] || FATAL 2 "Can't read the list of known good jobs ('${OldGoodJobList}')"
  readarray -t KnownGoodJobs < "$OldGoodJobList"
  echo "${#KnownGoodJobs[@]} jobs already checked as good."
fi

declare -a KnownBadJobs
if isFlagSet SkipBadJobs ; then
  [[ -n "$OldBadJobList" ]] || OldBadJobList="${FinalFailureXMLlist:-${DefaultFailureXMLlist}}"
  [[ -r "$OldBadJobList" ]] || FATAL 2 "Can't read the list of known bad jobs ('${OldBadJobList}')"
  readarray -t KnownBadJobs < "$OldBadJobList"
  echo "${#KnownBadJobs[@]} jobs already checked as bad."
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
  if [[ "$MaxProcessedJobs" -gt 0 ]] && [[ "$nXML" -gt "$MaxProcessedJobs" ]]; then
    echo "Reached the limit of ${MaxProcessedJobs} jobs."
    nXML="$MaxProcessedJobs"
    break
  fi
  ProcessXMLfile "$XMLfile" || let ++nErrors
done < "$XMLfileList"

if [[ $nErrors == 0 ]]; then
  echo "All ${nXML} job sets terminated without error."
else
  STDERR "${nErrors}/${nXML} job sets included jobs with errors (note: stderr stream contains only the error messages)."
fi


# renaming and print the name of the file lists (if available)
if [[ -n "$FinalSuccessXMLlist" ]]; then
  rm -f "$FinalSuccessXMLlist"
  if [[ -s "$SuccessXMLlist" ]]; then
    cp "$SuccessXMLlist" "$FinalSuccessXMLlist"
    echo "Good job configuration list written into '${FinalSuccessXMLlist}' ($(wc -l < "$FinalSuccessXMLlist") entries)"
  else
    echo "Good job configuration list would be empty: skipped."
  fi
fi
if [[ -n "$FinalFailureXMLlist" ]]; then
  rm -f "$FinalFailureXMLlist"
  if [[ -s "$FailureXMLlist" ]]; then
    cp "$FailureXMLlist" "$FinalFailureXMLlist"
    echo "Bad job configuration list written into '${FinalFailureXMLlist}' ($(wc -l < "$FinalFailureXMLlist") entries)"
  else
    echo "Bad job configuration list would be empty: skipped."
  fi
fi
if [[ -n "$FinalOutputList" ]]; then
  rm -f "$FinalOutputList"
  if [[ -s "$OutputFileList" ]]; then
    cp "$OutputFileList" "$FinalOutputList"
    echo "Output file list written into '${FinalOutputList}' ($(wc -l < "$FinalOutputList") entries)"
  else
    echo "Output file list would be empty: skipped."
  fi
fi


[[ $nErrors == 0 ]] # exit code
# ------------------------------------------------------------------------------
