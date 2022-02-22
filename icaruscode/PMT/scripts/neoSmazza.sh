#!/usr/bin/env bash
#
# This script creates in the specified directory a series of configuration files
# (FHiCL and XML respectively for art and project.py) to be submitted to build
# an entire photon visibility library.
# All parameters are hard-coded into the "Script settings" section;
# the script only takes as optional argument an alternative output directory.
#
# This script produces:
# * a lot of FHiCL files (in their own subdirectory), one per job
# * a lot of XML files (in their own subdirectory), one per job
# * one XML file list with absolute paths
# * one `project.py` script with all XML files
#
#
# Job details
# ------------
#
# Each job fills one voxel per event
# (i.e. a job with 10 events is going to fill 10 jobs).
# 
# 
# Run instructions
# -----------------
# 
# The script can and should run without arguments. For debugging purposes,
# the following arguments are accepted:
# 
#   neoSmazza.sh  [BaseOutputDir] [CampaignTag] [ScriptOutputDir]
# 
# * `BaseOutputDir` is the directory where job output will be copied;
# * `CampaignTag` is a reference tag for the production (default: today date);
# * `ScriptOutputDir` is the directory where scripts and configuration are
#     stored (same as `BaseOutputDir` by default)
#

SCRIPTNAME="$(basename "$0")"
SCRIPTVERSION="1.7"


################################################################################
###  Script settings
################################################################################

#
# physics configuration: geometry and number of voxels (and voxels per job);
# current implementation uses bash to do the math: round to integers!
# 
# The current code (v09_00_00) can autodetect the size of the volume,
# but it would make precise numbers; here we prefer to have round numbers and
# potentially miss a bit of the volume close to the cryostat border;
# these numbers describe the (rounded) volume of the default ICARUS geometry
# in `icaruscode` `v09_00_00`.
# These figures are from `icarus_v3` geometry.
#
declare XMin="-395" # cm
declare XMax=" -25" # cm
declare YMin="-215" # cm
declare YMax=" 170" # cm
declare ZMin="-985" # cm
declare ZMax=" 985" # cm

declare Step="5" # cm
declare XStep="$Step"
declare YStep="$Step"
declare ZStep="$Step"

declare -i PhotonsPerVoxel=1000000

#
# Job configuration for the generation in a few voxels (template):
# icaruscode version and FHiCL name
#
declare -r ProductionVersion="${ICARUS_VERSION:-v09_37_01}"
declare -r ReferenceConfiguration='photonlibrary_builder_icarus.fcl'

declare -ir VoxelsPerJob=814 # estimate: 70s/ 1M photons

declare -r DefaultUserOutputDir="/pnfs/icarus/scratch/users/${USER}/jobOutput"
declare -r DefaultCampaignTag="$(date '+%Y%m%d')" # current date in format YYYYMMDD

#
# Technical details that may need update because world goes on
#
declare -r Qualifiers='e20:prof' # GCC 9.3.0
declare -r ExecutionNodeOS='SL7' # Scientific Linux [Fermi] 7
declare    ExpectedJobTime='48h'
declare -r ExpectedMemoryUsage='2000'
declare -r GeneratorLabel='generator'

#
# TEST settings
#
if [[ "${THISISATEST:-0}" != "0" ]]; then
  if [[ "$THISISATEST" == 1 ]]; then
    PhotonsPerVoxel=$((${PhotonsPerVoxel:-1000000} / 100)) # 1% of the regular job
  else
    PhotonsPerVoxel="$THISISATEST"
  fi
  ExpectedJobTime='8h'
fi

declare -a ExtraJobsubOptions=(
  # this magic asks for more modern CPU:
#   "--append_condor_requirements='(TARGET.CpuFamily=?=6 &amp;&amp; TARGET.CpuModelNumber=?=85)'" 
)

# be careful about single quotes inside the amenment
declare ConfigurationAmend=''

################################################################################
# internal variables; can change if you really want to
# 
ReferenceBaseName="$(basename "${ReferenceConfiguration%.fcl}")"
XMLlistName="${ReferenceBaseName}-xml.list"
SubmitScriptName="${ReferenceBaseName}-submit.sh"

CampaignTag="${2:-"${DefaultCampaignTag}"}"
BaseOutputDir="${1:-"${DefaultUserOutputDir}/${ReferenceBaseName}/${CampaignTag}"}"
ScriptOutputDir="${3:-${BaseOutputDir}}"
BaseJobOutputDir="${BaseOutputDir%/}/output"
FHiCLdir="${ScriptOutputDir%/}/fcl"
XMLdir="${ScriptOutputDir%/}/xml"
XMLlistPath="${ScriptOutputDir%/}/${XMLlistName}"
SubmitScriptPath="${ScriptOutputDir%/}/${SubmitScriptName}"

StageName='LibraryBuild'

UserName="$USER"
HostName="$(hostname)"


#
# do the math
#
declare XRange="$(( XMax - XMin ))"
declare YRange="$(( YMax - YMin ))"
declare ZRange="$(( ZMax - ZMin ))"
declare -i NStepsX="$(( XRange / XStep ))"
declare -i NStepsY="$(( YRange / YStep ))"
declare -i NStepsZ="$(( ZRange / ZStep ))"
declare -i TotalVoxels="$(( NStepsX * NStepsY * NStepsZ ))"

declare -i NJobs="$(( TotalVoxels / VoxelsPerJob ))"
declare -i LastJobVoxels="$(( TotalVoxels - VoxelsPerJob * NJobs ))"
if [[ $LastJobVoxels -gt 0 ]]; then
  let ++NJobs
else
  LastJobVoxels="$VoxelsPerJob" # just as all others
fi
LastJob=$((NJobs - 1))


################################################################################
###  Support functions
################################################################################
function STDERR() { echo "$@" >&2 ; }
function FATAL() {
  local -i Code="$1"
  export KEEPTEMP=1
  shift
  STDERR "FATAL (${Code}): $*"
  if [[ -r "$ErrorLog" ]]; then
    echo "FATAL (${Code}): $*" >> "$ErrorLog"
    STDERR "  (log file available at: '${ErrorLog}')"
  fi
  exit "$Code"
} # FATAL()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function ClearTemporary() {
  [[ -d "$MyTempDir" ]] || return
  if [[ "${KEEPTEMP:-0}" != 0 ]]; then
    # if the temporary directory is empty, we remove it anyway.
    rmdir -p "$MyTempDir" >& /dev/null || STDERR "WARNING: temporary files in '${MyTempDir}' have been left behind."
  else
    rm -Rf "$MyTempDir"
  fi
} # ClearTemporary()
trap ClearTemporary EXIT
declare -r MyTempDir="$(mktemp --directory --tmpdir "${SCRIPTNAME%.sh}-${USER}-tmpXXXXXX" )"


function ScheduleForTransfer() {
  # adds the specified file in the list of files to transfer to ScriptOutputDir;
  # if the path is relative to our temporary directory, its relative path is preserved
  
  [[ -n "$MyTempDir" ]] || FATAL 1 "Logic error: temporary directory not defined!"
  [[ -n "$TransferListTempFile" ]] || FATAL 1 "Logic error: transfer file not defined!"
  [[ -n "$ScriptOutputDir" ]] || FATAL 1 "Logic error: script output directory not defined!"
  
  local SourceFile RelSourcePath DestPath
  for SourceFile in "$@" ; do
    RelSourcePath="${SourceFile#${MyTempDir}/}"
    [[ "${RelSourcePath:0:1}" == '/' ]] && RelSourcePath="$(basename "$SourceFile")"
    DestPath="${ScriptOutputDir%/}/${RelSourcePath}"
  
    echo "${SourceFile} ${DestPath}" >> "$TransferListTempFile"
  done
  
} # ScheduleForTransfer()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function CounterPadding() {
  declare -i Counts="$1"
  declare -i LastCount="$((Counts - 1))"
  echo "${#LastCount}"
} # CounterPadding()

function PadCounter() {
  declare -i Value="$1"
  declare -i Padding="$2"
  printf "%0*d" "$Padding" "$Value"
} # PadCounter()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function Pager() {
  local -i Iter="$1"
  local -i NIters="$2"
  local -i NSteps="${3:-10}"
  
  if [[ "$Iter" -ge "$NIters" ]]; then
    STDERR " 100%"
    return
  fi
  
  local -i Step="$(( NIters / NSteps ))"
  [[ "$NIters" -gt "$(( NSteps * Step ))" ]] && let ++Step
  
  local -i Page="$(( Iter / Step ))"
  if [[ "$(( Page * Step ))" == "$Iter" ]]; then
    STDERR -n " $(( 100 * Page / NSteps ))% "
  fi
  STDERR -n '.'
  
} # Pager()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function CreateDirectoryFor() {
  local Dest
  for Dest in "$@" ; do
    
    if [[ "${Dest: -1}" == '/' ]]; then
      Dir="$Dest"
    else
      Dir="$(dirname "$Dest")"
    fi
    
    [[ -d "$Dir" ]] && continue
    mkdir "$Dir" || FATAL $? "Can't create a directory for '${Dest}'."
  done
} # CreateDirectoryFor()


function FHiCLtoXMLname() {
  local TemplateName="$1"
  echo "${TemplateName%.fcl}.xml"
} # FHiCLtoXMLname()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function CreateFHiCL() {
  local TemplateFHiCL="$1"
  local FirstVoxel="$2"
  local NVoxels="$3"
  local ConfigurationFileName="$4"
  local JobTag="$5"
  
  # ----------------------------------------------------------------------------
  # ---  BEGIN FCL HERE
  # ----------------------------------------------------------------------------
  local ConfigurationFileBaseName
  if [[ -n "$ConfigurationFileName" ]]; then
    ConfigurationFileBaseName="$(basename "${ConfigurationFileName%.fcl}")"
    cat <<EOH
#
# File: ${ConfigurationFileName}
EOH
  fi
  
  cat <<EOF
#
#  Configuration file for ${NVoxels} voxels starting from ${FirstVoxel}${JobTag:+" (${JobTag}${NJobs:+" / ${NJobs} jobs"})"}
#
#  Template configuration file: '${TemplateFHiCL}'
#  
#  Created on $(date) for ${ReferenceBaseName} (${CampaignTag}) by ${UserName} on ${HostName}
#  with ${SCRIPTNAME} version ${SCRIPTVERSION}
#  

#include "${TemplateFHiCL}"

${ConfigurationAmend:+"# ------------------------------------------------------------------------------"}
${ConfigurationAmend:+"# Additional configuration "}
${ConfigurationAmend:+"# -------------------------"}
${ConfigurationAmend}

# ------------------------------------------------------------------------------

#
# set the library range
#
services.PhotonVisibilityService: {
  @table::services.PhotonVisibilityService
  
  UseCryoBoundary:  false
  XMin: ${XMin} # cm
  XMax: ${XMax} # cm
  YMin: ${YMin} # cm
  YMax: ${YMax} # cm
  ZMin: ${ZMin} # cm
  ZMax: ${ZMax} # cm
  
  NX: ${NStepsX}  # ${XStep} cm voxels filling ${XRange} cm
  NY: ${NStepsY}  # ${YStep} cm voxels filling ${YRange} cm
  NZ: ${NStepsZ}  # ${ZStep} cm voxels filling ${ZRange} cm

} # services.PhotonVisibilityService


#
# set the range of voxels processed by this job
#

physics.producers.${GeneratorLabel}: {
  @table::physics.producers.${GeneratorLabel}
  
  N:          ${PhotonsPerVoxel}
  FirstVoxel: ${FirstVoxel}
  LastVoxel: -1 # keeps increasing at each event (limited by \`maxEvents\` below)
  
} # physics.producers.${GeneratorLabel}

source.maxEvents: ${NVoxels}
${ConfigurationFileBaseName:+"services.TFileService.fileName: '${ConfigurationFileBaseName}-PhotonLibraryData.root'"}

# ------------------------------------------------------------------------------

EOF
  # ----------------------------------------------------------------------------
  # ---  END FCL
  # ----------------------------------------------------------------------------
  
} # CreateFHiCL()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function CreateXML() {
  
  local -r JobOutputDir="$1"
  local -r JobConfiguration="$2"
  local -i NVoxels="$3"
  local JobTag="$4"
  
  local -r JobConfigurationName="$(basename "$JobConfiguration")"
  local -r JobConfigurationDir="$(dirname "$JobConfiguration")"
  
  
  # ----------------------------------------------------------------------------
  # ---  BEGIN XML HERE
  # ----------------------------------------------------------------------------
  cat <<EOF
<?xml version="1.0"?>

<!-- Production Project -->

<!DOCTYPE project [
<!ENTITY release      "${ProductionVersion}" >
<!ENTITY file_type    "mc"        >
<!ENTITY run_type     "physics"   >
<!ENTITY name         "${ReferenceBaseName}${JobTag:+"-${JobTag}"}" >
<!ENTITY output_base_name "${JobOutputDir}" >
<!ENTITY qualifiers   "${Qualifiers}"  >
]>

<job>

  <project name="&name;">

    <!-- Project size -->
    <numevents>${NVoxels}</numevents>

    <!-- Operating System -->
    <os>${ExecutionNodeOS}</os>

    <!-- Batch resources -->
    <resource>DEDICATED,OPPORTUNISTIC</resource>
    <cpu>1</cpu>
    <memory>${ExpectedMemoryUsage}</memory>
    <disk>5GB</disk>

    <!-- LArSoft information -->
    <larsoft>
      <tag>&release;</tag>
      <qual>&qualifiers;</qual>
    </larsoft>
    
    <fcldir>${JobConfigurationDir}</fcldir>

    <!-- Project stages -->

    <stage name="${StageName}">
      
      <fcl>${JobConfigurationName}</fcl>
      <outdir>&output_base_name;/out</outdir>
      <workdir>&output_base_name;/log</workdir>
      <numjobs>1</numjobs>
      <datatier>generated</datatier>
      
      <!-- AutoRelease+Grace: if expected life time is exceeded, ask for the job to be rescheduled with 1.5 more days -->
      <!-- jobsub>--expected-lifetime ${ExpectedJobTime} --lines '+FERMIHTC_AutoRelease=True' --lines '+FERMIHTC_GraceLifetime=129600'${ExtraJobsubOptions:+ "${ExtraJobsubOptions[@]}"}</jobsub -->
      <!-- AutoRelease+Grace: release, with the same time requirement -->
      <!-- jobsub>--expected-lifetime ${ExpectedJobTime} --lines '+FERMIHTC_AutoRelease=True' --lines '+FERMIHTC_GraceLifetime=0'${ExtraJobsubOptions:+ "${ExtraJobsubOptions[@]}"}</jobsub -->
      
      <!-- if the job has been held (JobStatus==5) for more than 10 minutes ((CurrentTime-EnteredCurrentStatus)>600) and has not restarted fewer than 4 times (NumJobStarts<4),
           then release it; this is Vito's magic -->
      <jobsub>--expected-lifetime ${ExpectedJobTime} --lines='+PeriodicRelease=(JobStatus==5)&amp;&amp;(NumJobStarts&lt;4)&amp;&amp;((CurrentTime-EnteredCurrentStatus)&gt;600)'</jobsub>
      
      <TFileName>&name;-PhotonLibraryData.root</TFileName>
      
      <!-- analysis job, i.e. do not look for a art ROOT file -->
      <ana>1</ana>
      <anadatatier>root-tuple</anadatatier>
      
    </stage>
    
    <!-- file type -->
    <filetype>&file_type;</filetype>

    <!-- run type -->
    <runtype>&run_type;</runtype>

  </project>

</job>
EOF
  # ----------------------------------------------------------------------------
  # ---  END XML
  # ----------------------------------------------------------------------------
} # CreateXML()


################################################################################
###  main program
################################################################################

cat <<EOM

=================================================
   ICARUS photon library configuration builder
=================================================

Physics configuration
-----------------------

Coverage:
  X:  ${XMin} -- ${XMax} cm in ${NStepsX} x ${XStep}-cm steps
  Y:  ${YMin} -- ${YMax} cm in ${NStepsY} x ${YStep}-cm steps
  Z:  ${ZMin} -- ${ZMax} cm in ${NStepsZ} x ${ZStep}-cm steps

Total voxels: ${TotalVoxels} in ${NJobs} jobs with ${VoxelsPerJob} voxels each (last: ${LastJobVoxels}).

Photons: ${PhotonsPerVoxel} / voxel


Job configuration
-------------------

Job name:          '${ReferenceBaseName}'
Campaign tag:      '${CampaignTag}'
Output directory:  '${BaseOutputDir}'
EOM

[[ "$BaseOutputDir" != "$ScriptOutputDir" ]] && echo "Script directory:  '${ScriptOutputDir}'"

cat <<EOM

XML file list:     '${XMLlistPath}'
Submission script: '${SubmitScriptPath}'

EOM


if [[ "${#ExtraJobsubOptions[@]}" -gt 0 ]]; then
  cat <<EOM

Extra jobsub options: ${ExtraJobsubOptions[@]}
EOM
fi


if [[ -n "${ConfigurationAmend// }" ]]; then
  cat <<EOM

Extra configuration:
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
${ConfigurationAmend}
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EOM
fi

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#
# preparation of directories and scripts
#

# start the log (will be deleted on successful termination);
# not everything is going to be logged...
declare ErrorLog="${MyTempDir}/${SCRIPTNAME%.sh}.log"
cat <<EOH > "$ErrorLog"
================================================================================
${SCRIPTNAME} log file started at $(date) by ${USER}

Working directory: $(pwd)

Command: $*

================================================================================
EOH


# test the transfer to destination:
# we don't want to find out we can't transfer only at the end!
mkdir -p "$ScriptOutputDir" || FATAL $? "Can't create output directory '${ScriptOutputDir}'."
TestFile='test.file'
rm -f "$TestFile" "${ScriptOutputDir}/${TestFile}"
echo "Test file ($(date))" > "$TestFile"
declare -a Cmd=( ifdh cp '-D' "$TestFile" "$ScriptOutputDir" )
echo -n "Testing IFDH file transfer..."
cat <<EOM >> "$ErrorLog"
$(date) - testing IFDH file transfer...
CMD> ${Cmd[@]}
EOM
"${Cmd[@]}" >> "$ErrorLog" 2>&1 || FATAL $? "Test transfer failed."
[[ -r "${ScriptOutputDir}/${TestFile}" ]] || FATAL 2 "Test transfer completed, but test file is not available at '${ScriptOutputDir}/${TestFile}'."
rm -f "$TestFile" "${ScriptOutputDir}/${TestFile}"
echo " done."
cat <<EOM >> "$ErrorLog"
$(date) - IFDH file transfer test successful
--------------------------------------------------------------------------------
CMD> ${Cmd[@]}
EOM

# this file holds the temporary files and their final destination;
# IFDH will take care of the transfer at the end of the script.
TransferListTempFile="${MyTempDir}/TransferFileList.txt"


XMLlistTempFile="${MyTempDir}/${XMLlistName}"
SubmitTempFile="${MyTempDir}/${SubmitScriptName}"
CreateDirectoryFor "${MyTempDir}/fcl/" "${MyTempDir}/xml/" || FATAL 2 "Can't create temporary directories!"
ScheduleForTransfer "$XMLlistTempFile" "$SubmitTempFile"

CreateDirectoryFor "${FHiCLdir}/" "${XMLdir}/" "$XMLlistPath" "$SubmitScriptPath"
AbsFHiCLdir="$(readlink -f "$FHiCLdir")"


cat <<EOH > "$SubmitTempFile"
#!/usr/bin/env bash
#
# Run with no argument to execute '--submit', or with the arguments \`to project.py\`.
#

declare -ar Arguments=( "\${@:-"--submit"}" )
declare -i nErrors=0

cat <<EOM
Running project.py command: \${Arguments[@]}
  on ${NJobs} jobs from ${XMLdir}

EOM

EOH

#
# job loop
#
echo -n "$(date) - Creating ${NJobs} job XML and FHiCL files:"
declare -ri JobPadding="$(CounterPadding "$NJobs")"

declare -i FirstVoxel=0
declare -i iJob=0
while [[ $FirstVoxel -lt $TotalVoxels ]]; do
  
  JobTag="job$(PadCounter "$iJob" "$JobPadding")"
  
  if [[ $iJob == $LastJob ]]; then
    JobVoxels="$LastJobVoxels"
  else
    JobVoxels="$VoxelsPerJob"
  fi
  
  #
  # create the configuration file
  #
  JobConfigurationFileName="${ReferenceBaseName}-${JobTag}.fcl"
  JobConfigurationTempPath="${MyTempDir}/fcl/${JobConfigurationFileName}"
  JobConfigurationPath="${ScriptOutputDir%/}/fcl/${JobConfigurationFileName}"
  JobConfigurationXMLTempPath="${MyTempDir}/xml/$(FHiCLtoXMLname "$JobConfigurationFileName")"
  JobConfigurationXML="${ScriptOutputDir%/}/xml/$(FHiCLtoXMLname "$JobConfigurationFileName")"
  JobOutputDir="${BaseJobOutputDir%/}/${JobTag}"
  
  # for fast debug of the loop:
#   echo "Job ${iJob}: ${JobVoxels} voxels from ${FirstVoxel} => ${JobConfigurationXML}"
  
  # dCache also does not support overwriting with redirection...
  rm -f "$JobConfigurationTempPath" "$JobConfigurationXMLTempPath"
  CreateFHiCL  "$(basename "$ReferenceConfiguration")" "$FirstVoxel" "$JobVoxels" "$JobConfigurationFileName" "$JobTag" > "$JobConfigurationTempPath"
  CreateXML "$JobOutputDir" "$JobConfigurationPath" "$JobVoxels" "$JobTag" > "$JobConfigurationXMLTempPath"
  
  ScheduleForTransfer "$JobConfigurationTempPath" "$JobConfigurationXMLTempPath"
  echo "$JobConfigurationXML" >> "$XMLlistTempFile"
  
  echo "echo '[${JobTag}] ${JobConfigurationXML}' ; project.py --xml \"${JobConfigurationXML}\" --stage \"${StageName}\" \"\${Arguments[@]}\" || let ++nErrors" >> "$SubmitTempFile"
  
  # print on screen a dot for each iteration
  Pager "$iJob" "$NJobs" '20'
  
  let ++iJob
  let FirstVoxel+=JobVoxels
  
#   [[ $iJob -ge 4 ]] && break # another debug shortcut
done

# some sanity checks
[[ $FirstVoxel == $TotalVoxels ]] || FATAL 1 "Stop everything!! there should be ${TotalVoxels} in total, but we covered ${FirstVoxel}!!"
[[ $iJob == $NJobs ]] || FATAL 1 "Stop everything!! there should be ${NJobs} in total, but we created ${iJob}!!"

Pager "$iJob" "$NJobs" '20' # complete paging output

#
# finish
#
cat <<EOF >> "$SubmitTempFile"

echo "Execution of ${NJobs} commands completed (\${nErrors} errors)."
[[ \$nErrors == 0 ]]
EOF

#
# transfer files to the final destination;
# in principle if IFDH is not available we could work around it with NFS/POSIX,
# but we are trying to avoid that.
# 
# while read Source Dest ; do
#   cp -a "$Source" "$Dest"
# done < "$TransferListTempFile"
#
#

# remove the old versions of the output files, if any
# (except for the XML and FHiCL)
rm -f "$SubmitScriptPath" "$XMLlistPath"


echo -n "$(date) - Transferring the files to the final destination (it may take a while)..."
Cmd=( ifdh cp '-f' "$TransferListTempFile" )
cat <<EOM >> "$ErrorLog"
$(date) - final file transfer
CMD> ${Cmd[@]}
EOM
"${Cmd[@]}" >> "$ErrorLog" 2>&1 || FATAL $? "Error transferring the files to '${ScriptOutputDir}'."
chmod a+x "$SubmitScriptPath"
echo " done."

cat <<EOM >> "$ErrorLog"
$(date) - final file transfer completed.
--------------------------------------------------------------------------------
EOM

echo "$(date) - all done."


################################################################################
