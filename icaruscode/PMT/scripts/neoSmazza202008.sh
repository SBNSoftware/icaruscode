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
SCRIPTVERSION="1.2"


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
#
declare XMin="-405" # cm
declare XMax=" -35" # cm
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
declare -r ProductionVersion="${ICARUSCODE_VERSION:-"v09_00_00"}"
declare -r ReferenceConfiguration='photonlibrary_builder_icarus.fcl'

declare -ir VoxelsPerJob=1850 # estimate: 1'/ 1M photons

declare -r DefaultUserOutputDir="/pnfs/icarus/scratch/users/${USER}/jobOutput"
declare -r DefaultCampaignTag="$(date '+%Y%m%d')" # current date in format YYYYMMDD

#
# Technical details that may need update because world goes on
#
declare -r Qualifiers="${MRB_QUALS:-"e19:prof"}" # default: GCC 8.2.0
declare -r ExecutionNodeOS='SL7' # Scientific Linux [Fermi] 7
declare    ExpectedJobTime='72h'
declare -r ExpectedMemoryUsage='2000'
declare -r GeneratorLabel='generator'

#
# TEST settings
#
if [[ "${THISISATEST:-0}" != "0" ]]; then
  PhotonsPerVoxel=10000
  ExpectedJobTime='8h'
fi


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
  shift
  STDERR "FATAL (${Code}): $*"
  exit "$Code"
} # FATAL()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
declare -a TempFiles
function ClearTemporary() {
  local TempFile
  for TempFile in "${TempFiles[@]}" ; do
    rm -f "$TempFile"
  done
} # ClearTemporary()

function CreateTemporary() {
  
  local Pattern="$1"
  
  [[ "${#TempFiles[@]}" == 0 ]] && trap ClearTemporary EXIT
  
  local TempFile="$(mktemp --tmpdir ${Pattern} )"
  TempFiles+=( "$TempFile" )
  echo "$TempFile"
  
} # CreateTemporary()


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
function FHiCLtoXMLname() {
  local TemplateName="$1"
  echo "${TemplateName%.fcl}.xml"
} # FHiCLtoXMLname()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function CreateFHiCL() {
  
  local TemplateFHiCL="$1"
  local FirstVoxel="$2"
  local NVoxels="$3"
  local JobTag="$4"
  
  cat <<EOF
#
#  Configuration file for ${NVoxels} voxels starting from ${FirstVoxel}${JobTag:+" (${JobTag})"}
#
#  Template configuration file: '${TemplateFHiCL}'
#  
#  Created on $(date) for ${ReferenceBaseName} (${CampaignTag}) by ${UserName} on ${HostName}
#  with ${SCRIPTNAME} version ${SCRIPTVERSION}
#  

#include "${TemplateFHiCL}"

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


# ------------------------------------------------------------------------------

EOF
  
} # CreateFHiCL()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function CreateXML() {
  
  local -r JobOutputDir="$1"
  local -r JobConfiguration="$2"
  local -i NVoxels="$3"
  local JobTag="$4"
  
  local -r JobConfigurationName="$(basename "$JobConfiguration")"
  local -r JobConfigurationDir="$(dirname "$JobConfiguration")"
  
  
  # INPUT DATACARD
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
      
      <jobsub>--expected-lifetime ${ExpectedJobTime}</jobsub>
      
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#
# preparation of directories and scripts
#

# dCache does not support appending output to a file: we build the list locally and then we move it later
XMLlistTempFile="$(mktemp --tmpdir "${SCRIPTNAME%.sh}-${XMLlistName}.tmpXXXXXX")"
SubmitTempFile="$(mktemp --tmpdir "${SCRIPTNAME%.sh}-${SubmitScriptName}.tmpXXXXXX")"

if [[ ! -d "$FHiCLdir" ]]; then
  mkdir -p "$FHiCLdir" || FATAL 3 "Can't create FHiCL output directory '${FHiCLdir}'"
fi
AbsFHiCLdir="$(readlink -f "$FHiCLdir")"

if [[ ! -d "$XMLdir" ]]; then
  mkdir -p "$XMLdir" || FATAL 3 "Can't create XML output directory '${XMLdir}'"
fi
rm -f "$XMLlistPath"
touch "$XMLlistPath" || FATAL 2 "Can't create XML file list '${XMLlistPath}'!"

rm -f "$SubmitScriptPath"
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
  JobConfigurationFile="${AbsFHiCLdir}/${JobConfigurationFileName}"
  JobConfigurationXML="${XMLdir}/$(FHiCLtoXMLname "$JobConfigurationFileName")"
  JobOutputDir="${BaseJobOutputDir%/}/${JobTag}"
  
  # for fast debug of the loop:
#   echo "Job ${iJob}: ${JobVoxels} voxels from ${FirstVoxel} => ${JobConfigurationXML}"
  
  # dCache also does not support overwriting with redirection...
  rm -f "$JobConfigurationFile" "$JobConfigurationXML"
  CreateFHiCL "$(basename "$ReferenceConfiguration")" "$FirstVoxel" "$JobVoxels" "$JobTag" > "$JobConfigurationFile"
  CreateXML "$JobOutputDir" "$JobConfigurationFile" "$JobVoxels" "$JobTag" > "$JobConfigurationXML"
  
  echo "$JobConfigurationXML" >> "$XMLlistTempFile"
  
  echo "echo '[${JobTag}] ${JobConfigurationXML}' ; project.py --xml \"${JobConfigurationXML}\" --stage \"${StageName}\" \"\${Arguments[@]}\" || let ++nErrors" >> "$SubmitTempFile"
  
  # print on screen a dot for each iteration
  Pager "$iJob" "$NJobs" '20'
  
  let ++iJob
  let FirstVoxel+=JobVoxels
  
#   [[ $iJob -ge 4 ]] && break # another debug shortcut
done

#
# finish
#
cat <<EOF >> "$SubmitTempFile"

echo "Execution of ${NJobs} commands completed (\${nErrors} errors)."
[[ \$nErrors == 0 ]]
EOF

mkdir -p "$(dirname "$XMLlistPath")"
mv "$XMLlistTempFile" "$XMLlistPath"
mkdir -p "$(dirname "$SubmitScriptPath")"
mv "$SubmitTempFile" "$SubmitScriptPath"
chmod a+x "$SubmitScriptPath"

# some sanity checks
[[ $FirstVoxel == $TotalVoxels ]] || FATAL 1 "Stop everything!! there should be ${TotalVoxels} in total, but we covered ${FirstVoxel}!!"
[[ $iJob == $NJobs ]] || FATAL 1 "Stop everything!! there should be ${NJobs} in total, but we created ${iJob}!!"

Pager "$iJob" "$NJobs" '20'


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

exit
echo -n "Let's start the empty analisys  "
     Init
     QUANTI=-1
     PASSO=784
     while [ $QUANTI -lt 2340 ] ; do
	 echo $REC0 
         QUANTI=`expr $QUANTI + 1`
         #FirstVoxel=`expr $QUANTI * $PASSO `
         FirstVoxel=$(( QUANTI * PASSO ))
         RIQUANTI=`expr $QUANTI + 1`
         LastVoxel=$(( RIQUANTI * PASSO))
         echo $QUANTI $FirstVoxel $LastVoxel
	 RunNames
	 echo processing event for $OUT0 $OUT1
         cp -pvi icarus_prodsingle_buildopticallibrary.fcl test_job.fcl
         echo "physics.producers.generator.FirstVoxel: $FirstVoxel" >> test_job.fcl
         echo "physics.producers.generator.LastVoxel: $LastVoxel" >> test_job.fcl
	 mv -vf test_job.fcl $OUT1
	 CreateXML
	 mv -vf build_library.xml $OUT0 
         echo "project.py --xml " $OUT0 " --stage Lib --submit " >> lancio_processi.sh
         echo Now sleeping 1 second
         sleep $DELTAT
     done
 echo DONE



################################################################################
