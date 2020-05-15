#!/usr/bin/env bash

# note that the default test is ignored by this definition:
# its name would be an empty string, and it is not enclosed in quote marks
declare -ar DefaultTests=(
  {_{single,split}_induction,}{_{no,}overburden,}
)

declare -r ConfigurationTag='_CONFIGURATION'
declare -r TemplateSuffix='.template'
declare -rA TestConfigurationFiles=(
  ["geometry test"]="test_geometry${ConfigurationTag}_icarus.fcl${TemplateSuffix}"
  ["geometry dump"]="dump${ConfigurationTag}_icarus_geometry.fcl${TemplateSuffix}"
  ["channel mapping dump"]="dump${ConfigurationTag}_icarus_channelmap.fcl${TemplateSuffix}"
)

# ------------------------------------------------------------------------------
function STDERR() { echo "$*" >&2 ; }
function WARNING() { STDERR "WARNING: $*" ; }
function FATAL() {
  local Code="$1"
  shift
  STDERR "FATAL (${Code}): $*"
  exit $Code
} # FATAL()
function LASTFATAL() {
  local res="$?"
  [[ "$res" == 0 ]] || FATAL "$res" "$@"
} # LASTFATAL()

function Exec() {
  local -a Cmd=( "$@" )
  if [[ -n "${FAKE//0}" ]]; then
    STDERR "FAKE> ${Cmd[@]}"
  else
    "${Cmd[@]}"
  fi
} # Exec

# ------------------------------------------------------------------------------
function CMakeTestCommandHeader() {
  local TestName="$1"
  local TestType="$2"
  local ConfigurationFile="$3" # unused
  echo "# ${TestType} for geometry ${TestName}"
} # CMakeTestCommandHeader()


function CMakeTestCommand() {
  local TestName="$1"
  local TestType="$2"
  local ConfigurationFile="$3"
  cat <<EOC

$(CMakeTestCommandHeader "$TestName" "$TestType" "$ConfigurationFile")
cet_test(geometry${TestName}_${TestType// /_}_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ${ConfigurationFile}
)
EOC
} # CMakeTestCommand()


function makeJobConfigurationName() {
  local -r TestName="$1"
  local -r TemplatePath="$2"
  local TempPath="$TemplatePath"

  # keep basename
  TempPath="$(basename "${TempPath}")"

  # remove suffix
  TempPath="${TempPath%${TemplateSuffix}}"

  # replace configuration tag
  TempPath="${TempPath//${ConfigurationTag}/${TestName}}"

  # all done
  echo "$TempPath"
} # makeJobConfigurationName()


function addTest() {
  local TestName="$1"

  echo "Test: '${TestName}'"

  local TestType
  for TestType in "${!TestConfigurationFiles[@]}" ; do
    TemplatePath="${TestConfigurationFiles["$TestType"]}"
    local JobConfigurationFile="$(makeJobConfigurationName "$TestName" "$TemplatePath")"
    echo " - ${TestType}: creating configuration file: '${JobConfigurationFile}'"
    Exec sed -e "s/${ConfigurationTag}/${TestName}/g" "$TemplatePath" > "$JobConfigurationFile"
    LASTFATAL "Failed to create job file."

    if [[ -w "$CMakeListsFile" ]]; then
      if grep -q "^$(CMakeTestCommandHeader "$TestName" "$TestType" "$JobConfigurationFile")$" "$CMakeListsFile" ; then
        echo " - ${TestType}: test already present in '${CMakeListsFile}'"
      else
        echo " - ${TestType}: adding the test to '${CMakeListsFile}'"
        Exec CMakeTestCommand "$TestName" "$TestType" "$JobConfigurationFile" >> "$CMakeListsFile"
      fi
    fi
  done

} # addTest()



# ------------------------------------------------------------------------------
for TemplatePath in "${TestConfigurationFiles[@]}" ; do
  [[ -r "$TemplatePath" ]] || FATAL 2 "Could not find the template test job configuration '${TemplatePath}'!"
done

: ${CMakeListsDir:="$(dirname "$0")"}
declare -r CMakeListsFile="${CMakeListsDir%/}/CMakeLists.txt"

if [[ ! -r "$CMakeListsFile" ]]; then
  WARNING "File '${CMakeListsFile}' not found: tests will not be added."
elif [[ ! -w "$CMakeListsFile" ]]; then
  WARNING "File '${CMakeListsFile}' not writeable: tests will not be added."
fi

declare -a Tests=( "$@" )
[[ "${#Tests[@]}" == 0 ]] && Tests=( "${DefaultTests[@]}" )

for TestName in "${Tests[@]}" ; do
  addTest "$TestName"
done
