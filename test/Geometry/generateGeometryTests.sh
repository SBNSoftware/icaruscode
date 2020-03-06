#!/usr/bin/env bash

# note that the default test is ignored by this definition:
# its name would be an empty string, and it is not enclosed in quote marks
declare -ar DefaultTests=(
  {_{single,split}_induction,}{_{no,}overburden,}
)

declare -r ConfigurationTag='_CONFIGURATION'
declare -r TemplateSuffix='.template'
declare -r TestConfigurationFile="test_geometry${ConfigurationTag}_icarus.fcl${TemplateSuffix}"

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
  local ConfigurationFile="$2" # unused
  echo "# test for geometry ${TestName}"
} # CMakeTestCommandHeader()


function CMakeTestCommand() {
  local TestName="$1"
  local ConfigurationFile="$2"
  cat <<EOC

$(CMakeTestCommandHeader "$TestName" "$ConfigurationFile")
cet_test(geometry${TestName}_icarus HANDBUILT
  TEST_EXEC lar
  TEST_ARGS --rethrow-all --config ${ConfigurationFile}
)
EOC
} # CMakeTestCommand()


function makeJobConfigurationName() {
  local -r TestName="$1"
  local TempPath="$TestConfigurationFile"

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

  local -r JobConfigurationFile="$(makeJobConfigurationName "$TestName")"
  echo " - creating configuration file: '${JobConfigurationFile}'"
  Exec sed -e "s/${ConfigurationTag}/${TestName}/g" "$TestConfigurationFile" > "$JobConfigurationFile"
  LASTFATAL "Failed to create job file."

  if [[ -w "$CMakeListsFile" ]]; then
    if grep -q "^$(CMakeTestCommandHeader "$TestName" "$JobConfigurationFile")$" "$CMakeListsFile" ; then
      echo " - test already present in '${CMakeListsFile}'"
    else
      echo " - adding the test to '${CMakeListsFile}'"
      Exec CMakeTestCommand "$TestName" "$JobConfigurationFile" >> "$CMakeListsFile"
    fi
  fi

} # addTest()



# ------------------------------------------------------------------------------
[[ -r "$TestConfigurationFile" ]] || FATAL 2 "Could not find the template test job configuration '${TestConfigurationFile}'!"

: ${CMakeListsDir:="$(dirname "$TestConfigurationFile")"}
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
