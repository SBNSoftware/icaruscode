#!/bin/bash

unset PYTHONHOME
unset PYTHONPATH
rm -f tritonserver_nugraph-v0.log
export APPTAINER_BIND=/etc/hosts,/tmp,/cvmfs
/cvmfs/oasis.opensciencegrid.org/mis/apptainer/1.3.2/bin/apptainer run --pid --ipc --home ~/:${HOME} --pwd ${PWD} /cvmfs/icarus.opensciencegrid.org/containers/tritonserver/nugraph-v0/ >& tritonserver_nugraph-v0.log &
grep -q "Started" <(tail -f tritonserver_nugraph-v0.log)
