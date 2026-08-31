#!/bin/bash
#
unset PYTHONHOME
unset PYTHONPATH
rm -f tritonserver_nugraph-v0.log
export APPTAINER_BIND=/etc/hosts,/tmp,/cvmfs
#
export BASEPORT=8000
export NPORTS=3
#
export HTTPPORT=$BASEPORT
export GRPCPORT=$((BASEPORT+1))
export METRPORT=$((BASEPORT+2))
#
if [ -z "$FCL" ]; then 
    export FCL=$1
fi
#
while 2>/dev/null >"/dev/tcp/0.0.0.0/$BASEPORT"; do
    BASEPORT=$((BASEPORT+NPORTS))
    export HTTPPORT=$BASEPORT
    export GRPCPORT=$((BASEPORT+1))
    export METRPORT=$((BASEPORT+2))
done
#
echo "physics.producers.NuGraphCryoE.TritonConfig.serverURL: 'localhost:${GRPCPORT}'" >> "$FCL"
echo "physics.producers.NuGraphCryoW.TritonConfig.serverURL: 'localhost:${GRPCPORT}'" >> "$FCL"
#
/cvmfs/oasis.opensciencegrid.org/mis/apptainer/1.3.2/bin/apptainer exec --pid --ipc --home ~/:${HOME} --pwd ${PWD} /cvmfs/icarus.opensciencegrid.org/containers/tritonserver/nugraph-v0/ tritonserver --model-repository /triton-server-config/models --http-port=${HTTPPORT} --grpc-port=${GRPCPORT} --metrics-port=${METRPORT} >& tritonserver_nugraph-v0.log &
#
grep -q "Started" <(tail -f tritonserver_nugraph-v0.log)
