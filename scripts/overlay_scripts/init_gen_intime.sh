#!/bin/sh
############################################################
# Process the input options.                               #
############################################################


:'
FCLFILE=prodcorsika_proton_intime_filter_bnb.fcl
NEVTS=100
GENFILE=genfile.root.local
OUTPUTFILE=outputfile.root
INPUTDATAFILE=inputdata.root


# Get the options
while getopts ":h:f:n:g:o:s:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      f) # Enter a stage name
         FCLFILE=$OPTARG;;
      n) # Enter software version
         NEVTS=$OPTARG;;
      g) #Enter start run number
	    GENFILE=$OPTARG;;
      o) #Enter end run number
	    OUTPUTFILE=$OPTARG;;
      s) #make defs?
	    INPUTDATAFILE=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done
'

FCLFILE=${1}
NEVTS=${2}
GENFILE=genfile.root.local
OUTPUTFILE=${3}
INPUTDATAFILE=${4}

MY_OUT_LOG=larInitGen.out
MY_ERR_LOG=larInitGen.err

echo "#include \"$FCLFILE\"" > local_gen.fcl
#echo "physics.producers.generator.FluxCopyMethod: \"IFDH\"" >> local_gen.fcl
#echo "physics.producers.generator.MaxFluxFileMB: 500" >> local_gen.fcl
echo "services.IFDH: {}" >> local_gen.fcl
cat local_gen.fcl

echo "#include \"$FCLFILE\"" > local_gen_include.fcl
#echo "physics.producers.generator.FluxCopyMethod: \"IFDH\"" >> local_gen_include.fcl
#echo "physics.producers.generator.MaxFluxFileMB: 500" >> local_gen_include.fcl
echo "services.IFDH: {}" >> local_gen_include.fcl
echo "gen_detail: { physics: {@table::physics} services: {@table::services} outputs: {@table::outputs} source: {@table::source} process_name: @local::process_name }" >> local_gen_include.fcl



if [ -f $GENFILE ]; then
    echo "File $GENFILE exists, so not generating again."
else
    echo "File $GENFILE does not exist, so running generation of it..."
    echo "Running command 'lar -c local_gen.fcl -T ./genfile_hist.root -o $GENFILE -n $NEVTS > $MY_OUT_LOG 2> $MY_ERR_LOG' "
    lar -c local_gen.fcl -T ./genfile_hist.root -o $GENFILE -n $NEVTS > $MY_OUT_LOG 2> $MY_ERR_LOG
fi
echo "Running command 'lar -c intime_gen_overlay_SimInfoMixer1.fcl -s $INPUTDATAFILE -T ./genfile_pot_hist.root -n -1' "
echo "Running command 'lar -c intime_gen_overlay_SimInfoMixer2.fcl -s *_simmxd1.root -T ./genfile_pot_hist.root -n -1' "
lar -c intime_gen_overlay_SimInfoMixer1.fcl -s $INPUTDATAFILE -T ./genfile_pot_hist.root -n -1
lar -c intime_gen_overlay_SimInfoMixer2.fcl -s *_simmxd1.root -o $OUTPUTFILE -T ./genfile_pot_hist.root -n -1 
rm *_simmxd1.root
