#!/bin/sh

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "To run run_prodcorsika_proton_intime_filter_bnb.sh"
   echo
   echo "Syntax: sh run_prodcorsika_proton_intime_filter_bnb.sh [-s -n -g -o -f | -h ]"
   echo "options:"
   echo "s     Input raw data file for overlays"
   echo "n     Number of data events to loop over and use for overlays"
   echo "g     Number of generated simulation events. Use a really big number to overcome the filter"
   echo "o     Name of output file from the sim info mixer"
   echo "f     Name of the output MC file resulting from the intime particle generator"
   echo "h     Print this Help menu."
   echo
}

############################################################
# Process the input options.                               #
############################################################

FCLFILE="prodcorsika_proton_intime_filter_bnb.fcl"
NEVT=1
GNEVT=100
GENFILE="genfile.root.local"
OUTPUTFILE="gen_output_file.root"
OUTGENFILE="mcgen.root"
INPUTDATAFILE="/path/to/file.root"

# Get the options
while getopts ":s:n:g:o:c:f:" option; do
   case $option in
      s) # Enter input data file
         INPUTDATAFILE=$OPTARG;;
      n) # Enter number of events to process
         NEVT=$OPTARG;;
      g) #Enter number of events to generate
	 GNEVT=$OPTARG;;
      o) #Enter the output file name
	 OUTPUTFILE=$OPTARG;;
      c) #Name of input fcl file
	 FCLFILE=$OPTARG;;
      f) #Name of output genfile
	 OUTGENFILE=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option"
	 Help
         exit;;
   esac
done

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
    echo "Running command 'lar -c local_gen.fcl -T ./genfile_hist.root -o $GENFILE -n $GNEVT > $MY_OUT_LOG 2> $MY_ERR_LOG' "
    lar -c local_gen.fcl -T ./genfile_hist.root -o $GENFILE -n $GNEVT > $MY_OUT_LOG 2> $MY_ERR_LOG
fi
echo "fcl_file: $FCLFILE; GNEVT: $GNEVT; NEVT: $NEVT; OUTPUTFILE: $OUTPUTFILE; INPUTDATAFILE: $INPUTDATAFILE"

echo "Running command 'lar -c prodcorsika_proton_intime_filter_bnb_siminfomixer.fcl -s $INPUTDATAFILE -T ./genfile_pot_hist.root -n $NEVT -o ${OUTPUTFILE}' "
lar -c prodcorsika_proton_intime_filter_bnb_siminfomixer.fcl -s $INPUTDATAFILE -T ./genfile_pot_hist.root -n $NEVT -o ${OUTPUTFILE}

cp $GENFILE $OUTGENFILE
