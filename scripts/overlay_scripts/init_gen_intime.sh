#!/bin/sh

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "To run keepup_def_making.sh"
   echo
   echo "Syntax: sh keepup_def_making.sh [-s -v -f -l -m | -h ]"
   echo "options:"
   echo "s     Enter stage name. Default is 'raw'."
   echo "v     Enter software version. Default is 'v09_56_00_01'"
   echo "f     Enter starting run number for range. Default is '9301'"
   echo "l     Enter ending run number for range. Deafult is '9500'"
   echo "m     Enter 'y' to make definitions. Default is 'n'"
   echo "h     Print this Help menu."
   echo
}

############################################################
# Process the input options.                               #
############################################################


FCLFILE="prodcorsika_proton_intime_icarus_bnb_sce_on.fcl"
NEVT=50
GNEVT=2500
GENFILE="genfile.root.local"
OUTPUTFILE="gen_1234.root"
OUTGENFILE="mcgen.root"
INPUTDATAFILE="/path/to/file.root"

MY_OUT_LOG=larInitGen.out
MY_ERR_LOG=larInitGen.err

# Get the options
while getopts ":h:s:n:g:o:s:c:f:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      s) # Enter input data file
         INPUTDATAFILE=$OPTARG;;
      n) # Enter number of events to process
         NEVT=$OPTARG;;
      g) #Enter number of events to generate
	    GNEVT=$OPTARG;;
      o) #Enter the output file name
	    OUTPUTFILE=$OPTARG;;
      s) #name of input data file
	    INPUTDATAFILE=$OPTARG;;
      c) #Name of input fcl file
	    FCLFILE=$OPTARG;;
      f) #Name of output genfile
	    OUTGENFILE=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done



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

echo "Running command 'lar -c intime_gen_overlay_SimInfoMixer1.fcl -s $INPUTDATAFILE -T ./genfile_pot_hist.root -n -1' "
echo "Running command 'lar -c intime_gen_overlay_SimInfoMixer2.fcl -s *_simmxd1.root -o $OUTPUTFILE -T ./genfile_pot_hist.root -n -1' "
lar -c intime_gen_overlay_SimInfoMixer1.fcl -s $INPUTDATAFILE -T ./genfile_pot_hist.root -n $NEVT
lar -c intime_gen_overlay_SimInfoMixer2.fcl -s *_simmxd1.root -o $OUTPUTFILE -T ./genfile_pot_hist.root -n $NEVT 
cp $GENFILE $OUTGENFILE
rm *_simmxd1.root
