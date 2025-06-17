#!/bin/bash

# set -x  # Enable debug mode (prints executed commands)
# set -e  # Exit on error

# env | grep prox
# unset http_proxy https_proxy
# env | grep prox

date

#01- token check
alien-token-info &> tokenStatus.log
if [[ $(tail -1 tokenStatus.log) == "File >>><<< not found" ]];then  
   echo "ERROR :: Token Not Found, Create the token" 
   echo "$(date) ERROR :: Token Not Found, Create the token" >> tokenError.log
   exit 1  #stop execution if token not found
else
   echo "Token Found, Processing further ..."
fi

starttime=$(date)
rm -rf AnalysisResults.root
rm -rf AO2D_Derived.root

OutFile2="trkQADerive_output-02.log"
OutFile3="trkQADerive_output-03.log"

# filePath="InputFilePath"
# filePath="/alice/data/2024/LHC24ar/559903/apass1/0430/o2_ctf_run00559903_orbit0248688992_tf0001310015_epn191/009/AO2D.root"
# ROOT_File="Run_55990n_AO2D.root"
ROOT_File="/home/rahul/DataFiles/newData/Skimmed-PbPb_LHC24ar_559903_apass1_0430_o2_ctf_run00559903_orbit0248688992_tf0001310015_epn191_009_AO2D.root"
# echo "FILE DOWNLOAD START"
# time alien.py cp $filePath"@ALICE::CERN::EOS" file::$ROOT_File
# echo "FILE DOWNLOAD OVER"
echo "PROCESSING"
 INPUT_JSON="trackQADeriveInput_2025.05.23.json"
WRITER_JSON="trackQADeriveOutput_2025.05.23.json"

export OPTIONS=" -b --configuration json://$INPUT_JSON  --resources-monitoring 2  --aod-memory-rate-limit 1000000000  --shm-segment-size 7500000000"

o2-analysis-centrality-table      ${OPTIONS} |\
o2-analysis-trackselection        ${OPTIONS} |\
o2-analysis-track-to-collision-associator ${OPTIONS} |\
o2-analysis-hf-track-index-skim-creator ${OPTIONS} |\
o2-analysis-hf-candidate-creator-2prong ${OPTIONS} |\
o2-analysis-hf-pid-creator ${OPTIONS} |\
o2-analysis-hf-candidate-selector-d0 ${OPTIONS} |\
o2-analysis-lf-lambdakzerobuilder ${OPTIONS} |\
o2-analysis-trackqa-converter-003 ${OPTIONS} |\
o2-analysis-pid-tof-full          ${OPTIONS} |\
o2-analysis-pid-tpc-base          ${OPTIONS} |\
o2-analysis-pid-tpc               ${OPTIONS} |\
o2-analysis-pid-tof-beta          ${OPTIONS} |\
o2-analysis-pid-tof-base          ${OPTIONS} |\
o2-analysis-ft0-corrected-table   ${OPTIONS} |\
o2-analysis-event-selection       ${OPTIONS} |\
o2-analysis-multiplicity-table    ${OPTIONS} |\
o2-analysis-timestamp             ${OPTIONS} |\
o2-analysis-track-propagation     ${OPTIONS} |\
o2-analysis-occ-table-producer    ${OPTIONS} |\
o2-analysis-track-qa-derivedata   ${OPTIONS} --aod-file $ROOT_File --aod-writer-json $WRITER_JSON \
&> $OutFile2

# valgrind \
#   --tool=callgrind \
#   --callgrind-out-file=callgrind.trackqa.out \

# valgrind \
#   --tool=memcheck \
#   --leak-check=full \
#   --show-leak-kinds=all \
#   --track-origins=yes \
#   --log-file=valgrind-trackqa.log \


# rm -rf $ROOT_File
# echo "FILE DELETED"

# set +x  # Disable debug mode

echo $starttime >> $OutFile2
date >> $OutFile2
code $OutFile2

echo "Printing the Errors and Warnings " &> $OutFile3 
grep -e "IST" -e "\\[ERROR\\]" -e "\\[FATAL\\]" -e "segmentation" -e "Segmentation" -e "SEGMENTATION" -e "command not found" -e "Error:" -e "Error in " -e "\\[WARN\\]" -e "DEBUG" $OutFile2  >> $OutFile3
code $OutFile3
date

for i in $(seq 1 10); do echo -en "\a" ; sleep 0.25; done

echo "PROCESSING STATUS::COMPLETED"
