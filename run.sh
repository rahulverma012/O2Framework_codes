ROOT_File="Skimmed-Run_55990n_AO2D.root"
export OPTIONS=" --configuration json://myConfig.json --resources-monitoring 2 --aod-memory-rate-limit 1000000000  --shm-segment-size 7500000000"
o2-analysis-event-selection       ${OPTIONS} |\
o2-analysis-multiplicity-table    ${OPTIONS} |\
o2-analysis-occ-table-producer    ${OPTIONS} |\
o2-analysis-timestamp             ${OPTIONS} |\
o2-analysis-track-propagation     ${OPTIONS} |\
o2-analysis-occupancy-table-consumer   ${OPTIONS} --aod-file $ROOT_File
