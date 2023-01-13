# Tests

Overview of tests that should be feasible/benefitial to include.
These'll be organised per 'module'. Note that capitalisation of the 'Test' class is required to have it picked up.


# Modules

## demux - test_demux.py

 - [ ] misMatcher
 - [x] detMask
 - [x] hamming2Mismatch
 - [ ] writeDemuxSheet (use filecmp for this ?)
 - [x] readDemuxSheet

## drHouse - test_drHouse.py

 - [ ] matchOptdupsReqs

## misc - test_misc.py

- [x] getConf
- [ ] getNewFlowCell
- [x] parseRunInfo
- [x] hamming
- [x] joinLis
- [x] lenMask
- [x] P5Seriesret
- [ ] krakenfqs
- [x] moveOptDup
- [x] retBCstr
- [x] retIxtype
- [x] retMean_perc_Q
- [x] formatSeqRecipe
- [x] formatMisMatches
- [ ] fetchLatestSeqDir
- [x] umlautDestroyer
- [x] validateFqEnds
- [ ] matchingSheets

## postmux - test_postmux.py

 - [ ] matchIDtoName
 - [ ] renamefq
 - [ ] renameProject
