ln -nsf bkg_grp.root o2sim_grp.root
ln -nsf bkg_geometry.root o2sim_geometry.root

o2-tpc-reco-workflow --shm-segment-size 10000000000 --input-type digits --output-type clusters,tracks,send-clusters-per-sector  --configKeyValues "HBFUtils.orbitFirstSampled=0;HBFUtils.nHBFPerTF=256;GPU_global.continuousMaxTimeBin=100000;GPU_proc.ompThreads=4" -b
o2-its-reco-workflow --trackerCA --tracking-mode async --configKeyValues 'ITSVertexerParam.phiCut=0.5; ITSVertexerParam.clusterContributorsCut=3; ITSVertexerParam.tanLambdaCut=0.2' -b
o2-tpcits-match-workflow --shm-segment-size 10000000000 -b --use-ft0
o2-trd-tracklet-transformer --shm-segment-size 10000000000 -b
o2-trd-global-tracking --shm-segment-size 10000000000 -b
o2-tof-reco-workflow --shm-segment-size 10000000000 -b
o2-tof-matcher-workflow --shm-segment-size 10000000000 -b
o2-ft0-reco-workflow --shm-segment-size 10000000000 -b
o2-primary-vertexing-workflow --vertex-track-matching-sources ITS,TPC,ITS-TPC,TPC-TOF,ITS-TPC-TOF -b
o2-secondary-vertexing-workflow --vertexing-sources ITS,TPC,ITS-TPC,TPC-TOF,ITS-TPC-TOF --shm-segment-size 10000000000 -b
