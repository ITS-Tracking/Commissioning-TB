ConfPar='ITSVertexerParam.phiCut=0.5;ITSVertexerParam.clusterContributorsCut=3;ITSVertexerParam.tanLambdaCut=0.2;'
# ITSCATrackerParam.useDiamond=true; ITSCATrackerParam.diamondPos[2]=35.; ITSCATrackerParam.pvRes=10.'
if [ ! -f done ]; then
	if [ ! -f o2clus_its_.root ]; then
		o2-ctf-reader-workflow --onlyDet ITS  --delay 0.01 --ctf-input data.lst --shm-segment-size 24000000000  --max-cached-files 10  | o2-its-reco-workflow --clusters-from-upstream --trackerCA --tracking-mode sync_misaligned --disable-mc --configKeyValues "$ConfPar"  --autosave 10 --run -b && mv o2clus_its.root o2clus_its_.root
	else
		o2-its-cluster-reader-workflow --its-cluster-infile o2clus_its_.root | o2-its-reco-workflow --clusters-from-upstream --trackerCA --tracking-mode sync_misaligned --disable-mc --configKeyValues "$ConfPar"  --autosave 10 -b 
	fi
	touch done					
else
	echo "Clusters and tracks reconstruction done here."
fi

if [ ! -f pv_done ]; then
	o2-primary-vertexing-workflow --vertex-track-matching-sources ITS --vertexing-sources ITS -b --disable-mc --hbfutils-config none --configKeyValues "HBFUtils.orbitFirstSampled=295044346" --run
	touch pv_done
else
	echo "Primary vertex reconstructed"
fi
