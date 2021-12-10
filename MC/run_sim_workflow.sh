#!/usr/bin/env bash

ln -nsf bkg_grp.root o2sim_grp.root
ln -nsf bkg_geometry.root o2sim_geometry.root

export O2DPG_ROOT=/home/fmazzasc/alice/O2DPG



# # ----------- LOAD UTILITY FUNCTIONS --------------------------
. ${O2_ROOT}/share/scripts/jobutils.sh 

RNDSEED=42
NSIGEVENTS=100
NBKGEVENTS=5
NWORKERS=8
NTIMEFRAMES=1

${O2DPG_ROOT}/MC/bin/o2dpg_sim_workflow.py -e TGeant3 -eCM 13000 -gen boxgen -j ${NWORKERS} -ns ${NSIGEVENTS} -tf ${NTIMEFRAMES} -mod "--skipModules ZDC" \
	-confKey "BoxGun.pdg=1010010030;BoxGun.eta[0]=-0.9;BoxGun.eta[1]=0.9;BoxGun.prange[0]=0.5;BoxGun.prange[1]=10."  \
        -genBkg pythia8 -procBkg inel -colBkg pp --embedding -nb ${NBKGEVENTS}


# # run workflow
${O2DPG_ROOT}/MC/bin/o2_dpg_workflow_runner.py -f workflow.json