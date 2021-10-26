#!/usr/bin/env bash

set -x

MODULES="PIPE ITS TPC TRD TOF"
BKGEVENTS=200
SIGEVENTS=1000
NWORKERS=12

export O2DPG_ROOT=$PWD

# generate background

o2-sim -j ${NWORKERS} -n ${BKGEVENTS} -g pythia8pp -m ${MODULES} -o bkg \
       --configFile ${O2DPG_ROOT}/MC/config/common/ini/basic.ini \
       > logbkg 2>&1

# embed signal into background

o2-sim -j ${NWORKERS} -n ${SIGEVENTS} -g boxgen -m ${MODULES} -o sgn \
       --configKeyValues 'BoxGun.pdg=1010010030;BoxGun.eta[0]=-0.9;BoxGun.eta[1]=0.9;BoxGun.prange[0]=0.5;BoxGun.prange[1]=10.' \
       --embedIntoFile bkg_Kine.root \
       > logsgn 2>&1
o2-sim-digitizer-workflow --sims sgn -b
o2-trd-trap-sim -b
