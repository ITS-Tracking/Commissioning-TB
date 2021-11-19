#!/usr/bin/env bash

set -x

MODULES="PIPE ITS TPC TRD TOF FT0 FV0"
BKGEVENTS=1000
SIGEVENTS=1000
NWORKERS=$(grep 'cpu[0-9]' /proc/stat | wc -l)

export O2DPG_ROOT=/home/fmazzasc/alice/O2DPG

generate background

o2-sim -j ${NWORKERS} -n ${BKGEVENTS} -g pythia8pp -m ${MODULES} -o bkg \
       --configFile $O2DPG_ROOT/MC/config/common/ini/basic.ini \
       | tee logbkg 2>&1    

# embed signal into background

o2-sim -j ${NWORKERS} -n ${SIGEVENTS} -g boxgen -m ${MODULES} -o sgn \
       --configKeyValues 'BoxGun.pdg=1010010030;BoxGun.eta[0]=-0.9;BoxGun.eta[1]=0.9;BoxGun.prange[0]=0.5;BoxGun.prange[1]=10.' \
       --embedIntoFile bkg_Kine.root \
       | tee logsgn 2>&1
o2-sim-digitizer-workflow --sims bkg,sgn -b 
o2-trd-trap-sim -b
