
set -x

BKGEVENTS=10
SIGEVENTS=1000
INTERACTIONRATE=500000
NWORKERS=$(grep 'cpu[0-9]' /proc/stat | wc -l)

export O2DPG_ROOT=/home/fmazzasc/alice/O2DPG

# #generate background

o2-sim -j ${NWORKERS} -n ${BKGEVENTS} -g pythia8pp --skipModules MCH,MFT,ZDC -o bkg \
       --configFile $O2DPG_ROOT/MC/config/common/ini/basic.ini \
       | tee logbkg 2>&1    

# embed signal into background

o2-sim -j ${NWORKERS} -n ${SIGEVENTS} -g boxgen --skipModules MCH,MFT,ZDC -o sgn \
       --configKeyValues 'BoxGun.pdg=1010010030;BoxGun.eta[0]=-0.9;BoxGun.eta[1]=0.9;BoxGun.prange[0]=0.5;BoxGun.prange[1]=10.' \
       --embedIntoFile bkg_Kine.root \
       # | tee logsgn 2>&1


o2-steer-colcontexttool -i bkg,${INTERACTIONRATE:-400000},${SIGEVENTS}:${BKGEVENTS} sgn,@0:e1 --show-context
o2-sim-digitizer-workflow -b --sims bkg,sgn --incontext collisioncontext.root
o2-trd-trap-sim -b