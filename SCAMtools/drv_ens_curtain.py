import glob
import h0scam as h0

dir='/scratch/cluster/juliob/NCPL96_izumi_intel_nuopc_ENS'
fl = sorted( glob.glob( dir +'/x_E*') )
nf = len( fl )

#fl=fl[0:10]
ird=0
for f in fl:
    fsp=f.split('/')
    xp=fsp[-1]
    ensdir=f+'/run'
    bens=fsp[-2]
    bsp=bens.split('_ENS')
    base=bsp[-2]
    print(base)
