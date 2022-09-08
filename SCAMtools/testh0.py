import h0scam as h0

#----------------------------------
# Run from python prompt like this:
#  exec(open("./testh0.py").read())

#dir='/scratch/cluster/juliob/nCTOPb3_L58_080.0E_30.0N_2010-07-01_ENS/x_E03_L58_081.0E_33.0N_2010-07-01/run/'
#xp='x_E03_L58_081.0E_33.0N_2010-07-01'

#dir='/scratch/cluster/juliob/NCPL96_izumi_intel_nuopc_ENS/x_E01_L58_085.0E_32.0N_2010-04-01/run/'
#xp='x_E01_L58_085.0E_32.0N_2010-04-01'

#dir='/scratch/cluster/juliob/nCTOPb2_L58_080.0E_30.0N_2010-07-01/run/'
dir='/project/amp/juliob/scam/archive/nCTOPb2_L58_080.0E_30.0N_2010-07-01/atm/hist/'
xp='nCTOPb2_L58_080.0E_30.0N_2010-07-01'


xp='TEST_02_L58_180.0E_00.0_2010-01-01'
dir='/project/amp/juliob/scam/archive/' + xp + '/atm/hist/'

x = h0.h0scam(xp=xp,dir=dir)
x.basecase=xp

cu=x.curtain('U')
