#!/usr/bin/env python

import getopt as go
import scam_case as scm
import sys
import os

argv=sys.argv

case=scm.scam_case()
user=os.environ['USER']
spawncase=False
IOPcase=False

print(case.scmlon)

try:
   opts, args = go.getopt( argv[1:], "i:j:y:m:d:t:l:x:n:q:c:S:M:NI:", 
                           ["lon=","lat=","year=","month=","day=","tag=","nlev=","coupler=","nsteps="
                            ,"atm-ncpl=","compiler=","spawn=","machine=","NameByBuild=","IOP="] )
except:
    print( "something is wrong")
    exit()



##./scam_drv.py -i 180. -j 0 -t TEST_07 -n 69120 -m 1 -y 2010

date=case.startdate
print(date)

for opt, arg in opts:
    if opt in ("-i","--lon"):
        case.scmlon = float(arg)
    elif opt in ("-j","--lat"):
        case.scmlat = float(arg)
    elif opt in ("-y","--year"):
        case.startdate[0] = int(arg)
    elif opt in ("-m","--month"):
        case.startdate[1]  = int(arg)
    elif opt in ("-d","--day"):
        case.startdate[2]  = int(arg)
    elif opt in ("-l","--nlev"):
        case.nlev = int(arg)
    elif opt in ("-t","--tag"):
        case.tag = arg
    elif opt in ("-x","--coupler"):
        case.coupler = arg
    elif opt in ("-n","--nsteps"):
        case.nsteps = int(arg)
    elif opt in ("-q","--atm-ncpl"):
        case.atm_ncpl = int(arg)
    elif opt in ("-c","--compiler"):
        case.compiler = arg
    elif opt in ("-M","--machine"):
        case.machine = arg
    elif opt in ("-S","--spawn"):
        basecase = arg
        spawncase = True
    elif opt in ("-N","--NameByBuild"):
        case.NameByBuild=True
    elif opt in ("-I","--IOP"):
        case.IOP = arg
        IOPcase = True


date=case.startdate
print(date)

if (spawncase == False) and (IOPcase == False):
   case.base_case()
elif (spawncase == True) and (IOPcase == False):
   case.spawn_case(basecase)
elif (spawncase == False) and (IOPcase == True):
   case.IOP_case()

#fname = '../../cases/'+case_tag+'/'+'CaseInst.pysave'
#with open( fname, 'wb') as fob:
#   pickle.dump( case, fob )
