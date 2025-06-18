#!/usr/bin/env python
# Import packages 
import sys
workdir_ = '/glade/work/juliob'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} ")


import argparse as arg
import time
import numpy as np

# import modules in other directories
from PyRegridding.Regridder import GenRegrid as GnR
from PyRegridding.Regridder import Initialize as Init
from PyRegridding.Regridder import WriteDST as Wrt
from PyRegridding.Regridder import ReadInSrc as Rd 

import importlib
importlib.reload( Rd )

#Rdry = Con.Rdry() # 



def main(year, month, day, hour, Dst, DstVgrid, Src, IC_for_pg, RegridMethod = 'CONSERVE'):
    import calendar
    
    tic_total = time.perf_counter()
    days_in_month = calendar.monthrange(year,month)[1]

    print( f"About to process {year:04d}-{month:02d}-{day:02d}")

    #RegridMethod = 'CONSERVE'
    lnPS=False
    if(lnPS==True):
        ver='lnPS'
    else:
        ver=''

    # Override this when it can't be true
    if (Dst not in ('ne480np4','ne240np4','ne120np4','ne30np4')):
        IC_for_pg = False
        print(f'  Setting IC_for_pg={IC_for_pg} because cannot be otherwise for {Dst} ' )
    else: 
        print(f'  {IC_for_pg}: This is making an IC file for {Dst} ' )


    ret1 = Init.prep(Dst=Dst, DstVgrid=DstVgrid ,Src=Src, WOsrf=True, RegridMethod=RegridMethod , IC_for_pg=IC_for_pg )
    sys.stdout.flush()
    if (day==99):
        for iday in np.arange( days_in_month):
            ret2 = Rd.get_Src( year=year ,month=month ,day=iday+1 , hour0=99 )
            sys.stdout.flush()
            ret3 = GnR.xRegrid(HorzInterpLnPs=lnPS )
            sys.stdout.flush()
            ret4 = Wrt.write_netcdf(version=ver) #+'Test01')

    else:
        ret2 = Rd.get_Src( year=year ,month=month ,day=day , hour0=hour )
        sys.stdout.flush()
        ret3 = GnR.xRegrid(HorzInterpLnPs=lnPS )
        sys.stdout.flush()
        ret4 = Wrt.write_netcdf(version=ver) #+'Test01')
        
    code = 1
    toc_total = time.perf_counter()

    pTime = f"Total processing time was  {toc_total - tic_total:0.4f} seconds"
    print(pTime)

if __name__ == "__main__":
    
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--month",    type=int, default=1)
    my_parser.add_argument("--year",     type=int, default=2010)
    args = my_parser.parse_args()
    main(args.year, args.month )
