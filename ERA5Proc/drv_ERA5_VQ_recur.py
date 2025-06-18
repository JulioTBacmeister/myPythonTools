#!/usr/bin/env python
# Import packages 
import os
import argparse as arg

from mpi4py import MPI


import subprocess as sp
import update_config as uc
import make_ERA5_VQ as mVQ

######################################################################
# This function is called by PyBatch_ERA5regrid.csh, and
# then may also resubmit PyBatch_ERA5regrid.csh after incrementing
# month and decrementing Resubmit in config_ERA5regrid file.
#####################################################################

def main():

    file_path = './config_ERA5_VQ.yaml'  # Specify the path to your config file

    config = uc.read_config_yaml( file_path )
    print( config , flush=True )

    modes = config['modes']
    
    if ('hourly' in modes):
        mVQ.main( year=config['year'] , month=config['month'] )
    if ('daily' in modes ):
        mVQ.daily( year=config['year'] , month=config['month'] )
    if ('monthly' in modes):
        mVQ.monthly( year=config['year'] , month=config['month'] )
        
    
    #------------------------------
    config = uc.increment_month( config ) #, NoLeapYear=True )

    config = uc.decrement_Resubmit( config )
    print( config )
    uc.write_config_yaml(file_path, config)
   
    if ( config['Resubmit']>=0 ):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        if (rank == 0):
            print(f" Resubmitting myself through PyBatch ... .csh  ")
            
            sp.run(f"qsub PyBatch_ERA5_VQ.csh", 
                   shell=True )
            print(f"PyBatch ... " )
        
if __name__ == "__main__":
    main()
