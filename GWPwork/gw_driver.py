#!/usr/bin/env python
# Import packages 
import sys
import argparse as arg

import importlib
import glob
import copy
#import time
import os 
import subprocess as sp

from datetime import datetime, timedelta
import calendar

import numpy as np

#From in here
import gw_intr as GWi

importlib.reload( GWi )

def day_of_year_to_date(year, day_of_year):
    # Create a date object for January 1st of the given year
    jan_first = datetime(year, 1, 1)
    # Add the day_of_year to this date (subtract 1 because day_of_year is 1-based)
    target_date = jan_first + timedelta(days=day_of_year - 1)
    return target_date


def Driver(year,month,SourceMethod):

    days_in_month = calendar.monthrange(year, month)[1]

    for day in np.arange( start=1,stop=days_in_month+1):
        for hour in np.arange( start=0,stop=24 , step=6):
            rcode = GWi.tau_prof(date=[year,month,day,hour],SourceMethod=SourceMethod )


if __name__ == "__main__":
    
    #####################
    # Src     = 'ne30pg3'
    # Dst     = 'fv0.9x1.25'
    # hsPat = 'cam.h0'
    


    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--year", type=int, required=True, help="Year as a four-digit integer")
    my_parser.add_argument("--month", type=int, required=True, help="Month as a two-digit integer")
    my_parser.add_argument("--SourceMethod", type=str, default="vort500")
    args = my_parser.parse_args()

    print( args.year , args.month )
    Driver( year=args.year , month=args.month , SourceMethod=args.SourceMethod )
