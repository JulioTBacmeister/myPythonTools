import h0data as h0

#----------------------------------
# Run from python prompt like this:
#  exec(open("./testh0.py").read())


year=30
fld='TS'

dir='/glade/scratch/hannay/archive/'
xp='b.cesm3_cam058_mom_c.B1850WscMOM.ne30_L58_t061.009'

x = h0.h0(xp=xp,dir=dir)

ts = x.one_year_cube( fld=fld,year=year )

print("Out in testh0.py ..... ")
print(ts.dims)
