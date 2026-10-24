fields,lons,lats,clevs = [], [], [], []
title,cmaps=[], []

nobjs=len(Al)
Arrangement=[ 3, nobjs ]

"""
if (fld=='objprec'):
    flev=np.linspace( 0.,0.4,num=41 )
    cmapN='Blues'
    fldT=f"{A['objName']} Object Prec."
elif (fld=='prec'):
    flev=np.linspace( 0.,1.6,num=17 )
    cmapN='Blues'
    fldT=f"Total Prec."
elif (fld=='objfreq'):
    flev=np.linspace( 0.,0.25,num=6 )
    cmapN='Reds'
    fldT=f"{A['objName']} Object frequency"
elif (fld=='objIprec'):
    flev=np.linspace( 0.,5,num=51 )
    cmapN='Blues'
    fldT=f"{A['objName']} In-Object Prec."
elif (fld=='objFprec'):
    flev=np.linspace( 0.,1.,num=6 )
    cmapN='Reds'
    fldT=f"{A['objName']} Object Prec. Fraction"
"""

for A in Al:


    if (fld=='objprec'):
        flev=np.linspace( 0.,0.4,num=41 )
        cmapN='Blues'
        fldT=f"{A['objName']} Object Prec."
    elif (fld=='prec'):
        flev=np.linspace( 0.,1.6,num=17 )
        cmapN='Blues'
        fldT=f"Total Prec."
    elif (fld=='objfreq'):
        flev=np.linspace( 0.,0.25,num=6 )
        cmapN='Reds'
        fldT=f"{A['objName']} Object frequency"
    elif (fld=='objIprec'):
        flev=np.linspace( 0.,5,num=51 )
        cmapN='Blues'
        fldT=f"{A['objName']} In-Object Prec."
    elif (fld=='objFprec'):
        flev=np.linspace( 0.,1.,num=6 )
        cmapN='Reds'
        fldT=f"{A['objName']} Object Prec. Fraction"
    




    
    lonL,latL = A['lonL'], A['latL']
    lono,lato = A['lono'], A['lato']
    lonh,lath = A['lonh'], A['lath']

    field_L=A[f'{fld}_L_tv']
    fields.append(field_L) #, objprec_o_tv]
    lons.append(lonL)
    lats.append(latL)
    clevs.append(flev)
    title.append( f"{fldT} CAM LR" )
    cmaps.append(cmapN)
    
    field_h=A[f'{fld}_h_tv']
    fields.append(field_h) #, objprec_o_tv]
    lons.append(lonh)
    lats.append(lath)
    clevs.append(flev)
    title.append( f"{fldT} CAM HR" )
    cmaps.append(cmapN)
    
    field_o=A[f'{fld}_o_tv']
    fields.append(field_o) #, objprec_o_tv]
    lons.append(lono)
    lats.append(lato)
    clevs.append(flev)
    title.append( f"{fldT} Obs" )
    cmaps.append(cmapN)


LL.Maps( fields=fields ,lons=lons ,lats=lats ,clevs=clevs, cmaps=cmaps, titles=title , Arrangement=Arrangement  )