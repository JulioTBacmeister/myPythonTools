load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
    SUBREGION = False   ; Whether to regrid to a subregion of the FV grid

    if(SUBREGION) then
      wgt_file = "CAMSE_2_FV_subregion.nc"
      minlat   =  20       ; Select your subregion here
      maxlat   =  60
      minlon   = 235
      maxlon   = 295

      print("------------------------------------------------------------")
      print("Will be regridding to a subregion of the finite volume grid")
      print("------------------------------------------------------------")
    else
      wgt_file = "/glade/derecho/scratch/fudan/CAM7/CAM7_SE_NE120_to_latlon_721x1440.bilinear.nc"
      print("------------------------------------------------------------")
      print("Will be regridding to the full finite volume grid")
      print("------------------------------------------------------------")
    end if 

fili = systemfunc("ls f*h9i.2000*00.nc")
nfili = dimsizes(fili)
do nf=0,nfili-1
print("Processing " + fili(nf))

;---Data file containing source grid
  sfile = addfile(fili(nf),"r")
  afile = addfile(str_sub_str(fili(nf),"h9i","h9a"),"r")
;--- Read time dimension
  time = sfile->time
  tlen = dimsizes(time)
  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)
  utc_date = cd_calendar(time, 0)
  year   = tointeger(utc_date(:,0))    ; Convert to integer for
  month  = tointeger(utc_date(:,1))    ; use sprinti 
  day    = tointeger(utc_date(:,2))
  hour   = tointeger(utc_date(:,3))
  date_str = sprinti("%0.4i", year) + sprinti("%0.2i", month) + \
             sprinti("%0.2i", day)
  
;---Get variable to regrid
    src_lat  = sfile->lat
    src_lon  = sfile->lon

    
    U850_in   = sfile->U850 ;(0,:)              ; (time,ncol)
    V850_in   = sfile->V850 ;(0,:)              ; (time,ncol)
    U200_in   = sfile->U200 ;(0,:)              ; (time,ncol)
    V200_in   = sfile->V200 ;(0,:)              ; (time,ncol)
    T850_in   = sfile->T850 ;(0,:)              ; (time,ncol)
    Q850_in   = sfile->Q850 ;(0,:)              ; (time,ncol)
    Z500_in   = sfile->Z500 ;(0,:)              ; (time,ncol)
    PSL_in   = sfile->PSL ;(0,:)              ; (time,ncol)
    FLUT_in   = sfile->FLUT ;(0,:)              ; (time,ncol)
    uIVT_in   = sfile->uIVT ;(0,:)              ; (time,ncol)
    vIVT_in   = sfile->vIVT ;(0,:)              ; (time,ncol)
    PRECT_in   = afile->PRECT ;(0,:)              ; (time,ncol)
    PRECT_in   = PRECT_in*1000*3600 ; unit: mm/hr

;---Call the regridding function
    longitude = ESMF_regrid_with_weights(src_lon,wgt_file,False)
    mlon = 1440
    dlon = 360. / mlon
    lon = longitude(0,:)
    lon  = fspan (0, (mlon - 1) * dlon, mlon)    ;
    latitude = ESMF_regrid_with_weights(src_lat,wgt_file,False)
    mlat = 720
    dlat = 180. / mlat
    lat = latitude(:,0)
    lat  = fspan (-90, 90, mlat+1)    ;

    
    U850_rg = ESMF_regrid_with_weights(U850_in,wgt_file,False)
    V850_rg = ESMF_regrid_with_weights(V850_in,wgt_file,False)
    U200_rg = ESMF_regrid_with_weights(U200_in,wgt_file,False)
    V200_rg = ESMF_regrid_with_weights(V200_in,wgt_file,False)
    T850_rg = ESMF_regrid_with_weights(T850_in,wgt_file,False)
    Q850_rg = ESMF_regrid_with_weights(Q850_in,wgt_file,False)
    Z500_rg = ESMF_regrid_with_weights(Z500_in,wgt_file,False)
    PRECT_rg = ESMF_regrid_with_weights(PRECT_in,wgt_file,False)
    FLUT_rg = ESMF_regrid_with_weights(FLUT_in,wgt_file,False)
    PSL_rg = ESMF_regrid_with_weights(PSL_in,wgt_file,False)
    uIVT_rg = ESMF_regrid_with_weights(uIVT_in,wgt_file,False)
    vIVT_rg = ESMF_regrid_with_weights(vIVT_in,wgt_file,False)  

    PRECT_rg@units = "mm/hr"

    sigma=5.67*1E-8; % Stefan–Boltzmann constant; 5.67 E-8 W m-2 K-4;
    a=1.228;
    b=-1.106*1E-3; % -1.106*1E-3 K-1
    Tf=(FLUT_rg/sigma)^(0.25);
    Tb=(-a+sqrt(a^2+4*b*Tf))/(2*b); % flux equivalent brightness temperature in K
    Tb@long_name = "flux equivalent brightness temperature in K"
    Tb@units = "K"
    copy_VarCoords(FLUT_rg,Tb)  

    ;shift_positions = 180
    ;Tb = shift_data(Tb, shift_positions, 1) ;

 
    Y = conform_dims((/dimsizes(lat), dimsizes(lon)/), lat, 0)  
    Y!0 = "lat"
    Y!1 = "lon"
    X = conform_dims((/dimsizes(lat), dimsizes(lon)/), lon, 1)
    X!0 = "lat"
    X!1 = "lon"
    ;X = where(X.ge.180.0, X-360, X) 
    ;copy_VarCoords(longitude,X)
    ;copy_VarCoords(latitude,Y) 
    ;X=decimalPlaces(longitude ,2,True)
    ;Y=decimalPlaces(latitude ,2,True)

     PSL_rg!1 = "y"
    PSL_rg!2 = "x"
    V850_rg!1 = "y"
    V850_rg!2 = "x"
    U850_rg!1 = "y"
    U850_rg!2 = "x"
    V200_rg!1 = "y"
    V200_rg!2 = "x"
    U200_rg!1 = "y"
    U200_rg!2 = "x"
    vIVT_rg!1 = "y"
    vIVT_rg!2 = "x"
    uIVT_rg!1 = "y"
    uIVT_rg!2 = "x"
    T850_rg!1 = "y"
    T850_rg!2 = "x"
    Q850_rg!1 = "y"
    Q850_rg!2 = "x"
    PRECT_rg!1 = "y"
    PRECT_rg!2 = "x"
    FLUT_rg!1 = "y"
    FLUT_rg!2 = "x"
    Z500_rg!1 = "y"
    Z500_rg!2 = "x"

    lat!0 = "y"
    lon!0 = "x"
 
;---creating an output file for the variable data in new-grid format
  strso=(/"./CAM7_L58_HR_hourly.",date_str(0),".nc"/)
  ncdf     = addfile(str_concat(strso),"c")
;----define dimensions  
  dimNames = (/"time", "y", "x"/)
  dimSizes = (/ -1 ,  721,  1440 /)
  dimUnlim = (/ True , False, False/)
  filedimdef(ncdf, dimNames, dimSizes, dimUnlim)
;-----define variable and its attributes
  filevardef(ncdf,"lon", typeof(lon), getvardims(lon))
  filevarattdef(ncdf,"lon",src_lon)
  ncdf->lon= (/lon/)

  ;filevardef(ncdf,"X", typeof(X), getvardims(X))
  ;filevarattdef(ncdf,"X",src_lon)
  ncdf->X=(/X/)  

  filevardef(ncdf,"lat", typeof(lat), getvardims(lat))
  filevarattdef(ncdf,"lat",src_lat)
  ncdf->lat= (/lat/)

  ;filevardef(ncdf,"Y", typeof(Y), getvardims(Y))
  ;filevarattdef(ncdf,"Y",src_lat) 
  ncdf->Y=(/Y/)  

  filevardef(ncdf,"time", typeof(time), getvardims(time))
  filevarattdef(ncdf,"time",time)
  ncdf->time= (/time/) 
  

  filevardef(ncdf,"PSL", typeof(PSL_rg), getvardims(PSL_rg))
  filevarattdef(ncdf,"PSL",PSL_in)
  ncdf->PSL= (/PSL_rg/)
  
  filevardef(ncdf,"Tb", typeof(Tb), getvardims(PSL_rg))
  filevarattdef(ncdf,"Tb",Tb)
  ncdf->Tb= (/Tb/)

  filevardef(ncdf,"uIVT", typeof(uIVT_in), getvardims(PSL_rg))
  filevarattdef(ncdf,"uIVT",uIVT_in)
  ncdf->uIVT= (/uIVT_rg/)

  filevardef(ncdf,"vIVT", typeof(vIVT_in), getvardims(PSL_rg))
  filevarattdef(ncdf,"vIVT",vIVT_in)
  ncdf->vIVT= (/vIVT_rg/)

  filevardef(ncdf,"U200", typeof(U200_in), getvardims(PSL_rg))
  filevarattdef(ncdf,"U200",U200_in)
  ncdf->U200= (/U200_rg/)

  filevardef(ncdf,"V200", typeof(V200_in), getvardims(PSL_rg))
  filevarattdef(ncdf,"V200",V200_in)
  ncdf->V200= (/V200_rg/)

  filevardef(ncdf,"U850", typeof(U850_in), getvardims(PSL_rg))
  filevarattdef(ncdf,"U850",U850_in)
  ncdf->U850= (/U850_rg/)

  filevardef(ncdf,"V850", typeof(V850_in), getvardims(PSL_rg))
  filevarattdef(ncdf,"V850",V850_in)
  ncdf->V850= (/V850_rg/)

  filevardef(ncdf,"T850", typeof(T850_in), getvardims(PSL_rg))
  filevarattdef(ncdf,"T850",T850_in)
  ncdf->T850= (/T850_rg/)

  filevardef(ncdf,"Q850", typeof(Q850_in), getvardims(PSL_rg))
  filevarattdef(ncdf,"Q850",Q850_in)
  ncdf->Q850= (/Q850_rg/)

  filevardef(ncdf,"Z500", typeof(Z500_in), getvardims(PSL_rg))
  filevarattdef(ncdf,"Z500",Z500_in)
  ncdf->Z500= (/Z500_rg/)

  filevardef(ncdf,"PRECT", typeof(PRECT_in), getvardims(PSL_rg))
  filevarattdef(ncdf,"PRECT",PRECT_rg)
  ncdf->PRECT= (/PRECT_rg/)
  
end do
end   ; if subregion


