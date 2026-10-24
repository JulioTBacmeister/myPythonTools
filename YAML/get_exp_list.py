def Ctl_c6477():
    exp_file = 'Utility.yaml'
    #xl= ['control','NOmovmtn', 'NOfronts', 'movmtn_t01','movmtn_t02','movmtn_t03','movmtn_t04','fronts_t01', 'CMIP6_WACCM' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    xl= ['control','control_b' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= ['control','movmtn_t01','movmtn_t02','movmtn_t03','movmtn_t04', 'CMIP6_WACCM' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def ControlRuns():
    exp_file = 'keyfig_for_paper.yaml'
    xl= ['control','ctl_c64101' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    xl= ['control','control_2000','NOmovmtn','NOmovmtn_2000','CMIP6_WACCM' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def fGW_tests():
    exp_file = 'fGW_tests.yaml'
    #xl= ['control_c6477' ,'control_c64116' ,'fGWfix01_c64116' ,'fGWhack01_c64116' ,'fGWfix01_NOmmGW_c64116' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= ['NOmovmtn_c6477' ,'control_c64116' ,'fGWfix01_c64116' ,'fGWhack01_c64116' ,'fGWfix01_NOmmGW_c64116' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= ['control_c64116' , 'fGWfix01_c64116' ,'fGWhack01_c64116' ,]
    #xl= ['control_c64116' , 'fGWfix01_c64116' ,'fGWhack01_c64116' ,'NOmovmtn_c6477' , 'fGWfix01_NOmmGW_c64116' , ]
    #xl= ['control_c64116' , 'fGWfix01_c64116' ,'fGWfix04_c64116' , ]
    #xl= ['control_c64116' ,'fGWfix04_c64116' ,'fGWfix04_divdamp_c64116' , ]
    #xl= ['control_c64116' ,'fGWfix04_divdamp_c64116' ,'fGWfix04_divdamp_momcon_c64116' , 'control_c64126' ,]
    #xl= ['control_c64126' ,'NoFGW_c64126' ,'Sponge_c64126' ,]
    xl= ['CTL_c64144', 'NoMMGW_c64144', 'NoFGW_c64144' ,]
    return exp_file,xl

def tausrfRuns():
    #exp_file = 'tausrfRuns.yaml'
    #xl= ['control_c6477','tausrf200_c6477' ] 
    exp_file = 'OfficialRuns.yaml'
    xl= ['B1850_x234', 'B1850_x244', ] 
    return exp_file,xl

def DevelControlRuns():
    exp_file = 'ControlRuns.yaml'
    xl= ['control_c6477','control_c64116' ] 
    return exp_file,xl

def ModHeat():
    exp_file = 'Utility.yaml'
    xl= ['ctl_c64101', 'ModHeat_F03', ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def Utility():
    exp_file = 'Utility.yaml'
    #xl= ['control','NOmovmtn', 'NOfronts', 'movmtn_t01','movmtn_t02','movmtn_t03','movmtn_t04','fronts_t01', 'CMIP6_WACCM' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    xl= ['control','NOfronts','fronts_t01', 'CMIP6_WACCM' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= ['control','movmtn_t01','movmtn_t02','movmtn_t03','movmtn_t04', 'CMIP6_WACCM' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def Fig1Paper():
    exp_file = 'fig1_for_paper.yaml'
    #xl= [ 'NOmovmtn', 'NOmovmtn_2000', 'ne120_ctl', 'ne240_ctl' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= [ 'NOmovmtn', 'ne120_ctl', 'ne240_ctl' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    xl= [ 'NOmovmtn_2004_06_01', 'ne120_ctl', 'ne240_ctl' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= [ 'NOmovmtn_2000', 'NOmovmtn_2004_06_01','ne120_ctl', 'ne240_ctl' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= [ 'NOmovmtn_2000','ne120_ctl', 'ne240_ctl' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= [ 'NOmovmtn', 'ne120_ctl', 'ne240_2007' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def Fig_for_NSC2026():
    exp_file = 'fig1_for_paper.yaml'
    #xl= [ 'NOmovmtn', 'NOmovmtn_2000', 'ne120_ctl', 'ne240_ctl' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= [ 'NOmovmtn', 'ne120_ctl', 'ne240_ctl' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    xl= [ 'NOmovmtn', 'ne120_ctl', 'ne240_ctl' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= [ 'NOmovmtn_2000', 'NOmovmtn_2004_06_01','ne120_ctl', 'ne240_ctl' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= [ 'NOmovmtn_2000','ne120_ctl', 'ne240_ctl' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= [ 'NOmovmtn', 'ne120_ctl', 'ne240_2007' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def KeyFig_Paper():
    exp_file = 'keyfig_for_paper.yaml'
    xl= [ 'control','NOmovmtn', 'CMIP6_WACCM' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def JointGWTune():
    exp_file = 'Utility.yaml'
    #xl= ['control', 'movmtn_t02', 'fronts_t01',    'cfrntmm_t01' , 'cfrntmm_t02' ,  'cfrntmm_t03' , 'cfrntmm_t04' ,  'cfrntmm_t05' , 'cfrormm_t01' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    xl= ['control', 'NOmovmtn', 'NOfronts' ,'cfrntmm_t03' , 'cfrntmm_t04re' , 'cfrntmm_t04' , 'cfrntmm_t06' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def RadTest():
    exp_file = 'radtests.yaml'
    #xl= ['control', 'radtest_t04' , 'radtest_t04B' ,'radtest_t06' , ] 
    #xl= ['control', 'NOmovmtn','radtest_Q08_mpas120'  ,'radtest_W01_mpas120'  , ] 
    #xl= ['NOmovmtn','radtest_F02_mpas120'  ,'radtest_F04_mpas120'  , ] # 'radtest_F01_ne30pg3'  , ] 
    #xl= ['radtest_flt_J01_ne30pg3','radtest_flt_CTL_ne30pg3', ] # 'radtest_F01_ne30pg3'  , ] 
    #xl= ['NOmovmtn','radtest_J02_ne30pg3','radtest_J03_ne30pg3', ] 
    #xl= ['control','radtest_J03_ne30pg3', 'radtest_x01_ne30pg3',] 
    #xl= ['control','radtest_x01_ne30pg3','radtest_x02_ne30pg3','radtest_x03_ne30pg3','radtest_x04_ne30pg3',] 
    #xl= ['control','radtest_x03_ne30pg3','radtest_x04_ne30pg3','radtest_x04b_ne30pg3','radtest_x05_ne30pg3',] 
    #xl= ['control','radtest_S01_ne30pg3b','radtest_S02_ne30pg3b','radtest_T03_ne30pg3b',] 
    #xl= ['control','radtest_S01_ne30pg3b','radtest_T01_ne30pg3b','radtest_T02_ne30pg3b','radtest_T03_ne30pg3b',] 
    #xl= ['radtest_T03_ne30pg3b','radtest_U01_ne30pg3b',] 
    #xl= ['radtest_S01_ne30pg3b','radtest_U02_ne30pg3b',] 
    #xl= ['control','radtest_U02_ne30pg3b','radtest_U01_ne30pg3b','radtest_U03_ne30pg3b','radtest_U04_ne30pg3b','radtest_U05_ne30pg3b',] 
    xl= ['radtest_U02_ne30pg3b','radtest_S02_ne30pg3b','radtest_U01_ne30pg3b','radtest_U03_ne30pg3b','radtest_U04_ne30pg3b','radtest_U05_ne30pg3b',] 
    xl= ['control','radtest_U03_ne30pg3c',] 
    xl= ['radtest_U03_ne30pg3c','radtest_U03_ne30pg3d',] 
    #xl= ['radtest_flt_J01_mpas120','radtest_flt_CTL_mpas120', ] # 'radtest_F01_ne30pg3'  , ] 
    #xl= ['radtest_F01_ne30pg3','radtest_J01_ne30pg3', ] # 'radtest_F01_ne30pg3'  , ] 
    #xl= ['radtest_F01_ne30pg3','radtest_AICN03_mpas120'  , 'radtest_AICN04_mpas120'  , ] # 'radtest_F01_ne30pg3'  , ] 
    #xl= ['radtest_F02_mpas120'  ,'radtest_F04_mpas120'  , ] # 'radtest_F01_ne30pg3'  , ] 
    #xl= ['control', 'NOmovmtn','radtest_W01_mpas120'  , ] 
    #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= ['radtest_t01' , 'radtest_t02' ,  'radtest_t03' , 'radtest_t04' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def MPAS_Tests():
    exp_file = 'Utility.yaml'
    #xl= ['control', 'radtest_t04' , 'radtest_t04B' ,'radtest_t06' , ] 
    #xl= ['control', 'NOmovmtn','radtest_Q08_mpas120'  ,'radtest_W01_mpas120'  , ] 
    xl= ['NOmovmtn','radtest_W07_mpas120'  , 'radtest_F01_mpas120'  , ] 
    xl= ['NOmovmtn', 'radtest_mpasa120_c64109' , ] 
    #xl= ['control', 'NOmovmtn','radtest_W01_mpas120'  , ] 
    #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= ['radtest_t01' , 'radtest_t02' ,  'radtest_t03' , 'radtest_t04' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def MovmtnTune():
    exp_file = 'Utility.yaml'
    xl= ['control','NOmovmtn', 'movmtn_t01','movmtn_t02','movmtn_t03','movmtn_t04','movmtn_t05', ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def FrontalTune():
    exp_file = 'Utility.yaml'
    xl= ['control','NOfronts','fronts_t01', 'CMIP6_WACCM' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def Randall():
    exp_file='ForDaveRandal.yaml'
    xl= ['NOgwATall','noMM','ne240_ctl','tun00',] #'ne240x2_QxQsst']
    return exp_file,xl

def Xtreme():
    exp_file='GW_Xtreme_ppe.yaml'
    xl= ['control','nofgw','NOogw','NOgwATall',] #'ne240x2_QxQsst']
    #xl= ['control','NOgwATall',] #'ne240x2_QxQsst']
    return exp_file,xl

def Resolution():
    exp_file='Res_ensemble.yaml'
    xl=['ne30_ctl','ne120_ctl_0','ne240_ctl',]
    #xl=['ne30_ctl',]  
    return exp_file,xl

def Official():
    exp_file='OfficialRuns.yaml'
    #xl=['FMT_control','FMT_nu_div',]
    xl=['B1850_x194','B1850_x195',]
    xl=['B1850_x234','B1850_x245','B1850_x249','B1850_x254',]
    #xl=['ne30_ctl',]  
    return exp_file,xl


