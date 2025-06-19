def Ctl_c6477():
    exp_file = 'Utility.yaml'
    #xl= ['control','NOmovmtn', 'NOfronts', 'movmtn_t01','movmtn_t02','movmtn_t03','movmtn_t04','fronts_t01', 'CMIP6_WACCM' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    xl= ['control','control_b' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= ['control','movmtn_t01','movmtn_t02','movmtn_t03','movmtn_t04', 'CMIP6_WACCM' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def Utility():
    exp_file = 'Utility.yaml'
    #xl= ['control','NOmovmtn', 'NOfronts', 'movmtn_t01','movmtn_t02','movmtn_t03','movmtn_t04','fronts_t01', 'CMIP6_WACCM' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    xl= ['control','NOfronts','fronts_t01', 'CMIP6_WACCM' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    #xl= ['control','movmtn_t01','movmtn_t02','movmtn_t03','movmtn_t04', 'CMIP6_WACCM' ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def Fig1Paper():
    exp_file = 'fig1_for_paper.yaml'
    xl= [ 'NOmovmtn', 'ne120_ctl', 'ne240_ctl' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def KeyFig_Paper():
    exp_file = 'Utility.yaml'
    xl= [ 'control','NOmovmtn', 'CMIP6_WACCM' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def JointGWTune():
    exp_file = 'Utility.yaml'
    #xl= ['control', 'movmtn_t02', 'fronts_t01',    'cfrntmm_t01' , 'cfrntmm_t02' ,  'cfrntmm_t03' , 'cfrntmm_t04' ,  'cfrntmm_t05' , 'cfrormm_t01' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    xl= ['control', 'NOmovmtn', 'NOfronts' ,'cfrntmm_t03' , 'cfrntmm_t04re' , 'cfrntmm_t04' , 'cfrntmm_t06' , ] #'fmt_tau_amplify_150'] #,'tau_reduce_50','tau_amplify_200'] #'ne240x2_QxQsst']
    return exp_file,xl

def RadTest():
    exp_file = 'Utility.yaml'
    xl= ['control', 'radtest_t04' , 'radtest_t04B' ,'radtest_t06' , ] 
    xl= ['control', 'radtest_t06','radtest_Q05_mpas120'  ,'radtest_Q07_mpas120'  , ] 
    
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
    xl=['Coriolis_control','Non_trad_Coriolis',]
    #xl=['ne30_ctl',]  
    return exp_file,xl


