import pandas as pd
import numpy as np
import utilities as util
import writeup
import matplotlib
matplotlib.use('Agg')
#maplotlib.rc('text', usetex=True)
font = {'family':'sans-serif','sans-serif':['Helvetica'],
#        'weight' : 'normal',
        'size'   : 40}
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('axes', edgecolor = 'k')
matplotlib.rc('font', **font)
#maplotlib.rcParams['text.latex.preamble'] = [
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.linalg import pinv
from math import ceil

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) >= 0)

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

isfloat = np.vectorize(isfloat)

truepar2 = np.loadtxt("../Results/Estim/par_quad_toni_rr3_fixed.csv")
truepar3 = np.loadtxt("../Results/Estim/par_quad_drs_toni_rr3.csv")

nobj = 3
nres = nobj
nres_end = nres + 7
nmom     = nres_end
nmom_end = nmom + 8
nmom_kt  = nmom_end
nmom_kt_end = nmom_kt + 18
nmom_hp  = nmom_kt_end
nmom_hp_end    = nmom_hp + 18
nepfq  = nmom_hp_end
nepfq_end = nepfq + 18
nepfq_wvar = nepfq_end
nepfq_wvar_end    = nepfq_wvar + 21
nepfq_wvarex = nepfq_wvar_end
nepfq_wvarex_end    = nepfq_wvarex + 21

njac1 = nepfq_wvarex_end
njac1_end = njac1 + 8 * 7
njac2 = njac1_end
njac2_end = njac2 + 18 * 7

truepar = np.loadtxt("../Results/Estim/par_quad_toni_rr3.csv")
truepar[4] = truepar[6]
nobj = 3
nres = nobj
nres_end = nres + 7
nmom     = nres_end
nmom_end = nmom + 8
nmom_kt  = nmom_end
nmom_kt_end = nmom_kt + 18
nmom_hp  = nmom_kt_end
nmom_hp_end    = nmom_hp + 18
nepfq  = nmom_hp_end
nepfq_end = nepfq + 18
nepfq_wvar = nepfq_end
nepfq_wvar_end    = nepfq_wvar + 21
nepfq_wvarex = nepfq_wvar_end
nepfq_wvarex_end    = nepfq_wvarex + 21

njac1 = nepfq_wvarex_end
njac1_end = njac1 + 8 * 7
njac2 = njac1_end
njac2_end = njac2 + 18 * 7

# Moments
dictd_mom = {
                    'filename' : 'mom_rr10',
                    'dataname' : 'trial_data_rr_mom.csv',
                    'dataname2': 'trial_data_rr_quad.csv',
                    'wdataname': 'trial_data_rr_wmom.csv',
                    'vdataname': 'trial_data_rr_vmom.csv',
                    'filev2'   : 'trial_data_rr_vquad.csv',
                    'filev12'  : 'trial_data_rr_vmom_quad.csv',
                  }

dictv_mom = {}


# EPF
dictd_epfq = {
                    'filename' : 'quad_rr10',
                    'dataname' : 'trial_data_rr_quad.csv',
                    'dataname2': 'trial_data_rr_mom.csv',
                    'wdataname': 'trial_data_rr_wquad.csv',
                    'vdataname': 'trial_data_rr_vquad.csv',
                    'filev2'   : 'trial_data_rr_vmom.csv',
                    'filev12'  : 'trial_data_rr_vmom_quad.csv',
                  }

dictv_epfq = {
                    'nimom'     : nepfq,
                    'nimom_end' : nepfq_end,
                    'namom'     : nmom,
                    'namom_end' : nmom_end,
                    'njac'      : nepfq_wvarex_end + 8*7,
                    'njacend'   : nepfq_wvarex_end + 18*7 + 8*7,
                    'njaca'     : nepfq_wvarex_end,
                    'njacaend'  : nepfq_wvarex_end + 8*7,
                    'npar_est'  : 4,
                    'nresi'     : nres,
                    'nresi_end' : nres_end,
                    'transposer': False,
    
                  }


# EPF misspecified
dictd_epfq_mis = {
                    'filename' : 'quad_rr10_mis1',
                    'dataname' : 'trial_data_rr_mis_quad.csv',
                    'dataname2': 'trial_data_rr_mis_mom.csv',
                    'wdataname': 'trial_data_rr_mis_wquad.csv',
                    'vdataname': 'trial_data_rr_mis_vquad.csv',
                    'filev2'   : 'trial_data_rr_mis_vmom.csv',
                    'filev12'  : 'trial_data_rr_mis_vmom_quad.csv',
                  }

dictv_epfq_mis = {
                    'nimom'     : nepfq,
                    'nimom_end' : nepfq_end,
                    'namom'     : nmom,
                    'namom_end' : nmom_end,
                    'njac'      : nepfq_wvarex_end + 8*7,
                    'njacend'   : nepfq_wvarex_end + 18*7 + 8*7,
                    'njaca'     : nepfq_wvarex_end,
                    'njacaend'  : nepfq_wvarex_end + 8*7,
                    'npar_est'  : 4,
                    'nresi'     : nres,
                    'nresi_end' : nres_end,
                    'transposer': False,
                    'fixed'    : 0
                  }

# Moments misspecified
dictd_mom_mis = {
                    'filename' : 'mom_rr10_mis1',
                    'dataname' : 'trial_data_rr_mis_mom.csv',
                    'dataname2': 'trial_data_rr_mis_quad.csv',
                    'wdataname': 'trial_data_rr_mis_wmom.csv',
                    'vdataname': 'trial_data_rr_mis_vmom.csv',
                    'filev2'   : 'trial_data_rr_mis_vquad.csv',
                    'filev12'  : 'trial_data_rr_mis_vmom_quad.csv',
                  }

dictv_mom_mis = {
                    'fixed'    : 0
                  }

dict5 = [[dictd_epfq, dictv_epfq]]
dict7 = [[dictd_mom, dictv_mom]]

# Generate out of sample variances

for ddi in dict7:
    print(ddi[0]['filename'])
    util.generate_var_outofsample(ddi[0], ddi[1], 4, upper=100)
    #epf_par_out(truepar, ddi[0], ddi[1])


for ddi in dict5:
    print(ddi[0]['filename'])
    util.generate_var_outofsample(ddi[0], ddi[1], 4, upper=100)
    #epf_par_out(truepar, ddi[0], ddi[1])

print(dictd_epfq_mis['filename'])
util.generate_var_outofsample(dictd_epfq_mis, dictv_epfq_mis, 4, upper=100)
print(dictd_mom_mis['filename'])
util.generate_var_outofsample(dictd_mom_mis, dictv_mom_mis, 4, upper=100)

# Do analysis of Monte Carlo results

res5 = [util.epf_manipulate_alt(truepar, dict5[0][0], dict5[0][1], rhigh=100, rlow=0)]

res7 = [util.epf_manipulate_alt(truepar, dict7[0][0], dict7[0][1], rhigh=100, rlow=0)]

res_epf_mis = []
res_mom_mis = []

for i in range(0,5):
    res_epf_mis += [util.epf_manipulate_alt(truepar, dictd_epfq_mis, dictv_epfq_mis, rhigh=100, rlow=0, getlow = 1000*i, gethigh=1000*(i+1))]
    res_mom_mis += [util.epf_manipulate_alt(truepar, dictd_mom_mis, dictv_mom_mis, rhigh=100, rlow=0, getlow = 1000*i, gethigh= 1000*(i+1))]

# Create misspecification plots

writeup.mis_test_plot(res5 + res_epf_mis, 'quad_rr10')

writeup.mis_test_plot(res7 + res_mom_mis, 'mom_rr10')
