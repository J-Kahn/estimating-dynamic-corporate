# coding: utf-8

# In[570]:

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


# In[571]:

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) >= 0)


# In[572]:

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
dictd_mom_small = {
                    'filename' : 'mom_rr10_small',
                    'dataname' : 'trial_data_rr_small_mom.csv',
                    'dataname2': 'trial_data_rr_small_quad.csv',
                    'wdataname': 'trial_data_rr_small_wmom.csv',
                    'vdataname': 'trial_data_rr_small_vmom.csv',
                    'filev2'   : 'trial_data_rr_small_vquad.csv',
                    'filev12'  : 'trial_data_rr_small_vmom_quad.csv',
                  }

dictv_mom_small = {}

# Diagonal moments
dictd_mom_diag_small = {
                    'filename' : 'mom_rr10_small_diag',
                    'filename2': 'mom_rr10_small',
                    'dataname' : 'trial_data_rr_small_mom.csv',
                    'dataname2': 'trial_data_rr_small_quad.csv',
                    'wdataname': 'identity',
                    'vdataname': 'trial_data_rr_small_vmom.csv',
                    'filev2'   : 'trial_data_rr_small_vquad.csv',
                    'filev12'  : 'trial_data_rr_small_vmom_quad.csv',

                  }

dictv_mom_diag_small = {}

# EPF
dictd_epfq_small = {
                    'filename' : 'quad_rr10_small',
                    'dataname' : 'trial_data_rr_small_quad.csv',
                    'dataname2': 'trial_data_rr_small_mom.csv',
                    'wdataname': 'trial_data_rr_small_wquad.csv',
                    'vdataname': 'trial_data_rr_small_vquad.csv',
                    'filev2'   : 'trial_data_rr_small_vmom.csv',
                    'filev12'  : 'trial_data_rr7_small_vmom_quad.csv',
                  }

dictv_epfq_small = {
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

# Diagonal EPF
dictd_epfq_diag_small = {
                    'filename' : 'quad_rr10_small_diag',
                    'filename2': 'quad_rr10_small',
                    'dataname' : 'trial_data_rr_small_quad.csv',
                    'dataname2': 'trial_data_rr_small_mom.csv',
                    'wdataname': 'identity',
                    'vdataname': 'trial_data_rr_small_vquad.csv',
                    'filev2'   : 'trial_data_rr_small_vmom.csv',
                    'filev12'  : 'trial_data_rr7_small_vmom_quad.csv',
                  }

dictv_epfq_diag_small = {
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

dict6 = [[dictd_epfq_small, dictv_epfq_small], [dictd_epfq_diag_small, dictv_epfq_diag_small]]
dict8 = [[dictd_mom_small, dictv_mom_small], [dictd_mom_diag_small, dictv_mom_diag_small]]


# Generate out of sample variances

for ddi in dict8:
    print(ddi[0]['filename'])
    util.generate_var_outofsample(ddi[0], ddi[1], 4, upper=100)
    #epf_par_out(truepar, ddi[0], ddi[1])

for ddi in dict6:
    print(ddi[0]['filename'])
    util.generate_var_outofsample(ddi[0], ddi[1], 4, upper=100, nt = 1000)
    #epf_par_out(truepar, ddi[0], ddi[1])

# Do analysis of Monte Carlo results

res8 = []
res8 = [util.epf_ineff_manipulate_alt(truepar, dict8[1][0], dict8[1][1], rhigh=100, rlow=0, n = 1000),
        util.epf_manipulate_alt(truepar, dict8[0][0], dict8[0][1], rhigh=100, rlow=0, n=1000)]


res6 = []
res6 = [util.epf_ineff_manipulate_alt(truepar, dict6[1][0], dict6[1][1], rhigh=100, rlow=0, n = 1000),
        util.epf_manipulate_alt(truepar, dict6[0][0], dict6[0][1], rhigh=100, rlow=0, n = 1000)]


top1 = ''
top1 += '\\begin{table}[hhh]'
top1 += '  \\caption{Monte Carlo comparison of simulation estimators: Small sample size\\label{tab:montecarlo_small_sample}}'
top1 += '  \\begin{center}'
top1 += '  \\begin{tabular*}{1\\textwidth}{r@{\\extracolsep{\\fill}}rrrrr}'
top1 += '    \\hline'
top1 += '            &  \\multicolumn{2}{c}{Moments-based} & &\\multicolumn{2}{c}{EPF-based} \\\\'
top1 += '\\multicolumn{1}{l}{Parameter}            &  Identity  & Clustered &$\\qquad$ & Identity & Clustered \\\\'

bottom1 = ''
bottom1 += ' \\hline \n'
bottom1 += '\\end{tabular*} \n'
bottom1 += ' \n'
bottom1 += '\\parbox[c]{6.5in}{\\footnotesize Indicated expectations and probabilities \n'
bottom1 += 'are estimates based on 1,000 Monte Carlo samples of size 1,000.  The \n'
bottom1 += 'samples are generated from the model in Section \\ref{sec:model}.  We \n'
bottom1 += 'consider two different benchmarks:  traditional moments and empirical \n'
bottom1 += 'policy functions.  Both estimators use a clustered weight matrix.  For \n'
bottom1 += 'each parameter, we report three statistics.  Bias is expressed as a \n'
bottom1 += 'fraction of the true coefficient value.  RMSE indicates root mean \n'
bottom1 += 'squared error and is also expressed as a fraction of the true \n'
bottom1 += 'coefficient.  $\\Pr(t)$ is the fraction of the time we observe a nominal \n'
bottom1 += '5\\% rejection of the null hypothesis that a parameter equals its true \n'
bottom1 += 'value using a {\\em t}-test.  } \n'
bottom1 += ' \n'
bottom1 += '\\end{center} \n'
bottom1 += '\\end{table} \n'

des = ['$\\delta$ (depreciation rate)', '$\\lambda$ (equity issuance cost)', '$\\xi$ (collateral parameter)', '$\\gamma$ (investment adjustment cost)']

print('Small sample size')
sep1 = [' & ',' & ', ' & ', ' & ']
mid1 = ""
for ids, ds in enumerate(des):
  mid1 += writeup.print_var_alt(ds, ids, res8 + res6, sep1, factor = 100, lenm = writeup.find_max_len(des))

mid1 += '\\\\[5pt]'

tab01_pct = top1 +  mid1 + bottom1



f = open('../WR/small_table.tex', 'w')
f.write(tab01_pct)
f.close()
