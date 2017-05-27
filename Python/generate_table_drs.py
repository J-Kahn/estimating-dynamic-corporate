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

## Dictionary for epf with fixed cost
dictd_epfq_drs = {
                    'filename' : 'quad_drs_rr10',
                    'dataname' : 'trial_data_rr_drs_quad.csv',
                    'dataname2': 'trial_data_rr_drs_mom.csv',
                    'wdataname': 'trial_data_rr_drs_wquad.csv',
                    'vdataname': 'trial_data_rr_drs_vquad.csv',
                    'filev2'   : 'trial_data_rr_drs_vmom.csv',
                    'filev12'  : 'trial_data_rr_drs_vmom_quad.csv'
                  }

dictv_epfq_drs = {
                    'nimom'     : nepfq_wvarex_end + 18*7 + 8*7,
                    'nimom_end' : nepfq_wvarex_end + 18*7 + 8*7 + 18*1,
                    'namom'     : nmom,
                    'namom_end' : nmom_end,
                    'njac'      : nepfq_wvarex_end + 8*7,
                    'njacend'   : nepfq_wvarex_end + 18*7 + 8*7,
                    'njaca'     : nepfq_wvarex_end,
                    'njacaend'  : nepfq_wvarex_end + 8*7,
                    'npar'      : 7,
                    'npar_est'  : 4,
                    'nresi'     : nres,
                    'nresi_end' : nres_end,
                    'transposer': False
                  }

## Dictionary for epf with fixed cost
dictd_mom_drs = {
                    'filename' : 'mom_drs_rr10',
                    'dataname2' : 'trial_data_rr_drs_quad.csv',
                    'dataname': 'trial_data_rr_drs_mom.csv',
                    'wdataname': 'trial_data_rr_drs_wmom.csv',
                    'filev2': 'trial_data_rr_drs_vquad.csv',
                    'vdataname'   : 'trial_data_rr_drs_vmom.csv',
                    'filev12'  : 'trial_data_rr_drs_vmom_quad.csv'
                  }

dictv_mom_drs = {
                    'namom'     : nepfq,
                    'namom_end' : nepfq_end,
                    'nimom'     : nmom,
                    'nimom_end' : nmom_end,
                    'njaca'      : nepfq_wvarex_end + 8*7,
                    'njacaend'   : nepfq_wvarex_end + 18*7 + 8*7,
                    'njac'     : nepfq_wvarex_end,
                    'njacend'  : nepfq_wvarex_end + 8*7,
                    'npar'      : 7,
                    'npar_est'  : 4,
                    'nresi'     : nres,
                    'nresi_end' : nres_end,
}

dict4 = [[dictd_mom_drs, dictv_mom_drs], [dictd_epfq_drs, dictv_epfq_drs]]

# Generate out of sample variances

for ddi in dict4:
    print(ddi[0]['filename'])
    util.generate_var_outofsample(ddi[0], ddi[1], 4, upper=100)
    #epf_par_out(truepar, ddi[0], ddi[1])

# Do analysis of Monte Carlo results

res4 = []
for ddi in dict4:
    print(ddi[0]['filename'])
    res4 += [util.epf_manipulate_alt(truepar3, ddi[0], ddi[1], rhigh=100, rlow=0, n=75000)]
#res3 += [util.epf_manipulate_alt(truepar2, dict3[1][0], dict3[1][1], rhigh=100, rlow=0, n=75000)]


top1 = "\\begin{table}[hhh]\n  \\caption{Monte Carlo comparison of simulation estimators: Decreasing returns to scale model}\n  \\begin{center}\n  \\begin{tabular*}{1\\textwidth}{r@{\\extracolsep{\\fill}}rrrrr}\n    \\hline\n            &     \\multicolumn{1}{c}{Moments-based} & \\multicolumn{1}{c}{EPF-based} \\\\\n\\hline\n"

bottom1 = " \\hline \n \\end{tabular*} \n \\parbox[c]{7in}{\\footnotesize Indicated expectations and probabilities are estimates \n based on 1,000 Monte Carlo samples of size 75,000.  The samples are \n generated from the model in Section \\ref{sec:model}, with two modifications. We allow for decreasing returns to scale, with \n the profit function curvature given by $\\alpha$. We also omit investment adjustment costs. We consider two different benchmarks: traditiona moments and empirical \n policy functions. \n Both estimators use a clustered weight matrix.  For each parameter we report three statistics. Bias is expressed \n as a percentage of the true coefficient value.  RMSE indicates root mean \n squared error and is also expressed as a percentage of the true \n coefficient.  $\\Pr(t)$ is the fraction of the time we observe a nominal \n 5\\% rejection of the null hypothesis that a parameter equals its true \n value using a {\\em t}-test. \n We report the fraction of trials that \n produce a rejection of three additional tests.  The first is a nominal 5\\% test of \n the model overidentifying restrictions.  The second, out-of-sample test \n has two varieties.  For the EPF-based estimator, it is a chi-squared \n test of the null hypothesis that the moments equal their true values. \n For the moments-based estimator, it is a chi-squared test of the null \n hypothesis that the policy function slopes equal their true values.  For \n the {\\em t}-tests of individual moment conditions, we report the \n highest, median, and lowest rejection rates.  The third test is a test \n of the null hypothesis that features of the data that are not used in \n the estimation equal their true values. \n The moments-based estimator minimizes the \n distance between simulated and data moments.  The EPF-based estimator \n minimizes the distance between policy functions estimated from simulated \n data and those estimated from real data. \n }\n \n \\end{center} \n\\end{table}"

dest = ['$\\delta$ (depreciation rate)', '$\\lambda$ (equity issuance cost)', '$\\xi$ (collateral parameter)', '$\\psi$ (adjustment cost)']

mid1 = ''

sep3 = [' & ', ' & ']

print('Robustness DRS')


for ids, ds in enumerate(dest):
  mid1 += writeup.print_var_alt(ds, ids, res4, sep3, factor = 100)

mid1 += '\\\\[5pt]'

mid2 = writeup.mis_test_alt(75000, res4, sep3)

tab4_pct = top1 +  mid1 + mid2 + bottom1

f = open('../WR/drs_tables.tex', 'w')
f.write(tab4_pct)
f.close()
