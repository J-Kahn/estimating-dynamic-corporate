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

## Dictionary for epf transformed
dictd_epfq_trans = {
                    'filename' : 'quad_trans_rr10',
                    'dataname' : 'trial_data_rr_quad_trans.csv',
                    'dataname2': 'trial_data_rr_mom.csv',
                    'wdataname': 'trial_data_rr_wquad_trans.csv',
                    'vdataname': 'trial_data_rr_vquad_trans.csv',
                    'filev2'   : 'trial_data_rr_vmom.csv',
                    'filev12'  : 'trial_data_rr_vmom_quad_trans.csv',
                  }



dictv_epfq_trans = {
                    'nimom'     : nepfq_wvarex_end + 18*7 + 8*7,
                    'nimom_end' : nepfq_wvarex_end + 18*7 + 8*7 +18,
                    'namom'     : nmom,
                    'namom_end' : nmom_end,
                    'njac'      : nepfq_wvarex_end + 8*7,
                    'njacend'   : nepfq_wvarex_end + 18*7 + 8*7,
                    'njaca'     : nepfq_wvarex_end,
                    'njacaend'  : nepfq_wvarex_end + 8*7,
                    'transposer': False,
                    'fixed'   : 0
                  }

## Dictionary for epf cube
dictd_epfc = {
                    'filename' : 'cube_rr10',
                    'dataname' : 'trial_data_rr_cube.csv',
                    'dataname2': 'trial_data_rr_mom.csv',
                    'wdataname': 'trial_data_rr_wcube.csv',
                    'vdataname': 'trial_data_rr_vcube.csv',
                    'filev2'   : 'trial_data_rr_vmom.csv',
                    'filev12'  : 'trial_data_rr_vmom_cube.csv',
                  }

dictv_epfc = {
                    'nimom'     : nepfq_wvarex_end + 30*7 + 8*7,
                    'nimom_end' : nepfq_wvarex_end + 30*7 + 8*7 + 30,
                    'namom'     : nmom,
                    'namom_end' : nmom_end,
                    'njac'      : nepfq_wvarex_end + 8*7,
                    'njacend'   : nepfq_wvarex_end + 30*7 + 8*7,
                    'njaca'     : nepfq_wvarex_end,
                    'njacaend'  : nepfq_wvarex_end + 8*7,
                    'transposer': False,
                    'fixed'   : 0
                  }

## Dictionary for moments with covariances
dictd_mom_cv = {
                    'filename' : 'mom_cv_rr10',
                    'dataname' : 'trial_data_rr_mom_cv.csv',
                    'dataname2': 'trial_data_rr_quad.csv',
                    'wdataname': 'trial_data_rr_wmom_cv.csv',
                    'vdataname': 'trial_data_rr_vmom_cv.csv',
                    'filev2'   : 'trial_data_rr_vquad.csv',
                    'filev12'  : 'trial_data_rr_vmom_cv_quad.csv',
                  }

dictv_mom_cv = {
                    'nimom'     : nepfq_wvarex_end + 18*7 + 11*7,
                    'nimom_end' : nepfq_wvarex_end + 18*7 + 11*7 +11,
                    'namom'     : nepfq,
                    'namom_end' : nepfq_end,
                    'njaca'     : nepfq_wvarex_end + 11*7,
                    'njacaend'  : nepfq_wvarex_end + 11*7 + 18*7,
                    'njac'      : nepfq_wvarex_end,
                    'njacend'   : nepfq_wvarex_end + 11*7,
                    'fixed'    : 0
                  }


dict1 = [[dictd_mom_cv, dictv_mom_cv], [dictd_epfc, dictv_epfc], [dictd_epfq_trans, dictv_epfq_trans]]

# Generate out of sample variances

for ddi in dict1:
    print(ddi[0]['filename'])
    util.generate_var_outofsample(ddi[0], ddi[1], 4, upper=100)
    #epf_par_out(truepar, ddi[0], ddi[1])

# Do analysis of Monte Carlo results

res = []
for ddi in dict1:
    #for ddi in dict1[1:3]:
    print(ddi[0]['filename'])
    res  += [util.epf_manipulate_alt(truepar, ddi[0], ddi[1], rhigh=100, rlow=0)]

top1 = "\\begin{table}[hhh] \n  \\caption{Monte Carlo comparison of simulation estimators: EPFs, lage sample size} \n  \\begin{center} \n  \\hspace*{-0.5in} \n    \\begin{tabular}{rrrrrrrrr}\n    \\hline\n            &  \\multicolumn{1}{c}{Extra} \\\\\n            &  \\multicolumn{1}{c}{Covariance} & &\\multicolumn{1}{c}{High-order} & &\\multicolumn{1}{c}{Transformed} \\\\\n            &  \\multicolumn{1}{c}{Moments} & &\\multicolumn{1}{c}{EPF-based} & &\\multicolumn{1}{c}{EPF-based} \\\\\n\\hline\n"

bottom1 = " \\hline \n \\end{tabular} \n \\parbox[c]{7in}{\\footnotesize Indicated expectations and probabilities are estimates based on 1,000 Monte Carlo samples of size 75,000.  The samples are generated from the model in Section \\ref{sec:model}.  For each parameter we report three statistics. Bias is expressed as a percentage of the true coefficient value.  RMSE indicates root mean squared error and is also expressed as a percentage of the true coefficient.  $\\Pr(t)$ is the fraction of the time we observe a nominal 5\\% rejection of the null hypothesis that a parameter equals its true value using a {\\em t}-test. We report the fraction of trials that produce a rejection of three additional tests.  The first is a nominal 5\% test of the model overidentifying restrictions.  The second, out-of-sample test has two varieties.  For the EPF-based estimator, it is a chi-squared test of the null hypothesis that the moments equal their true values. For the moments-based estimator, it is a chi-squared test of the null hypothesis that the policy function slopes equal their true values.  For the {\\em t}-tests of individual moment conditions, we report the highest, median, and lowest rejection rates.  The third test is a test of the null hypothesis that features of the data that are not used in the estimation equal their true values. The moments-based estimator minimizes the distance between simulated and data moments.  The EPF based estimator minimizes the distance between policy functions estimated from simulated data and those estimated from real data.}\n \\end{center} \n \\end{table} "

des = ['$\\delta$ (depreciation rate)', '$\\lambda$ (equity issuance cost)', '$\\xi$ (collateral parameter)', '$\\gamma$ (investment adjustment cost)']

sep1 = [' & ',' && ', ' && ']

mid1 = ''

print('Robustness')
for ids, ds in enumerate(des):
  mid1 += writeup.print_var_alt(ds, ids, res, sep1, factor = 100,  lenm = writeup.find_max_len(des))

mid1 += '\\\\[5pt]'

mid2 = writeup.mis_test_alt(75000, res, sep1, lenm = writeup.find_max_len(des))

tab1_pct = top1 +  mid1  + mid2 + bottom1

f = open('../WR/robust_table.tex', 'w')
f.write(tab1_pct)
f.close()
