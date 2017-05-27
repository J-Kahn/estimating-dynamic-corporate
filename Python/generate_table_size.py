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


# Diagonal moments
dictd_mom_diag = {
                    'filename' : 'mom_diag_rr10',
                    'filename2': 'mom_rr10',
                    'dataname' : 'trial_data_rr_mom.csv',
                    'dataname2': 'trial_data_rr_quad.csv',
                    'wdataname': 'identity',
                    'vdataname': 'trial_data_rr_vmom.csv',
                    'filev2'   : 'trial_data_rr_vquad.csv',
                    'filev12'  : 'trial_data_rr_vmom_quad.csv',

                  }

dictv_mom_diag = {}

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

# Diagonal EPF
dictd_epfq_diag = {
                    'filename' : 'quad_diag_rr10',
                    'filename2': 'quad_rr10',
                    'dataname' : 'trial_data_rr_quad.csv',
                    'dataname2': 'trial_data_rr_mom.csv',
                    'wdataname': 'identity',
                    'vdataname': 'trial_data_rr_vquad.csv',
                    'filev2'   : 'trial_data_rr_vmom.csv',
                    'filev12'  : 'trial_data_rr_vmom_quad.csv',

                  }

dictv_epfq_diag = {
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

## Small sample
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

dict5 = [[dictd_epfq, dictv_epfq], [dictd_epfq_diag, dictv_epfq_diag]]
dict6 = [[dictd_epfq_small, dictv_epfq_small], [dictd_epfq_diag_small, dictv_epfq_diag_small]]
dict7 = [[dictd_mom, dictv_mom], [dictd_mom_diag, dictv_mom_diag]]
dict8 = [[dictd_mom_small, dictv_mom_small], [dictd_mom_diag_small, dictv_mom_diag_small]]


# Generate out of sample variance

for ddi in dict8:
    print(ddi[0]['filename'])
    util.generate_var_outofsample(ddi[0], ddi[1], 4, upper=100)
    #epf_par_out(truepar, ddi[0], ddi[1])


for ddi in dict7:
    print(ddi[0]['filename'])
    util.generate_var_outofsample(ddi[0], ddi[1], 4, upper=100)
    #epf_par_out(truepar, ddi[0], ddi[1])

for ddi in dict5:
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

res5 = []
res5 = [util.epf_ineff_manipulate_alt(truepar, dict5[1][0], dict5[1][1], rhigh=100, rlow=0),
        util.epf_manipulate_alt(truepar, dict5[0][0], dict5[0][1], rhigh=100, rlow=0)]

res7 = []
res7 = [util.epf_ineff_manipulate_alt(truepar, dict7[1][0], dict7[1][1], rhigh=100, rlow=0),
        util.epf_manipulate_alt(truepar, dict7[0][0], dict7[0][1], rhigh=100, rlow=0)]




mid2 = writeup.mis_test_alt(75000, res7 + res5, sep1)

mid21 = writeup.mis_test_alt(1000, res8 + res6, sep1)

top1 = ''
top1 += '\\begin{table} \n'
top1 += '   \\caption{Monte Carlo comparison of specification tests}\\label{tab:Jstats} \n'
top1 += ' \\begin{center} \n'
top1 += '  {\\small \n'
top1 += ' \\begin{tabular*}{1\\textwidth}{r@{\\extracolsep{\\fill}}ccccc} \n'
top1 += '%   \\begin{tabular}{lccccc} \n'
top1 += '             &  \\multicolumn{2}{c}{Moments-based} & \n'
top1 += ' &\\multicolumn{2}{c}{EPF-based} \\\\ \n'
top1 += '             &  Identity  & Clustered &$\\qquad$ & Identity & Clustered \\\\ \n'
top1 += '     \\hline \n'
top1 += '%\\multicolumn{1}{l}{Panel A: Minimized objective function} \\\\  \\\\ \n'
top1 += ' Sample size = 75,000 \\\\ \n'

bottom1 = ''
bottom1 += ' \\hline \n '
bottom1 += '   \\end{tabular*}} \n '
bottom1 += '   \\medskip \n '
bottom1 += ' \n '
bottom1 += '\\parbox{6.5in}{\\footnotesize Indicated expectations and probabilities \n '
bottom1 += 'are estimates based on 1,000 Monte Carlo samples of sizes 75,000 and \n '
bottom1 += '1,000.  The samples are generated from the model in Section \n '
bottom1 += '\\ref{sec:model}.  The moments-based estimator minimizes the distance \n '
bottom1 += 'between simulated and data moments.  The EPF-based estimator minimizes \n '
bottom1 += 'the distance between policy functions estimated from simulated data and \n '
bottom1 += 'those estimated from real data. We report the fraction of trials that \n '
bottom1 += 'produce a nominal 5\\% rejection of three additional tests.  The first is \n '
bottom1 += 'the test of the model overidentifying restrictions.  The second our \n '
bottom1 += 'out-of-sample test, which has two varieties.  For the EPF-based \n '
bottom1 += 'estimator, it is a chi-squared test of the null hypothesis that the \n '
bottom1 += 'moments equal their true values.  For the moments-based estimator, it is \n '
bottom1 += 'a chi-squared test of the null hypothesis that the policy function \n '
bottom1 += 'slopes equal their true values.  The third is a {\\em t}-test on \n '
bottom1 += 'individual moment conditions.  For these tests, we report the highest, \n '
bottom1 += 'median, and lowest rejection rates.} \n '
bottom1 += ' \n '
bottom1 += ' \\end{center} \n '
bottom1 += ' \\end{table} \n '

tab03 = top1 + mid2 + '\n Sample size = 1,000 \\\\ \n' + mid21 + bottom1


f = open('../WR/size_table.tex', 'w')
f.write(tab03)
f.close()
