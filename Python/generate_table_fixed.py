import pandas as pd
import numpy as np
import utilities as util
import writeup
import matplotlib
matplotlib.use('Agg')
matplotlib.rc('text', usetex=True)
font = {'family':'sans-serif','sans-serif':['Helvetica'],
#        'weight' : 'normal',
        'size'   : 40}
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('axes', edgecolor = 'k')
matplotlib.rc('font', **font)
matplotlib.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]
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
dictd_epfq_fixed = {
                    'filename' : 'quad_fixed_rr10',
                    'dataname' : 'trial_data_rr_fixed_quad.csv',
                    'dataname2': 'trial_data_rr_fixed_mom.csv',
                    'wdataname': 'trial_data_rr_fixed_wquad.csv',
                    'vdataname': 'trial_data_rr_fixed_vquad.csv',
                    'filev2'   : 'trial_data_rr_fixed_vmom.csv',
                    'filev12'  : 'trial_data_rr_fixed_vmom_quad.csv',
                  }

dictv_epfq_fixed = {
                    'nimom'     : nepfq_wvarex_end + 18*8 + 8*8+1,
                    'nimom_end' : nepfq_wvarex_end + 18*8 + 8*8 + 18*1 + 1,
                    'namom'     : nmom+1,
                    'namom_end' : nmom_end+1,
                    'njac'      : nepfq_wvarex_end + 8*8+1,
                    'njacend'   : nepfq_wvarex_end + 18*8 + 8*8+1,
                    'njaca'     : nepfq_wvarex_end+1,
                    'njacaend'  : nepfq_wvarex_end + 8*8+1,
                    'npar'      : 8,
                    'npar_est'  : 5,
                    'nresi'     : nres,
                    'nresi_end' : nres_end + 1,
                    'transposer': False,
                    'fixed'    : 0
                  }

## Dictionary for cubic epf with fixed cost
dictd_epfc_fixed = {
                    'filename' : 'cube_fixed_rr10',
                    'dataname' : 'trial_data_rr_fixed_cube.csv',
                    'dataname2': 'trial_data_rr_fixed_mom.csv',
                    'wdataname': 'trial_data_rr_fixed_wcube.csv',
                    'vdataname': 'trial_data_rr_fixed_vcube.csv',
                    'filev2'   : 'trial_data_rr_fixed_vmom.csv',
                    'filev12'  : 'trial_data_rr_fixed_vmom_cube.csv',
                  }

dictv_epfc_fixed = {
                    'nimom'     : nepfq_wvarex_end + 30*8 + 8*8+1,
                    'nimom_end' : nepfq_wvarex_end + 30*8 + 8*8 + 30+1,
                    'namom'     : nmom+1,
                    'namom_end' : nmom_end+1,
                    'njac'      : nepfq_wvarex_end + 8*8+1,
                    'njacend'   : nepfq_wvarex_end + 30*8 + 8*8+1,
                    'njaca'     : nepfq_wvarex_end+1,
                    'njacaend'  : nepfq_wvarex_end + 8*8+1,
                    'npar_est'  : 5,
                    'npar'      : 8,
                    'nresi'     : nres,
                    'nresi_end' : nres_end + 1,
                    'transposer': False,
                    'fixed'    : 0
                  }

dict3 = [[dictd_epfq_fixed, dictv_epfq_fixed], [dictd_epfc_fixed, dictv_epfc_fixed]]

# Generate out of sample variances

for ddi in dict3:
    print(ddi[0]['filename'])
    util.generate_var_outofsample(ddi[0], ddi[1], 4, upper=100)
    #epf_par_out(truepar2, ddi[0], ddi[1])


# Do analysis of Monte Carlo results
res3 = []
for ddi in dict3:
    print(ddi[0]['filename'])
    res3 += [util.epf_manipulate_alt(truepar2, ddi[0], ddi[1], rhigh=100, rlow=0, n=75000)]
#res3 += [util.epf_manipulate_alt(truepar2, dict3[1][0], dict3[1][1], rhigh=100, rlow=0, n=75000)]

sep3 = [' & ',' & ']

top1 = "\\begin{table}[hhh]\n  \\caption{Monte Carlo comparison of simulation estimators: Model with fixed adjustment costs}\n  \\begin{center}\n  \\begin{tabular*}{1\\textwidth}{r@{\\extracolsep{\\fill}}rrrrr}\n    \\hline\n            &     \\multicolumn{1}{c}{Second-order EPF-based}&     \\multicolumn{1}{c}{Third-order EPF-based} \\\\\n\\hline\n"

bottom1 = " \\hline \n \\end{tabular*} \n \\parbox[c]{7in}{\\footnotesize Indicated expectations and probabilities are estimates \n based on 1,000 Monte Carlo samples of size 75,000.  The samples are \n generated from the model in Section \\ref{sec:model}, with one addition. We allow for the presence of fixed costs of capital adjustment \n that take the form of $\\theta \\mathcal{I}(i\\ne0)$.  We consider two different benchmarks: a second-order approximation to the empirical policy function and \n  a third-order approximation. Both estimators use a clustered weight matrix. \n For each parameter we report three statistics. Bias is expressed \n as a percentage of the true coefficient value.  RMSE indicates root mean \n squared error and is also expressed as a percentage of the true \n coefficient.  $\\Pr(t)$ is the fraction of the time we observe a nominal \n 5\\% rejection of the null hypothesis that a parameter equals its true \n value using a {\\em t}-test. \n We report the fraction of trials that \n produce a rejection of three additional tests.  The first is a nominal 5\\% test of \n the model overidentifying restrictions.  The second, out-of-sample test \n has two varieties.  For the EPF-based estimator, it is a chi-squared \n test of the null hypothesis that the moments equal their true values. \n For the moments-based estimator, it is a chi-squared test of the null \n hypothesis that the policy function slopes equal their true values.  For \n the {\\em t}-tests of individual moment conditions, we report the \n highest, median, and lowest rejection rates.  The third test is a test \n of the null hypothesis that features of the data that are not used in \n the estimation equal their true values. \n The moments-based estimator minimizes the \n distance between simulated and data moments.  The EPF-based estimator \n minimizes the distance between policy functions estimated from simulated \n data and those estimated from real data.\n }\n \n \\end{center} \n\\end{table}"

dest = ['$\\delta$ (depreciation rate)', '$\\lambda$ (equity issuance cost)', '$\\xi$ (collateral parameter)', '$\\gamma$ (investment adjustment cost)', '$\\theta$ (fixed adjustment cost)']

mid1 = ''
print('Robustness fixed')

for ids, ds in enumerate(dest):
  mid1 += writeup.print_var_alt(ds, ids, res3, sep3, factor = 100,  lenm = writeup.find_max_len(dest))

mid1 += '\\\\[5pt]'

mid2 = writeup.mis_test_alt(75000, res3, sep3, lenm = writeup.find_max_len(dest))

tab3_pct = top1 +  mid1 + mid2 + bottom1


f = open('../WR/fixed_table.tex', 'w')
f.write(tab3_pct)
f.close()
