
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')
#maplotlib.rc('text', usetex=True)
matplotlib.rc('axes', edgecolor='k')
font = {'family':'sans-serif','sans-serif':['Helvetica'],
#        'weight' : 'normal',
        'size'   : 35}
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

matplotlib.rc('font', **font)
#maplotlib.rcParams['text.latex.preamble'] = [
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.linalg import pinv
from math import ceil

def tuftestyle(ax):
    """Styles an axes object to have minimalist graph features"""
    ax.grid(True, 'major', color='w', linestyle='-', linewidth=1.5,alpha=0.6)
    ax.patch.set_facecolor('white')

    #ax.set_axisbelow(True)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    #ax.yaxis.set_major_locator(MultipleLocator( (ax.get_yticks()[-1]-ax.get_yticks()[0]) / 0.1 ))
    #ax.get_xaxis().tick_bottom()
    #ax.get_yaxis().tick_left()

    #restyle the tick lines
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(5)
        line.set_color("black")
        line.set_markeredgewidth(0)

    for line in ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True):
        line.set_markersize(0)

    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.direction'] = 'in'
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

def mis_test_plot(d, name):
  mmed = [np.median(np.nanmean(np.abs(d1['tsn'][np.prod(np.isnan(d1['tsn'])==False,1)==1, :]) > 1.96,0)) for d1 in d]
  mmax = [np.max(np.nanmean(np.abs(d1['tsn'][np.prod(np.isnan(d1['tsn'])==False,1)==1, :]) > 1.96,0)) for d1 in d]
  mmin = [np.min(np.nanmean(np.abs(d1['tsn'][np.prod(np.isnan(d1['tsn'])==False,1)==1, :]) > 1.96,0)) for d1 in d]

  j = [np.nanmean(d1['pjstat'][np.prod(np.isnan(d1['tsn'])==False,1)==1]> 0.95) for d1 in d]
  jo = [np.nanmean(d1['pjstato'][np.prod(np.isnan(d1['tsn'])==False,1)==1] > 0.95) for d1 in d]

  debtc = [0.0, 0.0025, 0.005, 0.01, 0.02, 0.04]

  fig = plt.figure(figsize=(20,10))
  ax = fig.add_subplot(111)
  ax.plot(debtc, mmed, '-',lw=4)
  ax.plot(debtc, mmin, '--',lw=4)
  ax.plot(debtc, mmax, '-.',lw=4)
  plt.ylabel(r'Rejection rate')
  plt.xlabel(r'Debt cost')
  plt.tick_params(axis='both', which='major', labelsize=30, pad = 15)
  plt.tick_params(axis='both', which='minor', labelsize=30, pad = 15)
  ax.legend([r'Median',r'Minimum',r'Maximum'],loc=4,frameon=False,fontsize=30)
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.set_ylim([0,1.05])
  ax.set_xlim([0,0.04])
  plt.tight_layout()
  plt.savefig("../WR/reject_tstat_" + name + ".png")
  plt.close()

  fig = plt.figure(figsize=(20,10))
  ax = fig.add_subplot(111)
  ax.plot(debtc, j, '-',lw=4)
  plt.ylabel(r'Rejection rate')
  plt.xlabel(r'Debt cost')
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  plt.tick_params(axis='both', which='major', labelsize=30, pad = 15)
  plt.tick_params(axis='both', which='minor', labelsize=30, pad = 15)
  ax.set_ylim([0,1.05])
  ax.set_xlim([0,0.04])
  plt.tight_layout()
  plt.savefig("../WR/reject_jstat_" + name + ".png")
  plt.close()

  fig = plt.figure(figsize=(20,10))
  ax = fig.add_subplot(111)
  ax.plot(debtc, jo, '-',lw=4)
  plt.ylabel(r'Rejection rate')
  plt.xlabel(r'Debt cost')
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  plt.tick_params(axis='both', which='major', labelsize=30, pad = 15)
  plt.tick_params(axis='both', which='minor', labelsize=30, pad = 15)
  ax.set_ylim([0,1.05])
  ax.set_xlim([0,0.04])
  #tuftestyle(ax)
  plt.tight_layout()
  plt.savefig("../WR/reject_jstato_" + name + ".png")
  plt.close()


def strrnd(num, x = 10):
  if num == ' -nan':
    red = ''
  else:
    red = format(round(num,3),'.3f')
  red = red.ljust(10)
  return red

def sdstrrnd(num, x = 10):
  if num == ' -nan':
    red = ''
  else:
    red = format(round(num,3),'.3f')
  red = "(" + red + ")"
  red.ljust(10)
  return red

def find_max_len(des):
    desi = [len("\multicolumn{1}{l}{" + d + "} \\\\ ") for d in des ]
    return np.max(desi)

def print_var_alt(description, i, d, sep, factor = 1, lenm = 50):
  str1 = "\multicolumn{1}{l}{" + description + "} \\\\ \n"
  str1 += " Average $\%$ Bias".ljust(lenm)
  for iid, di in enumerate(d):
    str1 += sep[iid] + strrnd(di['bias_pct'][i]*factor)
  str1 += " \\\\\n"
  str1 += " RMSE $\%$".ljust(lenm)
  for iid, di in enumerate(d):
    str1 += sep[iid] + strrnd(di['rmse_pct'][i]*factor)
  str1 += " \\\\\n"
  str1 += " Mean Asymptotic SE $\%$".ljust(lenm)
  for iid, di in enumerate(d):
    str1 += sep[iid] + strrnd(di['sd_pct'][i]*factor)
  str1 += " \\\\\n"
  str1 += " $\Pr(t)$ ".ljust(lenm)
  for iid, di in enumerate(d):
    str1 += sep[iid] + strrnd(di['prt'][i])
  str1 += " \\\\\n  [10pt]\n "
  return(str1)


def mis_test_alt(size, d, sep, lenm = 50):
  m = []
  for di in d:
    print('Dropped: ' + str(np.mean(np.prod(np.isnan(di['tsn'])==False,1)==1)))
    print('Error on jstat: ' + str(np.mean(di['jstat']<0)))
    print('Error on out of sample jstat: ' + str(np.mean(di['jstato']<0)))
    m = m + [np.nanmean(np.abs(di['tsn'][np.prod(np.isnan(di['tsn'])==False,1)==1, :]) > 1.96, 0)]
  s1 = ''
  s1 += " Overidentification test rejection rate ".ljust(lenm)
  for iid, di in enumerate(d):
    s1 += sep[iid] + strrnd(np.nanmean(di['pjstat'][np.prod(np.isnan(di['tsn'])==False,1)==1] > 0.95))
  s1 += "\\\\ \n"
  s1 += " Out of sample test rejection rate ".ljust(lenm)
  for iid, di in enumerate(d):
    s1 += sep[iid] + strrnd(np.nanmean(di['pjstato'][np.prod(np.isnan(di['tsn'])==False,1)==1] > 0.95))
  s1 += "   \\\\ \n"
  s1 += "Moment {\\em t}-statistics:                  \\\\ \n"
  s1 += " \\mbox{  } maximum rejection rate ".ljust(lenm)
  for iid, mi in enumerate(m):
    s1 += sep[iid] + strrnd(np.max(mi))
  s1 += " \\\\ \n"
  s1 += " \\mbox{  } mean  rejection rate ".ljust(lenm)
  for iid, mi in enumerate(m):
    s1 += sep[iid] + strrnd(np.median(mi))
  s1 += " \\\\ \n"
  s1 += " \\mbox{  } minimum rejection rate ".ljust(lenm)
  for iid, mi in enumerate(m):
    s1 += sep[iid] + strrnd(np.min(mi))
  s1 += " \\\\\\\\ \n "
  return(s1)
