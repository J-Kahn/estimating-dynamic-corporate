
# coding: utf-8

# In[570]:

import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')
#maplotlib.rc('text', usetex=True)
font = {'family':'sans-serif','sans-serif':['Helvetica'],
#        'weight' : 'normal',
        'size'   : 35}
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlib.rc('axes', edgecolor='k')
matplotlib.rc('axes', linewidth = 2)
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

def quad(x, W):
    return np.dot(np.dot(x, W), x.transpose())

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

def quad(x, W):
    return np.dot(np.dot(x, W), x.transpose())

def nu_outofsample(jac1, jac2, omega1, omega2, omega12, w):
  return np.dot(np.dot(np.dot(jac2.transpose(), np.linalg.inv(quad(jac1, w),rcond=1e-3)), jac1), w)
def var_mom_outofsample(jac1, jac2, omega1, omega2, omega12, w):
  nu = np.dot(np.dot(np.dot(jac2.transpose(), np.linalg.inv(quad(jac1, w))), jac1), w)
  return omega2 - np.dot(nu, omega12) - np.dot(omega12.transpose(), nu.transpose()) + quad(nu, omega1)


# In[573]:

def get_var_outofsample(dictf, dictv):
    filev1 = dictf['vdataname']
    filev2 = dictf['filev2']
    filev12 = dictf['filev12']
    filew = dictf['wdataname']
    name = dictf['filename']
    dataname = dictf['dataname']
    dataname2 = dictf['dataname2']
    V1 = np.loadtxt("../Results/MC/" +  filev1, delimiter=",")
    V2 = np.loadtxt("../Results/MC/" + filev2, delimiter=",")
    V12 = np.loadtxt("../Results/MC/" + filev12, delimiter=",")
    return [V1, V2, V12]


# In[ ]:

def generate_var_outofsample(dictf, dictv, n_par, upper=30, direct = "", jchoice = -1, nhigh = 1000, nt = 75000):
    filev1 = dictf['vdataname']
    filev2 = dictf['filev2']
    filev12 = dictf['filev12']
    filew = dictf['wdataname']
    name = dictf['filename']
    dataname = dictf['dataname']
    dataname2 = dictf['dataname2']
    V1 = np.loadtxt("../Results/MC/" +  filev1, delimiter=",")
    V2 = np.loadtxt("../Results/MC/" + filev2, delimiter=",")
    V12 = np.loadtxt("../Results/MC/" + filev12, delimiter=",")
    m1 = np.loadtxt("../Results/MC/" + dataname, delimiter=",")
    m2 = np.loadtxt("../Results/MC/" + dataname2, delimiter=",")
    m1 = m1[0:nhigh,:]
    m2 = m2[0:nhigh,:]
    #v12 = np.cov(m1.transpose(), m2.transpose())
    #v12 = v12[:np.shape(m1)[1], np.shape(m1)[1]:] * nt
    if filew == 'identity':
      W = np.identity(4)
    else:
      W = np.loadtxt("../Results/MC/" + filew, delimiter=",")
    if 'fixed' in dictv:
        name = name + "_fix"
    if 'njac' in dictv:
      njac = dictv['njac']
    else:
      njac = njac1

    if 'njaca' in dictv:
      njaca = dictv['njaca']
    else:
      njaca = njac2

    if 'njacend' in dictv:
      njacend = dictv['njacend']
    else:
      njacend  = njac1_end
    
    if 'njacaend' in dictv:
      njacaend = dictv['njacaend']
    else:
      njacaend  = njac2_end
    
    if 'nimom' in dictv:
      nimom = dictv['nimom']
    else:
      nimom = nmom

    if 'nimom_end' in dictv:
      nimom_end = dictv['nimom_end']
    else:
      nimom_end = nmom_end

    if 'namom' in dictv:
      namom = dictv['namom']
    else:
      namom = nepfq

    if 'namom_end' in dictv:
      namom_end = dictv['namom_end']
    else:
      namom_end = nepfq_end

    if 'transposer' in dictv:
      transposer = dictv['transposer']
    else:
      transposer = True
    
    if 'npar_est' in dictv:
      npar_est = dictv['npar_est']
    else:
      npar_est = 4
    
    if 'npar' in dictv:
      npar = dictv['npar']
    else:
      npar = 7
    
    if 'jacname' in dictf:
      jacname = dictf['jacname']
    else:
      jacname = 'jacobians_clean' + name + "_"
    
    n_mom2 = namom_end - namom
    n_mom  = nimom_end - nimom
    if jchoice == -3:
        Jaccests = []
        for i in range(int(0),int(1000)):
            try:
                Jaccests = Jaccests + [pd.read_csv("../Results/MC/Trials/jacobians_" + jacname + "_" + str(i),header=None)]
            except:
                continue
        Jaccests = np.array(pd.concat(Jaccests))
        #Jaccests[np.abs(Jaccest)<1e-8] = 0
    err = 0
    skip = 0
    starter = 0
    for i in range(0, upper):
        try:
            Summest = np.array(pd.read_csv("../Results/MC/Trials/" + direct + "summary_" + name + "_" + str(i),header=None))
            if jchoice >= 0:
                Jaccest = np.array(pd.read_csv("../Results/MC/Trials/" + direct + jacname + str(i),header=None))
                Summest = Summest[0:np.shape(Jaccest)[0],:]
            if jchoice == -2:
                Jaccest = np.array(pd.read_csv("../Results/MC/Trials/" + direct + "summary_" + jacname + "_" + str(i),header=None))
                Summest = Summest[0:np.shape(Jaccest)[0],:]
            if jchoice == -3:
                n = np.shape(Summest)[0]
                Jaccest = Jaccests[starter:starter+n,:]
                starter = starter + n
            n = np.shape(Summest)[0]
            varoos = np.zeros((n, n_mom2 * n_mom2))
            for j in range(0, n):
                try:
                    if jchoice == -1:
                        #print(Summest)
                        jac1 = Summest[j,njac:njacend].astype(float).reshape(n_mom, npar).transpose()
                        jac2 = Summest[j,njaca:njacaend].astype(float).reshape(n_mom2, npar).transpose()
                        #print(jac1)
                        jac1 = jac1[0:npar_est, :]
                        jac2 = jac2[0:npar_est, :]
                    elif jchoice == -2:
                        #print(Summest)
                        jac1 = Jaccest[j,njac:njacend].astype(float).reshape(n_mom, npar).transpose()
                        jac2 = Jaccest[j,njaca:njacaend].astype(float).reshape(n_mom2, npar).transpose()
                        #print(jac1)
                        jac1 = jac1[0:npar_est, :]
                        jac2 = jac2[0:npar_est, :]
                    elif jchoice == -3:
                        jac1  = Jaccest[j,  0*(njac2_end - njac1) +  njac - njac1: 0 * (njac2_end - njac1) + njacend  - njac1].reshape(n_mom, npar).transpose()
                        jac2  = Jaccest[j,  0*(njac2_end - njac1) + njaca - njac1: 0 * (njac2_end - njac1) + njacaend - njac1].reshape(n_mom2, npar).transpose()
                        jac1 = jac1[0:npar_est, :]
                        jac2 = jac2[0:npar_est, :]
                    else:
                        jac1  = Jaccest[j, jchoice * (njac2_end - njac1) +  njac - njac1: jchoice * (njac2_end - njac1) + njacend  - njac1].reshape(n_mom, npar).transpose()
                        jac2  = Jaccest[j, jchoice * (njac2_end - njac1) + njaca - njac1: jchoice * (njac2_end - njac1) + njacaend - njac1].reshape(n_mom2, npar).transpose()
                        jac1 = jac1[0:npar_est, :]
                        jac2 = jac2[0:npar_est, :]
                    v1 = V1[Summest[j,0], :].reshape(n_mom, n_mom)
                    v2 = V2[Summest[j,0], :].reshape(n_mom2, n_mom2)

                    if transposer:
                        v12 = V12[Summest[j,0], :].reshape(n_mom2, n_mom).transpose()
                    else:
                        v12 = V12[Summest[j,0], :].reshape(n_mom, n_mom2)
                    if filew == 'identity':
                        w = np.identity(n_mom)
                    else:
                        w = W[Summest[j,0], :].reshape(n_mom, n_mom)
                    try:
                      varoos[j, :] = var_mom_outofsample(jac1, jac2, v1, v2, v12, w).reshape(1,n_mom2*n_mom2)
                    except:
                        raise
                except:
                  varoos[j, :] = np.nan
                  err += 1
                  print(name + " errors: " + str(err))
                jacdes = ""
                #if jchoice >= 0:
                #    jacdes = '_jac_' + str(jchoice)
                np.savetxt('../Results/MC/Trials/' + direct + 'varoos_' + name + "_" + str(i) + jacdes, varoos, delimiter=',')
        except:
            print(name + ' skipped ' + str(i))
            skip +=1
    print('Total incomplete ' + str(err + skip))

    

        




def epf_manipulate_alt(truepar, dictd, dictv, n=75000, rlow=0, rhigh=500, offset=0, per = 1, getlow = -1, gethigh = -1, direct="", jchoice = -1):
  filename = dictd['filename']
  dataname = dictd['dataname']
  dataname2 = dictd['dataname2']
  wdataname = dictd['wdataname']
  vdataname = dictd['vdataname']
  if 'fixed' in dictv:
        filename = filename + "_fix"
  Summest = []
  for i in range(int(rlow/per + offset/per),int(ceil(rhigh/per + offset/per))):
      try:
          Summest = Summest + [pd.read_csv("../Results/MC/Trials/" + direct + "summary_" + filename + "_" + str(i),header=None)]
      except:
          continue
  jacdes = ""
  if jchoice >= 0:
      jacdes = '_jac_' + str(jchoice)
  Varoos = []
  for i in range(int(rlow/per + offset/per),int(ceil(rhigh/per + offset/per))):
      try:
            Varoos = Varoos +  [pd.read_csv("../Results/MC/Trials/" + direct + "varoos_" + filename + "_" + str(i),header=None)]
      except:
          continue
  Summest = np.array(pd.concat(Summest))
  Varoos = np.array(pd.concat(Varoos))
  Qdat = np.array(pd.read_csv("../Results/MC/" + dataname, header=None))
  Qdat2 = np.array(pd.read_csv("../Results/MC/" + dataname2, header=None))
  Wdat = np.array(pd.read_csv("../Results/MC/" + wdataname, header=None))
  Vdat = np.array(pd.read_csv("../Results/MC/" + vdataname, header=None))
  #older = Summest[:,1] > 0
  #Summest = Summest[older, :]
  #Varoos = Varoos[older, :]
  #jacw = np.prod(np.isnan(Summest[:,njac1:njac1_end].astype(float))==False,1)

  if 'njac' in dictv:
    njac = dictv['njac']
  else:
    njac = njac1

  if 'njaca' in dictv:
    njaca = dictv['njaca']
  else:
    njaca = njac2

  if 'njacend' in dictv:
    njacend = dictv['njacend']
  else:
    njacend  = njac1_end

  if 'njacaend' in dictv:
    njacaend = dictv['njacaend']
  else:
    njacaend  = njac2_end

  if 'nimom' in dictv:
    nimom = dictv['nimom']
  else:
    nimom = nmom

  if 'nimom_end' in dictv:
    nimom_end = dictv['nimom_end']
  else:
    nimom_end = nmom_end

  if 'namom' in dictv:
    namom = dictv['namom']
  else:
    namom = nepfq

  if 'namom_end' in dictv:
    namom_end = dictv['namom_end']
  else:
    namom_end = nepfq_end

  if 'npar_est' in dictv:
    npar_est = dictv['npar_est']
  else:
    npar_est = 4
  if 'nresi' in dictv:
    nresi = dictv['nresi']
  else:
    nresi = nres

  if 'nresi_end' in dictv:
    nresi_end = dictv['nresi_end']
  else:
    nresi_end = nres_end

  if 'jacname' in dictd:
    jacname = dictd['jacname']
  else:
    jacname = 'jacobians_' + filename

  n_mom = nimom_end - nimom
  print(n_mom)
  n_mom2 = namom_end - namom

  if jchoice == -1:
    jacs = Summest[:, njac:njacend]
    jacs2 = Summest[:, njaca:njacaend]
    #jacs = Summest[:, njac:njacend]
    #jacs2 = Summest[:, njaca:njacaend]
    jacs[jacs == ' '] = ' -nan'
    jacs2[jacs2 == ' '] = ' -nan'
    jacs = jacs.astype(float)
    jacs2 = jacs2.astype(float)
  elif jchoice == -2:
    Jaccest = []
    for i in range(int(rlow/per + offset/per),int(ceil(rhigh/per + offset/per))):
          try:
              Jaccest = Jaccest + [pd.read_csv("../Results/MC/Trials/" + direct + "summary_" + jacname + "_" + str(i),header=None)]
          except:
              continue
    Jaccest = np.array(pd.concat(Jaccest))
    if np.shape(Jaccest)[0] > np.shape(Summest)[0]:
        Jaccest = Jaccest[0:np.shape(Summest)[0],:]
    jacs = Jaccest[:, njac:njacend]
    jacs2 = Jaccest[:, njaca:njacaend]
    #jacs = Summest[:, njac:njacend]
    #jacs2 = Summest[:, njaca:njacaend]
    jacs[jacs == ' '] = ' -nan'
    jacs2[jacs2 == ' '] = ' -nan'
    jacs = jacs.astype(float)
    jacs2 = jacs2.astype(float)
  else:
    Jaccest = []
    for i in range(int(rlow/per + offset/per),int(ceil(rhigh/per + offset/per))):
          try:
              Jaccest = Jaccest + [pd.read_csv("../Results/MC/Trials/" + direct + jacname + "_" + str(i),header=None)]
          except:
              continue
    Jaccest = np.array(pd.concat(Jaccest))
    #Jaccest = Jaccest[ind,:]
    #Jaccest = Jaccest[selected,:]
    jacs  = Jaccest[:,jchoice * (njac2_end - njac1) +  njac - njac1: jchoice * (njac2_end - njac1) + njacend  - njac1].astype(float)
    jacs2 = Jaccest[:,jchoice * (njac2_end - njac1) + njaca - njac1: jchoice * (njac2_end - njac1) + njacaend - njac1].astype(float)
  jacw = np.prod(np.isnan(jacs)==False,1)==1
  Summest = Summest[0:np.shape(jacs)[0],:]
  Varoos = Varoos[0:np.shape(jacs)[0],:]
  Summest = Summest[jacw, :]
  Varoos = Varoos[jacw, :]
  jacs = jacs[jacw,:]
  jacs2 = jacs2[jacw,:]

  #if n_mom == 8:
  #  jacw = np.prod(jacs2!=0,1)==1
  #else:
  #  jacw = np.prod(jacs!=0,1)==1
  #print(jacw)
  #Summest = Summest[jacw, :]
  #Varoos = Varoos[jacw, :]
  #jacs = jacs[jacw,:]
  #jacs2 = jacs2[jacw,:]
  trash, ind = np.unique(Summest[:,0], return_index = True)
  Summest = Summest[ind,:]
  Varoos = Varoos[ind,:]
  jacs = jacs[ind,:]
  jacs2 = jacs2[ind,:]
  if getlow < 0:
    getlow = np.nanmin(Summest[:,0].astype(float))
  if gethigh < 0:
    gethigh = np.nanmax(Summest[:,0].astype(float)) + 1
  selected = np.array([(si >= getlow) and (si < gethigh) for si in Summest[:,0].astype(float)])
  Summest = Summest[selected,:]
  Varoos = Varoos[selected,:]
  jacs = jacs[selected,:]
  jacs2 = jacs2[selected,:]

  Mest = np.array(Summest[:,nresi:nresi_end]).astype(float)
  n_par = Mest.shape[1]
  print(n_par)

  err = Mest - truepar
  #print(err)
  bias = np.nanmean(err,0)
  mse = np.nanmean(err**2,0).astype(float)
  #print(mse)
  #rmse = np.sqrt(np.array(mse))
  rmse = 0
  #for i in range(0, 4):
  #  bias[i] = np.nanmean(err[:,i][np.abs(err[:,i])>0])
  #  mse[i] = np.nanmean(err[:,i][np.abs(err[:,i])>0]**2)
  #print(mse)
  rmse = np.sqrt(np.array(mse))

  sd = np.zeros((Mest.shape[0], npar_est))
  ub = np.zeros((Mest.shape[0], npar_est))
  lb = np.zeros((Mest.shape[0], npar_est))
  jacworked = np.prod(np.isnan(jacs)==False,1)
  njacworked = np.sum(jacworked)
  ts = np.zeros((njacworked, npar_est))
  inb = np.zeros((njacworked, Mest.shape[1]))
  jn = 0
  tsn = np.zeros((njacworked, n_mom))
  tso = np.zeros((njacworked, n_mom))
  jstato = np.zeros(njacworked)
  Diff = np.zeros((Summest.shape[0], n_mom))
  Diff2 = np.zeros((njacworked, n_mom))
  Diffo = np.zeros((njacworked, n_mom2))
  Vo = np.zeros((njacworked,n_mom2 * n_mom2))
  V = np.zeros((njacworked, n_mom * n_mom))
  VV = np.zeros((njacworked, npar_est * npar_est))
  J2 = np.zeros((njacworked, n_mom2 * npar_est))
  J = np.zeros((njacworked, n_mom * npar_est))
  for i in range(0, Summest.shape[0]):
    jac = jacs[i, :].reshape(n_mom, n_par).transpose().astype(float)
    jac  = jac[0:npar_est,:]
    jac[np.abs(jac)<1e-8] = 0
    #jac[np.abs(jac) > 50] = 0
    jac2 = jacs2[i, :].reshape(n_mom2, n_par).transpose().astype(float)
    jac2 = jac2[0:npar_est,:]
    jac2[np.abs(jac2)<1e-8] = 0
    #jac2[np.abs(jac2) > 50] = 0
    #if moment:
    #    jac = JM
    #    jac2 = JE
    #else:
    #    jac = JE
    #    jac2 = JM
    #if n_mom == 8:
    #    if np.sum(jacs2[i,:]==0)>0:
    #        jac2 = jac2 * np.nan
    #else:
    #    if np.sum(jacs[i,:]==0)>0:
    #        jac = jac * np.nan
    try:
        #print(Summest[i,nimom:nimom_end])
        #print(Qdat[Summest[i,0],:])
        diff = Summest[i, nimom:nimom_end].astype(float) - Qdat[Summest[i,0],:]
    except:
        diff = Summest[i, nimom:nimom_end]*np.nan    #for i in range(0, n_par):
    #    for j in range(0, n_mom):
    #        jac[i, j] = float(jac[i,j])
    Diff[i, :] = diff
    try:
        w = Wdat[Summest[i,0], :].reshape(n_mom, n_mom)
        v = Vdat[Summest[i,0], :].reshape(n_mom, n_mom)
        diff2 = Qdat2[Summest[i,0], :] - Summest[i, namom:namom_end]
    except:
        w = np.ones((n_mom,n_mom)) * np.nan
        v = w
        diff2 = Summest[i, namom:namom_end]*np.nan
    try:
        jwj = np.linalg.inv(quad(jac, w))
        #jwjj = np.dot(jwj, jac)
        jwjj = np.linalg.solve(quad(jac,w),jac)
    except:
        jwj = np.ones((npar_est, npar_est)) * np.nan
        jwjj = np.ones((npar_est,n_mom)) * np.nan
    avar = jwj/n
    sd[i, :] = np.sqrt(np.diag(avar)*(1.0+1.0/10.0))
    ub[i, :]  = Mest[i, 0:npar_est] + 1.96 * sd[i, :]
    lb[i, :]  = Mest[i, 0:npar_est] - 1.96 * sd[i, :]
    if jacworked[i] == 1:
        for j in range(0, npar_est):
          if (ub[i, j] > truepar[j]) and (lb[i,j] < truepar[j]):
            inb[jn,j] = 1.0
        cv = np.dot(np.dot(jac.transpose(), jwjj), np.dot(w, v))
        vv = (v - cv - cv.transpose() + quad(np.dot(jac.transpose(), jwjj), quad(w, v))) * (1 + 10/11)
        try:
            wnew2 = np.linalg.pinv(Varoos[i, :].reshape(namom_end - namom, namom_end - namom)) * n * 10/11
            #if n_mom == 8:
            #    if np.sum(jacs2[i,:]==0)>0:
            #        wnew2 = wnew2 * np.nan
        except:
            wnew2 = np.ones((namom_end-namom,namom_end-namom))*np.nan
        jstato[jn] =  quad(diff2,wnew2)
        #if n_mom == 8:
        #    if np.sum(jacs2[i,:]==0)>0:
        #        jstato[jn] = np.nan
        Diff2[jn, :] = diff
        Diffo[jn, :] = diff2
        J2[jn, :] = jac2.reshape(n_mom2 * npar_est)
        J[jn, :] = jac.reshape(n_mom * npar_est)
        Vo[jn,:] = Varoos[i,:]
        V[jn, :] = vv.reshape(n_mom * n_mom)
        VV[jn, :] =  avar.reshape(1,npar_est*npar_est)
        try:
            ts[jn, :] = err[i, 0:npar_est] / sd[i, :] # / np.sqrt(n)
        except:
            ts[jn, :] = np.nan
        tsn[jn, :] = diff / np.sqrt(np.diag(vv)/n)
        jn += 1
  jstat_true = np.abs(Summest[isfloat(Summest[:,2]),2].astype(float)*n*(10/11))
  jstat = np.abs(Summest[isfloat(Summest[:,1]),1].astype(float)*n*(10/11))
  prt = np.nanmean(np.abs(ts[np.prod(np.isnan(ts)==False,1)==1, :]) > 1.96,0)
  dist = stats.chi2(n_mom - npar_est)
  dist2 = stats.chi2(n_mom)
  dist3 = stats.chi2(n_mom2)
  pjstat = dist.cdf(jstat)
  pjstat2 = dist2.cdf(jstat_true)
  pjstato = dist3.cdf(jstato)
  dnorm = stats.norm()
  results = {'mest' : np.nanmean(Mest,0), 'mestf' : Mest, 'bias' : bias, 'bias_pct' : bias/truepar,
             'mse' : mse, 'rmse' : rmse, 'rmse_pct' : rmse/truepar, 'sd_pct' : np.mean(sd[np.prod(np.isnan(tsn)==False,1)==1,:],0) / truepar[0:npar_est], 'sd' : sd, 'inb' : inb, 'prt': prt, 'summ' : Summest, 'jstat' : jstat,
            'jstat_true' : jstat_true, 'pjstat' : pjstat, 'pjstat2' : pjstat2, 'err' : err, 'ts': ts, 'tsn' : tsn,
             'jstato' : jstato, 'pjstato' : pjstato, 'qdat' : Qdat, 'qdat2' : Qdat2, 'w' : Wdat, 'jacs' : jacs, 'jacs2' : jacs2,
             'diff' : Diff, 'diff2' : Diff2, 'v' : V, 'vv' : VV, 'diffo' : Diffo, 'vo' : Vo, 'ju2' : J2, 'ju1' : J, 'dat' : Qdat}

  plt.figure(figsize=(20,20))
  plt.plot(dist.cdf(np.nanpercentile(jstat, np.arange(0,101))),np.arange(0,101)/100.0,color='r',lw = 4)
  plt.plot(dist2.cdf(np.nanpercentile(jstat_true, np.arange(0,101))),np.arange(0,101)/100.0,'b--',lw = 4)
  plt.plot(np.arange(0,101)/100.0, np.arange(0,101)/100.0,'k',lw = 2)
  try:
      plt.plot(dist3.cdf(np.nanpercentile(jstato, np.arange(0,101))),np.arange(0,101)/100.0,'g-.',lw = 4)
  except:
      print("No out of sample moments")
  plt.xlabel(r'Theoretical percentile',fontsize=40)
  plt.ylabel(r'Actual percentile',fontsize=40)
  plt.legend([r'Estimated parameters',r'True parameters', r'Theoretical',r'Out-of-sample test'],loc=2,frameon=False,fontsize = 35)
  plt.savefig("../WR/chi2plot_" + filename + ".png")
  plt.close()
  try:
      plt.figure()
      for i in range(0,tsn.shape[1]):
        plt.plot(np.nanpercentile(dnorm.cdf(tsn[:,i]), np.arange(0,101)),np.arange(0,101)/100.0)
      plt.plot(np.arange(0,101)/100.0, np.arange(0,101)/100.0,'k--')
      plt.xlabel('Theoretical percentile')
      plt.ylabel('Actual percentile')
      plt.savefig("../WR/tplot_" + filename + ".png")
      plt.close()
  except:
      print("no plot")
  try:
      plt.figure()
      for i in range(0,ts.shape[1]):
        plt.plot(np.nanpercentile(dnorm.cdf(ts[:,i]), np.arange(0,101)),np.arange(0,101)/100.0)
      plt.plot(np.arange(0,101)/100.0, np.arange(0,101)/100.0,'k--')
      plt.xlabel('Theoretical percentile')
      plt.ylabel('Actual percentile')
      plt.savefig("../WR/tplot_par_" + filename + ".png")
      plt.close()
  except:
      print("no plot")
  return results



from scipy.interpolate import spline





def epf_ineff_manipulate_alt(truepar, dictd, dictv, n=75000, rlow=0, rhigh=500, offset=0, per = 1, getlow = -1, gethigh = -1, direct="", jchoice = -1):
  filename = dictd['filename']
  dataname = dictd['dataname']
  dataname2 = dictd['dataname2']
  wdataname = dictd['wdataname']
  vdataname = dictd['vdataname']
  if 'filename2' in dictd:
      filename2 = dictd['filename2']
  if 'fixed' in dictv:
        filename = filename + '_fix'
  Summest = []
  for i in range(int(rlow/per + offset/per),int(ceil(rhigh/per + offset/per))):
      try:
          Summest = Summest + [pd.read_csv("../Results/MC/Trials/" + direct + "summary_" + filename + "_" + str(i),header=None)]
      except:
          continue
  if 'filename2' in dictd:
      Summest2 = []
      for i in range(int(rlow/per + offset/per),int(ceil(rhigh/per + offset/per))):
          try:
              Summest2 = Summest2 + [pd.read_csv("../Results/MC/Trials/" + direct + "summary_" + filename2 + "_" + str(i),header=None)]
          except:
              continue
  else:
      Summest2 = Summest
  jacdes = ""
  if jchoice >= 0:
      jacdes = '_jac_' + str(jchoice)
  Varoos = []
  for i in range(int(rlow/per + offset/per),int(ceil(rhigh/per + offset/per))):
      try:
            Varoos = Varoos +  [pd.read_csv("../Results/MC/Trials/" + direct + "varoos_" + filename + "_" + str(i),header=None)]
      except:
          continue
  Summest = np.array(pd.concat(Summest))
  Summest2 = np.array(pd.concat(Summest2))
  Varoos = np.array(pd.concat(Varoos))
  Qdat = np.array(pd.read_csv("../Results/MC/" + dataname, header=None))
  Qdat2 = np.array(pd.read_csv("../Results/MC/" + dataname2, header=None))
  if wdataname == 'identity':
    Wdat = np.identity(8)
  else:
    Wdat = np.array(pd.read_csv("../Results/MC/" + wdataname, header=None))
  Vdat = np.array(pd.read_csv("../Results/MC/" + vdataname, header=None))
  #older = Summest[:,1] > 0
  #Summest = Summest[older, :]
  #Varoos = Varoos[older, :]
  #jacw = np.prod(np.isnan(Summest[:,njac1:njac1_end].astype(float))==False,1)

  if 'njac' in dictv:
    njac = dictv['njac']
  else:
    njac = njac1

  if 'njaca' in dictv:
    njaca = dictv['njaca']
  else:
    njaca = njac2

  if 'njacend' in dictv:
    njacend = dictv['njacend']
  else:
    njacend  = njac1_end

  if 'njacaend' in dictv:
    njacaend = dictv['njacaend']
  else:
    njacaend  = njac2_end

  if 'nimom' in dictv:
    nimom = dictv['nimom']
  else:
    nimom = nmom

  if 'nimom_end' in dictv:
    nimom_end = dictv['nimom_end']
  else:
    nimom_end = nmom_end

  if 'namom' in dictv:
    namom = dictv['namom']
  else:
    namom = nepfq

  if 'namom_end' in dictv:
    namom_end = dictv['namom_end']
  else:
    namom_end = nepfq_end

  if 'npar_est' in dictv:
    npar_est = dictv['npar_est']
  else:
    npar_est = 4

  if 'nresi' in dictv:
    nresi = dictv['nresi']
  else:
    nresi = nres

  if 'nresi_end' in dictv:
    nresi_end = dictv['nresi_end']
  else:
    nresi_end = nres_end

  if 'jacname' in dictd:
    jacname = dictd['jacname']
  else:
    jacname = 'jacobians_' + filename

  n_mom = nimom_end - nimom
  print(n_mom)
  n_mom2 = namom_end - namom

  if jchoice == -1:
    jacs = Summest[:, njac:njacend]
    jacs2 = Summest[:, njaca:njacaend]
    #jacs = Summest[:, njac:njacend]
    #jacs2 = Summest[:, njaca:njacaend]
    jacs[jacs == ' '] = ' -nan'
    jacs2[jacs2 == ' '] = ' -nan'
    jacs = jacs.astype(float)
    jacs2 = jacs2.astype(float)
  elif jchoice == -2:
    Jaccest = []
    for i in range(int(rlow/per + offset/per),int(ceil(rhigh/per + offset/per))):
          try:
              Jaccest = Jaccest + [pd.read_csv("../Results/MC/Trials/" + direct + "summary_" + jacname + "_" + str(i),header=None)]
          except:
              continue
    Jaccest = np.array(pd.concat(Jaccest))
    if np.shape(Jaccest)[0] > np.shape(Summest)[0]:
        Jaccest = Jaccest[0:np.shape(Summest)[0]]
    jacs = Jaccest[:, njac:njacend]
    jacs2 = Jaccest[:, njaca:njacaend]
    #jacs = Summest[:, njac:njacend]
    #jacs2 = Summest[:, njaca:njacaend]
    jacs[jacs == ' '] = ' -nan'
    jacs2[jacs2 == ' '] = ' -nan'
    jacs = jacs.astype(float)
    jacs2 = jacs2.astype(float)
  else:
    Jaccest = [pd.read_csv("../Results/MC/Trials/" + direct + jacname + "_" + str(i),header=None) for i in range(int(rlow/per + offset/per),int(ceil(rhigh/per + offset/per)))]
    Jaccest = np.array(pd.concat(Jaccest))
    Jaccest = Jaccest[ind,:]
    Jaccest = Jaccest[selected,:]
    jacs  = Jaccest[:,jchoice * (njac2_end - njac1) +  njac - njac1: jchoice * (njac2_end - njac1) + njacend  - njac1].astype(float)
    jacs2 = Jaccest[:,jchoice * (njac2_end - njac1) + njaca - njac1: jchoice * (njac2_end - njac1) + njacaend - njac1].astype(float)
  jacw = np.prod(np.isnan(jacs)==False,1)==1
  Summest = Summest[0:np.shape(jacs)[0],:]
  Varoos = Varoos[0:np.shape(jacs)[0],:]
  Summest = Summest[jacw, :]
  Varoos = Varoos[jacw, :]
  jacs = jacs[jacw,:]
  jacs2 = jacs2[jacw,:]

  trash, ind = np.unique(Summest[:,0], return_index = True)
  Summest = Summest[ind,:]
  Varoos = Varoos[ind,:]
  jacs = jacs[ind,:]
  jacs2 = jacs2[ind,:]

  if getlow < 0:
    getlow = np.nanmin(Summest[:,0].astype(float))
  if gethigh < 0:
    gethigh = np.nanmax(Summest[:,0].astype(float)) + 1
  selected = np.array([(si >= getlow) and (si < gethigh) for si in Summest[:,0].astype(float)])
  Summest = Summest[selected,:]
  Varoos = Varoos[selected,:]
  jacs = jacs[selected,:]
  jacs2 = jacs2[selected,:]

  Mest = np.array(Summest[:,nres:nres_end]).astype(float)
  n_par = Mest.shape[1]
  print(n_par)

  err = Mest - truepar
  #print(err)
  bias = np.nanmean(err,0)
  mse = np.nanmean(err**2,0).astype(float)
  #print(mse)
  #rmse = np.sqrt(np.array(mse))
  rmse = 0
  #for i in range(0, 4):
  #  bias[i] = np.nanmean(err[:,i][np.abs(err[:,i])>0])
  #  mse[i] = np.nanmean(err[:,i][np.abs(err[:,i])>0]**2)
  #print(mse)
  rmse = np.sqrt(np.array(mse))
  sd = np.zeros((Mest.shape[0], 4))
  ub = np.zeros((Mest.shape[0], 4))
  lb = np.zeros((Mest.shape[0], 4))

  jacworked = np.prod(np.isnan(jacs)==False,1)
  njacworked = np.sum(jacworked)
  inb = np.zeros((njacworked, Mest.shape[1]))
  ts = np.zeros((njacworked, 4))
  jn = 0
  jstat = np.zeros(njacworked)
  jstato = np.zeros(njacworked)
  tsn = np.zeros((njacworked, n_mom))
  Diff = np.zeros((Summest.shape[0], n_mom))
  Diff2 = np.zeros((njacworked, n_mom))
  Diffo = np.zeros((njacworked, n_mom2))
  J = np.zeros((njacworked, n_mom * 4))
  J2 = np.zeros((njacworked, n_mom2 * 4))
  JWJ = np.zeros((njacworked, 4 * 4))
  CV = np.zeros((njacworked, n_mom * 4))
  V = np.zeros((njacworked, n_mom * n_mom))
  VV = np.zeros((njacworked, n_mom * n_mom))
  VVV = np.zeros((njacworked, n_mom * n_mom))
  VO = np.zeros((njacworked, n_mom2 * n_mom2))
  for i in range(0, Summest.shape[0]):
    jac = jacs[i, :].reshape(n_mom, n_par).transpose().astype(float)
    jac  = jac[0:npar_est,:]
    jac[np.abs(jac)<1e-8] = 0
    #jac[np.abs(jac) > 50] = 0
    jac2 = jacs2[i, :].reshape(n_mom2, n_par).transpose().astype(float)
    jac2 = jac2[0:npar_est,:]
    jac2[np.abs(jac2)<1e-8] = 0
    #jac2[np.abs(jac2) > 50] = 0
    #if moment:
    #    jac2 = JE
    #    jac = JM
    #else:
    #    jac = JE
    #    jac2 = JM
    try:
        diff = Summest[i, nimom:nimom_end] - Qdat[Summest[i,0],:]
    except:
        diff = Summest[i, nimom:nimom_end]*np.nan
    #if n_mom == 8:
    #    if np.sum(jacs2[i,:]==0)>0:
    #        jac2 = jac2 * np.nan
    #else:
    #    if np.sum(jacs[i,:]==0)>0:
    #        jac = jac * np.nan
    #for i in range(0, n_par):
    #    for j in range(0, n_mom):
    #        jac[i, j] = float(jac[i,j])
    Diff[i, :] = diff
    try:
        if wdataname == 'identity':
          w = np.identity(n_mom)
        else:
          w = Wdat[Summest[i,0], :].reshape(n_mom, n_mom)
        v = Vdat[Summest[i,0], :].reshape(n_mom, n_mom)
        diff2 = Qdat2[Summest[i,0], :] - Summest[i, namom:namom_end]
    except:
        w = np.ones((n_mom,n_mom)) * np.nan
        v = w
        diff2 = Summest[i, namom:namom_end]*np.nan
    try:
        jwj = np.linalg.inv(quad(jac, w))
        #jwjj = np.dot(jwj,jac)
        jwjj = np.linalg.solve(quad(jac,w),jac)
    except:
        jwj = np.ones((4,4)) * np.nan
        jwjj = np.ones((4,n_mom)) * np.nan
    brd = quad(w, v)
    jbrd = quad(jac, brd)
    avar = quad(jwj, jbrd)/n
    sd[i, :] = np.sqrt(np.diag(avar)*(1+1/10))
    ub[i, :] = Mest[i, 0:npar_est] + 1.96 * sd[i, :]
    lb[i, :] = Mest[i, 0:npar_est] - 1.96 * sd[i, :]
    if jacworked[i] == 1:
        for j in range(0, 4):
          if (ub[i, j] > truepar[ j]) and (lb[i,j] < truepar[j]):
            inb[jn,j] = 1.0
        #bread = np.eye(n_mom) - np.dot(np.dot(jac.transpose(), jwjj), w)
        #cv = np.dot(np.dot(jac.transpose(), jwjj), np.dot(w, v))
        cv = np.dot(np.dot(jac.transpose(), jwjj), np.dot(w, v))
        vv = (v - cv - cv.transpose() + quad(np.dot(jac.transpose(), jwjj), quad(w, v)))
        try:
          wnew = np.linalg.pinv(vv) * n * 10 /11
          #wnew = np.linalg.pinv(quad(bread, v) * (1 + 1/10)) * n
          #wnew  = np.linalg.pinv((v - cv*(1-1/10) - cv.transpose()*(1-1/10) + quad(np.dot(jac.transpose(), jwjj), quad(w, v))*(1+1/10))) * n
        except:
          wnew  = np.ones((nimom_end-nimom,nimom_end-nimom))*np.nan
        try:
          wnew2 = np.linalg.pinv(Varoos[i, :].reshape(namom_end - namom, namom_end - namom)) * n * 10 / 11
        except:
          wnew2 = np.ones((namom_end-namom,namom_end-namom))*np.nan
        #print(quad(diff,wnew))
        jstat[jn] = quad(diff, wnew)
        jstato[jn] = quad(diff2, wnew2)
        #if n_mom == 8:
        #    if np.sum(jacs2[i,:]==0)>0:
        #        jstato[jn] = np.nan
        #else:
        #    if np.sum(jacs[i,:]==0)>0:
        #        jstat[jn] = np.nan
        Diff2[jn, :] = diff
        vv = v - cv*(1-1/10) - cv.transpose()*(1-1/10) + quad(np.dot(jac.transpose(), jwjj), quad(w, v))*(1+1/10)
        VV[jn, :] = (v - cv*(1-1/10) - cv.transpose()*(1-1/10) + quad(np.dot(jac.transpose(), jwjj), quad(w, v))*(1+1/10)).reshape(1, n_mom*n_mom)
        V[jn, :] = wnew.reshape(1, n_mom * n_mom)
        VVV[jn, :] = v.reshape(1, n_mom * n_mom)
        VO[jn, :] = Varoos[i, :]
        Diffo[jn, :] = diff2
        JWJ[jn, :] = jwj.reshape(1, 4 * 4)
        J2[jn, :] = jac2.reshape(n_mom2 * npar_est)
        J[jn, :] = jac.reshape(n_mom * npar_est)
        CV[jn, :] = np.dot(jwj, np.dot(jac, np.dot(w, v))).reshape(1, 4 * n_mom)
        try:
          ts[jn, :] = err[i, 0:npar_est] / sd[i, :]
        except:
          ts[jn, :] = np.nan
        tsn[jn, :] = diff / np.sqrt(np.diag(vv)/n)
        jn += 1

  jstat_true = np.abs(Summest2[isfloat(Summest2[:,2]),2].astype(float)*n*(10/11))
  prt = np.nanmean(np.abs(ts[np.prod(np.isnan(ts)==False,1)==1, :]) > 1.96,0)
  dist = stats.chi2(n_mom - 4)
  dist2 = stats.chi2(n_mom)
  dist3 = stats.chi2(n_mom2)
  dnorm = stats.norm()
  pjstat = dist.cdf(jstat)
  pjstat2 = dist2.cdf(jstat_true)
  pjstato = dist3.cdf(jstato)
  if wdataname == 'ientity':
    Wdat = np.identity(n_mom)


  results = {'mest' : np.nanmean(Mest,0), 'mestf' : Mest, 'bias' : bias, 'bias_pct' : bias/truepar,
             'mse' : mse, 'rmse' : rmse, 'rmse_pct' : rmse/truepar, 'sd_pct' : np.mean(sd[np.prod(np.isnan(tsn)==False,1)==1,:]
,0) / truepar[0:npar_est], 'sd' : sd, 'inb' : inb, 'prt': prt, 'summ' : Summest, 'jstat' : jstat,
            'jstat_true' : jstat_true, 'pjstat' : pjstat, 'pjstat2' : pjstat2, 'err' : err, 'ts' : ts, 'tsn' : tsn, 'jstato' : jstato, 'pjstato' : pjstato,
             'w' : Wdat, 'diff' : Diff, 'diff2' : Diff2, 'v' : V, 'vv' : VV, 'vvv' : VVV, 'cv' : CV, 'jwj' : JWJ, 'ju1' : J, 'ju2': J2, 'vfile' : Vdat, 'vo' : VO, 'diffo' : Diffo, 'dat' : Qdat}

  plt.figure(figsize=(20,20))
  try:
      plt.plot(dist.cdf(np.nanpercentile(jstat, np.arange(0,101))),np.arange(0,101)/100.0,color='r',lw=4)
      plt.plot(dist2.cdf(np.nanpercentile(jstat_true, np.arange(0,101))),np.arange(0,101)/100.0,'b--',lw = 4)
  except:
      print("No in sample moments")
  plt.plot(np.arange(0,101)/100.0, np.arange(0,101)/100.0,'k',lw=2)
  try:
      plt.plot(dist3.cdf(np.nanpercentile(jstato, np.arange(0,101))),np.arange(0,101)/100.0,'g-.',lw=4)
  except:
      print("No out of sample moments")
  plt.xlabel(r'Theoretical percentile',fontsize=40)
  plt.ylabel(r'Actual percentile',fontsize=40)
  plt.legend([r'Estimated parameters',r'True parameters',r'Theoretical',r'Out-of-sample'],loc=2,frameon=False,fontsize=35)
  plt.savefig("../WR/chi2plot_" + filename + ".png")
  plt.close()

  try:
      plt.figure()
      for i in range(0,tsn.shape[1]):
        plt.plot(np.nanpercentile(dnorm.cdf(tsn[:,i]), np.arange(0,101)),np.arange(0,101)/100.0)
      plt.plot(np.arange(0,101)/100.0, np.arange(0,101)/100.0,'k--')
      plt.xlabel('Theoretical percentile')
      plt.ylabel('Actual percentile')
      plt.savefig("../WR/tplot_" + filename + ".png")
      plt.close()
  except:
      print("no plot")

  try:
      plt.figure()
      for i in range(0,ts.shape[1]):
        plt.plot(np.nanpercentile(dnorm.cdf(ts[:,i]), np.arange(0,101)),np.arange(0,101)/100.0)
      plt.plot(np.arange(0,101)/100.0, np.arange(0,101)/100.0,'k--')
      plt.xlabel('Theoretical percentile')
      plt.ylabel('Actual percentile')
      plt.savefig("../WR/tplot_par_" + filename + ".png")
      plt.close()
  except:
      print("no plot")
  return(results)







def epf_par_out(truepar, dictd, dictv, n=75000, rlow=0, rhigh=500, offset=0, per = 1, getlow = -1, gethigh = -1, direct="", jchoice = -1):
  filename = dictd['filename']
  dataname = dictd['dataname']
  dataname2 = dictd['dataname2']
  wdataname = dictd['wdataname']
  vdataname = dictd['vdataname']


  if 'njac' in dictv:
    njac = dictv['njac']
  else:
    njac = njac1

  if 'njaca' in dictv:
    njaca = dictv['njaca']
  else:
    njaca = njac2

  if 'njacend' in dictv:
    njacend = dictv['njacend']
  else:
    njacend  = njac1_end

  if 'njacaend' in dictv:
    njacaend = dictv['njacaend']
  else:
    njacaend  = njac2_end

  if 'nimom' in dictv:
    nimom = dictv['nimom']
  else:
    nimom = nmom

  if 'nimom_end' in dictv:
    nimom_end = dictv['nimom_end']
  else:
    nimom_end = nmom_end

  if 'namom' in dictv:
    namom = dictv['namom']
  else:
    namom = nepfq

  if 'namom_end' in dictv:
    namom_end = dictv['namom_end']
  else:
    namom_end = nepfq_end

  if 'npar_est' in dictv:
    npar_est = dictv['npar_est']
  else:
    npar_est = 4
  if 'nresi' in dictv:
    nresi = dictv['nresi']
  else:
    nresi = nres

  if 'nresi_end' in dictv:
    nresi_end = dictv['nresi_end']
  else:
    nresi_end = nres_end

  if 'jacname' in dictd:
    jacname = dictd['jacname']
  else:
    jacname = 'jacobians_' + filename

  n_mom = nimom_end - nimom
  print(n_mom)
  n_mom2 = namom_end - namom

  Summest = []
  for i in range(int(rlow/per + offset/per),int(ceil(rhigh/per + offset/per))):
      try:
          Summest = Summest + [pd.read_csv("../Results/MC/Trials/" + direct + "summary_" + filename + "_" + str(i),header=None)]
      except:
          continue
  Summest = np.array(pd.concat(Summest))
  Mest = np.array(Summest[:,nresi:nresi_end]).astype(float)
  n_par = Mest.shape[1]
  print(n_par)
  print(filename)
  np.savetxt('../Results/Test/resul_' + filename + ".csv",Mest, delimiter=',')


# In[450]:

def fix_jac_old(dictf, dictv,n_par, upper=30, direct = "", direct2 = "", jchoice = 0, nhigh = 1000, nt = 75000, fix = 2):
    filev1 = dictf['vdataname']
    filev2 = dictf['filev2']
    filev12 = dictf['filev12']
    filew = dictf['wdataname']
    name = dictf['filename']
    name2 = name + '_jac_fixed_'
    dataname = dictf['dataname']
    dataname2 = dictf['dataname2']
    V1 = np.loadtxt("../Results/MC/" +  filev1, delimiter=",")
    V2 = np.loadtxt("../Results/MC/" + filev2, delimiter=",")
    V12 = np.loadtxt("../Results/MC/" + filev12, delimiter=",")
    m1 = np.loadtxt("../Results/MC/" + dataname, delimiter=",")
    m2 = np.loadtxt("../Results/MC/" + dataname2, delimiter=",")
    m1 = m1[0:nhigh,:]
    m2 = m2[0:nhigh,:]
    #v12 = np.cov(m1.transpose(), m2.transpose())
    #v12 = v12[:np.shape(m1)[1], np.shape(m1)[1]:] * nt
    if filew == 'identity':
      W = np.identity(4)
    else:
      W = np.loadtxt("../Results/MC/" + filew, delimiter=",")

    if 'njac' in dictv:
      njac = dictv['njac']
    else:
      njac = njac1

    if 'njaca' in dictv:
      njaca = dictv['njaca']
    else:
      njaca = njac2

    if 'njacend' in dictv:
      njacend = dictv['njacend']
    else:
      njacend  = njac1_end

    if 'njacaend' in dictv:
      njacaend = dictv['njacaend']
    else:
      njacaend  = njac2_end

    if 'nimom' in dictv:
      nimom = dictv['nimom']
    else:
      nimom = nmom

    if 'nimom_end' in dictv:
      nimom_end = dictv['nimom_end']
    else:
      nimom_end = nmom_end

    if 'namom' in dictv:
      namom = dictv['namom']
    else:
      namom = nepfq

    if 'namom_end' in dictv:
      namom_end = dictv['namom_end']
    else:
      namom_end = nepfq_end

    if 'transposer' in dictv:
      transposer = dictv['transposer']
    else:
      transposer = True

    if 'npar_est' in dictv:
      npar_est = dictv['npar_est']
    else:
      npar_est = 4

    if 'npar' in dictv:
      npar = dictv['npar']
    else:
      npar = 7

    if 'jacname' in dictf:
      jacname = dictf['jacname']
    else:
      jacname = 'jacobians_clean' + name + "_"



    n_mom2 = namom_end - namom
    n_mom  = nimom_end - nimom

    if transposer:
        nmomj1 = n_mom
        nmomj2 = n_mom2
        njacj1 = njac
        njacj2 = njaca
        njacjend1 = njacend
        njacjend2 = njacaend
    else:
        nmomj1 = n_mom2
        nmomj2 = n_mom
        njacj1 = njaca
        njacj2 = njac
        njacjend1 = njacaend
        njacjend2 = njacend

    njactot = n_mom * npar + n_mom2 * npar
    Jaccests = []
    for i in range(int(0),int(1000)):
        try:
            Jaccests = Jaccests + [pd.read_csv("../Results/MC/Trials/" + direct + "jacobians_" + jacname + "_" + str(i),header=None)]
        except:
            continue
    Parests = np.array(pd.read_csv("../Results/Test/resul_" + name + ".csv",header=None))
    Jaccests = np.array(pd.concat(Jaccests))
    print(np.shape(Jaccests))
    if fix ==1:
        jc1 = Jaccests[:,0:11*7].reshape(np.shape(Jaccests)[0],11,7)
        jc2 = Jaccests[:,11*7:11*7+18*7].reshape(np.shape(Jaccests)[0],18,7)
        jc1 = jc1[:,0:8,:]
        for i in range(0,1000):
            jc1[i,:,:] = jc1[i,:,:] / np.abs(pc[i,:])/np.abs(pc[i,:])/1e-6/1e-6
            jc2[i,:,:] = jc2[i,:,:] / np.abs(pc[i,:])/np.abs(pc[i,:])/1e-6/1e-6
    else:
        jc1 = Jaccests[:,jchoice*(njactot):jchoice*(njactot)+nmomj1*npar].reshape(np.shape(Jaccests)[0],nmomj1,npar).astype(float)
        jc2 = Jaccests[:,jchoice*(njactot)+nmomj1*npar:(jchoice+1)*(njactot)].reshape(np.shape(Jaccests)[0],nmomj2,npar).astype(float)
        #jc1[np.abs(jc1) < 1e-8] =0
        #jc2[np.abs(jc2) < 1e08] =0
    err = 0
    skip = 0
    starter = 0
    for i in range(0, upper):
        Summest = np.array(pd.read_csv("../Results/MC/Trials/" + direct + "summary_" + name + "_" + str(i),header=None))
        n = np.shape(Summest)[0]
        try:
            J1      = jc1[starter:starter+n,:,:].reshape(n, nmomj1*npar)
            J2      = jc2[starter:starter+n,:,:].reshape(n, nmomj2*npar)
            Summest[:,njacj1:njacjend1] = J1
            Summest[:,njacj2:njacjend2] = J2
        except:
            Summest[:,njacj1:njacjend1] = np.nan
            Summest[:,njacj2:njacjend2] = np.nan
        pd.DataFrame(Summest).to_csv('../Results/MC/Trials/' + direct + 'summary_' + name + "_fix_" + str(i), header=None, index=False)
        starter += n
    return [jc1, jc2]


# In[451]:
