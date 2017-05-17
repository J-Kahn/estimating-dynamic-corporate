import numpy as np
import sys

import matplotlib.pyplot as plt

auxs = ['mom', 'quad', 'mom_cv', 'quad_trans', 'cube']
prefixes = sys.argv[1:]
folder = "../Results/MC/"
#prefixes = ["data_rr_small_","data_rr_","data_rr_fixed_","data_rr_mis_","data_rr_drs_"]

for prefix in prefixes:
    print(prefix)
    for aux in auxs:
      try:
        mom    = np.loadtxt(folder + "trial_" + prefix + "" + aux + ".csv",delimiter=",")
        vmom   = np.loadtxt(folder + "trial_"  + prefix + "v" + aux + ".csv",delimiter=",")
        #vmom2  = np.loadtxt(folder + "trial_"  + prefix + "v" + aux + "_unclust.csv",delimiter=",")
  
        print(aux + " ... " )
  
        jstats = np.zeros(mom.shape[0])
  
        rorder = np.random.permutation(range(0, mom.shape[0]))[0:10]
        wmom = np.zeros(vmom.shape)
        wmom2 = np.zeros(vmom.shape)
        wmom3 = np.zeros(vmom.shape)
        wmom4 = np.zeros(vmom.shape)
  
        m = mom[rorder[0:10],:].mean(0)
  
        for i in range(0, mom.shape[0]):
            w = np.linalg.inv(vmom[i, :].reshape(mom.shape[1],mom.shape[1]))
            w3 = np.linalg.inv(np.diag(np.diag(vmom[i, :].reshape(mom.shape[1],mom.shape[1]))))
            #w4 = np.linalg.inv(np.diag(np.diag(vmom2[i, :].reshape(mom.shape[1],mom.shape[1]))) )
            #w2 = np.linalg.inv(vmom2[i, :].reshape(mom.shape[1],mom.shape[1]))
            wmom[i, :] = w.reshape(mom.shape[1]*mom.shape[1])
            wmom3[i, :] = w3.reshape(mom.shape[1]*mom.shape[1])
            #wmom4[i, :] = w4.reshape(mom.shape[1]*mom.shape[1])
            #wmom2[i, :] = w2.reshape(mom.shape[1]*mom.shape[1])
            diff = mom[i, :] - m
            jstats[i] = np.dot(np.dot(diff.transpose(), w), diff)*10/11*71880
  
        np.savetxt(folder + "trial_" + prefix + "w" + aux + ".csv", wmom, delimiter = ",")
        np.savetxt(folder + "trial_" + prefix + "w" + aux + "_unclust.csv", wmom2, delimiter = ",")
        np.savetxt(folder + "trial_" + prefix + "w" + aux + "_diag.csv", wmom3, delimiter = ",")
        np.savetxt(folder + "trial_" + prefix + "w" + aux + "_unclust_diag.csv", wmom4, delimiter = ",")
      except:
        print('Error! ' + prefix + " " + aux)
