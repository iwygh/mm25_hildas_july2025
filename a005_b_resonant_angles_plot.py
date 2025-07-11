#%%
import numpy as np
import pandas as pd
import time
from matplotlib import pyplot as plt
plt.ioff()
jd = 2460796.5 # May 1, 2025 00:00:00
time_yrs  = int(2e6)
tstep_yrs = int(2e2)
dt_yrs = 0.2
tyrs_str = '2e6yr'
tstepyrs_str = '2e2yr'
dtyrs_str = '0.2yr'
yrsstr = tyrs_str + '_' + tstepyrs_str + '_' + dtyrs_str + '_jd' + str(jd)
nt = int(time_yrs/tstep_yrs+1)
tyrsvec = np.loadtxt('b002_tyrsvec_'+yrsstr+'.csv',delimiter=',')
angle32deg_mat = np.loadtxt('b005_32angledeg_astorb_' + yrsstr + '.csv',delimiter=',')
df = pd.read_csv('b001_astorb_horizons.csv')
nobj = angle32deg_mat.shape[0]
t0 = time.time()
# nobj = 10
for iobj in range(nobj):
    idhere = df['idhere'][iobj]
    idhere = idhere.rstrip()
    fig = plt.figure()
    plt.plot(tyrsvec/1e6,angle32deg_mat[iobj,:])
    plt.plot([0,2],[0,0],color='red',linewidth=1)
    plt.plot([0,2],[-180,-180],color='red',linewidth=1)
    plt.plot([0,2],[180,180],color='red',linewidth=1)
    titlestr = 'b005_angle32deg_'+ str(iobj+1) + '_' + str(nobj) + '_' + idhere
    plt.title(titlestr)
    # plt.xlim([0,0.1])
    savestr = titlestr + '_' + yrsstr + '.pdf'
    plt.savefig(savestr)
    t1 = time.time()
    elapsed_time_seconds = t1 - t0
    elapsed_time_minutes = elapsed_time_seconds/60
    elapsed_time_hours = elapsed_time_minutes/60
    # print(iobj+1,nobj)
    print(iobj+1,nobj,np.round(elapsed_time_seconds,3),np.round(elapsed_time_minutes,3),np.round(elapsed_time_hours,3))
    plt.close(fig)