#%%
import numpy as np
jd = 2460796.5 # May 1, 2025 00:00:00
time_yrs  = int(2e6)
tstep_yrs = int(2e2)
dt_yrs = 0.2
tyrs_str = '2e6yr'
tstepyrs_str = '2e2yr'
dtyrs_str = '0.2yr'
yrsstr = tyrs_str + '_' + tstepyrs_str + '_' + dtyrs_str + '_jd' + str(jd)
nt = int(time_yrs/tstep_yrs+1)
wdeg_mat_planets = np.loadtxt('b002_planets_perideg_' + yrsstr + '.csv',delimiter=',')
Wdeg_mat_planets = np.loadtxt('b002_planets_nodedeg_' + yrsstr + '.csv',delimiter=',')
Mdeg_mat_planets = np.loadtxt('b002_planets_Mdeg_' + yrsstr + '.csv',delimiter=',')
lambda_deg_mat_planets = Wdeg_mat_planets + wdeg_mat_planets + Mdeg_mat_planets
lambda_deg_jupiter = lambda_deg_mat_planets[0,:]
pomega_deg_mat_planets = Wdeg_mat_planets + wdeg_mat_planets
pomega_deg_jupiter = pomega_deg_mat_planets[0,:]
wdeg_mat_astorb = np.loadtxt('b003_astorb_perideg_' + yrsstr + '.csv',delimiter=',')
Wdeg_mat_astorb = np.loadtxt('b003_astorb_nodedeg_' + yrsstr + '.csv',delimiter=',')
Mdeg_mat_astorb = np.loadtxt('b003_astorb_Mdeg_' + yrsstr + '.csv',delimiter=',')
lambda_deg_mat_astorb = Wdeg_mat_astorb + wdeg_mat_astorb + Mdeg_mat_astorb
pomega_deg_mat_astorb = Wdeg_mat_astorb + wdeg_mat_astorb
nobj = wdeg_mat_astorb.shape[0]
phi32deg_mat_astorb = np.zeros((nobj,nt))
for iobj in range(nobj):
    for it in range(nt):
        lambda_prime = lambda_deg_jupiter[it]
        lambda_plain = lambda_deg_mat_astorb[iobj,it]
        pomega_plain = pomega_deg_mat_astorb[iobj,it]
        phi32deg = 3*lambda_prime - 2*lambda_plain - pomega_plain
        phi32deg = np.mod(phi32deg,360)
        if phi32deg > 180:
            phi32deg = -(360-phi32deg)
        phi32deg_mat_astorb[iobj,it] = phi32deg
np.savetxt('b005_32angledeg_astorb_' + yrsstr + '.csv',phi32deg_mat_astorb,delimiter=',')