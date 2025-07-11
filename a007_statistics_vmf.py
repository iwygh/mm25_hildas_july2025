#%%
def vmf_fun(xvec,yvec,zvec,probmass):
    import numpy as np
    nobj = len(xvec)
    Sx = np.sum(xvec)
    Sy = np.sum(yvec)
    Sz = np.sum(zvec)
    R = np.sqrt(Sx**2+Sy**2+Sz**2)
    xcc = Sx/R
    ycc = Sy/R
    zcc = Sz/R
    Rbar = R/nobj
    dsum = 0
    for iobj in range(nobj):
        dsum = dsum + (xvec[iobj]*xcc+yvec[iobj]*ycc+zvec[iobj]*zcc)**2
    d = 1 - 1/nobj * dsum
    sigmahat = np.sqrt(d/(nobj*Rbar**2))
    A = 1 - probmass
    angledeg = np.degrees(np.arcsin(sigmahat*np.sqrt(-np.log(A))))
    return xcc,ycc,zcc,angledeg
#%%
import numpy as np
import pandas as pd
jd = 2460796.5 # May 1, 2025 00:00:00
time_yrs  = int(2e6)
tstep_yrs = int(2e2)
dt_yrs = 0.2
tyrs_str = '2e6yr'
tstepyrs_str = '2e2yr'
dtyrs_str = '0.2yr'
yrsstr = tyrs_str + '_' + tstepyrs_str + '_' + dtyrs_str + '_jd' + str(jd)
tyrsvec = np.loadtxt('b002_tyrsvec_' + yrsstr + '.csv',delimiter=',')
nt = len(tyrsvec)
dfa = pd.read_csv('b006_astorb_labels.csv')
n_astorb = dfa.shape[0]
probmasses = [0.68,0.95,0.997]
Hmax = 15.7
ideg_mat = np.loadtxt('b003_astorb_ideg_'+yrsstr+'.csv',delimiter=',')
Wdeg_mat = np.loadtxt('b003_astorb_nodedeg_'+yrsstr+'.csv',delimiter=',')
# 7 rows of matrix are:
    # xcc,ycc,icc,Wcc,a1,a2,a3
    # all angles in degrees
    # a1,a2,a3 - vmf 68,95,99.7% uncertainty angles for xcc,ycc
    # no need to have separate uncertainty angles for icc,Wcc because those are
    # equal to uncertainty angles for mean pole
output_labels = ['P_L',      'P_LM',      'P_LMC',\
                 'H_L',      'H_LM',      'H_LMC',\
                 'S_L',      'S_LM',      'S_LMC',\
                 'B_L',      'B_LM',      'B_LMC',\
                 'BFG_L',    'BFG_LM',    'BFG_LMC',\
                 'BPHS_L',   'BPHS_LM',   'BPHS_LMC',\
                 'BPHSFG_L', 'BPHSFG_LM', 'BPHSFG_LMC']
input_conditions = [[['Potomac'],['librating']],[['Potomac'],['librating','maybe']],[['Potomac'],['librating','maybe','circulating']],\
    [['Hilda'],['librating']],[['Hilda'],['librating','maybe']],[['Hilda'],['librating','maybe','circulating']],\
    [['Schubart'],['librating']],[['Schubart'],['librating','maybe']],[['Schubart'],['librating','maybe','circulating']],\
    [['Background'],['librating']],[['Background'],['librating','maybe']],[['Background'],['librating','maybe','circulating']],\
    [['Background','Francette','Guinevere'],['librating']],[['Background','Francette','Guinevere'],['librating','maybe']],[['Background','Francette','Guinevere'],['librating','maybe','circulating']],\
    [['Background','Potomac','Hilda','Schubart'],['librating']],[['Background','Potomac','Hilda','Schubart'],['librating','maybe']],[['Background','Potomac','Hilda','Schubart'],['librating','maybe','circulating']],\
    [['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating']],[['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating','maybe']],[['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating','maybe','circulating']],\
    ]               
ncons = len(input_conditions)
for icon in range(ncons):
    output_label = output_labels[icon]
    input_groups = input_conditions[icon][0]
    input_labels = input_conditions[icon][1]
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) and (dfa['label'][i_astorb] in input_labels):
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    print(output_label,nobj)
    outmat = np.zeros((13,nt))
    irad_mat = np.radians(ideg_mat[astorb_indices,:])
    Wrad_mat = np.radians(Wdeg_mat[astorb_indices,:])
    xmat = np.sin(irad_mat)*np.cos(Wrad_mat)
    ymat = np.sin(irad_mat)*np.sin(Wrad_mat)
    zmat = np.cos(irad_mat)
    for it in range(nt):
        xvec = xmat[:,it]
        yvec = ymat[:,it]
        zvec = zmat[:,it]
        xcc,ycc,zcc,angle1deg = vmf_fun(xvec,yvec,zvec,probmasses[0])
        xcc,ycc,zcc,angle2deg = vmf_fun(xvec,yvec,zvec,probmasses[1])
        xcc,ycc,zcc,angle3deg = vmf_fun(xvec,yvec,zvec,probmasses[2])
        outmat[0,it] = xcc
        outmat[1,it] = ycc
        outmat[2,it] = angle1deg
        outmat[3,it] = angle2deg
        outmat[4,it] = angle3deg
    dictionary = {'xcc':outmat[0,:],\
                  'ycc':outmat[1,:],\
                  'angle_68_deg':outmat[2,:],\
                  'angle_95_deg':outmat[3,:],\
                  'angle_997_deg':outmat[4,:]}
    dfd = pd.DataFrame.from_dict(dictionary)
    dfd.to_csv('b007_'+output_label+'_statistics_vmf_'+yrsstr+'.csv',index=None)
#%%
# P_L 200
# P_LM 200
# P_LMC 200
# H_L 485
# H_LM 485
# H_LMC 485
# S_L 439
# S_LM 439
# S_LMC 439
# B_L 1056
# B_LM 1085
# B_LMC 1118
# BFG_L 1071
# BFG_LM 1100
# BFG_LMC 1133
# BPHS_L 2180
# BPHS_LM 2209
# BPHS_LMC 2242
# BPHSFG_L 2195
# BPHSFG_LM 2224
# BPHSFG_LMC 2257