import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import astropy.stats.circstats as acstats
jd = 2460796.5 # May 1, 2025 00:00:00
time_yrs  = int(2e6)
tstep_yrs = int(2e2)
dt_yrs = 0.2
tyrs_str = '2e6yr'
tstepyrs_str = '2e2yr'
dtyrs_str = '0.2yr'
yrsstr = tyrs_str + '_' + tstepyrs_str + '_' + dtyrs_str + '_jd' + str(jd)
output_labels = ['P_L',      'H_L',       'S_L',\
                 'B_L',      'B_LM',      'B_LMC',\
                 'BFG_L',    'BFG_LM',    'BFG_LMC',\
                 'BPHS_L',   'BPHS_LM',   'BPHS_LMC',\
                 'BPHSFG_L', 'BPHSFG_LM', 'BPHSFG_LMC']
nlabels = len(output_labels)
input_conditions = [[['Potomac'],['librating']],\
    [['Hilda'],['librating']],\
    [['Schubart'],['librating']],\
    [['Background'],['librating']],[['Background'],['librating','maybe']],[['Background'],['librating','maybe','circulating']],\
    [['Background','Francette','Guinevere'],['librating']],[['Background','Francette','Guinevere'],['librating','maybe']],[['Background','Francette','Guinevere'],['librating','maybe','circulating']],\
    [['Background','Potomac','Hilda','Schubart'],['librating']],[['Background','Potomac','Hilda','Schubart'],['librating','maybe']],[['Background','Potomac','Hilda','Schubart'],['librating','maybe','circulating']],\
    [['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating']],[['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating','maybe']],[['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating','maybe','circulating']],\
    ] 
distinct_groups = ['P_L',     'H_L',     'S_L',\
                   'B_L',     'B_M',     'B_C',\
                   'BFG_L',   'BFG_M',   'BFG_C',\
                   'BPHS_L',  'BPHS_M',  'BPHS_C',\
                   'BPHSFG_L','BPHSFG_M','BPHSFG_C']    
conditions_for_distinct_groups = [[['Potomac'],['librating']],\
    [['Hilda'],['librating']],\
    [['Schubart'],['librating']],\
    [['Background'],['librating']],[['Background'],['maybe']],[['Background'],['circulating']],\
    [['Background','Francette','Guinevere'],['librating']],[['Background','Francette','Guinevere'],['maybe']],[['Background','Francette','Guinevere'],['circulating']],\
    [['Background','Potomac','Hilda','Schubart'],['librating']],[['Background','Potomac','Hilda','Schubart'],['maybe']],[['Background','Potomac','Hilda','Schubart'],['circulating']],\
    [['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating']],[['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['maybe']],[['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['circulating']],\
    ] 
groups_in_conditions = [['P_L'],     ['H_L'],     ['S_L'],\
     ['B_L'],     ['B_L','B_M'],          ['B_L','B_M','B_C'],\
     ['BFG_L'],   ['BFG_L','BFG_M'],      ['BFG_L','BFG_M','BFG_C'],\
     ['BPHS_L'],  ['BPHS_L','BPHS_M'],    ['BPHS_L','BPHS_M','BPHS_C'],\
     ['BPHSFG_L'],['BPHSFG_L','BPHSFG_M'],['BPHSFG_L','BPHSFG_M','BPHSFG_C']] 
aau_mat_planets = np.loadtxt('b002_planets_aau_'+yrsstr+'.csv',delimiter=',')
xinvar = np.loadtxt('b002_planets_xinvar_'+yrsstr+'.csv',delimiter=',')
yinvar = np.loadtxt('b002_planets_yinvar_'+yrsstr+'.csv',delimiter=',')
it = 1
xI = xinvar[it]
yI = yinvar[it]
idegplanets = np.loadtxt('b002_planets_ideg_'+yrsstr+'.csv',delimiter=',')
Wdegplanets = np.loadtxt('b002_planets_nodedeg_'+yrsstr+'.csv',delimiter=',')
sini = np.sin(np.radians(idegplanets))
cosW = np.cos(np.radians(Wdegplanets))
sinW = np.sin(np.radians(Wdegplanets))
qplanets = sini * cosW
pplanets = sini * sinW
xJ = qplanets[0,it]
yJ = pplanets[0,it]
dfa = pd.read_csv('b006_astorb_labels.csv')
n_astorb = dfa.shape[0]
Hmax = 15.7
aau_astorb = np.loadtxt('b003_astorb_aau_'+yrsstr+'.csv',delimiter=',')
e_astorb = np.loadtxt('b003_astorb_e_'+yrsstr+'.csv',delimiter=',')
ideg_astorb = np.loadtxt('b003_astorb_ideg_'+yrsstr+'.csv',delimiter=',')
Wdeg_astorb = np.loadtxt('b003_astorb_nodedeg_'+yrsstr+'.csv',delimiter=',')
sini = np.sin(np.radians(ideg_astorb))
cosW = np.cos(np.radians(Wdeg_astorb))
sinW = np.sin(np.radians(Wdeg_astorb))
q_astorb = sini*cosW
p_astorb = sini*sinW
for ilabel in range(nlabels):
    label = output_labels[ilabel]
    df = pd.read_csv('b007_'+label+'_statistics_vmf_'+yrsstr+'.csv')
    xcc = df['xcc'][it]
    ycc = df['ycc'][it]
    df = pd.read_csv('b010_laplaceplane_'+label+'_'+yrsstr+'.csv')
    xl = df['laplacex'][it]
    yl = df['laplacey'][it]
    input_groups = input_conditions[ilabel][0]
    input_labels = input_conditions[ilabel][1]
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) and (dfa['label'][i_astorb] in input_labels):
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    print(label,nobj)
    it = 0
    aau_icon = aau_astorb[astorb_indices,it]
    e_icon = e_astorb[astorb_indices,it]
    ideg_icon = ideg_astorb[astorb_indices,it]
    q_icon = q_astorb[astorb_indices,it]
    p_icon = p_astorb[astorb_indices,it]
    q = q_icon
    p = p_icon
    dx = q - xl
    dy = p - yl
    dxsum = np.sum(dx)
    dysum = np.sum(dy)
    print('dxsum_'+label+'_wrt_laplace = ',dxsum)
    print('dysum_'+label+'_wrt_laplace = ',dysum)
    dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
    dsini = np.sin(np.radians(dideg))
    dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
    plt.figure()
    plt.hist(dWdeg)
    title = 'b018_Wfree_deg_'+label+'_wrt_laplace'
    plt.title(title)
    plt.savefig(title+'.png')
    plt.show()
    pval = acstats.rayleightest(np.radians(dWdeg))
    print('p_W_'+label+'_laplace = ',pval)
    dx = q - xJ
    dy = p - yJ
    dxsum = np.sum(dx)
    dysum = np.sum(dy)
    print('dxsum_'+label+'_wrt_xJ = ',dxsum)
    print('dysum_'+label+'_wrt_yJ = ',dysum)
    dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
    dsini = np.sin(np.radians(dideg))
    dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
    plt.figure()
    plt.hist(dWdeg)
    title = 'b018_Wfree_deg_'+label+'_wrt_xJ'
    plt.title(title)
    plt.savefig(title+'.png')
    plt.show()
    pval = acstats.rayleightest(np.radians(dWdeg))
    print('p_W_'+label+'_xJ = ',pval)
    dx = q - xcc
    dy = p - ycc
    dxsum = np.sum(dx)
    dysum = np.sum(dy)
    print('dxsum_'+label+'_wrt_xcc = ',dxsum)
    print('dysum_'+label+'_wrt_ycc = ',dysum)
    dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
    dsini = np.sin(np.radians(dideg))
    dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
    plt.figure()
    plt.hist(dWdeg)
    title = 'b018_Wfree_deg_'+label+'_wrt_xcc'
    plt.title(title)
    plt.savefig(title+'.png')
    plt.show()
    pval = acstats.rayleightest(np.radians(dWdeg))
    print('p_W_'+label+'_xcc = ',pval)
    dx = q - xI
    dy = p - yI
    dxsum = np.sum(dx)
    dysum = np.sum(dy)
    print('dxsum_'+label+'_wrt_xI = ',dxsum)
    print('dysum_'+label+'_wrt_xI = ',dysum)
    dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
    dsini = np.sin(np.radians(dideg))
    dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
    plt.figure()
    plt.hist(dWdeg)
    title = 'b018_Wfree_deg_'+label+'_wrt_xI'
    plt.title(title)
    plt.savefig(title+'.png')
    plt.show()
    pval = acstats.rayleightest(np.radians(dWdeg))
    print('p_W_'+label+'_xI = ',pval)
    print('')
    print('')


















# import numpy as np
# import pandas as pd
# from matplotlib import pyplot as plt
# import astropy.stats.circstats as acstats
# #%%
# for famstr in['Potomac','Hilda','Schubart','background','astorb','bphs']:
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_laplacexvec_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_laplaceyvec_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
#     dsini = np.sin(np.radians(dideg))
#     # plt.figure()
#     # plt.hist(dsini)
#     # plt.title('sin(ifree) '+famstr+'_wrt_laplace (librating)')
#     # plt.show()
#     dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
#     plt.figure()
#     plt.hist(dWdeg)
#     title = 'Wfree(deg) '+famstr+'_wrt_laplace (librating)'
#     plt.title(title)
#     plt.savefig(title+'.png')
#     plt.show()
#     pval = acstats.rayleightest(np.radians(dWdeg))
#     print('p_W_'+famstr+'_laplace_librating = ',pval)
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_laplacexvec_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_laplaceyvec_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
#     dsini = np.sin(np.radians(dideg))
#     # plt.figure()
#     # plt.hist(dsini)
#     # plt.title('sin(ifree) '+famstr+'_wrt_laplace (wc)')
#     # plt.show()
#     dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
#     plt.figure()
#     plt.hist(dWdeg)
#     title = 'Wfree(deg) '+famstr+'_wrt_laplace (wc)'
#     plt.title(title)
#     plt.savefig(title+'.png')
#     plt.show()
#     pval = acstats.rayleightest(np.radians(dWdeg))
#     print('p_W_'+famstr+'_laplace = ',pval)
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_xccvec_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_yccvec_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
#     dsini = np.sin(np.radians(dideg))
#     # plt.figure()
#     # plt.hist(dsini)
#     # plt.title('sin(ifree) '+famstr+'_wrt_xcc (librating)')
#     # plt.show()
#     dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
#     plt.figure()
#     plt.hist(dWdeg)
#     title = 'Wfree(deg) '+famstr+'_wrt_xcc (librating)'
#     plt.title(title)
#     plt.savefig(title+'.png')
#     plt.show()
#     pval = acstats.rayleightest(np.radians(dWdeg))
#     print('p_W_'+famstr+'_xcc_librating = ',pval)
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_xccvec_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_yccvec_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
#     dsini = np.sin(np.radians(dideg))
#     # plt.figure()
#     # plt.hist(dsini)
#     # plt.title('sin(ifree) '+famstr+'_wrt_xcc (wc)')
#     # plt.show()
#     dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
#     plt.figure()
#     plt.hist(dWdeg)
#     title = 'Wfree(deg) '+famstr+'_wrt_xcc (wc)'
#     plt.title(title)
#     plt.savefig(title+'.png')
#     plt.show()
#     pval = acstats.rayleightest(np.radians(dWdeg))
#     print('p_W_'+famstr+'_xcc = ',pval)
#     print('')
# #%%
# for famstr in['Potomac','Hilda','Schubart','background','astorb','bphs']:
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     print('nobj = ',idegmat.shape[0])
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_xccvec_Hlt15.7_librating_arithmeticaverage_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_yccvec_Hlt15.7_librating_arithmeticaverage_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
#     dsini = np.sin(np.radians(dideg))
#     # plt.figure()
#     # plt.hist(dsini)
#     # plt.title('sin(ifree) '+famstr+'_wrt_xcc (librating) (ar_avg)')
#     # plt.show()
#     dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
#     plt.figure()
#     plt.hist(dWdeg)
#     title = 'Wfree(deg) '+famstr+'_wrt_xcc (librating) (ar_avg)'
#     plt.title(title)
#     plt.savefig(title+'.png')
#     plt.show()
#     pval = acstats.rayleightest(np.radians(dWdeg))
#     print('p_W_'+famstr+'_xcc_librating_arithmeticaverage = ',pval)
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     print('nobj (wc) = ',idegmat.shape[0])
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_xccvec_Hlt15.7_librating_arithmeticaverage_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_yccvec_Hlt15.7_librating_arithmeticaverage_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
#     dsini = np.sin(np.radians(dideg))
#     # plt.figure()
#     # plt.hist(dsini)
#     # plt.title('sin(ifree) '+famstr+'_wrt_xcc (wc) (ar_avg)')
#     # plt.show()
#     dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
#     plt.figure()
#     plt.hist(dWdeg)
#     title = 'Wfree(deg) '+famstr+'_wrt_xcc (wc) (ar_avg)'
#     plt.title(title)
#     plt.savefig(title+'.png')
#     plt.show()
#     pval = acstats.rayleightest(np.radians(dWdeg))
#     print('p_W_'+famstr+'_xcc_arithmeticaverage = ',pval)
#     print('')
# #%%
# for famstr in['background','astorb','bphs']:
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_xccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_yccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
#     dsini = np.sin(np.radians(dideg))
#     # plt.figure()
#     # plt.hist(dsini)
#     # plt.title('sin(ifree) '+famstr+'_wrt_xcc (librating) (vmf)')
#     # plt.show()
#     dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
#     plt.figure()
#     plt.hist(dWdeg)
#     title = 'Wfree(deg) '+famstr+'_wrt_xcc (librating) (vmf)'
#     plt.title(title)
#     plt.savefig(title+'.png')
#     plt.show()
#     pval = acstats.rayleightest(np.radians(dWdeg))
#     print('p_W_'+famstr+'_xcc_librating_vmf = ',pval)
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_xccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_yccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
#     dsini = np.sin(np.radians(dideg))
#     # plt.figure()
#     # plt.hist(dsini,bins=60)
#     # plt.title('sin(ifree) '+famstr+'_wrt_xcc (wc) (vmf)')
#     # plt.show()
#     dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
#     plt.figure()
#     plt.hist(dWdeg)
#     title = 'Wfree(deg) '+famstr+'_wrt_xcc (wc) (vmf)'
#     plt.title(title)
#     plt.savefig(title+'.png')
#     plt.show()
#     pval = acstats.rayleightest(np.radians(dWdeg))
#     print('p_W_'+famstr+'_xcc_vmf = ',pval)
#     print('')
# #%%
# for famstr in['Potomac','Hilda','Schubart','background','astorb','bphs']:
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_xccvec_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_yccvec_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dxsum = np.sum(dx)
#     dysum = np.sum(dy)
#     print('dxsum_'+famstr+'_xcc_librating = ',dxsum)
#     print('dysum_'+famstr+'_xcc_librating = ',dysum)
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_xccvec_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_yccvec_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dxsum = np.sum(dx)
#     dysum = np.sum(dy)
#     print('dxsum_'+famstr+'_xcc (wc) = ',dxsum)
#     print('dysum_'+famstr+'_xcc (wc) = ',dysum)
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_xccvec_Hlt15.7_librating_arithmeticaverage_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_yccvec_Hlt15.7_librating_arithmeticaverage_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dxsum = np.sum(dx)
#     dysum = np.sum(dy)
#     print('dxsum_'+famstr+'_xcc_librating_arithmeticaverage = ',dxsum)
#     print('dysum_'+famstr+'_xcc_librating_arithmeticaverage = ',dysum)
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_xccvec_Hlt15.7_arithmeticaverage_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_yccvec_Hlt15.7_arithmeticaverage_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dxsum = np.sum(dx)
#     dysum = np.sum(dy)
#     print('dxsum_'+famstr+'_xcc_arithmeticaverage (wc) = ',dxsum)
#     print('dysum_'+famstr+'_xcc_arithmeticaverage (wc) = ',dysum)
#     print('')
# #%%
# for famstr in['background','astorb','bphs']:
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_librating_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_xccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_yccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dxsum = np.sum(dx)
#     dysum = np.sum(dy)
#     print('dxsum_'+famstr+'_xcc_librating = ',dxsum)
#     print('dysum_'+famstr+'_xcc_librating = ',dysum)
#     idegmat = np.loadtxt(famstr+'_ideg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ideg = idegmat[:,0]
#     Wdegmat = np.loadtxt(famstr+'_nodedeg_Hlt15.7_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     Wdeg = Wdegmat[:,0]
#     irad = np.radians(ideg)
#     Wrad = np.radians(Wdeg)
#     q = np.sin(irad)*np.cos(Wrad)
#     p = np.sin(irad)*np.sin(Wrad)
#     xlvec = np.loadtxt(famstr+'_xccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     ylvec = np.loadtxt(famstr+'_yccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
#     xl = xlvec[0]
#     yl = ylvec[0]
#     dx = q - xl
#     dy = p - yl
#     dxsum = np.sum(dx)
#     dysum = np.sum(dy)
#     print('dxsum_'+famstr+'_xcc (wc) (vmf) = ',dxsum)
#     print('dysum_'+famstr+'_xcc (wc) (vmf) = ',dysum)
#     print('')
    