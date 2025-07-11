#%%
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
jd = 2460796.5 # May 1, 2025 00:00:00
time_yrs  = int(2e6)
tstep_yrs = int(2e2)
dt_yrs = 0.2
tyrs_str = '2e6yr'
tstepyrs_str = '2e2yr'
dtyrs_str = '0.2yr'
yrsstr = tyrs_str + '_' + tstepyrs_str + '_' + dtyrs_str + '_jd' + str(jd)
aau_mat_planets = np.loadtxt('b002_planets_aau_'+yrsstr+'.csv',delimiter=',')
tyrsvec = np.loadtxt('b002_tyrsvec_' + yrsstr + '.csv',delimiter=',')
nt = len(tyrsvec)
xinvar = np.loadtxt('b002_planets_xinvar_'+yrsstr+'.csv',delimiter=',')
yinvar = np.loadtxt('b002_planets_yinvar_'+yrsstr+'.csv',delimiter=',')
idegplanets = np.loadtxt('b002_planets_ideg_'+yrsstr+'.csv',delimiter=',')
Wdegplanets = np.loadtxt('b002_planets_nodedeg_'+yrsstr+'.csv',delimiter=',')
sini = np.sin(np.radians(idegplanets))
cosW = np.cos(np.radians(Wdegplanets))
sinW = np.sin(np.radians(Wdegplanets))
qplanets = sini * cosW
pplanets = sini * sinW
xJ = qplanets[0,:]
yJ = pplanets[0,:]
dfa = pd.read_csv('b006_astorb_labels.csv')
n_astorb = dfa.shape[0]
aau_astorb = np.loadtxt('b003_astorb_aau_'+yrsstr+'.csv',delimiter=',')
e_astorb = np.loadtxt('b003_astorb_e_'+yrsstr+'.csv',delimiter=',')
ideg_astorb = np.loadtxt('b003_astorb_ideg_'+yrsstr+'.csv',delimiter=',')
Wdeg_astorb = np.loadtxt('b003_astorb_nodedeg_'+yrsstr+'.csv',delimiter=',')
sini = np.sin(np.radians(ideg_astorb))
cosW = np.cos(np.radians(Wdeg_astorb))
sinW = np.sin(np.radians(Wdeg_astorb))
q_astorb = sini*cosW
p_astorb = sini*sinW
output_labels = ['P_L','H_L','S_L','BFG_LMC','BPHSFG_LMC']
input_conditions = [[['Potomac'],['librating']],\
                    [['Hilda'],['librating']],\
                    [['Schubart'],['librating']],\
                    [['Background','Francette','Guinevere'],['librating','maybe','circulating']],\
                    [['Background','Francette','Guinevere','Potomac','Hilda','Schubart'],['librating','maybe','circulating']],\
                    ] 
clone_labels = ['Potomac','Hilda','Schubart','Ismene']
fig = plt.figure(figsize=(5,8)) # (width,height)
nrows = 5
ncols = 2
ax = [fig.add_subplot(nrows,ncols,i+1) for i in range(nrows*ncols)]
plt.rcParams.update({'font.size':6})
for iax in [0,1,2,3,4,5,6,7,8]:
    a = ax[iax]
    # a.set_xlabel('Time, Myr')
    a.tick_params(axis='x',direction='in')
    a.tick_params(axis='y',direction='in')
fig.subplots_adjust(wspace=0.25, hspace=0.1)
# subplot titles go across each row, then on to the next row
subplottitles = ['Potomacs Δ (deg)',   'Potomac clones Δ (deg)',\
                 'Hildas Δ (deg)',     'Hilda clones Δ (deg)',\
                 'Schubarts Δ (deg)',  'Schubart clones Δ (deg)',\
                 'Backgrounds Δ (deg)','Ismene clones Δ (deg)',\
                 'All Δ (deg)']
for iax in [0,1,2,3,4,5,6,7,8]:
    ax[iax].text(0.06,0.9,subplottitles[iax],horizontalalignment='left',verticalalignment='top',\
               transform=ax[iax].transAxes,bbox=dict(facecolor='white',alpha=1))
    ax[iax].grid('on')
for iax in [0,2,4,6,8]:
    averaging_window = 150
    output_label = output_labels[int(iax/2)]
    input_condition = input_conditions[int(iax/2)]
    dfvmf = pd.read_csv('b007_'+output_label+'_statistics_vmf_'+yrsstr+'.csv')
    xcc = np.array(dfvmf['xcc'].to_list())
    ycc = np.array(dfvmf['ycc'].to_list())
    dfl_icon = pd.read_csv('b010_laplaceplane_'+output_label+'_'+yrsstr+'.csv')
    dflx = dfl_icon['laplacex'].to_list()
    dfly = dfl_icon['laplacey'].to_list()
    deltadeg_laplace_here = np.degrees(np.arcsin(np.sqrt((xcc-dflx)**2+(ycc-dfly)**2)))
    deltadeg_jupiter_here = np.degrees(np.arcsin(np.sqrt((xcc-xJ)**2+(ycc-yJ)**2)))
    deltadeg_invariable_here = np.degrees(np.arcsin(np.sqrt((xcc-xinvar)**2+(ycc-yinvar)**2)))
    deltadeg_laplace_here_av = []
    deltadeg_jupiter_here_av = []
    deltadeg_invariable_here_av = []
    for it in range(nt-averaging_window):
        ddlh = np.mean(deltadeg_laplace_here[it:it+averaging_window])
        ddjh = np.mean(deltadeg_jupiter_here[it:it+averaging_window])
        ddih = np.mean(deltadeg_invariable_here[it:it+averaging_window])
        deltadeg_laplace_here_av.append(ddlh)
        deltadeg_jupiter_here_av.append(ddjh)
        deltadeg_invariable_here_av.append(ddih)
    # ax[iax].semilogy(tyrsvec/1e6,deltadeg_laplace_here,'r-',lw=0.5)
    # ax[iax].semilogy(tyrsvec/1e6,deltadeg_jupiter_here,'b-',lw=0.5)
    # ax[iax].semilogy(tyrsvec/1e6,deltadeg_invariable_here,'k-',lw=0.5)
    ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,deltadeg_laplace_here_av,'r-',lw=0.25)
    ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,deltadeg_jupiter_here_av,'b-',lw=0.25)
    ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,deltadeg_invariable_here_av,'k-',lw=0.25)
    # ax[iax].set_ylim([0.01,10])
    # ax[iax].set_yticks([0.01,0.1,1])
    ax[iax].grid('on')
    ax[iax].set_xlim([0,2])
    # ax[iax].set_xticks([0,0.5,1,1.5])
    # ax[iax].set_xticklabels([])
for iax in [1,3,5,7]:
    averaging_window = 1
    clone_label = clone_labels[int((iax-1)/2)]
    output_label = output_labels[int((iax-1)/2)]
    input_condition = input_conditions[int((iax-1)/2)]
    dfcf = pd.read_csv('b008_'+clone_label+'_clones_statistics_cf_'+yrsstr+'.csv')
    xcc = np.array(dfcf['xcc'].to_list())
    ycc = np.array(dfcf['ycc'].to_list())
    dfl_icon = pd.read_csv('b010_laplaceplane_'+output_label+'_'+yrsstr+'.csv')
    dflx = dfl_icon['laplacex'].to_list()
    dfly = dfl_icon['laplacey'].to_list()
    deltadeg_laplace_here = np.degrees(np.arcsin(np.sqrt((xcc-dflx)**2+(ycc-dfly)**2)))
    deltadeg_jupiter_here = np.degrees(np.arcsin(np.sqrt((xcc-xJ)**2+(ycc-yJ)**2)))
    deltadeg_invariable_here = np.degrees(np.arcsin(np.sqrt((xcc-xinvar)**2+(ycc-yinvar)**2)))
    # deltadeg_laplace_here_av = []
    # deltadeg_jupiter_here_av = []
    # deltadeg_invariable_here_av = []
    # for it in range(nt-averaging_window):
    #     ddlh = np.mean(deltadeg_laplace_here[it:it+averaging_window])
    #     ddjh = np.mean(deltadeg_jupiter_here[it:it+averaging_window])
    #     ddih = np.mean(deltadeg_invariable_here[it:it+averaging_window])
    #     deltadeg_laplace_here_av.append(ddlh)
    #     deltadeg_jupiter_here_av.append(ddjh)
    #     deltadeg_invariable_here_av.append(ddih)
    # ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,deltadeg_laplace_here_av,'r-',lw=0.25)
    # ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,deltadeg_jupiter_here_av,'b-',lw=0.25)
    # ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,deltadeg_invariable_here_av,'k-',lw=0.25)
    ax[iax].plot(tyrsvec/1e6,deltadeg_laplace_here,'r-',lw=0.25)
    ax[iax].plot(tyrsvec/1e6,deltadeg_jupiter_here,'b-',lw=0.25)
    ax[iax].plot(tyrsvec/1e6,deltadeg_invariable_here,'k-',lw=0.25)
    ax[iax].grid('on')
    ax[iax].set_xlim([0,2])
ax[0].set_ylim([0,2])
ax[2].set_ylim([0,1])
ax[4].set_ylim([0,1])
ax[6].set_ylim([0,1])
ax[8].set_ylim([0,1])
ax[1].set_ylim([0,1])
ax[3].set_ylim([0,1])
ax[5].set_ylim([0,1])
ax[7].set_ylim([0,1])
ax[9].set_xticks([0,0.5,1,1.5])
ax[1].set_yticks([0,0.25,0.5,0.75])
ax[3].set_yticks([0,0.25,0.5,0.75])
ax[5].set_yticks([0,0.25,0.5,0.75])
ax[7].set_yticks([0,0.25,0.5,0.75])
ax[9].set_yticks([])
ax[0].set_yticks([0,0.5,1,1.5])
ax[2].set_yticks([0,0.25,0.5,0.75])
ax[4].set_yticks([0,0.25,0.5,0.75])
ax[6].set_yticks([0,0.25,0.5,0.75])
ax[8].set_yticks([0,0.25,0.5,0.75])
for iax in range(8):
    ax[iax].set_xticklabels([])
ax[8].set_xlabel('Time, Myr')
ax[9].set_xlabel('Time, Myr')
outfile_plot = 'b013_deltadegplots_vmfobj_cfclones.png'
plt.savefig(outfile_plot,dpi=400)