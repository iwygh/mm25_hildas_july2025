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
for iax in range(nrows*ncols-1):
    a = ax[iax]
    a.tick_params(axis='x',direction='in')
    a.tick_params(axis='y',direction='in')
    a.set_xticks([0,0.5,1,1.5])
    a.set_xlim([0,2])
fig.subplots_adjust(wspace=0.05, hspace=0.05)
# subplot titles go across each row, then on to the next row
# subplottitles = ['Potomacs','Hildas','Schubarts','Backgrounds','All',\
                 # 'Potomac clones','Hilda clones','Schubart clones','Ismene clones']
subplottitles = ['Potomacs','Potomac clones',\
                 'Hildas','Hilda clones',\
                 'Schubarts','Schubart clones',\
                 'Backgrounds','Ismene clones',\
                 'All']
# for iax in [0,1,2,3,4,  5,6,7,8,9,  10,11,12,13,  15,16,17,18]:
for iax in [0,1,  2,3,  4,5,  6,7,  8]:
    ax[iax].text(0.06,0.9,subplottitles[iax],horizontalalignment='left',verticalalignment='top',\
               transform=ax[iax].transAxes,bbox=dict(facecolor='white',alpha=1))
    ax[iax].grid('on')
for iax in [8,9]:
    ax[iax].set_xlabel('Time, Myr')
    ax[iax].set_xlim([0,2])
    ax[iax].set_xticks([0,0.5,1,1.5])
for iax in [0,1,  2,3,  4,5,  6,7]:
    ax[iax].set_xticklabels([])
for iax in [9]:
    ax[iax].set_yticks([])
for iax in [1,3,5,7]:
    ax[iax].set_yticklabels([])
for iax in [0,2,4,6,8]:
    output_label = output_labels[int(iax/2)]
    input_condition = input_conditions[int(iax/2)]
    dfvmf = pd.read_csv('b007_'+output_label+'_statistics_vmf_'+yrsstr+'.csv')
    xcc = np.array(dfvmf['xcc'].to_list())
    ycc = np.array(dfvmf['ycc'].to_list())
    angle_95_deg = np.array(dfvmf['angle_95_deg'].to_list())
    gap95 = np.sin(np.radians(angle_95_deg))
    dfl_icon = pd.read_csv('b010_laplaceplane_'+output_label+'_'+yrsstr+'.csv')
    dflx = dfl_icon['laplacex'].to_list()
    dfly = dfl_icon['laplacey'].to_list()
    ax[iax].scatter(tyrsvec/1e6,xcc-dflx,color='blue',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
    # ax[iax].scatter(tyrsvec/1e6,xcc-xJ,color='red',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
    ax[iax].scatter(tyrsvec/1e6,xcc-xinvar,color='green',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
    ax[iax].scatter(tyrsvec/1e6,+gap95,color='black',alpha=1,edgecolor='none',s=1,marker='.',linewidth=1)
    ax[iax].scatter(tyrsvec/1e6,-gap95,color='black',alpha=1,edgecolor='none',s=1,marker='.',linewidth=1)
    ax[iax].set_ylim([-0.03,0.03])
    ax[iax].set_yticks([-0.015,0,0.015])
for iax in [1,3,5,7]:
    output_label = output_labels[int(iax/2-1/2)]
    input_condition = input_conditions[int(iax/2-1/2)]
    clone_label = clone_labels[int(iax/2-1/2)]
    dfcf = pd.read_csv('b008_'+clone_label+'_clones_statistics_cf_'+yrsstr+'.csv')
    xcc = np.array(dfcf['xcc'].to_list())
    ycc = np.array(dfcf['ycc'].to_list())
    angle_95_deg = np.array(dfcf['sigx_95_deg'].to_list())
    gap95a = np.sin(np.radians(angle_95_deg))
    dfl_icon = pd.read_csv('b010_laplaceplane_'+output_label+'_'+yrsstr+'.csv')
    dflx = dfl_icon['laplacex'].to_list()
    dfly = dfl_icon['laplacey'].to_list()
    ax[iax].scatter(tyrsvec/1e6,xcc-dflx,color='blue',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
    # ax[iax].scatter(tyrsvec/1e6,xcc-xJ,color='red',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
    ax[iax].scatter(tyrsvec/1e6,xcc-xinvar,color='green',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
    ax[iax].scatter(tyrsvec/1e6,+gap95a,color='black',alpha=1,edgecolor='none',s=1,marker='.',linewidth=1)
    ax[iax].scatter(tyrsvec/1e6,-gap95a,color='black',alpha=1,edgecolor='none',s=1,marker='.',linewidth=1)
    ax[iax].set_ylim([-0.03,0.03])
    ax[iax].set_yticks([-0.015,0,0.015])
for iax in range(nrows*ncols):
    ax[iax].set_xlim([0,2])
outfile_plot = 'b014_10plots_q_xcc_vmfobj_cfclones.png'
plt.savefig(outfile_plot,dpi=400)


# fig = plt.figure(figsize=(9,6)) # (width,height)
# nrows = 4
# ncols = 5
# ax = [fig.add_subplot(nrows,ncols,i+1) for i in range(nrows*ncols)]
# plt.rcParams.update({'font.size':6})
# for iax in range(nrows*ncols-1):
#     a = ax[iax]
#     a.tick_params(axis='x',direction='in')
#     a.tick_params(axis='y',direction='in')
#     a.set_xticks([0,0.5,1,1.5])
#     a.set_xlim([0,2])
# fig.subplots_adjust(wspace=0.05, hspace=0.05)
# # subplot titles go across each row, then on to the next row
# subplottitles = ['Potomacs center, q','Hildas center, q','Schubarts center, q','Backgrounds center, q','All center, q',\
#                  'Relative to Laplace, q','Relative to Laplace, q','Relative to Laplace, q','Relative to Laplace, q','Relative to Laplace, q',\
#                  'Potomac clones center, q','Hilda clones center, q','Schubart clones center, q','Ismene clones center, q','',\
#                  'Relative to Laplace, q','Relative to Laplace, q','Relative to Laplace, q','Relative to Laplace, q','']
# for iax in [0,1,2,3,4,  5,6,7,8,9,  10,11,12,13,  15,16,17,18]:
#     ax[iax].text(0.1,0.9,subplottitles[iax],horizontalalignment='left',verticalalignment='top',\
#                transform=ax[iax].transAxes,bbox=dict(facecolor='white',alpha=1))
#     ax[iax].grid('on')
# for iax in [15,16,17,18,19]:
#     ax[iax].set_xlabel('Time, Myr')
# for iax in range(15):
#     ax[iax].set_xticklabels([])
# for iax in [14,19]:
#     ax[iax].set_yticks([])
# for iax in [1,2,3,4,  6,7,8,9,  11,12,13,  16,17,18]:
#     ax[iax].set_yticklabels([])
# for iax in [0,1,2,3,4]:
#     output_label = output_labels[iax]
#     input_condition = input_conditions[iax]
#     dfvmf = pd.read_csv('b007_'+output_label+'_statistics_vmf_'+yrsstr+'.csv')
#     xcc = np.array(dfvmf['xcc'].to_list())
#     ycc = np.array(dfvmf['ycc'].to_list())
#     angle_95_deg = np.array(dfvmf['angle_95_deg'].to_list())
#     gap95 = np.sin(np.radians(angle_95_deg))
#     dfl_icon = pd.read_csv('b010_laplaceplane_'+output_label+'_'+yrsstr+'.csv')
#     dflx = dfl_icon['laplacex'].to_list()
#     dfly = dfl_icon['laplacey'].to_list()
#     ax[iax].scatter(tyrsvec/1e6,dflx,color='red',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax].scatter(tyrsvec/1e6,xcc,color='blue',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax].scatter(tyrsvec/1e6,xcc+gap95,color='cadetblue',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax].scatter(tyrsvec/1e6,xcc-gap95,color='cadetblue',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax+5].scatter(tyrsvec/1e6,xcc-dflx,color='blue',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax+5].scatter(tyrsvec/1e6,xcc-xJ,color='red',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax+5].scatter(tyrsvec/1e6,xcc-xinvar,color='green',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax+5].scatter(tyrsvec/1e6,+gap95,color='black',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax+5].scatter(tyrsvec/1e6,-gap95,color='black',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax].set_ylim([-0.05,0.05])
#     ax[iax].set_yticks([-0.03,0,0.03])
#     ax[iax+5].set_ylim([-0.03,0.03])
#     ax[iax+5].set_yticks([-0.015,0,0.015])
# for iax in [10,11,12,13]:
#     output_label = output_labels[iax-10]
#     input_condition = input_conditions[iax-10]
#     clone_label = clone_labels[iax-10]
#     dfcf = pd.read_csv('b008_'+clone_label+'_clones_statistics_cf_'+yrsstr+'.csv')
#     xcc = np.array(dfcf['xcc'].to_list())
#     ycc = np.array(dfcf['ycc'].to_list())
#     angle_95_deg = np.array(dfcf['sigx_95_deg'].to_list())
#     gap95 = np.sin(np.radians(angle_95_deg))
#     dfl_icon = pd.read_csv('b010_laplaceplane_'+output_label+'_'+yrsstr+'.csv')
#     dflx = dfl_icon['laplacex'].to_list()
#     dfly = dfl_icon['laplacey'].to_list()
#     ax[iax].scatter(tyrsvec/1e6,dflx,color='red',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax].scatter(tyrsvec/1e6,xcc,color='blue',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax].scatter(tyrsvec/1e6,xcc+gap95,color='cadetblue',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax].scatter(tyrsvec/1e6,xcc-gap95,color='cadetblue',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax+5].scatter(tyrsvec/1e6,xcc-dflx,color='blue',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax+5].scatter(tyrsvec/1e6,xcc-xJ,color='red',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax+5].scatter(tyrsvec/1e6,xcc-xinvar,color='green',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax+5].scatter(tyrsvec/1e6,+gap95,color='black',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax+5].scatter(tyrsvec/1e6,-gap95,color='black',alpha=1,edgecolor='none',s=1,marker='.',linewidth=0)
#     ax[iax].set_ylim([-0.05,0.05])
#     ax[iax].set_yticks([-0.03,0,0.03])
#     ax[iax+5].set_ylim([-0.03,0.03])
#     ax[iax+5].set_yticks([-0.015,0,0.015])