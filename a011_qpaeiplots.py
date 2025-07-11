#%%
def plot_lines_and_spokes(ax):
    import numpy as np
    th = np.linspace(start=0,stop=2*np.pi,num=100,endpoint=True)
    costh = np.cos(th)
    sinth=  np.sin(th)
    for angle in [90,80,70,60,50,40,30,20,10,0]: # 10-deg latitude lines
        ax.plot(np.cos(np.radians(angle))*costh,np.cos(np.radians(angle))*sinth,color='gray',\
                linestyle='-',linewidth=0.25)
    for integer in [0,1,2,3,4,5,6,7,8,9,10,11]: # 30-deg longitude lines
        ax.plot([0,np.cos(np.radians(integer*30))],[0,np.sin(np.radians(integer*30))],color='gray',\
                linestyle='-',linewidth=0.25)
    # x and y axis lines
    ax.axhline(color='black',linestyle='-',linewidth=1)
    ax.axvline(color='black',linestyle='-',linewidth=1)
    return ax
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
dfa = pd.read_csv('b006_astorb_labels.csv')
n_astorb = dfa.shape[0]
probmasses = [0.68,0.95,0.997]
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
output_labels = ['P_L',      'H_L',       'S_L',\
                 'B_L',      'B_LM',      'B_LMC',\
                 'BFG_L',    'BFG_LM',    'BFG_LMC',\
                 'BPHS_L',   'BPHS_LM',   'BPHS_LMC',\
                 'BPHSFG_L', 'BPHSFG_LM', 'BPHSFG_LMC']
input_conditions = [[['Potomac'],['librating']],\
    [['Hilda'],['librating']],\
    [['Schubart'],['librating']],\
    [['Background'],['librating']],[['Background'],['librating','maybe']],[['Background'],['librating','maybe','circulating']],\
    [['Background','Francette','Guinevere'],['librating']],[['Background','Francette','Guinevere'],['librating','maybe']],[['Background','Francette','Guinevere'],['librating','maybe','circulating']],\
    [['Background','Potomac','Hilda','Schubart'],['librating']],[['Background','Potomac','Hilda','Schubart'],['librating','maybe']],[['Background','Potomac','Hilda','Schubart'],['librating','maybe','circulating']],\
    [['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating']],[['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating','maybe']],[['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating','maybe','circulating']],\
    ] 
colors = ['goldenrod','red','blue',\
          'green','green','green',\
          'green','green','green',\
          'black','black','black',\
          'black','black','black'] 
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
sizes = [5,20,40]
alphas = [0.1,0.3,1]
markers = ['.','s','o']
ncons = len(input_conditions)
#%%
outfile_plot = 'b011_qp_all'+yrsstr+'.png'
it = 0
# labels_to_use = ['P_L','H_L','S_L','BFG_LMC','BPHSFG_LMC']
labels_to_use = ['BFG_LMC','P_L','H_L','S_L']
nuse = len(labels_to_use)
xlim = [-0.25,0.25]
ylim = [-0.25,0.25]
fig = plt.figure(figsize=(7,7))
plt.rcParams['font.size'] = 18
ax = fig.add_subplot(111)
ax = plot_lines_and_spokes(ax)
for iuse in range(nuse):
    label_use = labels_to_use[iuse]
    indx = output_labels.index(label_use)
    input_groups = input_conditions[indx][0]
    input_labels = input_conditions[indx][1]
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) and (dfa['label'][i_astorb] in input_labels):
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    print(label_use,nobj)
    it = 0
    aau_icon = aau_astorb[astorb_indices,it]
    e_icon = e_astorb[astorb_indices,it]
    ideg_icon = ideg_astorb[astorb_indices,it]
    q_icon = q_astorb[astorb_indices,it]
    p_icon = p_astorb[astorb_indices,it]
    x = q_icon
    y = p_icon
    ax.scatter(x,y,color=colors[indx],s=20,alpha=1,edgecolor='none',linewidth=0,marker='.')
plt.axis('equal')
ax.set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
ax.set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_box_aspect(1)
plt.tight_layout()
plt.savefig(outfile_plot,dpi=300)
plt.show()    
#%%
# for icon in range(ncons):
#     output_label = output_labels[icon]
#     input_groups = input_conditions[icon][0]
#     input_labels = input_conditions[icon][1]
#     astorb_indices = []
#     for i_astorb in range(n_astorb):
#         if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) and (dfa['label'][i_astorb] in input_labels):
#             astorb_indices.append(i_astorb)
#     nobj = len(astorb_indices)
#     print(output_label,nobj)
#     it = 0
#     aau_icon = aau_astorb[astorb_indices,it]
#     e_icon = e_astorb[astorb_indices,it]
#     ideg_icon = ideg_astorb[astorb_indices,it]
#     q_icon = q_astorb[astorb_indices,it]
#     p_icon = p_astorb[astorb_indices,it]
#     outfile_plot = 'b011_qp_' + output_label + '_all_'+yrsstr+'.png'
#     xlim = [-0.35,0.35]
#     ylim = [-0.35,0.35]
#     fig = plt.figure(figsize=(7,7))
#     plt.rcParams['font.size'] = 18
#     ax = fig.add_subplot(111)
#     ax = plot_lines_and_spokes(ax)
#     x = q_icon
#     y = p_icon
#     ax.scatter(x,y,color=colors[icon],s=10,alpha=1,edgecolor='none',linewidth=0,marker='.')
#     plt.axis('equal')
#     ax.set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
#     ax.set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
#     ax.set_xlim(xlim)
#     ax.set_ylim(ylim)
#     ax.set_box_aspect(1)
#     plt.title(output_label+' '+str(nobj))
#     plt.tight_layout()
#     plt.savefig(outfile_plot,dpi=300)
#     plt.show()    
#     outfile_plot = 'b011_qp_' + output_label + '_distinctLMC_'+yrsstr+'.png'
#     it = 0
#     xlim = [-0.35,0.35]
#     ylim = [-0.35,0.35]
#     fig = plt.figure(figsize=(7,7))
#     plt.rcParams['font.size'] = 18
#     ax = fig.add_subplot(111)
#     ax = plot_lines_and_spokes(ax)
#     groups_in_this_condition = groups_in_conditions[icon]
#     ngroups_con = len(groups_in_this_condition)
#     for igroup in range(ngroups_con):
#         individual_group = groups_in_this_condition[igroup]
#         individual_group_index = distinct_groups.index(individual_group)
#         input_groups = conditions_for_distinct_groups[individual_group_index][0]
#         input_labels = conditions_for_distinct_groups[individual_group_index][1]
#         astorb_indices = []
#         for i_astorb in range(n_astorb):
#             if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) and (dfa['label'][i_astorb] in input_labels):
#                 astorb_indices.append(i_astorb)
#         nobj = len(astorb_indices)
#         aau_icon = aau_astorb[astorb_indices,it]
#         e_icon = e_astorb[astorb_indices,it]
#         ideg_icon = ideg_astorb[astorb_indices,it]
#         q_icon = q_astorb[astorb_indices,it]
#         p_icon = p_astorb[astorb_indices,it]
#         x = q_icon
#         y = p_icon
#         ax.scatter(x,y,color=colors[icon],s=sizes[igroup],alpha=1,edgecolor='none',linewidth=0,marker=markers[igroup])
#     plt.axis('equal')
#     ax.set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
#     ax.set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
#     ax.set_xlim(xlim)
#     ax.set_ylim(ylim)
#     ax.set_box_aspect(1)
#     plt.title(output_label+' distinct LMC')
#     plt.tight_layout()
#     plt.savefig(outfile_plot,dpi=300)
#     plt.show()
#%%
outfile_plot = 'b011_aei_all'+yrsstr+'.png'
it = 0
# labels_to_use = ['P_L','H_L','S_L','BFG_LMC','BPHSFG_LMC']
labels_to_use = ['BFG_LMC','P_L','H_L','S_L']
nuse = len(labels_to_use)
fig = plt.figure(figsize=(7,7))
plt.rcParams['font.size'] = 18
ax = fig.add_subplot(223)
for iuse in range(nuse):
    label_use = labels_to_use[iuse]
    indx = output_labels.index(label_use)
    input_groups = input_conditions[indx][0]
    input_labels = input_conditions[indx][1]
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) and (dfa['label'][i_astorb] in input_labels):
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    print(label_use,nobj)
    it = 0
    aau_icon = aau_astorb[astorb_indices,it]
    e_icon = e_astorb[astorb_indices,it]
    ideg_icon = ideg_astorb[astorb_indices,it]
    x = aau_icon
    y = e_icon
    ax.scatter(x,y,color=colors[indx],s=10,alpha=1,edgecolor='none',linewidth=0,marker='.')
ax.set_xlabel('a')
ax.set_ylabel('e')
ax = fig.add_subplot(221)
for iuse in range(nuse):
    label_use = labels_to_use[iuse]
    indx = output_labels.index(label_use)
    input_groups = input_conditions[indx][0]
    input_labels = input_conditions[indx][1]
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) and (dfa['label'][i_astorb] in input_labels):
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    print(label_use,nobj)
    it = 0
    aau_icon = aau_astorb[astorb_indices,it]
    e_icon = e_astorb[astorb_indices,it]
    ideg_icon = ideg_astorb[astorb_indices,it]
    x = aau_icon
    y = ideg_icon
    ax.scatter(x,y,color=colors[indx],s=10,alpha=1,edgecolor='none',linewidth=0,marker='.')
ax.set_xlabel('a')
ax.set_ylabel('i (deg)')
ax = fig.add_subplot(224)
for iuse in range(nuse):
    label_use = labels_to_use[iuse]
    indx = output_labels.index(label_use)
    input_groups = input_conditions[indx][0]
    input_labels = input_conditions[indx][1]
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) and (dfa['label'][i_astorb] in input_labels):
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    print(label_use,nobj)
    it = 0
    aau_icon = aau_astorb[astorb_indices,it]
    e_icon = e_astorb[astorb_indices,it]
    ideg_icon = ideg_astorb[astorb_indices,it]
    x = e_icon
    y = ideg_icon
    ax.scatter(x,y,color=colors[indx],s=10,alpha=1,edgecolor='none',linewidth=0,marker='.')
ax.set_xlabel('e')
ax.set_ylabel('i (deg)')
plt.tight_layout()
plt.savefig(outfile_plot,dpi=300)
plt.show()    
#%%
# for icon in range(ncons):
#     output_label = output_labels[icon]
#     input_groups = input_conditions[icon][0]
#     input_labels = input_conditions[icon][1]
#     astorb_indices = []
#     for i_astorb in range(n_astorb):
#         if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) and (dfa['label'][i_astorb] in input_labels):
#             astorb_indices.append(i_astorb)
#     nobj = len(astorb_indices)
#     print(output_label,nobj)
#     it = 0
#     aau_icon = aau_astorb[astorb_indices,it]
#     e_icon = e_astorb[astorb_indices,it]
#     ideg_icon = ideg_astorb[astorb_indices,it]
#     outfile_plot = 'b011_aei_' + output_label + '_all_'+yrsstr+'.png'
#     fig = plt.figure(figsize=(7,7))
#     plt.rcParams['font.size'] = 18
#     ax = fig.add_subplot(223)
#     ax.scatter(aau_icon,e_icon,color=colors[icon],s=10,alpha=1,edgecolor='none',linewidth=0,marker='.')
#     ax.set_xlabel('a')
#     ax.set_ylabel('e')
#     ax = fig.add_subplot(221)
#     ax.scatter(aau_icon,ideg_icon,color=colors[icon],s=10,alpha=1,edgecolor='none',linewidth=0,marker='.')
#     ax.set_xlabel('a')
#     ax.set_ylabel('i (deg)')
#     ax = fig.add_subplot(224)
#     ax.scatter(e_icon,ideg_icon,color=colors[icon],s=10,alpha=1,edgecolor='none',linewidth=0,marker='.')
#     ax.set_xlabel('e')
#     ax.set_ylabel('i (deg)')
#     plt.suptitle(output_label+' a, e, i')
#     plt.tight_layout()
#     plt.savefig(outfile_plot,dpi=300)
#     plt.show()    
#     outfile_plot = 'b011_aei_' + output_label + '_distinctLMC_'+yrsstr+'.png'
#     it = 0
#     fig = plt.figure(figsize=(7,7))
#     plt.rcParams['font.size'] = 18
#     ax = fig.add_subplot(223)
#     ax2 = fig.add_subplot(221)
#     ax3 = fig.add_subplot(224)
#     groups_in_this_condition = groups_in_conditions[icon]
#     ngroups_con = len(groups_in_this_condition)
#     for igroup in range(ngroups_con):
#         individual_group = groups_in_this_condition[igroup]
#         individual_group_index = distinct_groups.index(individual_group)
#         input_groups = conditions_for_distinct_groups[individual_group_index][0]
#         input_labels = conditions_for_distinct_groups[individual_group_index][1]
#         astorb_indices = []
#         for i_astorb in range(n_astorb):
#             if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) and (dfa['label'][i_astorb] in input_labels):
#                 astorb_indices.append(i_astorb)
#         nobj = len(astorb_indices)
#         aau_icon = aau_astorb[astorb_indices,it]
#         e_icon = e_astorb[astorb_indices,it]
#         ideg_icon = ideg_astorb[astorb_indices,it]       
#         ax.scatter(aau_icon,e_icon,color=colors[icon],s=sizes[igroup],alpha=1,edgecolor='none',linewidth=0,marker=markers[igroup])
#         ax2.scatter(aau_icon,ideg_icon,color=colors[icon],s=sizes[igroup],alpha=1,edgecolor='none',linewidth=0,marker=markers[igroup])
#         ax3.scatter(e_icon,ideg_icon,color=colors[icon],s=sizes[igroup],alpha=1,edgecolor='none',linewidth=0,marker=markers[igroup])
#     ax.set_xlabel('a')
#     ax.set_ylabel('e')
#     ax2.set_xlabel('a')
#     ax2.set_ylabel('i (deg)')
#     ax3.set_xlabel('e')
#     ax3.set_ylabel('i (deg)')
#     plt.suptitle(output_label+' distinct LMC a, e, i')
#     plt.tight_layout()
#     plt.savefig(outfile_plot,dpi=300)
#     plt.show()