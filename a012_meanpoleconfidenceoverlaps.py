#%%
def ellipse_points_2(a0,b0,siga,sigb,rhoab,probmass):
    import numpy as np
    from scipy.stats import chi2 as sschi2
    cov = np.array([[siga**2,rhoab*siga*sigb],[rhoab*siga*sigb,sigb**2]])
    lenth = 1001
    th = np.linspace(start=0,stop=2*np.pi,num=lenth,endpoint=True)
    costh = np.cos(th)
    sinth=  np.sin(th)
    chi2val = sschi2.ppf(probmass,2)
    # print(chi2val)
    eigenvalues,eigenvectors = np.linalg.eig(cov)
    eigmin = np.min(eigenvalues)
    eigmax = np.max(eigenvalues)
    eigmaxindex = np.argmax(eigenvalues)
    eigmaxvec = eigenvectors[:,eigmaxindex]
    rotation_angle = np.arctan2(eigmaxvec[1],eigmaxvec[0])
    rotation_matrix = np.array([[np.cos(rotation_angle),-np.sin(rotation_angle)],\
                                [np.sin(rotation_angle),np.cos(rotation_angle)]])
    semimajoraxis = np.sqrt(chi2val*eigmax)
    semiminoraxis = np.sqrt(chi2val*eigmin)
    ellipse_x_vec = costh*semimajoraxis
    ellipse_y_vec = sinth*semiminoraxis
    for i in range(lenth):
        xhere = np.array([ellipse_x_vec[i],ellipse_y_vec[i]])
        xhere2 = np.matmul(rotation_matrix.reshape(2,2),xhere.reshape(2,1))
        xherex = xhere2[0][0]
        xherey = xhere2[1][0]
        ellipse_x_vec[i] = xherex
        ellipse_y_vec[i] = xherey
    ellipse_x_vec = ellipse_x_vec + a0
    ellipse_y_vec = ellipse_y_vec + b0
    return ellipse_x_vec,ellipse_y_vec
#%%
def plot_fitted_circles_covariance_3(ax,xcc,ycc,siga,sigb,rhoab,probmass,edgecolor,linestyle_in,lw):
    ellipse_x,ellipse_y = ellipse_points_2(xcc,ycc,siga,sigb,rhoab,probmass)
    ax.plot(ellipse_x,ellipse_y,color=edgecolor,linestyle=linestyle_in,linewidth=lw,alpha=0.5)
    return ax
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
tyrsvec = np.loadtxt('b002_tyrsvec_' + yrsstr + '.csv',delimiter=',')
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
plot_labels = ['P','H','S','B','B','B','B','B','B','A','A','A','A','A','A']
input_conditions = [[['Potomac'],['librating']],\
    [['Hilda'],['librating']],\
    [['Schubart'],['librating']],\
    [['Background'],['librating']],[['Background'],['librating','maybe']],[['Background'],['librating','maybe','circulating']],\
    [['Background','Francette','Guinevere'],['librating']],[['Background','Francette','Guinevere'],['librating','maybe']],[['Background','Francette','Guinevere'],['librating','maybe','circulating']],\
    [['Background','Potomac','Hilda','Schubart'],['librating']],[['Background','Potomac','Hilda','Schubart'],['librating','maybe']],[['Background','Potomac','Hilda','Schubart'],['librating','maybe','circulating']],\
    [['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating']],[['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating','maybe']],[['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating','maybe','circulating']],\
    ] 
colors = ['goldenrod','red','blue',\
          'lightgreen','limegreen','darkgreen',\
          'magenta','deeppink','purple',\
          'cyan','cadetblue','steelblue',\
          'black','dimgray','darkgray'] 
ncons = len(input_conditions)
xlim = [-0.015,0.01]
xspan = xlim[1]-xlim[0]
ymin = 0.005
ylim = [ymin,ymin+xspan]
outfile_plot = 'b012_meanpoleconfidenceoverlaps_selected_'+yrsstr+'.png'
fig = plt.figure(figsize=(7,7))
plt.rcParams['font.size'] = 10
fig,ax = plt.subplots(1,1)
thetavec = np.linspace(start=0,stop=2*np.pi,num=1001,endpoint=False)
lw_here = 0.5
alpha_here = 0.5
# for icon in range(ncons):
# for icon in [0,1,2, 3,4,5, 6,7,8, 9,10,11, 12,13,14]:
for icon in [0,1,2, 8,   14]:
    output_label = output_labels[icon]
    plot_label = plot_labels[icon]
    dfvmf = pd.read_csv('b007_'+output_label+'_statistics_vmf_'+yrsstr+'.csv')
    xcc = np.array(dfvmf['xcc'].to_list())
    ycc = np.array(dfvmf['ycc'].to_list())
    angle1rad = np.radians(np.array(dfvmf['angle_68_deg'].to_list()))
    angle2rad = np.radians(np.array(dfvmf['angle_95_deg'].to_list()))
    angle3rad = np.radians(np.array(dfvmf['angle_997_deg'].to_list()))
    x1_vec = xcc[it] + np.sin(angle1rad[it])*np.cos(thetavec)
    x2_vec = xcc[it] + np.sin(angle2rad[it])*np.cos(thetavec)
    x3_vec = xcc[it] + np.sin(angle3rad[it])*np.cos(thetavec)
    y1_vec = ycc[it] + np.sin(angle1rad[it])*np.sin(thetavec)
    y2_vec = ycc[it] + np.sin(angle2rad[it])*np.sin(thetavec)
    y3_vec = ycc[it] + np.sin(angle3rad[it])*np.sin(thetavec)
    ax.plot(x1_vec,y1_vec,color=colors[icon],linestyle='solid',lw=lw_here,alpha=alpha_here)
    ax.plot(x2_vec,y2_vec,color=colors[icon],linestyle='solid',lw=lw_here,alpha=alpha_here)
    ax.plot(x3_vec,y3_vec,color=colors[icon],linestyle='solid',lw=lw_here,alpha=alpha_here)
    ax.scatter(xcc[it],ycc[it],color=colors[icon],s=12)
    ax.text(xcc[it],ycc[it],plot_label)
    dfl_icon = pd.read_csv('b010_laplaceplane_'+output_label+'_'+yrsstr+'.csv')
    dflx = dfl_icon['laplacex'].to_list()
    dfly = dfl_icon['laplacey'].to_list()
    ax.scatter(dflx[it],dfly[it],color=colors[icon],linestyle='solid',lw=0.2,alpha=0.1)
ax.scatter(xJ,yJ,color='mediumvioletred',s=12)
ax.scatter(xI,yI,color='magenta',s=12)
ax.text(xJ,yJ,'J')
ax.text(xI,yI,'I')
ax.text(dflx[it],dfly[it],'L')
ax.set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
ax.set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
ax.set_xlim(xlim)
ax.set_ylim(ylim)
# title = 'Mean pole confidence overlaps '+str(tyrsvec[it])
title = 'Mean pole confidence overlaps'
# ax.set_title(title)
ax.grid(True)
ax.set_box_aspect(1)
plt.tight_layout()
plt.savefig(outfile_plot,dpi=400)