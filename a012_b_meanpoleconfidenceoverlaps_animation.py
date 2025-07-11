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
import matplotlib.animation as animation
import os
# ffmpeg_path = "/path/to/ffmpeg"  
ffmpeg_path = '/usr/local/bin/ffmpeg'
os.environ['PATH'] += f':{os.path.dirname(ffmpeg_path)}'
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
it = 5000
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
          'black','dimgray','darkgray',\
          'cyan','cadetblue','steelblue',\
          'magenta','deeppink','purple'] 
plot_labels = ['P','H','S','B','B','B','B','B','B','A','A','A','A','A','A']
ncons = len(input_conditions)
xlim = [-0.05,0.025]
xspan = xlim[1]-xlim[0]
ymin = -0.01
ylim = [ymin,ymin+xspan]
outfile_plot = 'b012_b_meanpoleconfidenceoverlaps_animation_PHB'+yrsstr+'.mp4'
fig = plt.figure(figsize=(7,7))
plt.rcParams['font.size'] = 10
fig,ax = plt.subplots(1,1)
thetavec = np.linspace(start=0,stop=2*np.pi,num=1001,endpoint=False)
lw_here = 0.5
alpha_here = 0.5
def update(frame):
    print(frame)
    ax.clear()
    # for icon in range(ncons):
    # for icon in [0,1,2, 3,4,5, 6,7,8, 9,10,11, 12,13,14]:
    for icon in [0,1, 8    ]:
    # for icon in [6,12]:
        output_label = output_labels[icon]
        plot_label = plot_labels[icon]
        dfvmf = pd.read_csv('b007_'+output_label+'_statistics_vmf_'+yrsstr+'.csv')
        xcc = np.array(dfvmf['xcc'].to_list())
        ycc = np.array(dfvmf['ycc'].to_list())
        angle1rad = np.radians(np.array(dfvmf['angle_68_deg'].to_list()))
        angle2rad = np.radians(np.array(dfvmf['angle_95_deg'].to_list()))
        angle3rad = np.radians(np.array(dfvmf['angle_997_deg'].to_list()))
        x1_vec = xcc[frame] + np.sin(angle1rad[frame])*np.cos(thetavec)
        x2_vec = xcc[frame] + np.sin(angle2rad[frame])*np.cos(thetavec)
        x3_vec = xcc[frame] + np.sin(angle3rad[frame])*np.cos(thetavec)
        y1_vec = ycc[frame] + np.sin(angle1rad[frame])*np.sin(thetavec)
        y2_vec = ycc[frame] + np.sin(angle2rad[frame])*np.sin(thetavec)
        y3_vec = ycc[frame] + np.sin(angle3rad[frame])*np.sin(thetavec)
        # ax.scatter(xcc[:frame],ycc[:frame],color=colors[icon],s=1,alpha=0.1)
        ax.plot(x1_vec,y1_vec,color=colors[icon],linestyle='solid',lw=lw_here,alpha=alpha_here)
        ax.plot(x2_vec,y2_vec,color=colors[icon],linestyle='solid',lw=lw_here,alpha=alpha_here)
        ax.plot(x3_vec,y3_vec,color=colors[icon],linestyle='solid',lw=lw_here,alpha=alpha_here)
        ax.scatter(xcc[frame],ycc[frame],color=colors[icon],s=12)
        # ax.text(xcc[frame],ycc[frame],output_label)
        ax.text(xcc[frame],ycc[frame],plot_label)
        dfl_icon = pd.read_csv('b010_laplaceplane_'+output_label+'_'+yrsstr+'.csv')
        dflx = dfl_icon['laplacex'].to_list()
        dfly = dfl_icon['laplacey'].to_list()
        # ax.scatter(dflx[:frame],dfly[:frame],color=colors[icon],s=1,alpha=0.1)
        ax.scatter(dflx[frame],dfly[frame],color=colors[icon],s=12,alpha=1)
    # ax.scatter(qplanets[0,:frame],pplanets[0,:frame],color='mediumvioletred',s=1,alpha=0.1)
    # ax.scatter(xinvar[:frame],yinvar[:frame],color='magenta',s=1,alpha=0.1)
    ax.scatter(qplanets[0,frame],pplanets[0,frame],color='mediumvioletred',s=12,alpha=1)
    ax.scatter(xinvar[frame],yinvar[frame],color='magenta',s=12,alpha=1)
    ax.text(qplanets[0,frame],pplanets[0,frame],'J')
    ax.text(xinvar[frame],yinvar[frame],'I')
    ax.set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
    ax.set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    title = 'Mean pole confidence overlaps '+str(tyrsvec[frame])
    ax.set_title(title)
    ax.grid(True)
    ax.set_box_aspect(1)
    # plt.tight_layout()
    # plt.savefig(outfile_plot,dpi=400)
fps = 12
nt = len(tyrsvec)
framevec = np.arange(start=0,stop=nt,step=5)
ani = animation.FuncAnimation(fig, update, frames=framevec, interval=1000/fps)
writer = animation.FFMpegWriter(fps=fps)
# outfile = 'animation_laplaceJupiter' + yrsstr + '_2Myrplot.mp4'
ani.save(outfile_plot,writer=writer)
# import numpy as np
# import pandas as pd
# from matplotlib import pyplot as plt
# # laplace indices for 3.8026 to 4.1608 au
# laplaceindex1 = 134
# laplaceindex2 = 186
# laplacemid = int((laplaceindex1+laplaceindex2)/2)
# xL = np.loadtxt('laplaceplane_2_6_q.csv',delimiter=',')
# yL = np.loadtxt('laplaceplane_2_6_p.csv',delimiter=',')
# # xL = xL[laplaceindex1:laplaceindex2]
# # yL = yL[laplaceindex1:laplaceindex2]
# xL = xL[laplacemid]
# yL = yL[laplacemid]
# xI = np.loadtxt('005_a_planets_xinvar_2e6yr_1e2yr_0.2yr_jd2460796.5.csv')
# yI = np.loadtxt('005_a_planets_yinvar_2e6yr_1e2yr_0.2yr_jd2460796.5.csv')
# xI = xI[100] # arbitrary timepoint
# yI = yI[100] # arbitrary timepoint
# probmasses = [0.68,0.95,0.997]
# jd = 2460796.5 # May 1, 2025 00:00:00
# time_yrs  = int(2e6)
# tstep_yrs = int(2e2)
# dt_yrs = 0.2
# tyrs_str = '2e6yr'
# tstepyrs_str = '1e2yr'
# dtyrs_str = '0.2yr'
# yrsstr = tyrs_str + '_' + tstepyrs_str + '_' + dtyrs_str + '_jd' + str(jd)
# aau_mat_planets = np.loadtxt('005_a_planets_aau_' + yrsstr + '.csv',delimiter=',')
# e_mat_planets = np.loadtxt('005_a_planets_e_' + yrsstr + '.csv',delimiter=',')
# irad_mat_planets = np.radians(np.loadtxt('005_a_planets_ideg_' + yrsstr + '.csv',delimiter=','))
# wrad_mat_planets = np.radians(np.loadtxt('005_a_planets_perideg_' + yrsstr + '.csv',delimiter=','))
# Wrad_mat_planets = np.radians(np.loadtxt('005_a_planets_nodedeg_' + yrsstr + '.csv',delimiter=','))
# Mrad_mat_planets = np.radians(np.loadtxt('005_a_planets_Mdeg_' + yrsstr + '.csv',delimiter=','))
# iJ = irad_mat_planets[0,0]
# WJ = Wrad_mat_planets[0,0]
# xJ = np.sin(iJ)*np.cos(WJ)
# yJ = np.sin(iJ)*np.sin(WJ)
#%%
# title = 'Confidence ellipses for vmf'
# outfile_plot = 'mean_pole_confidence_overlaps_start_may2025_Hlt15.7_librating_vmf.png'
# xlim = [-0.04,0.04]
# xspan = xlim[1]-xlim[0]
# ymin = -0.02
# ylim = [ymin,ymin+xspan]
# df = pd.read_csv('a002_Potomac_Hlt15.7_librating_family_horizons.csv')
# irad_Potomac = np.radians(np.array(df['incl_deg'].tolist()))
# Wrad_Potomac = np.radians(np.array(df['Omega_deg'].tolist()))
# x_Potomac = np.sin(irad_Potomac)*np.cos(Wrad_Potomac)
# y_Potomac = np.sin(irad_Potomac)*np.sin(Wrad_Potomac)
# df = pd.read_csv('a002_Hilda_Hlt15.7_librating_family_horizons.csv')
# irad_Hilda = np.radians(np.array(df['incl_deg'].tolist()))
# Wrad_Hilda = np.radians(np.array(df['Omega_deg'].tolist()))
# x_Hilda = np.sin(irad_Hilda)*np.cos(Wrad_Hilda)
# y_Hilda = np.sin(irad_Hilda)*np.sin(Wrad_Hilda)
# df = pd.read_csv('a002_Schubart_Hlt15.7_librating_family_horizons.csv')
# irad_Schubart = np.radians(np.array(df['incl_deg'].tolist()))
# Wrad_Schubart = np.radians(np.array(df['Omega_deg'].tolist()))
# x_Schubart = np.sin(irad_Schubart)*np.cos(Wrad_Schubart)
# y_Schubart = np.sin(irad_Schubart)*np.sin(Wrad_Schubart)
# df = pd.read_csv('a002_background_Hlt15.7_librating_family_horizons.csv')
# irad_background = np.radians(np.array(df['incl_deg'].tolist()))
# Wrad_background = np.radians(np.array(df['Omega_deg'].tolist()))
# x_background = np.sin(irad_background)*np.cos(Wrad_background)
# y_background = np.sin(irad_background)*np.sin(Wrad_background)
# df = pd.read_csv('a002_astorb_Hlt15.7_librating_family_horizons.csv')
# irad_astorb = np.radians(np.array(df['incl_deg'].tolist()))
# Wrad_astorb = np.radians(np.array(df['Omega_deg'].tolist()))
# x_astorb = np.sin(irad_astorb)*np.cos(Wrad_astorb)
# y_astorb = np.sin(irad_astorb)*np.sin(Wrad_astorb)
# xccvec_Potomac = np.loadtxt('Potomac_xccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# yccvec_Potomac = np.loadtxt('Potomac_yccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle1_deg_vec_Potomac = np.loadtxt('Potomac_angle1_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle2_deg_vec_Potomac = np.loadtxt('Potomac_angle2_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle3_deg_vec_Potomac = np.loadtxt('Potomac_angle3_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# xcc_Potomac = xccvec_Potomac[0]
# ycc_Potomac = yccvec_Potomac[0]
# angle1_deg_Potomac = angle1_deg_vec_Potomac[0]
# angle2_deg_Potomac = angle2_deg_vec_Potomac[0]
# angle3_deg_Potomac = angle3_deg_vec_Potomac[0]
# xccvec_Hilda = np.loadtxt('Hilda_xccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# yccvec_Hilda = np.loadtxt('Hilda_yccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle1_deg_vec_Hilda = np.loadtxt('Hilda_angle1_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle2_deg_vec_Hilda = np.loadtxt('Hilda_angle2_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle3_deg_vec_Hilda = np.loadtxt('Hilda_angle3_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# xcc_Hilda = xccvec_Hilda[0]
# ycc_Hilda = yccvec_Hilda[0]
# angle1_deg_Hilda = angle1_deg_vec_Hilda[0]
# angle2_deg_Hilda = angle2_deg_vec_Hilda[0]
# angle3_deg_Hilda = angle3_deg_vec_Hilda[0]
# xccvec_Schubart = np.loadtxt('Schubart_xccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# yccvec_Schubart = np.loadtxt('Schubart_yccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle1_deg_vec_Schubart = np.loadtxt('Schubart_angle1_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle2_deg_vec_Schubart = np.loadtxt('Schubart_angle2_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle3_deg_vec_Schubart = np.loadtxt('Schubart_angle3_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# xcc_Schubart = xccvec_Schubart[0]
# ycc_Schubart = yccvec_Schubart[0]
# angle1_deg_Schubart = angle1_deg_vec_Schubart[0]
# angle2_deg_Schubart = angle2_deg_vec_Schubart[0]
# angle3_deg_Schubart = angle3_deg_vec_Schubart[0]
# xccvec_background = np.loadtxt('background_xccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# yccvec_background = np.loadtxt('background_yccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle1_deg_vec_background = np.loadtxt('background_angle1_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle2_deg_vec_background = np.loadtxt('background_angle2_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle3_deg_vec_background = np.loadtxt('background_angle3_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# xcc_background = xccvec_background[0]
# ycc_background = yccvec_background[0]
# angle1_deg_background = angle1_deg_vec_background[0]
# angle2_deg_background = angle2_deg_vec_background[0]
# angle3_deg_background = angle3_deg_vec_background[0]
# xccvec_astorb = np.loadtxt('astorb_xccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# yccvec_astorb = np.loadtxt('astorb_yccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle1_deg_vec_astorb = np.loadtxt('astorb_angle1_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle2_deg_vec_astorb = np.loadtxt('astorb_angle2_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle3_deg_vec_astorb = np.loadtxt('astorb_angle3_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# xcc_astorb = xccvec_astorb[0]
# ycc_astorb = yccvec_astorb[0]
# angle1_deg_astorb = angle1_deg_vec_astorb[0]
# angle2_deg_astorb = angle2_deg_vec_astorb[0]
# angle3_deg_astorb = angle3_deg_vec_astorb[0]
# xccvec_bphs = np.loadtxt('bphs_xccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# yccvec_bphs = np.loadtxt('bphs_yccvec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle1_deg_vec_bphs = np.loadtxt('bphs_angle1_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle2_deg_vec_bphs = np.loadtxt('bphs_angle2_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle3_deg_vec_bphs = np.loadtxt('bphs_angle3_deg_vec_Hlt15.7_librating_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# xcc_bphs = xccvec_bphs[0]
# ycc_bphs = yccvec_bphs[0]
# angle1_deg_bphs = angle1_deg_vec_bphs[0]
# angle2_deg_bphs = angle2_deg_vec_bphs[0]
# angle3_deg_bphs = angle3_deg_vec_bphs[0]
# fig = plt.figure(figsize=(7,7))
# plt.rcParams['font.size'] = 10
# fig,ax = plt.subplots(1,1)
# lw_here = 0.5
# # for probmass in probmasses:
#     # ax = plot_fitted_circles_covariance_3(ax,xcc_Potomac,ycc_Potomac,siga_Potomac,sigb_Potomac,rhoab_Potomac,probmass,edgecolor='goldenrod',linestyle_in='dashdot',lw=lw_here)
#     # ax = plot_fitted_circles_covariance_3(ax,xcc_Hilda,ycc_Hilda,siga_Hilda,sigb_Hilda,rhoab_Hilda,probmass,edgecolor='red',linestyle_in='solid',lw=lw_here)
#     # ax = plot_fitted_circles_covariance_3(ax,xcc_Schubart,ycc_Schubart,siga_Schubart,sigb_Schubart,rhoab_Schubart,probmass,edgecolor='blue',linestyle_in='dotted',lw=lw_here)
# thetavec = np.linspace(start=0,stop=2*np.pi,num=1001,endpoint=False)
# x1_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle1_deg_Potomac))*np.cos(thetavec)
# x2_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle2_deg_Potomac))*np.cos(thetavec)
# x3_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle3_deg_Potomac))*np.cos(thetavec)
# y1_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle1_deg_Potomac))*np.sin(thetavec)
# y2_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle2_deg_Potomac))*np.sin(thetavec)
# y3_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle3_deg_Potomac))*np.sin(thetavec)
# x1_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle1_deg_Hilda))*np.cos(thetavec)
# x2_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle2_deg_Hilda))*np.cos(thetavec)
# x3_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle3_deg_Hilda))*np.cos(thetavec)
# y1_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle1_deg_Hilda))*np.sin(thetavec)
# y2_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle2_deg_Hilda))*np.sin(thetavec)
# y3_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle3_deg_Hilda))*np.sin(thetavec)
# x1_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle1_deg_Schubart))*np.cos(thetavec)
# x2_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle2_deg_Schubart))*np.cos(thetavec)
# x3_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle3_deg_Schubart))*np.cos(thetavec)
# y1_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle1_deg_Schubart))*np.sin(thetavec)
# y2_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle2_deg_Schubart))*np.sin(thetavec)
# y3_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle3_deg_Schubart))*np.sin(thetavec)
# x1_vec_background = xcc_background + np.sin(np.radians(angle1_deg_background))*np.cos(thetavec)
# x2_vec_background = xcc_background + np.sin(np.radians(angle2_deg_background))*np.cos(thetavec)
# x3_vec_background = xcc_background + np.sin(np.radians(angle3_deg_background))*np.cos(thetavec)
# y1_vec_background = ycc_background + np.sin(np.radians(angle1_deg_background))*np.sin(thetavec)
# y2_vec_background = ycc_background + np.sin(np.radians(angle2_deg_background))*np.sin(thetavec)
# y3_vec_background = ycc_background + np.sin(np.radians(angle3_deg_background))*np.sin(thetavec)
# x1_vec_astorb = xcc_astorb + np.sin(np.radians(angle1_deg_astorb))*np.cos(thetavec)
# x2_vec_astorb = xcc_astorb + np.sin(np.radians(angle2_deg_astorb))*np.cos(thetavec)
# x3_vec_astorb = xcc_astorb + np.sin(np.radians(angle3_deg_astorb))*np.cos(thetavec)
# y1_vec_astorb = ycc_astorb + np.sin(np.radians(angle1_deg_astorb))*np.sin(thetavec)
# y2_vec_astorb = ycc_astorb + np.sin(np.radians(angle2_deg_astorb))*np.sin(thetavec)
# y3_vec_astorb = ycc_astorb + np.sin(np.radians(angle3_deg_astorb))*np.sin(thetavec)
# x1_vec_bphs = xcc_bphs + np.sin(np.radians(angle1_deg_bphs))*np.cos(thetavec)
# x2_vec_bphs = xcc_bphs + np.sin(np.radians(angle2_deg_bphs))*np.cos(thetavec)
# x3_vec_bphs = xcc_bphs + np.sin(np.radians(angle3_deg_bphs))*np.cos(thetavec)
# y1_vec_bphs = ycc_bphs + np.sin(np.radians(angle1_deg_bphs))*np.sin(thetavec)
# y2_vec_bphs = ycc_bphs + np.sin(np.radians(angle2_deg_bphs))*np.sin(thetavec)
# y3_vec_bphs = ycc_bphs + np.sin(np.radians(angle3_deg_bphs))*np.sin(thetavec)
# ax.plot(x1_vec_Potomac,y1_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Potomac,y2_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Potomac,y3_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_Hilda,y1_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Hilda,y2_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Hilda,y3_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_Schubart,y1_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Schubart,y2_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Schubart,y3_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_background,y1_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_background,y2_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_background,y3_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_astorb,y1_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_astorb,y2_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_astorb,y3_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_bphs,y1_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_bphs,y2_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_bphs,y3_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.scatter(xL,yL,color='cyan',s=12)
# ax.scatter(xJ,yJ,color='cyan',s=12)
# ax.scatter(xI,yI,color='magenta',s=12)
# ax.scatter(xcc_Potomac,ycc_Potomac,color='goldenrod',s=12)
# ax.scatter(xcc_Hilda,ycc_Hilda,color='red',s=12)
# ax.scatter(xcc_Schubart,ycc_Schubart,color='blue',s=12)
# ax.scatter(xcc_background,ycc_background,color='black',s=12)
# ax.scatter(xcc_astorb,ycc_astorb,color='green',s=12)
# ax.scatter(xcc_bphs,ycc_bphs,color='magenta',s=12)
# ax.text(xL,yL,'L')
# ax.text(xJ,yJ,'J')
# ax.text(xI,yI,'I')
# ax.text(xcc_Potomac,ycc_Potomac,'P')
# ax.text(xcc_Hilda,ycc_Hilda,'H')
# # ax.text(xcc_Schubart,ycc_Schubart,'S')
# ax.text(xcc_background,ycc_background,'B')
# ax.text(xcc_astorb,ycc_astorb,'A')
# ax.text(xcc_bphs,ycc_bphs,'Z')
# ax.set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
# ax.set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
# ax.set_xlim(xlim)
# ax.set_ylim(ylim)
# ax.set_title(title)
# # ax.set_xticks([-0.015,-0.010,-0.005,0,0.005])
# # ax.set_yticks([0.015,0.020,0.025,0.030])
# ax.set_box_aspect(1)
# ax.grid()
# plt.savefig(outfile_plot,dpi=400)
# plt.tight_layout()
#%%
# outfile_plot = 'mean_pole_confidence_overlaps_start_may2025_Hlt15.7_librating_2panels_vmf.png'
# xlim = [-0.04,0.04]
# xspan = xlim[1]-xlim[0]
# ymin = -0.02
# ylim = [ymin,ymin+xspan]
# fig = plt.figure(figsize=(6,3))
# plt.rcParams['font.size'] = 10
# ax = fig.add_subplot(121)
# lw_here = 1.5
# # for probmass in probmasses:
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Potomac,ycc_Potomac,siga_Potomac,sigb_Potomac,rhoab_Potomac,probmass,edgecolor='goldenrod',linestyle_in='dashdot',lw=lw_here)
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Hilda,ycc_Hilda,siga_Hilda,sigb_Hilda,rhoab_Hilda,probmass,edgecolor='red',linestyle_in='solid',lw=lw_here)
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Schubart,ycc_Schubart,siga_Schubart,sigb_Schubart,rhoab_Schubart,probmass,edgecolor='blue',linestyle_in='dotted',lw=lw_here)
# thetavec = np.linspace(start=0,stop=2*np.pi,num=1001,endpoint=False)
# x1_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle1_deg_Potomac))*np.cos(thetavec)
# x2_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle2_deg_Potomac))*np.cos(thetavec)
# x3_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle3_deg_Potomac))*np.cos(thetavec)
# y1_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle1_deg_Potomac))*np.sin(thetavec)
# y2_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle2_deg_Potomac))*np.sin(thetavec)
# y3_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle3_deg_Potomac))*np.sin(thetavec)
# x1_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle1_deg_Hilda))*np.cos(thetavec)
# x2_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle2_deg_Hilda))*np.cos(thetavec)
# x3_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle3_deg_Hilda))*np.cos(thetavec)
# y1_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle1_deg_Hilda))*np.sin(thetavec)
# y2_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle2_deg_Hilda))*np.sin(thetavec)
# y3_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle3_deg_Hilda))*np.sin(thetavec)
# x1_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle1_deg_Schubart))*np.cos(thetavec)
# x2_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle2_deg_Schubart))*np.cos(thetavec)
# x3_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle3_deg_Schubart))*np.cos(thetavec)
# y1_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle1_deg_Schubart))*np.sin(thetavec)
# y2_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle2_deg_Schubart))*np.sin(thetavec)
# y3_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle3_deg_Schubart))*np.sin(thetavec)
# x1_vec_background = xcc_background + np.sin(np.radians(angle1_deg_background))*np.cos(thetavec)
# x2_vec_background = xcc_background + np.sin(np.radians(angle2_deg_background))*np.cos(thetavec)
# x3_vec_background = xcc_background + np.sin(np.radians(angle3_deg_background))*np.cos(thetavec)
# y1_vec_background = ycc_background + np.sin(np.radians(angle1_deg_background))*np.sin(thetavec)
# y2_vec_background = ycc_background + np.sin(np.radians(angle2_deg_background))*np.sin(thetavec)
# y3_vec_background = ycc_background + np.sin(np.radians(angle3_deg_background))*np.sin(thetavec)
# x1_vec_astorb = xcc_astorb + np.sin(np.radians(angle1_deg_astorb))*np.cos(thetavec)
# x2_vec_astorb = xcc_astorb + np.sin(np.radians(angle2_deg_astorb))*np.cos(thetavec)
# x3_vec_astorb = xcc_astorb + np.sin(np.radians(angle3_deg_astorb))*np.cos(thetavec)
# y1_vec_astorb = ycc_astorb + np.sin(np.radians(angle1_deg_astorb))*np.sin(thetavec)
# y2_vec_astorb = ycc_astorb + np.sin(np.radians(angle2_deg_astorb))*np.sin(thetavec)
# y3_vec_astorb = ycc_astorb + np.sin(np.radians(angle3_deg_astorb))*np.sin(thetavec)
# x1_vec_bphs = xcc_bphs + np.sin(np.radians(angle1_deg_bphs))*np.cos(thetavec)
# x2_vec_bphs = xcc_bphs + np.sin(np.radians(angle2_deg_bphs))*np.cos(thetavec)
# x3_vec_bphs = xcc_bphs + np.sin(np.radians(angle3_deg_bphs))*np.cos(thetavec)
# y1_vec_bphs = ycc_bphs + np.sin(np.radians(angle1_deg_bphs))*np.sin(thetavec)
# y2_vec_bphs = ycc_bphs + np.sin(np.radians(angle2_deg_bphs))*np.sin(thetavec)
# y3_vec_bphs = ycc_bphs + np.sin(np.radians(angle3_deg_bphs))*np.sin(thetavec)
# ax.plot(x1_vec_Potomac,y1_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Potomac,y2_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Potomac,y3_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_Hilda,y1_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Hilda,y2_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Hilda,y3_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_Schubart,y1_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Schubart,y2_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Schubart,y3_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_background,y1_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_background,y2_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_background,y3_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_astorb,y1_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_astorb,y2_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_astorb,y3_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_bphs,y1_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_bphs,y2_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_bphs,y3_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.scatter(xL,yL,color='cyan',s=12)
# ax.scatter(xJ,yJ,color='cyan',s=12)
# ax.scatter(xI,yI,color='magenta',s=12)
# ax.scatter(xcc_Potomac,ycc_Potomac,color='goldenrod',s=12)
# ax.scatter(xcc_Hilda,ycc_Hilda,color='red',s=12)
# ax.scatter(xcc_Schubart,ycc_Schubart,color='blue',s=12)
# ax.scatter(xcc_background,ycc_background,color='black',s=12)
# ax.scatter(xcc_astorb,ycc_astorb,color='green',s=12)
# ax.scatter(xcc_bphs,ycc_bphs,color='magenta',s=12)
# ax.text(xL,yL,'L')
# ax.text(xJ,yJ,'J')
# ax.text(xI,yI,'I')
# ax.text(xcc_Potomac,ycc_Potomac,'P')
# ax.text(xcc_Hilda,ycc_Hilda,'H')
# ax.text(xcc_Schubart,ycc_Schubart,'S')
# ax.text(xcc_background,ycc_background,'B')
# ax.text(xcc_astorb,ycc_astorb,'A')
# ax.text(xcc_bphs,ycc_bphs,'Z')
# ax.set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
# ax.set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
# ax.set_xlim(xlim)
# ax.set_ylim(ylim)
# # ax.set_xticks([-0.015,-0.010,-0.005,0,0.005])
# # ax.set_yticks([0.015,0.020,0.025,0.030])
# ax.set_box_aspect(1)
# ax.grid()
# xlim = [-0.016,0.007]
# xspan = xlim[1]-xlim[0]
# ymin = 0.006
# ylim = [ymin,ymin+xspan]
# ax = fig.add_subplot(122)
# lw_here = 2
# # for probmass in probmasses:
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Potomac,ycc_Potomac,siga_Potomac,sigb_Potomac,rhoab_Potomac,probmass,edgecolor='goldenrod',linestyle_in='dashdot',lw=lw_here)
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Hilda,ycc_Hilda,siga_Hilda,sigb_Hilda,rhoab_Hilda,probmass,edgecolor='red',linestyle_in='solid',lw=lw_here)
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Schubart,ycc_Schubart,siga_Schubart,sigb_Schubart,rhoab_Schubart,probmass,edgecolor='blue',linestyle_in='dotted',lw=lw_here)
# thetavec = np.linspace(start=0,stop=2*np.pi,num=1001,endpoint=False)
# x1_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle1_deg_Potomac))*np.cos(thetavec)
# x2_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle2_deg_Potomac))*np.cos(thetavec)
# x3_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle3_deg_Potomac))*np.cos(thetavec)
# y1_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle1_deg_Potomac))*np.sin(thetavec)
# y2_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle2_deg_Potomac))*np.sin(thetavec)
# y3_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle3_deg_Potomac))*np.sin(thetavec)
# x1_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle1_deg_Hilda))*np.cos(thetavec)
# x2_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle2_deg_Hilda))*np.cos(thetavec)
# x3_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle3_deg_Hilda))*np.cos(thetavec)
# y1_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle1_deg_Hilda))*np.sin(thetavec)
# y2_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle2_deg_Hilda))*np.sin(thetavec)
# y3_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle3_deg_Hilda))*np.sin(thetavec)
# x1_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle1_deg_Schubart))*np.cos(thetavec)
# x2_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle2_deg_Schubart))*np.cos(thetavec)
# x3_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle3_deg_Schubart))*np.cos(thetavec)
# y1_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle1_deg_Schubart))*np.sin(thetavec)
# y2_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle2_deg_Schubart))*np.sin(thetavec)
# y3_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle3_deg_Schubart))*np.sin(thetavec)
# x1_vec_background = xcc_background + np.sin(np.radians(angle1_deg_background))*np.cos(thetavec)
# x2_vec_background = xcc_background + np.sin(np.radians(angle2_deg_background))*np.cos(thetavec)
# x3_vec_background = xcc_background + np.sin(np.radians(angle3_deg_background))*np.cos(thetavec)
# y1_vec_background = ycc_background + np.sin(np.radians(angle1_deg_background))*np.sin(thetavec)
# y2_vec_background = ycc_background + np.sin(np.radians(angle2_deg_background))*np.sin(thetavec)
# y3_vec_background = ycc_background + np.sin(np.radians(angle3_deg_background))*np.sin(thetavec)
# x1_vec_astorb = xcc_astorb + np.sin(np.radians(angle1_deg_astorb))*np.cos(thetavec)
# x2_vec_astorb = xcc_astorb + np.sin(np.radians(angle2_deg_astorb))*np.cos(thetavec)
# x3_vec_astorb = xcc_astorb + np.sin(np.radians(angle3_deg_astorb))*np.cos(thetavec)
# y1_vec_astorb = ycc_astorb + np.sin(np.radians(angle1_deg_astorb))*np.sin(thetavec)
# y2_vec_astorb = ycc_astorb + np.sin(np.radians(angle2_deg_astorb))*np.sin(thetavec)
# y3_vec_astorb = ycc_astorb + np.sin(np.radians(angle3_deg_astorb))*np.sin(thetavec)
# x1_vec_bphs = xcc_bphs + np.sin(np.radians(angle1_deg_bphs))*np.cos(thetavec)
# x2_vec_bphs = xcc_bphs + np.sin(np.radians(angle2_deg_bphs))*np.cos(thetavec)
# x3_vec_bphs = xcc_bphs + np.sin(np.radians(angle3_deg_bphs))*np.cos(thetavec)
# y1_vec_bphs = ycc_bphs + np.sin(np.radians(angle1_deg_bphs))*np.sin(thetavec)
# y2_vec_bphs = ycc_bphs + np.sin(np.radians(angle2_deg_bphs))*np.sin(thetavec)
# y3_vec_bphs = ycc_bphs + np.sin(np.radians(angle3_deg_bphs))*np.sin(thetavec)
# ax.plot(x1_vec_Potomac,y1_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Potomac,y2_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Potomac,y3_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_Hilda,y1_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Hilda,y2_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Hilda,y3_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_Schubart,y1_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Schubart,y2_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Schubart,y3_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_background,y1_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_background,y2_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_background,y3_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_astorb,y1_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_astorb,y2_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_astorb,y3_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_bphs,y1_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_bphs,y2_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_bphs,y3_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.scatter(xL,yL,color='cyan',s=12)
# ax.scatter(xJ,yJ,color='cyan',s=12)
# ax.scatter(xI,yI,color='magenta',s=12)
# ax.scatter(xcc_Potomac,ycc_Potomac,color='goldenrod',s=12)
# ax.scatter(xcc_Hilda,ycc_Hilda,color='red',s=12)
# ax.scatter(xcc_Schubart,ycc_Schubart,color='blue',s=12)
# ax.scatter(xcc_background,ycc_background,color='black',s=12)
# ax.scatter(xcc_astorb,ycc_astorb,color='green',s=12)
# ax.scatter(xcc_bphs,ycc_bphs,color='magenta',s=12)
# ax.text(xL,yL,'L')
# ax.text(xJ,yJ,'J')
# ax.text(xI,yI,'I')
# ax.text(xcc_Potomac,ycc_Potomac,'P')
# ax.text(xcc_Hilda,ycc_Hilda,'H')
# ax.text(xcc_Schubart,ycc_Schubart,'S')
# ax.text(xcc_background,ycc_background,'B')
# ax.text(xcc_astorb,ycc_astorb,'A')
# ax.text(xcc_bphs,ycc_bphs,'Z')
# ax.set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
# # ax.set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
# ax.set_xlim(xlim)
# ax.set_ylim(ylim)
# # ax.set_xticks([-0.015,-0.010,-0.005,0,0.005])
# # ax.set_yticks([0.015,0.020,0.025,0.030])
# ax.set_box_aspect(1)
# ax.grid()
# plt.suptitle('Confidence ellipses for vmf')
# plt.tight_layout()
# plt.savefig(outfile_plot,dpi=300)
# plt.show()
# # fig.subplots_adjust(wspace=0, hspace=0)
# print('aei_all32s_Hlt15.7_librating')
#%%
# title = 'Confidence ellipses for vmf (wc)'
# outfile_plot = 'mean_pole_confidence_overlaps_start_may2025_Hlt15.7_vmf.png'
# xlim = [-0.04,0.04]
# xspan = xlim[1]-xlim[0]
# ymin = -0.02
# ylim = [ymin,ymin+xspan]
# df = pd.read_csv('a002_Potomac_Hlt15.7_family_horizons.csv')
# irad_Potomac = np.radians(np.array(df['incl_deg'].tolist()))
# Wrad_Potomac = np.radians(np.array(df['Omega_deg'].tolist()))
# x_Potomac = np.sin(irad_Potomac)*np.cos(Wrad_Potomac)
# y_Potomac = np.sin(irad_Potomac)*np.sin(Wrad_Potomac)
# df = pd.read_csv('a002_Hilda_Hlt15.7_family_horizons.csv')
# irad_Hilda = np.radians(np.array(df['incl_deg'].tolist()))
# Wrad_Hilda = np.radians(np.array(df['Omega_deg'].tolist()))
# x_Hilda = np.sin(irad_Hilda)*np.cos(Wrad_Hilda)
# y_Hilda = np.sin(irad_Hilda)*np.sin(Wrad_Hilda)
# df = pd.read_csv('a002_Schubart_Hlt15.7_family_horizons.csv')
# irad_Schubart = np.radians(np.array(df['incl_deg'].tolist()))
# Wrad_Schubart = np.radians(np.array(df['Omega_deg'].tolist()))
# x_Schubart = np.sin(irad_Schubart)*np.cos(Wrad_Schubart)
# y_Schubart = np.sin(irad_Schubart)*np.sin(Wrad_Schubart)
# df = pd.read_csv('a002_background_Hlt15.7_family_horizons.csv')
# irad_background = np.radians(np.array(df['incl_deg'].tolist()))
# Wrad_background = np.radians(np.array(df['Omega_deg'].tolist()))
# x_background = np.sin(irad_background)*np.cos(Wrad_background)
# y_background = np.sin(irad_background)*np.sin(Wrad_background)
# df = pd.read_csv('a002_astorb_Hlt15.7_family_horizons.csv')
# irad_astorb = np.radians(np.array(df['incl_deg'].tolist()))
# Wrad_astorb = np.radians(np.array(df['Omega_deg'].tolist()))
# x_astorb = np.sin(irad_astorb)*np.cos(Wrad_astorb)
# y_astorb = np.sin(irad_astorb)*np.sin(Wrad_astorb)
# xccvec_Potomac = np.loadtxt('Potomac_xccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# yccvec_Potomac = np.loadtxt('Potomac_yccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle1_deg_vec_Potomac = np.loadtxt('Potomac_angle1_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle2_deg_vec_Potomac = np.loadtxt('Potomac_angle2_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle3_deg_vec_Potomac = np.loadtxt('Potomac_angle3_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# xcc_Potomac = xccvec_Potomac[0]
# ycc_Potomac = yccvec_Potomac[0]
# angle1_deg_Potomac = angle1_deg_vec_Potomac[0]
# angle2_deg_Potomac = angle2_deg_vec_Potomac[0]
# angle3_deg_Potomac = angle3_deg_vec_Potomac[0]
# xccvec_Hilda = np.loadtxt('Hilda_xccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# yccvec_Hilda = np.loadtxt('Hilda_yccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle1_deg_vec_Hilda = np.loadtxt('Hilda_angle1_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle2_deg_vec_Hilda = np.loadtxt('Hilda_angle2_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle3_deg_vec_Hilda = np.loadtxt('Hilda_angle3_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# xcc_Hilda = xccvec_Hilda[0]
# ycc_Hilda = yccvec_Hilda[0]
# angle1_deg_Hilda = angle1_deg_vec_Hilda[0]
# angle2_deg_Hilda = angle2_deg_vec_Hilda[0]
# angle3_deg_Hilda = angle3_deg_vec_Hilda[0]
# xccvec_Schubart = np.loadtxt('Schubart_xccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# yccvec_Schubart = np.loadtxt('Schubart_yccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle1_deg_vec_Schubart = np.loadtxt('Schubart_angle1_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle2_deg_vec_Schubart = np.loadtxt('Schubart_angle2_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle3_deg_vec_Schubart = np.loadtxt('Schubart_angle3_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# xcc_Schubart = xccvec_Schubart[0]
# ycc_Schubart = yccvec_Schubart[0]
# angle1_deg_Schubart = angle1_deg_vec_Schubart[0]
# angle2_deg_Schubart = angle2_deg_vec_Schubart[0]
# angle3_deg_Schubart = angle3_deg_vec_Schubart[0]
# xccvec_background = np.loadtxt('background_xccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# yccvec_background = np.loadtxt('background_yccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle1_deg_vec_background = np.loadtxt('background_angle1_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle2_deg_vec_background = np.loadtxt('background_angle2_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle3_deg_vec_background = np.loadtxt('background_angle3_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# xcc_background = xccvec_background[0]
# ycc_background = yccvec_background[0]
# angle1_deg_background = angle1_deg_vec_background[0]
# angle2_deg_background = angle2_deg_vec_background[0]
# angle3_deg_background = angle3_deg_vec_background[0]
# xccvec_astorb = np.loadtxt('astorb_xccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# yccvec_astorb = np.loadtxt('astorb_yccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle1_deg_vec_astorb = np.loadtxt('astorb_angle1_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle2_deg_vec_astorb = np.loadtxt('astorb_angle2_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle3_deg_vec_astorb = np.loadtxt('astorb_angle3_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# xcc_astorb = xccvec_astorb[0]
# ycc_astorb = yccvec_astorb[0]
# angle1_deg_astorb = angle1_deg_vec_astorb[0]
# angle2_deg_astorb = angle2_deg_vec_astorb[0]
# angle3_deg_astorb = angle3_deg_vec_astorb[0]
# xccvec_bphs = np.loadtxt('bphs_xccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# yccvec_bphs = np.loadtxt('bphs_yccvec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle1_deg_vec_bphs = np.loadtxt('bphs_angle1_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle2_deg_vec_bphs = np.loadtxt('bphs_angle2_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# angle3_deg_vec_bphs = np.loadtxt('bphs_angle3_deg_vec_Hlt15.7_vmf_2e6yr_1e2yr_0.2yr_jd2460796.5.csv',delimiter=',')
# xcc_bphs = xccvec_bphs[0]
# ycc_bphs = yccvec_bphs[0]
# angle1_deg_bphs = angle1_deg_vec_bphs[0]
# angle2_deg_bphs = angle2_deg_vec_bphs[0]
# angle3_deg_bphs = angle3_deg_vec_bphs[0]
# fig = plt.figure(figsize=(7,7))
# plt.rcParams['font.size'] = 10
# fig,ax = plt.subplots(1,1)
# lw_here = 0.5
# # for probmass in probmasses:
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Potomac,ycc_Potomac,siga_Potomac,sigb_Potomac,rhoab_Potomac,probmass,edgecolor='goldenrod',linestyle_in='dashdot',lw=lw_here)
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Hilda,ycc_Hilda,siga_Hilda,sigb_Hilda,rhoab_Hilda,probmass,edgecolor='red',linestyle_in='solid',lw=lw_here)
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Schubart,ycc_Schubart,siga_Schubart,sigb_Schubart,rhoab_Schubart,probmass,edgecolor='blue',linestyle_in='dotted',lw=lw_here)
# thetavec = np.linspace(start=0,stop=2*np.pi,num=1001,endpoint=False)
# x1_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle1_deg_Potomac))*np.cos(thetavec)
# x2_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle2_deg_Potomac))*np.cos(thetavec)
# x3_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle3_deg_Potomac))*np.cos(thetavec)
# y1_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle1_deg_Potomac))*np.sin(thetavec)
# y2_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle2_deg_Potomac))*np.sin(thetavec)
# y3_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle3_deg_Potomac))*np.sin(thetavec)
# x1_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle1_deg_Hilda))*np.cos(thetavec)
# x2_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle2_deg_Hilda))*np.cos(thetavec)
# x3_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle3_deg_Hilda))*np.cos(thetavec)
# y1_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle1_deg_Hilda))*np.sin(thetavec)
# y2_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle2_deg_Hilda))*np.sin(thetavec)
# y3_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle3_deg_Hilda))*np.sin(thetavec)
# x1_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle1_deg_Schubart))*np.cos(thetavec)
# x2_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle2_deg_Schubart))*np.cos(thetavec)
# x3_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle3_deg_Schubart))*np.cos(thetavec)
# y1_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle1_deg_Schubart))*np.sin(thetavec)
# y2_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle2_deg_Schubart))*np.sin(thetavec)
# y3_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle3_deg_Schubart))*np.sin(thetavec)
# x1_vec_background = xcc_background + np.sin(np.radians(angle1_deg_background))*np.cos(thetavec)
# x2_vec_background = xcc_background + np.sin(np.radians(angle2_deg_background))*np.cos(thetavec)
# x3_vec_background = xcc_background + np.sin(np.radians(angle3_deg_background))*np.cos(thetavec)
# y1_vec_background = ycc_background + np.sin(np.radians(angle1_deg_background))*np.sin(thetavec)
# y2_vec_background = ycc_background + np.sin(np.radians(angle2_deg_background))*np.sin(thetavec)
# y3_vec_background = ycc_background + np.sin(np.radians(angle3_deg_background))*np.sin(thetavec)
# x1_vec_astorb = xcc_astorb + np.sin(np.radians(angle1_deg_astorb))*np.cos(thetavec)
# x2_vec_astorb = xcc_astorb + np.sin(np.radians(angle2_deg_astorb))*np.cos(thetavec)
# x3_vec_astorb = xcc_astorb + np.sin(np.radians(angle3_deg_astorb))*np.cos(thetavec)
# y1_vec_astorb = ycc_astorb + np.sin(np.radians(angle1_deg_astorb))*np.sin(thetavec)
# y2_vec_astorb = ycc_astorb + np.sin(np.radians(angle2_deg_astorb))*np.sin(thetavec)
# y3_vec_astorb = ycc_astorb + np.sin(np.radians(angle3_deg_astorb))*np.sin(thetavec)
# x1_vec_bphs = xcc_bphs + np.sin(np.radians(angle1_deg_bphs))*np.cos(thetavec)
# x2_vec_bphs = xcc_bphs + np.sin(np.radians(angle2_deg_bphs))*np.cos(thetavec)
# x3_vec_bphs = xcc_bphs + np.sin(np.radians(angle3_deg_bphs))*np.cos(thetavec)
# y1_vec_bphs = ycc_bphs + np.sin(np.radians(angle1_deg_bphs))*np.sin(thetavec)
# y2_vec_bphs = ycc_bphs + np.sin(np.radians(angle2_deg_bphs))*np.sin(thetavec)
# y3_vec_bphs = ycc_bphs + np.sin(np.radians(angle3_deg_bphs))*np.sin(thetavec)
# ax.plot(x1_vec_Potomac,y1_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Potomac,y2_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Potomac,y3_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_Hilda,y1_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Hilda,y2_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Hilda,y3_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_Schubart,y1_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Schubart,y2_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Schubart,y3_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_background,y1_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_background,y2_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_background,y3_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_astorb,y1_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_astorb,y2_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_astorb,y3_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_bphs,y1_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_bphs,y2_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_bphs,y3_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.scatter(xL,yL,color='cyan',s=12)
# ax.scatter(xJ,yJ,color='cyan',s=12)
# ax.scatter(xI,yI,color='magenta',s=12)
# ax.scatter(xcc_Potomac,ycc_Potomac,color='goldenrod',s=12)
# ax.scatter(xcc_Hilda,ycc_Hilda,color='red',s=12)
# ax.scatter(xcc_Schubart,ycc_Schubart,color='blue',s=12)
# ax.scatter(xcc_background,ycc_background,color='black',s=12)
# ax.scatter(xcc_astorb,ycc_astorb,color='green',s=12)
# ax.scatter(xcc_bphs,ycc_bphs,color='magenta',s=12)
# ax.text(xL,yL,'L')
# ax.text(xJ,yJ,'J')
# ax.text(xI,yI,'I')
# ax.text(xcc_Potomac,ycc_Potomac,'P')
# ax.text(xcc_Hilda,ycc_Hilda,'H')
# ax.text(xcc_Schubart,ycc_Schubart,'S')
# ax.text(xcc_background,ycc_background,'B')
# ax.text(xcc_astorb,ycc_astorb,'A')
# ax.text(xcc_bphs,ycc_bphs,'Z')
# ax.set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
# ax.set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
# ax.set_xlim(xlim)
# ax.set_ylim(ylim)
# ax.set_title(title)
# # ax.set_xticks([-0.015,-0.010,-0.005,0,0.005])
# # ax.set_yticks([0.015,0.020,0.025,0.030])
# ax.set_box_aspect(1)
# ax.grid()
# plt.savefig(outfile_plot,dpi=400)
# plt.tight_layout()
#%%
# outfile_plot = 'mean_pole_confidence_overlaps_start_may2025_Hlt15.7_2panels_vmf.png'
# xlim = [-0.04,0.04]
# xspan = xlim[1]-xlim[0]
# ymin = -0.02
# ylim = [ymin,ymin+xspan]
# fig = plt.figure(figsize=(6,3))
# plt.rcParams['font.size'] = 10
# ax = fig.add_subplot(121)
# lw_here = 1.5
# # for probmass in probmasses:
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Potomac,ycc_Potomac,siga_Potomac,sigb_Potomac,rhoab_Potomac,probmass,edgecolor='goldenrod',linestyle_in='dashdot',lw=lw_here)
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Hilda,ycc_Hilda,siga_Hilda,sigb_Hilda,rhoab_Hilda,probmass,edgecolor='red',linestyle_in='solid',lw=lw_here)
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Schubart,ycc_Schubart,siga_Schubart,sigb_Schubart,rhoab_Schubart,probmass,edgecolor='blue',linestyle_in='dotted',lw=lw_here)
# thetavec = np.linspace(start=0,stop=2*np.pi,num=1001,endpoint=False)
# x1_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle1_deg_Potomac))*np.cos(thetavec)
# x2_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle2_deg_Potomac))*np.cos(thetavec)
# x3_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle3_deg_Potomac))*np.cos(thetavec)
# y1_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle1_deg_Potomac))*np.sin(thetavec)
# y2_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle2_deg_Potomac))*np.sin(thetavec)
# y3_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle3_deg_Potomac))*np.sin(thetavec)
# x1_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle1_deg_Hilda))*np.cos(thetavec)
# x2_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle2_deg_Hilda))*np.cos(thetavec)
# x3_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle3_deg_Hilda))*np.cos(thetavec)
# y1_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle1_deg_Hilda))*np.sin(thetavec)
# y2_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle2_deg_Hilda))*np.sin(thetavec)
# y3_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle3_deg_Hilda))*np.sin(thetavec)
# x1_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle1_deg_Schubart))*np.cos(thetavec)
# x2_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle2_deg_Schubart))*np.cos(thetavec)
# x3_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle3_deg_Schubart))*np.cos(thetavec)
# y1_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle1_deg_Schubart))*np.sin(thetavec)
# y2_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle2_deg_Schubart))*np.sin(thetavec)
# y3_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle3_deg_Schubart))*np.sin(thetavec)
# x1_vec_background = xcc_background + np.sin(np.radians(angle1_deg_background))*np.cos(thetavec)
# x2_vec_background = xcc_background + np.sin(np.radians(angle2_deg_background))*np.cos(thetavec)
# x3_vec_background = xcc_background + np.sin(np.radians(angle3_deg_background))*np.cos(thetavec)
# y1_vec_background = ycc_background + np.sin(np.radians(angle1_deg_background))*np.sin(thetavec)
# y2_vec_background = ycc_background + np.sin(np.radians(angle2_deg_background))*np.sin(thetavec)
# y3_vec_background = ycc_background + np.sin(np.radians(angle3_deg_background))*np.sin(thetavec)
# x1_vec_astorb = xcc_astorb + np.sin(np.radians(angle1_deg_astorb))*np.cos(thetavec)
# x2_vec_astorb = xcc_astorb + np.sin(np.radians(angle2_deg_astorb))*np.cos(thetavec)
# x3_vec_astorb = xcc_astorb + np.sin(np.radians(angle3_deg_astorb))*np.cos(thetavec)
# y1_vec_astorb = ycc_astorb + np.sin(np.radians(angle1_deg_astorb))*np.sin(thetavec)
# y2_vec_astorb = ycc_astorb + np.sin(np.radians(angle2_deg_astorb))*np.sin(thetavec)
# y3_vec_astorb = ycc_astorb + np.sin(np.radians(angle3_deg_astorb))*np.sin(thetavec)
# x1_vec_bphs = xcc_bphs + np.sin(np.radians(angle1_deg_bphs))*np.cos(thetavec)
# x2_vec_bphs = xcc_bphs + np.sin(np.radians(angle2_deg_bphs))*np.cos(thetavec)
# x3_vec_bphs = xcc_bphs + np.sin(np.radians(angle3_deg_bphs))*np.cos(thetavec)
# y1_vec_bphs = ycc_bphs + np.sin(np.radians(angle1_deg_bphs))*np.sin(thetavec)
# y2_vec_bphs = ycc_bphs + np.sin(np.radians(angle2_deg_bphs))*np.sin(thetavec)
# y3_vec_bphs = ycc_bphs + np.sin(np.radians(angle3_deg_bphs))*np.sin(thetavec)
# ax.plot(x1_vec_Potomac,y1_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Potomac,y2_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Potomac,y3_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_Hilda,y1_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Hilda,y2_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Hilda,y3_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_Schubart,y1_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Schubart,y2_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Schubart,y3_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_background,y1_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_background,y2_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_background,y3_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_astorb,y1_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_astorb,y2_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_astorb,y3_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_bphs,y1_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_bphs,y2_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_bphs,y3_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.scatter(xL,yL,color='cyan',s=12)
# ax.scatter(xJ,yJ,color='cyan',s=12)
# ax.scatter(xI,yI,color='magenta',s=12)
# ax.scatter(xcc_Potomac,ycc_Potomac,color='goldenrod',s=12)
# ax.scatter(xcc_Hilda,ycc_Hilda,color='red',s=12)
# ax.scatter(xcc_Schubart,ycc_Schubart,color='blue',s=12)
# ax.scatter(xcc_background,ycc_background,color='black',s=12)
# ax.scatter(xcc_astorb,ycc_astorb,color='green',s=12)
# ax.scatter(xcc_bphs,ycc_bphs,color='magenta',s=12)
# ax.text(xL,yL,'L')
# ax.text(xJ,yJ,'J')
# ax.text(xI,yI,'I')
# ax.text(xcc_Potomac,ycc_Potomac,'P')
# ax.text(xcc_Hilda,ycc_Hilda,'H')
# ax.text(xcc_Schubart,ycc_Schubart,'S')
# ax.text(xcc_background,ycc_background,'B')
# ax.text(xcc_astorb,ycc_astorb,'A')
# ax.text(xcc_bphs,ycc_bphs,'Z')
# ax.set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
# ax.set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
# ax.set_xlim(xlim)
# ax.set_ylim(ylim)
# # ax.set_xticks([-0.015,-0.010,-0.005,0,0.005])
# # ax.set_yticks([0.015,0.020,0.025,0.030])
# ax.set_box_aspect(1)
# ax.grid()
# xlim = [-0.016,0.007]
# xspan = xlim[1]-xlim[0]
# ymin = 0.006
# ylim = [ymin,ymin+xspan]
# ax = fig.add_subplot(122)
# lw_here = 2
# # for probmass in probmasses:
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Potomac,ycc_Potomac,siga_Potomac,sigb_Potomac,rhoab_Potomac,probmass,edgecolor='goldenrod',linestyle_in='dashdot',lw=lw_here)
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Hilda,ycc_Hilda,siga_Hilda,sigb_Hilda,rhoab_Hilda,probmass,edgecolor='red',linestyle_in='solid',lw=lw_here)
# #     ax = plot_fitted_circles_covariance_3(ax,xcc_Schubart,ycc_Schubart,siga_Schubart,sigb_Schubart,rhoab_Schubart,probmass,edgecolor='blue',linestyle_in='dotted',lw=lw_here)
# thetavec = np.linspace(start=0,stop=2*np.pi,num=1001,endpoint=False)
# x1_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle1_deg_Potomac))*np.cos(thetavec)
# x2_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle2_deg_Potomac))*np.cos(thetavec)
# x3_vec_Potomac = xcc_Potomac + np.sin(np.radians(angle3_deg_Potomac))*np.cos(thetavec)
# y1_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle1_deg_Potomac))*np.sin(thetavec)
# y2_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle2_deg_Potomac))*np.sin(thetavec)
# y3_vec_Potomac = ycc_Potomac + np.sin(np.radians(angle3_deg_Potomac))*np.sin(thetavec)
# x1_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle1_deg_Hilda))*np.cos(thetavec)
# x2_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle2_deg_Hilda))*np.cos(thetavec)
# x3_vec_Hilda = xcc_Hilda + np.sin(np.radians(angle3_deg_Hilda))*np.cos(thetavec)
# y1_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle1_deg_Hilda))*np.sin(thetavec)
# y2_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle2_deg_Hilda))*np.sin(thetavec)
# y3_vec_Hilda = ycc_Hilda + np.sin(np.radians(angle3_deg_Hilda))*np.sin(thetavec)
# x1_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle1_deg_Schubart))*np.cos(thetavec)
# x2_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle2_deg_Schubart))*np.cos(thetavec)
# x3_vec_Schubart = xcc_Schubart + np.sin(np.radians(angle3_deg_Schubart))*np.cos(thetavec)
# y1_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle1_deg_Schubart))*np.sin(thetavec)
# y2_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle2_deg_Schubart))*np.sin(thetavec)
# y3_vec_Schubart = ycc_Schubart + np.sin(np.radians(angle3_deg_Schubart))*np.sin(thetavec)
# x1_vec_background = xcc_background + np.sin(np.radians(angle1_deg_background))*np.cos(thetavec)
# x2_vec_background = xcc_background + np.sin(np.radians(angle2_deg_background))*np.cos(thetavec)
# x3_vec_background = xcc_background + np.sin(np.radians(angle3_deg_background))*np.cos(thetavec)
# y1_vec_background = ycc_background + np.sin(np.radians(angle1_deg_background))*np.sin(thetavec)
# y2_vec_background = ycc_background + np.sin(np.radians(angle2_deg_background))*np.sin(thetavec)
# y3_vec_background = ycc_background + np.sin(np.radians(angle3_deg_background))*np.sin(thetavec)
# x1_vec_astorb = xcc_astorb + np.sin(np.radians(angle1_deg_astorb))*np.cos(thetavec)
# x2_vec_astorb = xcc_astorb + np.sin(np.radians(angle2_deg_astorb))*np.cos(thetavec)
# x3_vec_astorb = xcc_astorb + np.sin(np.radians(angle3_deg_astorb))*np.cos(thetavec)
# y1_vec_astorb = ycc_astorb + np.sin(np.radians(angle1_deg_astorb))*np.sin(thetavec)
# y2_vec_astorb = ycc_astorb + np.sin(np.radians(angle2_deg_astorb))*np.sin(thetavec)
# y3_vec_astorb = ycc_astorb + np.sin(np.radians(angle3_deg_astorb))*np.sin(thetavec)
# x1_vec_bphs = xcc_bphs + np.sin(np.radians(angle1_deg_bphs))*np.cos(thetavec)
# x2_vec_bphs = xcc_bphs + np.sin(np.radians(angle2_deg_bphs))*np.cos(thetavec)
# x3_vec_bphs = xcc_bphs + np.sin(np.radians(angle3_deg_bphs))*np.cos(thetavec)
# y1_vec_bphs = ycc_bphs + np.sin(np.radians(angle1_deg_bphs))*np.sin(thetavec)
# y2_vec_bphs = ycc_bphs + np.sin(np.radians(angle2_deg_bphs))*np.sin(thetavec)
# y3_vec_bphs = ycc_bphs + np.sin(np.radians(angle3_deg_bphs))*np.sin(thetavec)
# ax.plot(x1_vec_Potomac,y1_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Potomac,y2_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Potomac,y3_vec_Potomac,color='goldenrod',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_Hilda,y1_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Hilda,y2_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Hilda,y3_vec_Hilda,color='red',linestyle='solid',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_Schubart,y1_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_Schubart,y2_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_Schubart,y3_vec_Schubart,color='blue',linestyle='dotted',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_background,y1_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_background,y2_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_background,y3_vec_background,color='black',linestyle='dashed',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_astorb,y1_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_astorb,y2_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_astorb,y3_vec_astorb,color='green',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x1_vec_bphs,y1_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x2_vec_bphs,y2_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.plot(x3_vec_bphs,y3_vec_bphs,color='magenta',linestyle='dashdot',lw=lw_here,alpha=0.5)
# ax.scatter(xL,yL,color='cyan',s=12)
# ax.scatter(xJ,yJ,color='cyan',s=12)
# ax.scatter(xI,yI,color='magenta',s=12)
# ax.scatter(xcc_Potomac,ycc_Potomac,color='goldenrod',s=12)
# ax.scatter(xcc_Hilda,ycc_Hilda,color='red',s=12)
# ax.scatter(xcc_Schubart,ycc_Schubart,color='blue',s=12)
# ax.scatter(xcc_background,ycc_background,color='black',s=12)
# ax.scatter(xcc_astorb,ycc_astorb,color='green',s=12)
# ax.scatter(xcc_bphs,ycc_bphs,color='magenta',s=12)
# ax.text(xL,yL,'L')
# ax.text(xJ,yJ,'J')
# ax.text(xI,yI,'I')
# ax.text(xcc_Potomac,ycc_Potomac,'P')
# ax.text(xcc_Hilda,ycc_Hilda,'H')
# ax.text(xcc_Schubart,ycc_Schubart,'S')
# ax.text(xcc_background,ycc_background,'B')
# ax.text(xcc_astorb,ycc_astorb,'A')
# ax.text(xcc_bphs,ycc_bphs,'Z')
# ax.set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
# # ax.set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
# ax.set_xlim(xlim)
# ax.set_ylim(ylim)
# # ax.set_xticks([-0.015,-0.010,-0.005,0,0.005])
# # ax.set_yticks([0.015,0.020,0.025,0.030])
# ax.set_box_aspect(1)
# ax.grid()
# plt.suptitle('Confidence ellipses for vmf (wc)')
# plt.tight_layout()
# plt.savefig(outfile_plot,dpi=300)
# plt.show()
# # fig.subplots_adjust(wspace=0, hspace=0)
# print('aei_all32s_Hlt15.7')