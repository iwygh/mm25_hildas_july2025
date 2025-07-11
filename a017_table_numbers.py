import pandas as pd
import numpy as np
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
it = 0
for ilabel in range(nlabels):
    label = output_labels[ilabel]
    df = pd.read_csv('b007_'+label+'_statistics_vmf_'+yrsstr+'.csv')
    xcc = df['xcc'][it]
    ycc = df['ycc'][it]
    angle68 = df['angle_68_deg'][it]
    angle95 = df['angle_95_deg'][it]
    angle997 = df['angle_997_deg'][it]
    icc = np.degrees(np.arcsin(np.sqrt(xcc**2+ycc**2)))
    sini = np.sin(np.radians(icc))
    Wcc = np.degrees(np.arctan2(ycc/sini,xcc/sini))
    print(label)
    print('i0 deg = ',icc)
    print('W0 deg = ',Wcc)
    print('q0 = ',xcc)
    print('p0 = ',ycc)
    print('angle_68_deg = ',angle68)
    print('angle_95_deg = ',angle95)
    print('angle_997_deg = ',angle997)
    print('')