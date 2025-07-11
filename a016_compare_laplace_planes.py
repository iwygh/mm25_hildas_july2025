import numpy as np
import pandas as pd

df = pd.read_csv('b010_laplaceplane_P_L_2e6yr_2e2yr_0.2yr_jd2460796.5.csv')
laplacexvec_potomacs = np.array(df['laplacex'].to_list())
laplaceyvec_potomacs = np.array(df['laplacey'].to_list())
df = pd.read_csv('b010_laplaceplane_H_L_2e6yr_2e2yr_0.2yr_jd2460796.5.csv')
laplacexvec_hildacolls = np.array(df['laplacex'].to_list())
laplaceyvec_hildacolls = np.array(df['laplacey'].to_list())
df = pd.read_csv('b010_laplaceplane_S_L_2e6yr_2e2yr_0.2yr_jd2460796.5.csv')
laplacexvec_schubarts = np.array(df['laplacex'].to_list())
laplaceyvec_schubarts = np.array(df['laplacey'].to_list())
df = pd.read_csv('b010_laplaceplane_BFG_LMC_2e6yr_2e2yr_0.2yr_jd2460796.5.csv')
laplacexvec_hildabackgrounds = np.array(df['laplacex'].to_list())
laplaceyvec_hildabackgrounds = np.array(df['laplacey'].to_list())
df = pd.read_csv('b010_laplaceplane_Potomac_clones_2e6yr_2e2yr_0.2yr_jd2460796.5.csv')
laplacexvec_potomacclones = np.array(df['laplacex'].to_list())
laplaceyvec_potomacclones = np.array(df['laplacey'].to_list())
df = pd.read_csv('b010_laplaceplane_Hilda_clones_2e6yr_2e2yr_0.2yr_jd2460796.5.csv')
laplacexvec_hildaclones = np.array(df['laplacex'].to_list())
laplaceyvec_hildaclones = np.array(df['laplacey'].to_list())
df = pd.read_csv('b010_laplaceplane_Schubart_clones_2e6yr_2e2yr_0.2yr_jd2460796.5.csv')
laplacexvec_schubartclones = np.array(df['laplacex'].to_list())
laplaceyvec_schubartclones = np.array(df['laplacey'].to_list())
df = pd.read_csv('b010_laplaceplane_Ismene_clones_2e6yr_2e2yr_0.2yr_jd2460796.5.csv')
laplacexvec_ismeneclones = np.array(df['laplacex'].to_list())
laplaceyvec_ismeneclones = np.array(df['laplacey'].to_list())


potomacs_hildacolls = np.degrees(np.max(np.sqrt(\
    (laplacexvec_potomacs-laplacexvec_hildacolls)**2+\
    (laplaceyvec_potomacs-laplaceyvec_hildacolls)**2)))
potomacs_schubarts = np.degrees(np.max(np.sqrt(\
    (laplacexvec_potomacs-laplacexvec_schubarts)**2+\
    (laplaceyvec_potomacs-laplaceyvec_schubarts)**2)))
potomacs_hildabackgrounds = np.degrees(np.max(np.sqrt(\
    (laplacexvec_potomacs-laplacexvec_hildabackgrounds)**2+\
    (laplaceyvec_potomacs-laplaceyvec_hildabackgrounds)**2)))
potomacs_potomacclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_potomacs-laplacexvec_potomacclones)**2+\
    (laplaceyvec_potomacs-laplaceyvec_potomacclones)**2)))
potomacs_hildaclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_potomacs-laplacexvec_hildaclones)**2+\
    (laplaceyvec_potomacs-laplaceyvec_hildaclones)**2)))
potomacs_schubartclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_potomacs-laplacexvec_schubartclones)**2+\
    (laplaceyvec_potomacs-laplaceyvec_schubartclones)**2)))
potomacs_ismeneclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_potomacs-laplacexvec_ismeneclones)**2+\
    (laplaceyvec_potomacs-laplaceyvec_ismeneclones)**2)))

hildacolls_schubarts = np.degrees(np.max(np.sqrt(\
    (laplacexvec_hildacolls-laplacexvec_schubarts)**2+\
    (laplaceyvec_hildacolls-laplaceyvec_schubarts)**2)))
hildacolls_hildabackgrounds = np.degrees(np.max(np.sqrt(\
    (laplacexvec_hildacolls-laplacexvec_hildabackgrounds)**2+\
    (laplaceyvec_hildacolls-laplaceyvec_hildabackgrounds)**2)))
hildacolls_potomacclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_hildacolls-laplacexvec_potomacclones)**2+\
    (laplaceyvec_hildacolls-laplaceyvec_potomacclones)**2)))
hildacolls_hildaclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_hildacolls-laplacexvec_hildaclones)**2+\
    (laplaceyvec_hildacolls-laplaceyvec_hildaclones)**2)))
hildacolls_schubartclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_hildacolls-laplacexvec_schubartclones)**2+\
    (laplaceyvec_hildacolls-laplaceyvec_schubartclones)**2)))
hildacolls_ismeneclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_hildacolls-laplacexvec_ismeneclones)**2+\
    (laplaceyvec_hildacolls-laplaceyvec_ismeneclones)**2)))
    
schubarts_hildabackgrounds = np.degrees(np.max(np.sqrt(\
    (laplacexvec_schubarts-laplacexvec_hildabackgrounds)**2+\
    (laplaceyvec_schubarts-laplaceyvec_hildabackgrounds)**2)))
schubarts_potomacclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_schubarts-laplacexvec_potomacclones)**2+\
    (laplaceyvec_schubarts-laplaceyvec_potomacclones)**2)))
schubarts_hildaclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_schubarts-laplacexvec_hildaclones)**2+\
    (laplaceyvec_schubarts-laplaceyvec_hildaclones)**2)))
schubarts_schubartclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_schubarts-laplacexvec_schubartclones)**2+\
    (laplaceyvec_schubarts-laplaceyvec_schubartclones)**2)))
schubarts_ismeneclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_schubarts-laplacexvec_ismeneclones)**2+\
    (laplaceyvec_schubarts-laplaceyvec_ismeneclones)**2)))
    
hildabackgrounds_potomacclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_hildabackgrounds-laplacexvec_potomacclones)**2+\
    (laplaceyvec_hildabackgrounds-laplaceyvec_potomacclones)**2)))
hildabackgrounds_hildaclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_hildabackgrounds-laplacexvec_hildaclones)**2+\
    (laplaceyvec_hildabackgrounds-laplaceyvec_hildaclones)**2)))
hildabackgrounds_schubartclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_hildabackgrounds-laplacexvec_schubartclones)**2+\
    (laplaceyvec_hildabackgrounds-laplaceyvec_schubartclones)**2)))
hildabackgrounds_ismeneclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_hildabackgrounds-laplacexvec_ismeneclones)**2+\
    (laplaceyvec_hildabackgrounds-laplaceyvec_ismeneclones)**2)))
    
potomacclones_hildaclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_potomacclones-laplacexvec_hildaclones)**2+\
    (laplaceyvec_potomacclones-laplaceyvec_hildaclones)**2)))
potomacclones_schubartclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_potomacclones-laplacexvec_schubartclones)**2+\
    (laplaceyvec_potomacclones-laplaceyvec_schubartclones)**2)))
potomacclones_ismeneclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_potomacclones-laplacexvec_ismeneclones)**2+\
    (laplaceyvec_potomacclones-laplaceyvec_ismeneclones)**2)))
    
hildaclones_schubartclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_hildaclones-laplacexvec_schubartclones)**2+\
    (laplaceyvec_hildaclones-laplaceyvec_schubartclones)**2)))
hildaclones_ismeneclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_hildaclones-laplacexvec_ismeneclones)**2+\
    (laplaceyvec_hildaclones-laplaceyvec_ismeneclones)**2)))
    
schubartclones_ismeneclones = np.degrees(np.max(np.sqrt(\
    (laplacexvec_schubartclones-laplacexvec_ismeneclones)**2+\
    (laplaceyvec_schubartclones-laplaceyvec_ismeneclones)**2)))
    
    
distance_list = [potomacs_hildacolls,potomacs_schubarts,potomacs_hildabackgrounds,\
                     potomacs_potomacclones,potomacs_hildaclones,potomacs_schubartclones,potomacs_ismeneclones,\
                hildacolls_schubarts,hildacolls_hildabackgrounds,hildacolls_potomacclones,hildacolls_hildaclones,\
                     hildacolls_schubartclones,hildacolls_ismeneclones,\
                schubarts_hildabackgrounds,schubarts_potomacclones,schubarts_hildaclones,schubarts_schubartclones,\
                    schubarts_ismeneclones,\
                hildabackgrounds_potomacclones,hildabackgrounds_hildaclones,hildabackgrounds_schubartclones,\
                    hildabackgrounds_ismeneclones,\
                potomacclones_hildaclones,potomacclones_schubartclones,potomacclones_ismeneclones,\
                hildaclones_schubartclones,hildaclones_ismeneclones,\
                schubartclones_ismeneclones]
print(np.max(np.array(distance_list)))
