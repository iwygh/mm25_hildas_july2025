#%%
import numpy as np
import pandas as pd
cl = [   3,   67,  145,  129,  174,  235,  342,  416,  404,  650,\
       743,  739,  954, 1143, 1125, 1118, 1198, 1197, 1172, 1291,\
      1319, 1372, 1360, 1521, 1508, 1590, 1690, 1802, 1891, 1915,\
      2027, 2024, 2098, 2109, 2269, 2399, 2447, 2412, 2477, 2515,\
      2750, 2736, 2794, 2786, 2886, 2883, 2941, 2931, 2971, 2967,\
      2964, 3271, 3265, 3384, 3490, 3534, 3578, 3681, 3791, 3870,\
      3950, 4152, 4387, 4356, 4353, 4477, 4613, 4666, 4783, 4812,\
      4939, 4999, 5031, 5090, 5060, 5168, 5159, 5158, 5257, 5396,\
      5367, 5529, 5630, 5602, 5739, 5734, 5789, 5837, 5826, 5948,\
      5905, 6046, 6059, 6117, 6208, 6327, 6357,\
     ] # circulating objects - numbers start at 1 (like matlab), not 0 (pythonic)
ml = [ 112,  226,  210,  381,  359,  583,  555,  644,  676,  660,\
       849,  852,  944, 1095, 1300, 1272, 1324, 1301, 1586, 1556,\
      1665, 1821, 1861, 1860, 2040, 2037, 2028, 2019, 2017, 2361,\
      2369, 2514, 2620, 2973, 2958, 3455, 3875, 3855, 4083, 4112,\
      4317, 4632, 4880, 5148, 5147, 5146, 5165, 5443, 5498, 5531,\
      5516, 5593, 5628, 5698, 5723, 5860, 5926, 5962, 5952, 6049,\
      6008, 6387,\
     ] # maybe circulating objects - numbers start at 1 (like matlab), not 0 (pythonic)
horizons_file = 'b001_astorb_horizons.csv'
astorb_file = 'a000_nesvorny_astorb_astdys_wise_akari_sloan.dat_WITH_PROPER_ELMTS'
Hilda_file = 'a000_nesvorny_153_Hilda_family.list_090_WO_INTERLOPERS'
Schubart_file = 'a000_nesvorny_1911_Schubart_family.list_060_WO_INTERLOPERS'
Potomac_file = 'a000_nesvorny_1345_Potomac_family.list_140_WO_INTERLOPERS'
Francette_file = 'a000_nesvorny_1212_Francette_family.list_030_WO_INTERLOPERS'
Guinevere_file = 'a000_nesvorny_2483_Guinevere_family.list_040_WO_INTERLOPERS'
# a2008TG106_file = 'a000_nesvorny_269345_2008TG106_family.list_054_WO_INTERLOPERS'
# leave out 2008 TG106 family because it has only 17 members, 16 of which are also in Schubart family - assign remaining member to Background
family_files = [Hilda_file,Schubart_file,Potomac_file,Francette_file,Guinevere_file]
family_top_lines = [12,12,12,11,11]
family_names = ['Hilda','Schubart','Potomac','Francette','Guinevere']
ngroups = len(family_names)
idlists_families = []
idlist2 = []
for igroup in range(ngroups):
    infile = family_files[igroup]
    top_lines = family_top_lines[igroup]
    file1 = open(infile, 'r')
    Lines = file1.readlines()
    ids = []
    nobj = len(Lines) - top_lines
    iobj = top_lines
    while iobj < len(Lines):
        line = Lines[iobj]
        idhere = line[7:26]
        while '_' in idhere:
            iii = idhere.index('_')
            idhere = idhere[0:iii] + ' ' + idhere[iii+1:]
            idhere = idhere.lstrip()
            idhere = idhere.rstrip()
        ids.append(idhere)
        iobj = iobj + 1
    idlists_families.append(ids)
dfa = pd.read_csv(horizons_file)
nobj = dfa.shape[0]
idlist_astorb = dfa['idhere'].to_list()
Hlist_astorb = dfa['H_mag'].to_list()
family_list = []
for iobj in range(nobj):
    idhere = idlist_astorb[iobj]
    idhere2 = idhere.lstrip()
    idhere3 = idhere2.rstrip()
    check = 0
    append = 0
    appendages = []
    for igroup in range(ngroups):
        if idhere3 in idlists_families[igroup]:
            family_list.append(family_names[igroup])
            check = 1
            append = append + 1
            appendages.append(family_names[igroup])
    if check == 0:
        family_list.append('Background')
        append = append + 1
        appendages.append('Background')
    if append != 1:
        print(iobj,idhere3,appendages)
label_list = []
for iobj in range(nobj):
    if iobj+1 in cl:
        label_list.append('circulating')
    elif iobj+1 in ml:
        label_list.append('maybe')
    else:
        label_list.append('librating')
for iobj in range(nobj):
    idhere = idlist_astorb[iobj]
    idhere2 = idhere.rstrip()
    idhere3 = idhere2.lstrip()
    idlist2.append(idhere3)
dictionary = {'idhere':idlist2,\
              'H_mag':Hlist_astorb,\
              'family':family_list,\
              'label':label_list}
dfo = pd.DataFrame.from_dict(dictionary)
dfo.to_csv('b006_astorb_labels.csv',index=False)