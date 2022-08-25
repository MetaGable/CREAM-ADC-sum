#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 09:46:24 2022

Collect individual pickled gain table into a complete table
Gain columns
1. rib_ind
2. layer
3. rib
4. original gain
5. original threshold
6. data_source
7. which_plot
8. cut_type
9. gain
10. std of different fits
11. Mid ch thres
12. Low ch thres
13. remarks


@author: gene
"""
import csv
import pickle
import os
# import pandas as pd
import numpy as np
import glob as glob
import re
'''put the fitting information in your own raw_gain to compute the weighted
 average to put into the gain list'''
reconstructed_raw_gains = pickle.load(open("/Users/gene/gains/R_fit_info_raw",'rb'))
full_raw_gains = pickle.load(open("/Users/gene/gains/F_raw_gain_table",'rb'))
merge = pickle.load(open("/Users/gene/Desktop/code/merge",'rb'))

auto_tracked = pickle.load(open("/Users/gene/L0Tables/Auto_tracked",'rb'))
auto_full_raw = pickle.load(open("/Users/gene/L0Tables/Auto_full????",'rb'))
F_raw_gain_table = pickle.load(open("/Users/gene/L0Tables/F_raw_gain_table",'rb'))
manual_track = pickle.load(open("/Users/gene/L0Tables/manual_track_raw",'rb'))
raw_gain_tables= {}
for key in ["taB","taL", "taR"]:
    raw_gain_tables[key] = auto_tracked
for key in ["tmB", "tmL", "tmR"]:
    raw_gain_tables[key] = manual_track
for key in ["faB","faL", "faR"]:
    raw_gain_tables[key] = auto_full_raw
for key in ["fmB", "fmL", "fmR"]:
    raw_gain_tables[key] = F_raw_gain_table

file_dir = {"taB":"Tracked/track Both (801)",
           "taL": "Tracked/track L (7)",
           "taR": "Tracked/track R (26)", 
           "tmB": "Tracked/track manual (29)", 
           "tmL": "",
           "tmR": "", #Manual ribbon unfinished
           "tt" : "Tracked/total count (994)", 
           "ft" : "full/count", 
           "nf" : "Tracked/need full", 
           "pi" : "Tracked/potential improvement", 
           "faB": "full/both",
           "faL": "full/Full L",
           "faR": "full/Full R",
           "fmB": "full/Full manual",
           "fmL": "",#The manual ribbons aren't organized yet
           "fmR": "",
           "nfd": "../need_full_downloaded"}

gain_spec = {"taB":('tracked', 'Both', 'auto'),
            "taL": ('tracked', 'Left', 'auto'),
            "taR": ('tracked', 'Right', 'auto'),
            "tmB": ('tracked', 'Both', 'manual'),
            "tmL": ('tracked', 'Left', 'manual'),
            "tmR": ('tracked', 'Right', 'manual'),
            "faB": ('full', 'Both', 'auto'),
            "faL": ('full', 'Left', 'auto'),
            "faR": ('full', 'Right', 'auto'),
            "fmB": ('full', 'Both', 'manual'),
            "fmB": ('full', 'Left', 'manual'),
            "fmB": ('full', 'Right', 'manual')}
s
def ind_lr(ind):
    l = ind // 50 + 1
    r = ind % 50 + 1
    return l, r

def lr_ind(l,r):
    return l * 50 + r - 51

def load(path):
    '''load the name of the files in the current directory in to python
    used inside [filling]'''
    good_files = glob.glob(f'/Users/gene/gains/{file_dir[path]}/scatter*')
    good_ribbons = []
    for filename in good_files:
        layer = re.findall(r"er_\d+_",filename)[0][3:-1]
        ribbon = re.findall(r"_\d+.png",filename)[0][1:-4]
        good_ribbons.append((int(layer),int(ribbon)))
    good_ribbons_index = sorted([lr_ind(i[0],i[1]) for i in good_ribbons])
    return good_ribbons, good_ribbons_index


def flling(path, source = 'full', data_type = 'reconstructed'):
    '''flling in the info'''
    _ , good_ribbons_index=load(path)
    table = raw_gain_tables[path]
    std='NA'
    for rib_ind in good_ribbons_index:
        temp = sorted(table[rib_ind],key=lambda x: x[1])
        if 4>=len(temp)>=2:
            std=round(np.std(np.transpose(temp)[0]),2)
        elif len(temp)>=4:
            std=round(np.std(np.transpose(temp)[0][0:4]),2)
        low_x, low_y, remarks= 'NA', 'NA', 'NA'
        for i in temp:
            if i[2]>30000:
                low_x, low_y=round(i[2]), round(i[3])
                break
        layer, rib =  (rib_ind+50)//50, (rib_ind)%50+1
        gain = round(temp[0][0], 2)
        data_source, which_plot, cut_type = gain_spec[path]
        gain_list[rib_ind]=(rib_ind, layer, rib, merge[rib_ind][4],
                    merge[rib_ind][5], data_source, which_plot, cut_type,
                    gain, std, low_x, low_y, remarks)

if os.path.exists('/Users/gene/fitted_gain.csv'):
    '''adding to existing file'''
    fields = []
    gain_list = []
    with open('/Users/gene/fitted_gain.csv', 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        fields = next(csvreader)
        for row in csvreader:
            gain_list.append(row)
        print("Total no. of rows: %d"%(csvreader.line_num))
    gain_list.pop(0)
else:
    '''Starting a new file'''
    '''header will mess up the indicies'''
    gain_list = [(i, merge[i][0], merge[i][1], merge[i][4], merge[i][5],'',
                  '','','', '', '', '', 'Fit unavailable') for i in range(1000)]


flling("taB",'full', 'reconstructed')
flling("taR",'R','reconstructed')
# flling("tmL",'L','reconstructed')
# flling("faB",'full', 'all data')
# flling("faR",'R','all data')
# flling("faL",'L','all data')

with open('/Users/gene/fitted_gain3.csv', 'w') as csvfile:
    csvwriter = csv.writer(csvfile) 
    csvwriter.writerow(('rib_ind', 'layer', 'rib', 'original gain',
                'original threshold', 'data_source', 'which_plot', 'cut_type',
                'fitted gain', 'std of different fits', 'Mid ch thres', 
                'Low ch thres', 'remarks'))
    csvwriter.writerows(gain_list)














# """auto_tracked: working"""
# auto_tracked = pickle.load(open("/Users/gene/Auto_tracked",'rb'))
# # auto_tracked = [i for i in auto_tracked if i!=[]] #len 991

# """auto_full: not raw, need to regenerate the raw table"""
# auto_full_raw = pickle.load(open("/Users/gene/Auto_full???? ",'rb'))
# # auto_full_raw = [i for i in auto_full_raw if type(i)!=str] #len 102
# """only one ribbon"""
# auto_full = pickle.load(open("/Users/gene/Auto_full",'rb')) #1

# """full manual: working
# need to label time window"""
# manual_full = pickle.load(open("/Users/gene/Manual_full_newest",'rb'))
# # manual_full = [i for i in manual_full if i!=[]] #len 27
# F_raw_gain_table = pickle.load(open("/Users/gene/F_raw_gain_table",'rb'))
# # F_raw_gain_table = [i for i in F_raw_gain_table if i!=[]] #len 68

# """tracked manual: some """
# manual_track = pickle.load(open("/Users/gene/manual_track_raw",'rb'))
# # [i for i in manual_track if i !=[]]

# raw_gain_tables= {"taB": auto_tracked,"taL": auto_tracked,
#                   "tm": manual_track, "taR":auto_tracked,
#                   "tmL": manual_track,"tmR": manual_track,
#                   "faB": auto_full_raw,"faL": auto_full_raw,
#                   "fm": F_raw_gain_table, "faR":auto_full_raw,
#                   "fmL": F_raw_gain_table,"fmR": F_raw_gain_table}

# """Counting missing ribs"""
# at = pickle.load(open("/Users/gene/Auto_tracked",'rb'))
# at = [1 if i!=[] else 0 for i in at]
# af = pickle.load(open("/Users/gene/Auto_full????",'rb'))
# af = [1 if type(i)!=str else 0 for i in af]
# mf = pickle.load(open("/Users/gene/F_raw_gain_table",'rb'))
# mf = [1 if i!=[] else 0 for i in mf]
# mt = pickle.load(open("/Users/gene/manual_track_raw",'rb'))
# mt = [1 if i!=[] else 0 for i in mt]
# count=zip(at, af, mf, mt)


# missing = [76, 800, 804, 806, 984]
# missinglr = [ (2, 27), (17, 1), (17, 5), (17, 6), (20, 35)]
# """ Only 17_6 and 20_35 salvable, Use manual fitting later on. 
#     17_1 and 17_5 is a unfathomable big mess; 
#     2_27 lacks data"""
# for i1, i2 in missinglr:
#     print(f'get scatter_{i1}_{i2}.png')


"""filtering auto track"""


_=[np.transpose(sorted(auto_tracked[i],key = lambda x:x[1]))[0] for i in ribs]
rib_std=[i for i in zip(np.arange(798), ribs, std, _)]
temp=sorted(rib_std,key=lambda x:x[2])
rib_std2=[i for i in zip(ribs, _) if len(i[1])>=2]
temp2=sorted(rib_std2,key=lambda x:abs(x[1][1]-x[1][0]))
temp3 = [abs(i[1][1]-i[1][0]) for i in temp2]


missing_pi_ind=[i for i in pi_ind if i not in ftotal_ind]
missing_pi_ind=[i for i in pi_ind if i not in missing_nf_ind]



for i in missing_nf:
    print(f'get full_auto/scatter_{i[0]}_{i[1]}.png')

"""
Remarks
2_40 label time window. 
20_42 lack full plot
Full version, bad cut. Need to re-generate ped intact full data without cut
    19_49 ped intact need plot
    6_11 both need plot
    11_43 both need plot
    18_35 both need plot
    17_14 bad cut
    17_10 bad cut

"""
gain_list[lr_ind(2,40)][12]="removed 2017Dec13~14, 2018Sep11~13"


