#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 23:24:38 2022
Combined Plot_all_evt witg all_manual allow omitting specified cuts

@author: gene
"""


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np
import glob as glob
import time
import pickle
import csv
from heapq import nsmallest
# import datetime as dt

where = "remote"
# where = "local"
table_path={"local":"/Users/gene/gains/FManual_raw_gain_table",
          "remote":"/home/genec420/full_auto/FManual_raw_gain_table",}
merge_path={"local":"/Users/gene/Desktop/code/merge", 
            "remote":"/home/genec420/misc/merge"}
color_list=['tab:blue', 'tab:orange', 'tab:purple', 'tab:red', 'tab:green',
            'tab:pink', 'tab:cyan', 'tab:olive', 'tab:brown']
merge = pickle.load(open(merge_path[where],'rb'))
# data_type = input('track / full?/n')
data_type = 'full'
directory={'track':'gain2', 'full':'gain3', 'R':'no*', 'C':"C*"}

def load_bar(iteration, total, prefix='', suffix='', decimal=1, length=90):
    if total == -1:
        """for testing"""
        for i in np.arange(100):
            time.sleep(0.04)
            load_bar(i+1,100,'Display Testing', 'Complete',length = 40)
        return
    percent = ('{0:.' + str(decimal) + 'f}').format(100 * (iteration/float(total)))
    filled_length = int(length*iteration//total)
    bar = '>' * filled_length + '-' * (length-filled_length)
    print(f'\r {prefix} |{bar}| {percent}% {suffix}',end='\r')
    if iteration == total:
        print()

def ind_lr(ind):
    l = ind // 50 + 1
    r = ind % 50 + 1
    return l, r

def lr_ind(l,r):
    return int(l * 50 + r - 51)

def uni_load_scatter(ribbon_num,which = 'C'):  
    '''First data cut'''
    if data_type == 'track':
        file_prefix = 'T_'
        directory = 'gain2'
    elif data_type == 'full':
        file_prefix = 'F_'
        directory = 'gain3'
    if os.path.exists(file_prefix+str(ribbon_num)+str(which)):
        scatter = pickle.load(open(file_prefix+str(ribbon_num)+str(which),'rb'))
        # print('------ quick load complete------')
    else:
        scatter=[[],[]]
        if which =='R':
            input_type = '/home/genec420/'+directory+'/no*'
        else:
            input_type = '/home/genec420/'+directory+'/C*'
        input_paths = glob.glob(input_type)
        total = len(input_paths)
        for num, path in enumerate(input_paths):
            input_data = pickle.load(open(path, 'rb'))
            scatter[0]+=input_data[ribbon_num][1]
            scatter[1]+=input_data[ribbon_num][0]
            load_bar(num+1, total)
        scatter = np.array(scatter,dtype=float)    
        with open("F_"+str(ribbon_num)+str(which), 'wb') as fp:
            pickle.dump(scatter, fp)    
    if which =='C':
        mask1 = [(i[0]!=0 and i[1]!=0)and (i[0]<30000) for i in np.transpose(scatter)]
        # x, y = scatter[0][mask1], scatter[1][mask1]
    else:
        mask1 = [i[0]>0 and i[1]>0 for i in np.transpose(scatter)]
    x, y = scatter[0][mask1], scatter[1][mask1]
    Length0 = len(x)
    '''define coordinates for the slanted cut'''
    try:
        x_factor, y_factor = 7, 2
        low_y = np.median(y)-y_factor*np.std(y)
        mask = [i>low_y for i in y]
        x1,y1 = x[mask],y[mask]
        low_x = np.median(x1)-x_factor*np.std(x)
        mask = [i>low_x for i in x1]
        x2,y2 = x1[mask],y1[mask]
        mid_x, mid_y = np.max(x2)-(np.max(x2)-np.min(x2))/6, min(y2)
        mask3 = [merge[ribbon_num][4]*(i[0]-mid_x)<(i[1]-mid_y) for i in np.transpose([x2,y2])]
        x3,y3 = x2[mask3],y2[mask3]
        print(f"1st data cut: {Length0} -> {len(x)} events left")
    except:
        if len(x2)>80:
            x, y = x2, y2
        elif len(x1)>80:
            x, y = x1, y1
        print(f"No Cuts could be made, {Length0} events left")
        return x, y            
    return x3,y3

def data_cut(rib, x, y,x_fac, y_fac, fit, Cut_info, Slope_info):
    sufficient_data = True
    LCut, RCut, BCut, TCut = Cut_info
    auto_LCut, auto_BCut = np.median(x)-x_fac*np.std(x),np.median(y)-y_fac*np.std(y)   
    if not bool(LCut): LCut = auto_LCut
    if not bool(BCut): BCut = auto_BCut
    mask = [i[0]>LCut and i[1]>BCut for i in np.transpose([x,y])]
    x,y = x[mask],y[mask]
    if RCut!="":
        mask = [i[0]<RCut for i in np.transpose([x,y])]
        x,y = x[mask],y[mask]
    if TCut!="":
        mask = [i[1]<TCut for i in np.transpose([x,y])]    
        x,y = x[mask],y[mask]
    len2 = len(x)
    if len2<20:
        sufficient_data=False
        print(f"factors: X={x_factor:<4}, Y={y_factor:<4}; events > {len2:^4} Cut aborted #1")
        return x, y ,sufficient_data
    _x,_y = x, y
    Bmid_x_factor, B_slope, B_point_x, B_point_y = Slope_info[0:4]
    if B_slope!='0':
        if not bool(B_slope):B_slope=merge[rib][4]
        if not bool(Bmid_x_factor):Bmid_x_factor = 1/7 
        if not bool(B_point_x):
            B_point_x = np.min(x)+(np.max(x)-np.min(x))*Bmid_x_factor
            B_point_y = auto_BCut
        mask = [B_slope*(i[0]-B_point_x)<(i[1]-B_point_y) for i in np.transpose([x,y])]
        x,y = x[mask],y[mask]
        len3 = sum(mask)
        if len3<15:
            sufficient_data=True
            print(f"factors: X={x_factor:<4}, Y={y_factor:<4}; events > {len2:^4} > {len3:^3}, Fall back #2")
            fit['plot_cut'].append((LCut, RCut, BCut, TCut, min(_x),max(_x),min(_y),max(_y)))
            fit['plot_coord'].append(((0,0,0),(0,0,0)))
            return _x, _y ,sufficient_data
    _x,_y = x, y
    Tmid_x_factor, T_slope, T_point_x, T_point_y = Slope_info[4:8]
    if T_slope!='0':    
        if not bool(T_slope): T_slope=merge[rib][4]
        if not bool(Tmid_x_factor):Tmid_x_factor = 1/2
        if not bool(T_point_x):
            T_point_x = float(np.min(x)+(np.max(x)-np.min(x))*Tmid_x_factor)
            T_point_y = float(np.max(y))
        # print(type(x[0]))
        mask = [T_slope*(i[0]-T_point_x)>(i[1]-T_point_y) for i in np.transpose([x,y])]
        x,y = x[mask],y[mask]
        len4 = len(x)
        if len4<15:
            sufficient_data=True
            print(f"factors: X={x_factor:<4}, Y={y_factor:<4}; events: > {len2:^4} > {len3:^3} > {len4:^3}, Fall back #3")
            fit['plot_cut'].append((LCut, RCut, BCut, TCut, min(_x),max(_x),min(_y),max(_y)))
            fit['plot_coord'].append(((0,0,0),(B_slope, B_point_x, B_point_y)))
            return _x, _y ,sufficient_data
        print(f"factors: X={x_factor:<4}, Y={y_factor:<4}; events: > {len2:^4} > {len3:^3} > {len4:^3}")
    else:
        print(f"factors: X={x_factor:<4}, Y={y_factor:<4}; events: > {len2:^4} > {len3:^3}")
    fit['plot_cut'].append((LCut, RCut, BCut, TCut,np.min(x),max(x),min(y),max(y)))
    fit['plot_coord'].append(((T_slope,T_point_x,T_point_y),(B_slope, B_point_x, B_point_y)))
    return x, y, sufficient_data
    

def fit_plot(ax1,x,y,sufficient_data, fit):
    if sufficient_data:
        (m, b), [residual], _, _, _ = np.polyfit(x,y,1,cov=False,full=True)
        red_chi = round(residual/(len(x)-2))
        '''store the fit info for selection'''
        fit['parameters'].append((m, b))
        fit['reduced_chi'].append(red_chi)
        fit['evt_count'].append(len(x))
        fit['x_y_factor'].append((x_factor,y_factor))


def choose_fit(ax, x, y, plot_type, fit, total_plot_num, gain_info):
    _ = ax.scatter(x, y, alpha = 0.6,s=13)
    if fit['evt_count']==[]: 
        print('No fit to choose (all data cut away)')
        sufficient_data = False
        return sufficient_data, gain_info
    evt_num_Rrange = np.min(fit['evt_count']) + 0.1*np.std(fit['evt_count']) + 1
    fit_mask = [i<evt_num_Rrange for i in fit['evt_count']]
    reduced_chi = np.array(fit['reduced_chi'])
    minimals = nsmallest(total_plot_num, set(reduced_chi[fit_mask]))
    order = 0
    for minimal in minimals:
        ind = np.where(reduced_chi==minimal)[0][0]        
        (low_x, hi_x, low_y, hi_y, min_x,max_x,min_y,max_y) = fit['plot_cut'][ind]
        ((T_slope,T_point_x,T_point_y),(B_slope, B_point_x, B_point_y)) = fit['plot_coord'][ind]
        m, b = fit['parameters'][ind]
        plot_x = np.linspace(min_x,max_x,10)
        _ = ax.plot(plot_x,m*plot_x+b,alpha = 0.6, 
                    label=r'$fit{0}|{3}events|Gain:{1}|\chi^2_\nu:{2}$'.format(
                        order,round(m,1),reduced_chi[ind],fit['evt_count'][ind]))
        gain_info=(m, reduced_chi[ind], low_y, low_x, rib)
        if bool(hi_x):
            _ = ax.vlines(hi_x, min_y, max_y, alpha = 0.6, 
                color=color_list[order%9], linestyle='--')
        if bool(hi_y):
            _ = ax.hlines(hi_y, low_x, max_x, alpha = 0.6, 
                  color=color_list[order%9], linestyle='--')
        if bool(T_point_x):
            T_intersect_x = (max_y-T_point_y)/T_slope+T_point_x
            T_intersect_y = (low_x-T_point_x)*T_slope+T_point_y
            slope_range_x = np.linspace(low_x,T_intersect_x,10)
            _ = ax.plot(slope_range_x, T_slope*(slope_range_x-T_point_x)+T_point_y, 
                       alpha = 0.5, color=color_list[order%9],linestyle='--')
            _ = ax.vlines(low_x, min_y, T_intersect_y, alpha = 0.6, 
                        color=color_list[order%9], linestyle='--')
        else:
            _ = ax.vlines(low_x, min_y, max_y, alpha = 0.6, 
                    color=color_list[order%9], linestyle='--')
        if bool(B_slope):            
            B_point_min_x = (low_y-B_point_y)/B_slope+B_point_x
            if B_point_min_x < max_x:
                slope_range_x = np.linspace(B_point_min_x,max_x,10)
                _ = ax.plot(slope_range_x, B_slope*(slope_range_x-B_point_x)+B_point_y, 
                           alpha = 0.5, color=color_list[order%9],linestyle='--')
                _ = ax.hlines(low_y, low_x, B_point_min_x, alpha = 0.6,
                              color=color_list[order%9], linestyle='--')
            else:
                _ = ax.hlines(low_y, low_x, max_x, alpha = 0.6, 
                              color=color_list[order%9], linestyle='--')
        order += 1
    max_x += (max_x-min_x)/9
    max_y += (max_y-min_y)/8
    min_x -= (max_x-min_x)/7
    min_y -= (max_y-min_y)/7
    _=ax.set_ylim(bottom = min_y, top = max_y)
    _=ax.set_xlim(left = min_x, right = max_x)
    sufficient_data = True
    return sufficient_data, gain_info


def plot_labels(ax, plot_type='r', sufficient_data=True):
    _=ax.set_ylabel('low-ranged ch ADC',fontdict = {'fontsize' : 16})
    _=ax.set_xlabel('mid-ranged ch ADC',fontdict = {'fontsize' : 16})
    _=ax.grid(alpha=0.5)
    if plot_type=='r':
        ''' TO DO: make x & y factor as input'''
        _=ax.set_title("With pedestal intact",fontdict = {'fontsize' : 17})
        if not sufficient_data:
            ax.text(0.5, 0.5, 'not enough data to fit ')
            return
    elif plot_type=='c':
        _=ax.set_title("Subtracted pedestal with Cmean",fontdict = {'fontsize' : 17})
        if not sufficient_data:
            ax.text(0.5, 0.5, 'not enough data to fit ')
            return
    _=ax.legend(loc=2, fontsize='small',framealpha=0.4)


n = 2
# n = 1
if not os.path.exists(table_path[where]):
    with open(table_path[where], 'wb') as fp:
        gain_list = [[] for i in range(1000)]
        pickle.dump(gain_list, fp)
else:
    gain_list = pickle.load(open(table_path[where],'rb'))


"""For some ribs not showing events, Jun 10th"""
missing_ribs=[101, 105, 109, 114, 116, 119, 121, 122, 124, 125, 128, 130,
              184, 201, 211, 227, 229, 237, 241, 268, 300, 308, 364, 401,
              404, 615, 623, 668, 745, 825, 831, 878, 880, 881, 884, 928]
'''need to regenerate following ribbons on terminal to append gainlist'''
# missing_ribs = [881, 880, 831, 825, 668, 268, 241, 184, 130, 128, 125, 122, 121, 119,114]

parameter_list=[]
with open('full_manual_record2.csv', newline='') as csvfile:
# with open('/home/genec420/full_manual/full_manual_record.csv', newline='') as csvfile:
# with open('/home/genec420/full_manual/windowçed_full_manual_record.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    headers = next(reader)
    for row in reader:
        parameter_list.append([float(i) if i !='' else '' for i in row ])

for parameters in parameter_list:
    layer, ribbon = parameters[0:2]
    rib = lr_ind(layer, ribbon)
    x_factor_list=[2.2, 1.2, -1, -2.5]
    y_factor_list=[-2, -4, -6, -8,]
    print('\n------------- plotting raw ADC ', rib,'------------')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))    
    fig.set_dpi(300)
    _=fig.suptitle('''layer {0} ribbon {1}: ch{2} mid-range ADC over ch{3} low-range ADC
                   From gain table: Gain = {4}, threshold = {5}'''.format(
                   merge[rib][0],merge[rib][1],merge[rib][3],merge[rib][2], 
                   merge[rib][4],merge[rib][5]),fontsize=19)
    fit = {'ind':0, 'parameters':[], 'evt_count':[], 'reduced_chi':[],
           'x_y_factor':[], 'plot_cut':[], 'plot_coord':[]}
    rx,ry=uni_load_scatter(rib,'no ped')
    L_Cut_info, L_Slope_info = parameters[2:6],parameters[6:14]
    R_Cut_info, R_Slope_info = parameters[14:18],parameters[18:26]
    if bool(L_Cut_info[0]):
        x_factor_list=[1]
    if bool(L_Cut_info[2]):
        y_factor_list=[1]
    for x_factor in x_factor_list:
        for y_factor in y_factor_list:
            x,y,sufficient_data = data_cut(rib, rx,ry,x_factor,y_factor,fit, L_Cut_info, L_Slope_info)
            fit_plot(ax1,x,y,sufficient_data,fit)
            fit['ind']+=1
    sufficient_data, gain_info = choose_fit(ax1, rx, ry, 'r', fit, n, gain_list)
    plot_labels(ax1,'r',sufficient_data)
    del(x, y, rx, ry, sufficient_data, fit)
    print('--------- plotting subtracted  ADC ', rib,'--------')
    fit = {'ind':0, 'parameters':[], 'evt_count':[], 'reduced_chi':[],
           'x_y_factor':[], 'plot_cut':[], 'plot_coord':[]}
    cx,cy = uni_load_scatter(rib)
    x_factor_list=[2.2, 1.2, -1, -2.5]
    y_factor_list=[-2, -4, -6, -8,]
    if bool(R_Cut_info[0]):
        x_factor_list=[1]
    if bool(R_Cut_info[2]):
        y_factor_list=[1]
    for x_factor in x_factor_list:
        for y_factor in y_factor_list:
            x,y,sufficient_data = data_cut(rib, cx,cy,x_factor,y_factor,fit, R_Cut_info, R_Slope_info)
            fit_plot(ax2,x,y,sufficient_data,fit)
            fit['ind']+=1
    sufficient_data, gain_info = choose_fit(ax2, cx, cy, 'c', fit, n, gain_list)
    plot_labels(ax2,'c', sufficient_data)
    # ax2.scatter(cx, cy)
    del(x, y, cx, cy, sufficient_data, fit)
    fig.tight_layout()
    gain_list[rib]=gain_info
    with open(table_path[where], 'wb') as fp:
        pickle.dump(gain_list, fp)
    fn = "scatter_{}_{}F".format(merge[rib][0],merge[rib][1])
    plt.savefig(fn)
    plt.close('all')


ribs = [lr_ind(i[0],i[1]) for i in parameter_list]
for rib in ribs:
    fn = "get full_auto/scatter_{}_{}F.png".format(merge[rib][0],merge[rib][1])
    print(fn)

# temp=[]
# for parameters in parameter_list:
#     layer, ribbon = parameters[0:2]
#     rib = lr_ind(layer, ribbon)
#     if not os.path.exists('F_'+str(rib)+'no ped'):
#         # print(f'get F_{rib}C\nget F_{rib}no\ ped')
#         fn = f"get full_auto/F_{rib}no\ ped"
#         print(fn)
        # temp.append(rib)
# print('------ quick load complete------')# for rib in missing_ribs:
#     fn = f"get full_auto/F_{rib}C"
#     print(fn)
#     fn = f"get full_auto/F_{rib}no\ ped"
#     print(fn)
