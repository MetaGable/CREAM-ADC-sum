#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 8 07:53:54 2021
Last modified Jan 25th
Plot for all events (for files in /home/genec420/gain3)
@author: gene
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np
import glob as glob
# import datetime as dt
import time
import pickle
from heapq import nsmallest

where = "remote"
# where = "local"
table_path={"local":"/Users/gene/gains/F_raw_gain_table",
          "remote":"/home/genec420/full_auto/F_raw_gain_table",}
merge_path={"local":"/Users/gene/Desktop/code/merge", 
            "remote":"/home/genec420/misc/merge"}
color_list=['tab:blue', 'tab:orange', 'tab:purple', 'tab:red', 'tab:green',
            'tab:pink', 'tab:cyan', 'tab:olive', 'tab:brown']
merge = pickle.load(open(merge_path[where],'rb'))

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
    return l * 50 + r - 51


def uni_load_scatter(ribbon_num,which = 'C'):  
    '''First data cut'''
    if os.path.exists("F_"+str(ribbon_num)+str(which)):
        scatter = pickle.load(open("F_"+str(ribbon_num)+str(which),'rb'))
    else:
        scatter=[[],[]]
        input_type = '/home/genec420/gain3/C*'
        if which !='C': input_type = '/home/genec420/gain3/no*'
        input_paths = glob.glob(input_type)
        total = len(input_paths)
        for num, path in enumerate(input_paths):
            input_data = pickle.load(open(path, 'rb'))
            scatter[0]+=input_data[ribbon_num][1]
            scatter[1]+=input_data[ribbon_num][0]
            load_bar(num+1, total)
        scatter = np.array(scatter)    
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
    x_factor, y_factor = 10, 3
    low_x, low_y = np.median(x)-x_factor*np.std(x),np.median(y)-y_factor*np.std(y)
    # mask2 = [(i[0]>low_x and i[1]>low_y) for i in np.transpose([x,y])]
    # x,y = x[mask2],y[mask2]
    """Jun 13th experiment: separate vertical and horizontal cuts"""
    try:
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


def slant_data_cut(rib, rx, ry, x_factor, y_factor, fit):
    '''Don't plots the cuts'''
    low_x, low_y = np.median(rx)-x_factor*np.std(rx),np.median(ry)-y_factor*np.std(ry)
    sufficient_data = True
    '''retangular cut at bottom left'''
    mask = [(i[0]>low_x and i[1]>low_y) for i in np.transpose([rx,ry])]
    x,y = rx[mask],ry[mask]
    """Jun 13th experiment: separate vertical and horizontal cuts"""
    # mask = [i>low_y for i in ry]
    # x,y = rx[mask],ry[mask]
    # low_x = np.median(rx)-x_factor*np.std(rx)
    # mask = [i>low_x for i in rx]
    # x,y = rx[mask],ry[mask]
    if sum(mask)<30:
        sufficient_data=False
        print(f"factors: X={x_factor:<4}, Y={y_factor:<4}; not enough data")
        return x, y ,sufficient_data
    '''define coordinates for the slanted cut'''
    mid_x, mid_y = np.max(x)-(np.max(x)-np.min(x))*(6/7), min(y)
    mask = [merge[rib][4]*(i[0]-mid_x)<(i[1]-mid_y) for i in np.transpose([x,y])]
    x,y = x[mask],y[mask]
    if sum(mask)<20:
        sufficient_data=False
        print(f"factors: X={x_factor:<4}, Y={y_factor:<4}; not enough data")
        return x, y ,sufficient_data
    print(f"factors: X={x_factor:<4}, Y={y_factor:<4}; {len(x):>4} events left")
    intercept_x = (low_y-min(y))/merge[rib][4]+mid_x
    fit['plot_cut'].append((low_y,low_x, x_factor, y_factor))
    fit['plot_coord'].append((intercept_x, np.min(x), np.max(x), mid_x, np.min(y), merge[rib][4], np.max(y)))
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
    #     print(len(x)," points;    Chi^2 = ", red_chi, )
    # else:
    #     lenx = len(x)
    #     print(f"X factor={x_factor}, y factor={y_factor}, {lenx} events left")


def choose_fit(ax, x, y, plot_type, fit, total_plot_num, gain_info):
    _ = ax.scatter(x, y, alpha = 0.6,s=13)
    if fit['evt_count']==[]: 
        print('No fit to choose (all data cut away)')
        sufficient_data = False
        return sufficient_data, gain_info
    evt_num_Rrange = np.min(fit['evt_count']) + 0.1*np.std(fit['evt_count']) + 1
    # print("min = ", np.min(fit['evt_count']), "std = ", np.std(fit['evt_count']),
    #       'Evt Num cutoff: ', evt_num_Rrange)
    fit_mask = [i<evt_num_Rrange for i in fit['evt_count']]
    reduced_chi = np.array(fit['reduced_chi'])
    minimals = nsmallest(total_plot_num, set(reduced_chi[fit_mask]))
    for minimal in minimals:
        order = 0
        ind = np.where(reduced_chi==minimal)[0][0]
        (low_y, low_x, x_factor, y_factor) = fit['plot_cut'][ind]
        (intercept_x, min_x, max_x, mid_x, min_y, slope, max_y) = fit['plot_coord'][ind]
        m, b = fit['parameters'][ind]
        plot_x = np.linspace(min_x,max_x,10)
        _ = ax.plot(plot_x,m*plot_x+b,alpha = 0.6, 
                    label=r'$fit{0}|{3}events|Gain:{1}|\chi^2_\nu:{2}$'.format(
                        order,round(m,1),reduced_chi[ind],fit['evt_count'][ind]))
        # _ = ax.vlines(low_x, low_y, max_y, alpha = 0.6, label=x_factor, 
        #             color=color_list[0], linestyle='--')
        _ = ax.vlines(low_x, low_y, max_y, alpha = 0.6, 
                    color=color_list[0], linestyle='--')
        gain_info=(m, reduced_chi[ind], low_y, low_x, rib)
        order += 1
        if intercept_x < max_x:
            slope_range_x = np.linspace(intercept_x,max_x,10)
            _ = ax.plot(slope_range_x, slope*(slope_range_x-mid_x)+min_y, 
                       alpha = 0.5, color=color_list[order%9],linestyle='--')
            # _ = ax.hlines(low_y, low_x, intercept_x, alpha = 0.6, label=y_factor,
            #               color=color_list[order%9], linestyle='--')
            _ = ax.hlines(low_y, low_x, intercept_x, alpha = 0.6,
                          color=color_list[order%9], linestyle='--')
        else:
            # _ = ax.hlines(low_y, low_x, max_x, alpha = 0.6, label=y_factor, 
            #               color=color_list[order%9], linestyle='--')
            _ = ax.hlines(low_y, low_x, max_x, alpha = 0.6, 
                          color=color_list[order%9], linestyle='--')
    max_x += (max_x-min_x)/9
    max_y += (max_y-min_y)/9
    min_x -= (max_x-min_x)/7
    min_y -= (max_y-min_y)/7
    _=ax.set_ylim(bottom = min_y, top = max_y)
    _=ax.set_xlim(left = min_x, right = max_x)
    sufficient_data = True
    return sufficient_data, gain_info


def plot_labels(ax, plot_type='r', sufficient_data=True):
    _=ax.legend(loc=2, fontsize='small',framealpha=0.4)
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
        

y_factor_list=[-2, -4, -6, -8,]
x_factor_list=[2.2, 1.2, -1, -2.5]
# """Latest cut parameters for normal full ribbons:"""
# y_factor_list=[-2, -4, -8, -12]
# x_factor_list=[1.2, -0.5, -1, -1.5]
# """cut parameters for bad full ribbons:"""
# y_factor_list=[-10, -100]
# x_factor_list=[1.2, -0.5, -1, -1.5]
# """experiment with including all data"""
# x_factor_list = [10]
# y_factor_list = [10]
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


for rib in missing_ribs:
    print('\n------------- plotting raw ADC ', rib,'------------')
    rx,ry=uni_load_scatter(rib,'no ped')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    _=fig.suptitle('''layer {0} ribbon {1}: ch{2} mid-range ADC over ch{3} low-range ADC
                   From gain table: Gain = {4}, threshold = {5}'''.format(
                   merge[rib][0],merge[rib][1],merge[rib][3],merge[rib][2], 
                   merge[rib][4],merge[rib][5]),fontsize=19)
    fit = {'ind':0, 'parameters':[], 'evt_count':[], 'reduced_chi':[],
           'x_y_factor':[], 'plot_cut':[], 'plot_coord':[]}
    for x_factor in x_factor_list:
        for y_factor in y_factor_list:
            x,y,sufficient_data = slant_data_cut(rib, rx,ry,x_factor,y_factor,fit)
            # x,y,sufficient_data = slant_data_cut2(ax1, rx,ry,x_factor ,y_factor ,fit)
            fit_plot(ax1,x,y,sufficient_data,fit)
            fit['ind']+=1
    sufficient_data, gain_info = choose_fit(ax1, rx, ry, 'r', fit, n, gain_list)
    plot_labels(ax1,'r',sufficient_data)
    # ax1.scatter(rx,ry)
    del(x, y, rx, ry, sufficient_data, fit)
    print('--------- plotting subtracted  ADC ', rib,'--------')
    fit = {'ind':0, 'parameters':[], 'evt_count':[], 'reduced_chi':[],
           'x_y_factor':[], 'plot_cut':[], 'plot_coord':[]}
    cx,cy = uni_load_scatter(rib)
    for x_factor in x_factor_list:
        for y_factor in y_factor_list:
            x,y,sufficient_data = slant_data_cut(rib, cx,cy,x_factor,y_factor,fit)
            # x,y,sufficient_data = slant_data_cut2(ax2, cx,cy,x_factor,y_factor, fit)
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
    # print('###################################\n get full_auto/'+fn+".png \n###################################")
    plt.close('all')

for rib in missing_ribs:
    fn = "get full_auto/scatter_{}_{}F.png".format(merge[rib][0],merge[rib][1])
    print(fn)

for rib in missing_ribs:
    fn = f"get full_auto/F_{rib}C"
    print(fn)
    fn = f"get full_auto/F_{rib}no\ ped"
    print(fn)
