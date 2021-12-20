#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 8 07:53:54 2021
Last modified Dec 17th
Plot for all events (for files in /home/genec420/gain3)
@author: gene
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob as glob
# import datetime as dt
import pickle
from heapq import nsmallest
merge = pickle.load(open("/home/genec420/misc/merge",'rb'))
color_list=['tab:blue', 'tab:orange', 'tab:purple', 'tab:red', 'tab:green',
            'tab:pink', 'tab:cyan', 'tab:olive', 'tab:brown']

def lr(layer,ribbon):
    '''converts (layer, ribbon) into (the number of the ribbon from 0~999)'''
    return layer*50+ribbon-51


def uni_load_scatter(ribbon_num,which = 'C'):  
    '''First data cut'''
    x_factor, y_factor = 3, 3
    scatter=[[],[]]
    input_type = '/home/genec420/gain3/C*'
    if which !='C': input_type = '/home/genec420/gain3/no*'
    input_paths = glob.glob(input_type)
    total = len(input_paths)
    for num, path in enumerate(input_paths):
    # ''' Only loaded some of the data for faster debugging'''
    # for path in input_paths[::5]:
        input_data = pickle.load(open(path, 'rb'))
        scatter[0]+=input_data[ribbon_num][1]
        scatter[1]+=input_data[ribbon_num][0]
        print(num,'/',total)
    scatter = np.array(scatter)
    if which =='C':
        mask1 = [(i[0]!=0 and i[1]!=0)and (i[0]<30000) for i in np.transpose(scatter)]
    else:
        mask1 = [(i[0]>0 and i[1]>0) for i in np.transpose(scatter)]
    '''define coordinates for the slanted cut'''
    x, y = scatter[0][mask1], scatter[1][mask1]
    '''test'''
    low_x, low_y = np.median(x)-x_factor*np.std(x),np.median(y)-y_factor*np.std(y)
    mask2 = [(i[0]>low_x and i[1]>low_y) for i in np.transpose([x,y])]
    x,y = x[mask2],y[mask2]
    mid_x, mid_y = (np.max(x)+np.min(x))/2, min(y)
    mask3 = [8.4*(i[0]-mid_x)<(i[1]-mid_y) for i in np.transpose([x,y])]
    x,y = x[mask3],y[mask3]
    # print('{} events loaded'.format(len(x)))
    return x,y


def slant_data_cut(rib, rx, ry, x_factor, y_factor, fit):
    '''Don't plots the cuts'''
    low_x, low_y = np.median(rx)-x_factor*np.std(rx),np.median(ry)-y_factor*np.std(ry)
    sufficient_data = True
    '''retangular cut at bottom left'''
    mask = [(i[0]>low_x and i[1]>low_y) for i in np.transpose([rx,ry])]
    fit['plot_cut'].append((low_y,low_x))
    x,y = rx[mask],ry[mask]
    if sum(mask)<30:
        sufficient_data=False
        '''not enough data'''
        return x,y,sufficient_data
    '''define coordinates for the slanted cut'''
    mid_x, mid_y = (np.max(x)+np.min(x))/2, (np.max(y)+np.min(y))/2-(np.max(y)-np.min(y))/3
    mask = [merge[rib][4]*(i[0]-mid_x)<(i[1]-mid_y) for i in np.transpose([x,y])]
    if sum(mask)<30:
            sufficient_data=False
            '''not enough data'''
            return x,y,sufficient_data
    x,y = x[mask],y[mask]
    intercept_x = (low_y-min(y))/merge[rib][4]+mid_x
    fit['plot_coord'].append((intercept_x, np.min(x), np.max(x), mid_x, min(y), merge[rib][4]))
    return x, y, sufficient_data


def slant_data_cut2(ax1, rx, ry, x_factor, y_factor, fit):
    low_x, low_y = np.median(rx)-x_factor*np.std(rx),np.median(ry)-y_factor*np.std(ry)
    sufficient_data = True
    num = fit['ind']
    '''retangular cut at bottom left'''
    mask = [(i[0]>low_x and i[1]>low_y) for i in np.transpose([rx,ry])]
    fit['plot_cut'].append((low_y,low_x))
    x,y = rx[mask],ry[mask]
    if sum(mask)<15:
            sufficient_data=False
            '''not enough data'''
            return x,y,sufficient_data
    '''define coordinates for the slanted cut'''
    # mid_x, mid_y = (np.max(x)+np.min(x))/2, (np.max(y)+np.min(y))/2-(np.max(y)-np.min(y))/3
    mid_x, mid_y = (np.max(x)+np.min(x))/2, min(y)
    mask = [8.4*(i[0]-mid_x)<(i[1]-mid_y) for i in np.transpose([x,y])]
    if sum(mask)>10:
        x,y = x[mask],y[mask]
        min_x, max_x, min_y, max_y = np.min(x), np.max(x), np.min(y), np.max(y)
        intercept_x = (low_y-mid_y)/8.4+mid_x
        _=ax1.vlines(low_x, min_y, max_y, alpha = 0.6, color=color_list[num%9], linestyle='--') 
        if intercept_x < max_x:
            slope_range_x = np.linspace(intercept_x,max_x,10)
            _=ax1.plot(slope_range_x, 8.4*(slope_range_x-mid_x)+mid_y, 
                       alpha = 0.5, color=color_list[num%9],linestyle='--')
            _=ax1.hlines(low_y, min_x, intercept_x, alpha = 0.6, color=color_list[num%9], linestyle='--')
        else:
            print('intercept Beyond bound:',intercept_x, max_x)
            _=ax1.hlines(low_y, min_x, max_x, alpha = 0.6, color=color_list[num%9], linestyle='--')
    else:
        sufficient_data=False
        print('not enough data, #', num)
    return x,y,sufficient_data


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
    #     print(len(x)," points;    not enough data to fit")


def choose_fit(ax, plot_type, fit, total_plot_num, gain_list):
    if fit['evt_count']==[]: 
        print('No fit to choose (all data cut away)')
        sufficient_data = False
        return sufficient_data
    if plot_type == 'c':
        _ = ax.scatter(cx, cy, alpha = 0.6,s=13)
    else:
        _ = ax.scatter(rx, ry, alpha = 0.6,s=13)
    evt_num_Rrange = np.min(fit['evt_count']) + 0.1*np.std(fit['evt_count']) + 1
    # print("min = ", np.min(fit['evt_count']), "std = ", np.std(fit['evt_count']),
    #       'Evt Num cutoff: ', evt_num_Rrange)
    fit_mask = [i<evt_num_Rrange for i in fit['evt_count']]
    reduced_chi = np.array(fit['reduced_chi'])
    minimals = nsmallest(total_plot_num, set(reduced_chi[fit_mask]))
    for minimal in minimals:
        order = 0
        ind = np.where(reduced_chi==minimal)[0][0]
        (low_y, low_x) = fit['plot_cut'][ind]
        (intercept_x, min_x, max_x, mid_x, mid_y, slope) = fit['plot_coord'][ind]
        m, b = fit['parameters'][ind]
        plot_x = np.linspace(min_x,max_x,10)
        _ = ax.plot(plot_x,m*plot_x+b,alpha = 0.6, 
                    label=r'$fit{0}|{3}events|Gain:{1}|\chi^2_\nu:{2}$'.format(
                        order,round(m,1),reduced_chi[ind],fit['evt_count'][ind]))
        _ = ax.vlines(low_x, low_y, max(y), alpha = 0.6, 
                    color=color_list[0], linestyle='--')
        gain_list.append((m, reduced_chi[ind], low_y, low_x))
        order += 1
        if intercept_x < max_x:
            slope_range_x = np.linspace(intercept_x,max_x,10)
            _ = ax.plot(slope_range_x, slope*(slope_range_x-mid_x)+mid_y, 
                       alpha = 0.5, color=color_list[order%9],linestyle='--')
            _ = ax.hlines(low_y, min_x, intercept_x, alpha = 0.6, color=color_list[order%9], linestyle='--')
        else:
            _ = ax.hlines(low_y, min_x, max_x, alpha = 0.6, color=color_list[order%9], linestyle='--')
    sufficient_data = True
    return sufficient_data


def plot_labels(ax, plot_type='r', sufficient_data=True):
    x_factor, y_factor = 6, 6
    _=ax.legend(loc=2, fontsize='small',framealpha=0.4)
    _=ax.set_ylabel('low-ranged ch ADC',fontdict = {'fontsize' : 16})
    _=ax.set_xlabel('mid-ranged ch ADC',fontdict = {'fontsize' : 16})
    _=ax.grid(alpha=0.5)
    if plot_type=='r':
        ''' TO DO: make x & y factor as input'''
        low_x, low_y = np.median(rx)-x_factor*np.std(rx),np.median(ry)-y_factor*np.std(ry)
        if not sufficient_data:
            ax.text(0.5, 0.5, 'not enough data to fit ')
            return
        _=ax.set_title("With pedestal intact",fontdict = {'fontsize' : 17})
        _=ax.set_xlim(left=low_x)
        _=ax.set_ylim(bottom=low_y)
    elif plot_type=='c':
        low_x, low_y = np.median(cx)-x_factor*np.std(cx),np.median(cy)-y_factor*np.std(cy)
        _=ax.set_title("Subtracted pedestal with Cmean",fontdict = {'fontsize' : 17})
        if not sufficient_data:
            ax.text(0.5, 0.5, 'not enough data to fit ')
            return
        _=ax.set_ylim(bottom = low_y)
        _=ax.set_xlim(left = low_x)


gain_list = [[] for i in range(1000)]
y_factor_list=[-4, -8, -12]
# x_factor_list=[1.5, 1]
x_factor_list=[1.5, -0.5]
# n = len(x_factor_list)*len(y_factor_list)
n = 2
# for rib in [lr(1, 23), lr(2,50), lr(12, 9), lr(19,42)]:
for rib in np.arange(0,1000):
    print('----------------------------------plotting raw ADC ', 'rib','---------------------------------')
    rx,ry=uni_load_scatter(rib,'no ped')
    # print('X data range: ', np.min(rx),' ~ ', np.max(rx), "; ", len(rx), ' events')
    # print('Y data range: ', np.min(ry),' ~ ', np.max(ry), "; ", len(ry), ' events')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    _=fig.suptitle('''layer {0} ribbon {1}: ch{2} mid-range ADC over ch{3} low-range ADC
                   From gain table: Gain = {4}, threshold = {5}'''.format(
                   merge[rib][0],merge[rib][1],merge[rib][2],merge[rib][3], 
                   merge[rib][4],merge[rib][5]),fontsize=19)
    fit = {'ind':0, 'parameters':[], 'evt_count':[], 'reduced_chi':[],
           'x_y_factor':[], 'plot_cut':[], 'plot_coord':[]}
    for x_factor in x_factor_list:
        for y_factor in y_factor_list:
            x,y,sufficient_data = slant_data_cut(rib, rx,ry,x_factor,y_factor,fit)
            # x,y,sufficient_data = slant_data_cut2(ax1, rx,ry,x_factor ,y_factor ,fit)
            # print('X_low cut: {}, Y_low cut {}, x_factor: {}, y_factor: {}, evt_num: {}'.format(
            #     int(np.median(x)-x_factor*np.std(x)),int(np.median(y)-y_factor*np.std(y)), x_factor, y_factor, len(x)))
            fit_plot(ax1,x,y,sufficient_data,fit)
            fit['ind']+=1
    sufficient_data = choose_fit(ax1, 'r', fit, n, gain_list[rib])
    plot_labels(ax1,'r',sufficient_data)
    # ax1.scatter(rx,ry)
    del(x, y, rx, ry, sufficient_data, fit)
    print('------------------------------plotting subtracted  ADC------------------------------')
    fit = {'ind':0, 'parameters':[], 'evt_count':[], 'reduced_chi':[],
           'x_y_factor':[], 'plot_cut':[], 'plot_coord':[]}
    cx,cy = uni_load_scatter(rib)
    # print('X data range: ', np.min(cx),' ~ ', np.max(cx), "; ", len(cx), ' events')
    # print('Y data range: ', np.min(cy),' ~ ', np.max(cy), "; ", len(cy), ' events')
    for x_factor in x_factor_list:
        for y_factor in y_factor_list:
            x,y,sufficient_data = slant_data_cut(rib, cx,cy,x_factor,y_factor,fit)
            # x,y,sufficient_data = slant_data_cut2(ax2, cx,cy,x_factor,y_factor, fit)
            # print('X_low cut: {}, Y_low cut {}, x_factor: {}, y_factor: {}, evt_num: {}'.format(
            #     int(np.median(x)-x_factor*np.std(x)),int(np.median(y)-y_factor*np.std(y)), x_factor, y_factor, len(x)))
            fit_plot(ax2,x,y,sufficient_data,fit)
            fit['ind']+=1
    sufficient_data = choose_fit(ax2, 'c', fit, n, gain_list[rib])
    plot_labels(ax2,'c', sufficient_data)
    # ax2.scatter( cx, cy)
    del(x, y, cx, cy, sufficient_data, fit)
    fig.tight_layout()
    fn = "scatter_{}_{}".format(merge[rib][0],merge[rib][1])
    plt.savefig(fn)
    print('get gain3/plots/'+fn+".png")
    plt.close('all')


with open('gain_table_L0', 'wb') as fp:
    pickle.dump(gain_list, fp)

