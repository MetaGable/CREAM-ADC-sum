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
    percent = ('{0:.' + str(decimal) + 'f}').format(100 * (iteration/float(total)))
    filled_length = int(length*iteration//total)
    bar = '>' * filled_length + '-' * (length-filled_length)
    print(f'\r {prefix} |{bar}| {percent}% {suffix}',end='\r')
    if iteration == total:
        print()


def lr(layer,ribbon):
    '''converts (layer, ribbon) into (the number of the ribbon from 0~999)'''
    return layer*50+ribbon-51


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
        '''define coordinates for the slanted cut'''
    x, y = scatter[0][mask1], scatter[1][mask1]
    x_factor, y_factor = 3, 3
    low_x, low_y = np.median(x)-x_factor*np.std(x),np.median(y)-y_factor*np.std(y)
    mask2 = [(i[0]>low_x and i[1]>low_y) for i in np.transpose([x,y])]
    x,y = x[mask2],y[mask2]
    # mid_x, mid_y = (np.max(x)+np.min(x))/2, min(y)
    # mask3 = [merge[ribbon_num][4]*(i[0]-mid_x)<(i[1]-mid_y) for i in np.transpose([x,y])]
    # x,y = x[mask3],y[mask3]
    return x,y


def slant_data_cut(rib, rx, ry, x_factor, y_factor, fit):
    '''Don't plots the cuts'''
    low_x, low_y = np.median(rx)-x_factor*np.std(rx),np.median(ry)-y_factor*np.std(ry)
    sufficient_data = True
    '''retangular cut at bottom left'''
    mask = [(i[0]>low_x and i[1]>low_y) for i in np.transpose([rx,ry])]
    x,y = rx[mask],ry[mask]
    if sum(mask)<30:
        sufficient_data=False
        print('not enough data')
        return x,y,sufficient_data
    '''define coordinates for the slanted cut'''
    mid_x, mid_y = np.min(x)+(np.max(x)-np.min(x))/6, (np.max(y)+np.min(y))/2-(np.max(y)-np.min(y))/3
    mask = [merge[rib][4]*(i[0]-mid_x)<(i[1]-mid_y) for i in np.transpose([x,y])]
    if sum(mask)<30:
        sufficient_data=False
        print('not enough data')
        return x,y,sufficient_data
    x,y = x[mask],y[mask]
    # print(f"Data cut: {len(rx)} -> {len(x)} events left")
    intercept_x = (low_y-min(y))/merge[rib][4]+mid_x
    fit['plot_cut'].append((low_y,low_x, x_factor, y_factor))
    fit['plot_coord'].append((intercept_x, np.min(x), np.max(x), mid_x, np.min(y), merge[rib][4], np.max(y)))
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
            _=ax1.hlines(low_y, min_x, intercept_x, alpha = 0.6, color=color_list[num%9], linestyle='--', label = y_factor)
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
    else:
        print(len(x)," points;    not enough data to fit")


def choose_fit(ax, x, y, plot_type, fit, total_plot_num, gain_list):
    if fit['evt_count']==[]: 
        print('No fit to choose (all data cut away)')
        sufficient_data = False
        return sufficient_data, gain_list
    _ = ax.scatter(x, y, alpha = 0.6,s=13)
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
        gain_list[rib].append((m, reduced_chi[ind], low_y, low_x, rib))
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
    return sufficient_data, gain_list


def plot_labels(ax, plot_type='r', sufficient_data=True):
    _=ax.legend(loc=2, fontsize='small',framealpha=0.4)
    _=ax.set_ylabel('low-ranged ch ADC',fontdict = {'fontsize' : 16})
    _=ax.set_xlabel('mid-ranged ch ADC',fontdict = {'fontsize' : 16})
    _=ax.grid(alpha=0.5)
    if plot_type=='r':
        ''' TO DO: make x & y factor as input'''
        if not sufficient_data:
            ax.text(0.5, 0.5, 'not enough data to fit ')
            return
        _=ax.set_title("With pedestal intact",fontdict = {'fontsize' : 17})
    elif plot_type=='c':
        _=ax.set_title("Subtracted pedestal with Cmean",fontdict = {'fontsize' : 17})
        if not sufficient_data:
            ax.text(0.5, 0.5, 'not enough data to fit ')
            return



"""For potential improvement, Jun 6th"""
missing_ribs=[303, 320, 330, 350, 363, 400, 415, 418, 425, 442, 448, 449, 459, 461, 469, 479, 481, 482, 485, 493, 498, 499, 500,
    506, 508, 509, 512, 519, 531, 535, 559, 565, 568, 572, 580, 582, 585, 587, 599, 601, 606, 611, 620, 638, 645, 649,
    653, 671, 673, 680, 681, 682, 687, 688, 709, 711, 715, 718, 719, 720, 735, 738, 747, 765, 766, 785, 791, 801, 809,
    811, 813, 821, 822, 843, 844, 850, 854, 864, 868, 870, 872, 873, 874, 875, 879, 883, 889, 893, 896, 900, 910, 911,
    918, 920, 930, 935, 937, 945, 954, 955, 957, 958, 959, 962, 971, 977, 978, 981, 987, 988, 990, 991, 993, 994, 997]
# y_factor_list=[-4, -8, -12]
# x_factor_list=[1.2, -0.5]
y_factor_list=[-2, -4, -8, -12]
x_factor_list=[1.2, -0.5, -1, -1.5]
n = 2
if not os.path.exists(table_path[where]):
    with open(table_path[where], 'wb') as fp:
        gain_list = [[] for i in range(1000)]
        pickle.dump(gain_list, fp)
else:
    gain_list = pickle.load(open(table_path[where],'rb'))



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
    sufficient_data, gain_list = choose_fit(ax1, rx, ry, 'r', fit, n, gain_list)
    plot_labels(ax1,'r',sufficient_data)
    # ax1.scatter(rx,ry)
    del(x, y, rx, ry, sufficient_data, fit)
    print('--------- plotting subtracted  ADC ', rib,'-------')
    fit = {'ind':0, 'parameters':[], 'evt_count':[], 'reduced_chi':[],
           'x_y_factor':[], 'plot_cut':[], 'plot_coord':[]}
    cx,cy = uni_load_scatter(rib)
    for x_factor in x_factor_list:
        for y_factor in y_factor_list:
            x,y,sufficient_data = slant_data_cut(rib, cx,cy,x_factor,y_factor,fit)
            # x,y,sufficient_data = slant_data_cut2(ax2, cx,cy,x_factor,y_factor, fit)
            fit_plot(ax2,x,y,sufficient_data,fit)
            fit['ind']+=1
    sufficient_data, gain_list = choose_fit(ax2, cx, cy, 'c', fit, n, gain_list)
    plot_labels(ax2,'c', sufficient_data)
    # ax2.scatter(cx, cy)
    del(x, y, cx, cy, sufficient_data, fit)
    fig.tight_layout()
    saved_gain_list = pickle.load(open(table_path[where],'rb'))
    saved_gain_list[rib]=gain_list
    with open(table_path[where], 'wb') as fp:
        pickle.dump(gain_list, fp)
    fn = "scatter_{}_{}F".format(merge[rib][0],merge[rib][1])
    plt.savefig(fn)
    # print('###################################\n get full_auto/'+fn+".png \n###################################")
    plt.close('all')



# for i in missing_ribs:
#     print(f'get full_auto/scatter_{ind_lr(i)[0]}_{ind_lr(i)[1]}.png')

# saved_gain_list = pickle.load(open("/home/genec420/gain4/gain_table_L0",'rb'))


# with open('gain_table_L0', 'wb') as fp:
#     pickle.dump([], fp)

# with open('gain_table_L0', 'wb') as fp:
#     pickle.dump(saved_gain_list.append(gain_list), fp)

# with open('gain_table_L0', 'wb') as fp:
#     pickle.dump(gain_list, fp)
'''Mar14th (already done 7, 23, 602, 736, 840, 949)'''
# missing_ribs = [46, 135, 136, 141, 148, 149, 150, 164, 166, 168, 178, 180
#                  , 194, 203, 207, 231, 239, 247, 249, 256, 260, 296, 298, 303, 
#                  320, 330, 350, 363, 400, 415, 418, 425, 442, 448, 449, 459, 
#                   461, 469, 479, 481, 482, 485, 493, 498, 499, 500, 506, 508,
#                  509, 512, 519, 531, 535, 559, 565, 568, 572, 580, 582, 585,
#                  587, 599, 601, 606, 611, 620, 638, 645, 649, 653, 671, 673,
#                  680, 681, 682, 687, 688, 709, 715, 718, 719, 720, 735, 738,
#                  747, 765, 766, 785, 791, 801, 809, 811, 813, 821, 822, 843, 
#                  844, 850, 854, 868, 870, 872, 873, 874, 875, 879, 883, 889,
#                  893, 896, 900, 910, 911, 918, 920, 930, 935, 937, 945, 954,
#                  955, 957, 958, 959, 962, 971, 977, 978, 981, 987, 988, 990,
#                  991, 993, 994]
# missing_ribs = [131, 137, 138, 139, 140, 147, 153, 158, 160, 181]
# missing_ribs = [23, 52, 54, 57, 59, 61, 62, 63, 65,
#                 67, 68, 69, 70, 72, 74, 75, 77, 78, 79, 80, 81, 83, 85, 87, 89, 
#                 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104,
#                 106, 107, 108, 110, 111, 112, 113, 115, 117, 118, 120, 123, 126,
#                 127, 129, 130, 133, 145, 272, 636, 704, 736, 840, 852, 915, 921,
#                 931, 933, 949, 969, 985, 989]
# missing_ribs = [ 76, 800, 804, 806, 984]#May 15th
# missing_ribs = [131, 137, 138, 139, 140, 147, 153, 158, 160, 181,
#  184, 185, 186, 187, 195, 201, 209, 211, 215, 217, 227, 229, 233, 237, 241,
#  268, 281, 299, 300, 301, 308, 325, 344, 364, 367, 401, 402, 404, 407, 423,
#  429, 446, 453, 463, 487, 513, 518, 520, 542, 544, 546, 555, 593, 604, 608,
#  612, 615, 622, 623, 626, 655, 665, 668, 679, 689, 711, 731, 745, 746, 760,
#  764, 788, 790, 803, 812, 825, 826, 831, 840, 842, 853, 864, 878, 880, 881,
#  884, 888, 890, 891, 895, 897, 898, 923, 928, 932, 944, 948, 956, 996, 997]
