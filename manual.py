#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6th
Modified last on Feb 26th
Made to fit reconstructed data manually
@author: gene
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob as glob
import pickle
import os
merge = pickle.load(open("/home/genec420/misc/merge",'rb'))
def load_scatter(ribbon_num,which = 'C'):  
    '''First data cut'''
    x_factor, y_factor = 3, 3
    scatter=[[],[]]
    if which =='R':
        input_type = '/home/genec420/gain2/no*'
    else:
        input_type = '/home/genec420/gain2/C*'
    input_paths = glob.glob(input_type)
    for path in input_paths:
        input_data = pickle.load(open(path, 'rb'))
        scatter[0]+=input_data[ribbon_num][1]
        scatter[1]+=input_data[ribbon_num][0]
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
    # print('{} out of {} events loaded'.format(len(x),sum(mask1)))
    return x,y

def data_cut(rib, ax, x, y, fit, Cut_info, Slope_info):
    LCut, RCut, BCut, TCut = Cut_info
    lower_slope, upper_slope = Slope_info[0], Slope_info[3]
    mask = [(i[0]>LCut and i[0]<RCut and i[1]>BCut and i[1]<TCut
             ) for i in np.transpose([x,y])]
    fit['plot_cut'].append((BCut, LCut))
    x,y = x[mask],y[mask]
    print('{} data points left after first cut'.format(len(x)))
    _=ax.hlines([BCut,TCut], LCut, RCut, alpha = 0.5, linestyle='--')
    _=ax.vlines([LCut,RCut], BCut, TCut, alpha = 0.5, linestyle='--')
    line_x = np.linspace(np.min(x),np.max(x),10)
    if bool(lower_slope):
        lower_slope=float(lower_slope)
        lower_point_x, lower_point_y = Slope_info[1:3]
        upper_point_x, upper_point_y = Slope_info[4:]
        mask = [lower_slope*(i[0]-lower_point_x)<(i[1]-lower_point_y) for i in np.transpose([x,y])]
        x,y = x[mask],y[mask]
        low_line_y = [lower_slope*(i-lower_point_x)+lower_point_y for i in line_x]
        _ = ax.plot(line_x, low_line_y, alpha = 0.5,linestyle='--')
    if bool(upper_slope):
        upper_slope=float(upper_slope)
        mask = [upper_slope*(i[0]-upper_point_x)>(i[1]-upper_point_y) for i in np.transpose([x,y])]
        x,y = x[mask],y[mask]
        hi_line_y = [upper_slope*(i-upper_point_x)+upper_point_y for i in line_x]
        _ = ax.plot(line_x, hi_line_y, alpha = 0.5,linestyle='--')
    print('{} data points left finally'.format(len(x)))
    return x, y

def fit_plot(ax,x,y, fit, gain_list):
    (m, b), [residual], _, _, _ = np.polyfit(x,y,1,cov=False,full=True)
    red_chi = round(residual/(len(x)-2))
    plot_x = np.linspace(np.min(x),np.max(x),10)
    _ = ax.plot(plot_x,m*plot_x+b,alpha = 0.6, 
            label=r'${2}events|Gain:{0}|\chi^2_\nu:{1}$'.format(
            round(m,1),red_chi ,len(x)))
    BCut, LCut = fit['plot_cut'][0]
    gain_list[ribbon_ind] = (m, red_chi, LCut, BCut, ribbon_ind)
    '''store the fit info for selection'''

def plot_labels(ax, x, y, title, low_x='', xlim='', tightcuts=False):
    x_factor, y_factor = 1.5, 1
    _=ax.set_title(title,fontdict = {'fontsize' : 17})
    _=ax.legend(loc=2, fontsize='small',framealpha=0.4)
    _=ax.set_ylabel('low-ranged ch ADC',fontdict = {'fontsize' : 16})
    _=ax.set_xlabel('mid-ranged ch ADC',fontdict = {'fontsize' : 16})
    _=ax.grid(alpha=0.5)
    _=ax.scatter(x, y, alpha = 0.6, s=13)
    low_y = np.median(y)-y_factor*np.std(y)
    # if low_x !='':
    #     _=ax.set_xlim(left=float(low_x))
    # else:
    #     temp = np.median(x)-x_factor*np.std(x)
    #     _=ax.set_xlim(left = temp)
    # if xlim != '':
    #     _=ax.set_xlim(right=float(xlim))
    _=ax.set_ylim(bottom = low_y)
    if tightcuts!=False:
        LCut, RCut, BCut, TCut = tightcuts
        LCut-=15
        RCut+=25
        TCut+=75
        _=ax.set_xlim(LCut, RCut)
        _=ax.set_ylim(top = TCut)

def specify_input():
    LCut = float(input('Min x:\n'))
    RCut = float(input('Max x:\n'))
    BCut = float(input('Min y:\n'))
    TCut = float(input('Max y:\n'))
    lower_slope = input('Slope of the lower slant cut:\n (if no change needed, press enter) \n')
    if lower_slope != '':
        lower_point_x = float(input('x coordinate of the point of the lower slant cut:\n'))
        lower_point_y = float(input('y coordinate of the point of the lower slant cut:\n'))
    else:
        lower_point_x, lower_point_y = '', ''
    upper_slope = input('Slope of the upper slant cut:\n (if no change needed, press enter) \n')
    if upper_slope != '':
        upper_point_x = float(input('x coordinate of the point of the upper slant cut:\n'))
        upper_point_y = float(input('y coordinate of the point of the upper slant cut:\n'))
    else:
        upper_point_x, upper_point_y = '', ''
    Cut_info=[LCut, RCut, BCut, TCut]
    Slope_info = [lower_slope, lower_point_x, lower_point_y,
                  upper_slope, upper_point_x, upper_point_y]
    return Cut_info, Slope_info


layer = int(input('Enter layer number (1~20):\n'))
ribbon = int(input('Enter ribbon number (1~50):\n'))
which = input('Use which plots? (l/r/both)\n')   
plot1 = False
plot2 = False
if which == 'l' or which =='both':
    plot1 = True
    print('For the Left plot:')
    L_Cut_info, L_Slope_info = specify_input()
if which == 'r' or which =='both':
    plot2 = True    
    print('For the Right plot:')
    R_Cut_info, R_Slope_info = specify_input()
print('')
# Llow_x = input('Formatting the Left Plot Limit\n (if no change needed, press enter):\n minimum x:\n')
# Lxlim = input('maximum x:\n')
# Rlow_x = input('Formatting the Right Plot Limit\n (if no change needed, press enter):\n minimum x:\n')
# Rxlim = input('maximum x:\n')
ribbon_ind = ribbon + 50 * layer - 51
if os.path.exists('gain_table_L0'):
    gain_list = pickle.load(open("gain_table_L0",'rb'))
else:
    gain_list = [[] for i in range(1000)]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
_=fig.suptitle('''layer {0} ribbon {1}: ch{2} mid-range ADC over ch{3} low-range ADC
               From gain table: Gain = {4}, threshold = {5}'''.format(
               merge[ribbon_ind][0],merge[ribbon_ind][1],merge[ribbon_ind][2],
               merge[ribbon_ind][3],merge[ribbon_ind][4],merge[ribbon_ind][5]),fontsize=19)
rx, ry=load_scatter(ribbon_ind, 'R')
if plot1:
    fit = { 'parameters':[], 'reduced_chi':[], 
       'x_y_factor':[], 'plot_cut':[], 'plot_coord':[]}
    print('Left plot (Ped intact)')
    x, y = data_cut(ribbon_ind, ax1, rx, ry, fit, L_Cut_info, L_Slope_info)
    fit_plot(ax1,x,y,fit, gain_list)
    del(x,y,fit)
# plot_labels(ax1,rx, ry, "With pedestal intact", Llow_x, Lxlim)
plot_labels(ax1,rx, ry, "With pedestal intact", tightcuts=L_Cut_info)
del(rx, ry)
cx,cy=load_scatter(ribbon_ind, 'C')
if plot2:
    print('Right plot (Ped subtracted)')
    fit = {'ind':0, 'parameters':[], 'evt_count':[], 'reduced_chi':[],
       'x_y_factor':[], 'plot_cut':[], 'plot_coord':[]}
    x, y= data_cut(ribbon_ind, ax2, cx, cy, fit, R_Cut_info, R_Slope_info)
    fit_plot(ax2,x,y,fit,gain_list)
    del(x,y,fit)
# plot_labels(ax2, cx, cy, "Subtracted pedestal with Cmean", Rlow_x, Rxlim)
plot_labels(ax2, cx, cy, "Subtracted pedestal with Cmean", tightcuts=R_Cut_info)
del(cx, cy)
fig.tight_layout()
fn = "scatter_{}_{}".format(merge[ribbon_ind][0],merge[ribbon_ind][1])
plt.savefig(fn)
print('#####################################\nget '+fn+".png \n#####################################")
plt.close('all')
with open('gain_table_L0', 'wb') as fp:
    pickle.dump(gain_list, fp)



