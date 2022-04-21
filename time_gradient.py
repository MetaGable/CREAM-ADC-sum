#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 2022, ver1

Used to identify the time range for events in
low vs mid ADC scatter plots. 

@author: gene
"""


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import csv
import numpy as np
import glob as glob
import pickle
import datetime as dt
import matplotlib.dates as mdates
merge = pickle.load(open("/home/genec420/misc/merge",'rb'))
def lr(layer,ribbon):
    '''converts (layer, ribbon) into (the number of the ribbon from 0~999)'''
    return layer * 50 + ribbon - 51


def load_bar(iteration, total, prefix='', suffix='', decimal=1, length=90):
    percent = ('{0:.' + str(decimal) + 'f}').format(100 * (iteration/float(total)))
    filled_length = int(length*iteration//total)
    bar = '>' * filled_length + '-' * (length-filled_length)
    print(f'\r {prefix} |{bar}| {percent}% {suffix}',end='\r')
    # print(f'\r{prefix} |{bar}| {percent}% {suffix}',end='\r')
    if iteration == total:
        print()

data_type = 'full'
directory={'track':'gain2', 'full':'gain3', 'R':'no*', 'C':"C*"}
def load_scatter(ribbon_num, which = 'C'):  
    if os.path.exists(f"GradientV2_{ribbon_num:03d}{which}"):
        x, y, time, c = pickle.load(open(f"GradientV2_{ribbon_num:03d}{which}",'rb'))
        print('------- quick load complete -------')
    else:
        scatter, time, c = [[],[]], [], []
        input_paths = glob.glob(f'/home/genec420/{directory[data_type]}/{directory[which]}')
        input_paths = sorted(input_paths, key=lambda x: int("".join([i for i in x if i.isdigit()])))
        total = len(input_paths)
        for num, path in enumerate(input_paths):
            input_data = pickle.load(open(path, 'rb'))
            length = len(input_data[ribbon_num][0])
            scatter[0]+=input_data[ribbon_num][1]
            scatter[1]+=input_data[ribbon_num][0]
            time += [path[-15: ]]*length
            c += [num]*length
            load_bar(num+1, total)
        scatter = np.array(scatter)
        c = np.array(c)/(np.float(num)+1)
        if which =='C':
            mask1 = [(i[0]!=0 and i[1]!=0)and (i[0]<30000) for i in np.transpose(scatter)]
        else:
            mask1 = [(i[0]>0 and i[1]>0) for i in np.transpose(scatter)]
        scatter = np.array(scatter)
        x, y = scatter[0][mask1], scatter[1][mask1]
        time = np.array(time)[mask1]
        c=np.array(c)[mask1]
        with open("GradientV2_"+str(ribbon_num).zfill(3)+str(which), 'wb') as fp:
           pickle.dump([x, y, time, c], fp)
        print('------ quick load created-------')
    return x, y, time, c

def color_date(ax, x, y, c0, dates, windows):
    mask = [1]*len(dates)
    for (start_day, start_plus, end_day, end_plus) in windows:
        indicies_range, cache_timestamp, flag =[], 0, False
        for ind, date in enumerate(dates):
            if (date[0:8]==start_day) and start_plus!=-1:
                if (not flag) & (cache_timestamp!=date):
                    start_plus -= 1
                    cache_timestamp=date
                    if start_plus==0:
                        print(f'start at {date}, #{ind:7d}')
                        flag = True
                        start_plus = -1
                        # indicies_range.append((ind,date))
                        ind1 = ind
            if end_plus!=-1 and (date[0:8] == end_day):
                if flag & (cache_timestamp!=date):
                    end_plus -= 1
                    cache_timestamp=date
                    if flag & end_plus==0:
                        '''Need padding'''
                        print(f'end at   {date}, #{ind:7d}')
                        flag = False
                        end_plus = -1
                        # indicies_range.append((ind,date))
                        ind2 = ind
                        '''Plot events in the window'''
                        label = f'{start_day}~{end_day[-4:]}'
                        _ = ax.scatter(x[ind1:ind2],y[ind1:ind2],c=c0[in1:ind2], marker=2,s=3,label=label)
                        break
        mask[ind1:ind2]=[0]*(ind2-ind1)
    '''Plot events outside the window'''
    mask = [bool(i) for i in mask]
    _ = ax.scatter(x[mask], y[mask], alpha = 0.6, s = 11, c = c0[mask], cmap = 'gist_ncar', label='Rest of the data')    


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
        mask = [lower_slope*(i[0]-lower_point_x)<(i[1]-lower_point_y) for i in np.transpose([x,y])]
        x,y = x[mask],y[mask]
        low_line_y = [lower_slope*(i-lower_point_x)+lower_point_y for i in line_x]
        _ = ax.plot(line_x, low_line_y, alpha = 0.5,linestyle='--', color = 'tab:blue')
    if bool(upper_slope):
        upper_slope=float(upper_slope)
        upper_point_x, upper_point_y = Slope_info[4:]
        mask = [upper_slope*(i[0]-upper_point_x)>(i[1]-upper_point_y) for i in np.transpose([x,y])]
        x,y = x[mask],y[mask]
        hi_line_y = [upper_slope*(i-upper_point_x)+upper_point_y for i in line_x]
        _ = ax.plot(line_x, hi_line_y, alpha = 0.5,linestyle='--', color = 'tab:blue')
    print('{} data points left finally'.format(len(x)))
    return x, y

# 'gist_ncar', 'hsv', 
def plot_labels(ax, title, low_x='', xlim='', tightcuts=False):
    _=ax.set_title(title,fontdict = {'fontsize' : 17})
    _=ax.legend(loc=2, fontsize='small',framealpha=0.4)
    _=ax.set_ylabel('low-ranged ch ADC',fontdict = {'fontsize' : 16})
    _=ax.set_xlabel('mid-ranged ch ADC',fontdict = {'fontsize' : 16})
    _=ax.grid(alpha=0.5)
    # if tightcuts!=False:
    LCut, RCut, BCut, TCut = tightcuts
    LCut-=0.12*(RCut-LCut)
    RCut+=0.12*(RCut-LCut)
    TCut+=0.1*(TCut-BCut)
    BCut-=0.18*(TCut-BCut)
    _=ax.set_xlim(LCut, RCut)
    _=ax.set_ylim(BCut, TCut)


def specify_input():
    LCut = float(input('Min x:\n'))
    RCut = float(input('Max x:\n'))
    BCut = float(input('Min y:\n'))
    TCut = float(input('Max y:\n'))
    Cut_info=[LCut, RCut, BCut, TCut]
    return Cut_info

def auto():
    parameter_list=[]
    with open('/home/genec420/time_grad/TG.csv', newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            parameter_list.append(row)
    for raw_parameters in parameter_list:
        layer, ribbon = raw_parameters[0:2]
        print('-------------------------------------\nstart layer {}, ribbon {}:'.format(layer, ribbon))
        parameters = [float(i) if ((ind!=2) & (i!='')) else i for ind, i in enumerate(raw_parameters)]
        L_Cut_info, R_Cut_info = parameters[3:7], parameters[8:12]
        '''Plot individual time windows'''
        windows, ind_shift = [], 0
        while raw_parameters[13+ind_shift]!='' :
            start, end= raw_parameters[13+ind_shift:15+ind_shift]
            # if start_plus=='': start_plus=1
            # if end_plus=='': end_plus=1
            # start_plus, end_plus=int(start_plus), int(end_plus)
            start_plus, end_plus=1,1
            windows.append((start, start_plus, end, end_plus))
            length1 = len(glob.glob(f'/home/genec420/gain3/C{start}*'))
            length2 = len(glob.glob(f'/home/genec420/gain3/C{end}*'))
            print(f'{length1:3d} files in {start} and {length2:3d} files in {end}')
            ind_shift+=2
            print(len(raw_parameters)<14+ind_shift)
            if (len(raw_parameters)<13+ind_shift): break
        ribbon_ind = int(ribbon) + 50 * int(layer) - 51
        fig = plt.figure(dpi=200)
        fig.set_figheight(9)
        fig.set_figwidth(24)
        ax3 = plt.subplot2grid(shape=(12, 4), loc=(10, 0), rowspan=1, colspan=4)
        _=fig.suptitle('''{2} layer {0} ribbon {1};'''.format(
            merge[ribbon_ind][0],merge[ribbon_ind][1], data_type),fontsize=19)
        rx, ry, time, c0 = load_scatter(ribbon_ind, 'R')
        '''Color reference subplot'''
        dt_time = [dt.datetime.strptime(i,"%Y%m%d-%H%M%S") for i in time]
        ax3.scatter(dt_time,[0]*len(c0),marker='|',s=200,c=c0,cmap='gist_ncar')
        _=ax3.grid(alpha=0.5,axis='x')
        loc = mdates.DayLocator(bymonthday=(1, 16),interval=1)
        _ = ax3.get_xaxis().set_major_locator(loc)
        loc = mdates.DayLocator(bymonthday=(1, 9, 16, 24),interval=1)
        _ = ax3.get_xaxis().set_minor_locator(loc)
        ax3.tick_params(axis="x", rotation=90)
        # fmt = mdates.DateFormatter('%b\n%Y')
        # ax3.xaxis.set_major_formatter(fmt)
        '''Left Plot'''
        ax1 = plt.subplot2grid(shape=(12, 4), loc=(0, 0), rowspan=9, colspan=2)
        color_date(ax1, rx, ry, c0, time, windows)
        plot_labels(ax1, "With pedestal intact", tightcuts=L_Cut_info)
        del( rx, ry, time, c0)
        '''Right Plot'''
        cx,cy, time, c0 =load_scatter(ribbon_ind,'C')
        ax2 = plt.subplot2grid(shape=(12, 4), loc=(0, 2), rowspan=9,  colspan=2)
        color_date(ax2, cx, cy, c0, time, windows)
        plot_labels(ax2, "Subtracted pedestal with Cmean", tightcuts=R_Cut_info)
        del(cx, cy, time)
        fig.tight_layout()
        fn = "scatter_{:0>2}_{:0>2}".format(merge[ribbon_ind][0],merge[ribbon_ind][1])
        plt.savefig(fn)
        print('\nget time_grad/'+fn+".png \n")
        plt.close('all')



auto()


