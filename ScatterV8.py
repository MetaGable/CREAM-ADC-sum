#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 14:47:16 2021
Using tracking log
@author: gene
"""

import numpy as np
import glob as glob
import root_pandas as rpd
import datetime as dt
import pickle
import os

cache = -1
# np.array(pickle.load(open("normalC2017-08-22",'rb')))
merge = np.array(pickle.load(open("/home/genec420/misc/merge",'rb')),dtype='int')
all_ribbon = [(i[2], i[3], i[8]) for i in merge]
date_str,time_str,time_us = np.transpose(np.loadtxt("/home/genec420/misc/All_ntracked_with_evttime.log",dtype='str'))
track = np.loadtxt("/home/genec420/misc/All_ntracked_with_evttime.log",dtype='str')
layers = []
for i in range(0,1000,50):
    one_layer = []
    for j in range(i,i+50):
        one_layer.append(merge[:,2][j])
    layers.append(set(one_layer))


def get_path(date,file_type="spi"):
    Year = str(date.year)
    Month = str(date)[5:7]
    Day = str(date)[8:10]  
    if file_type=="spi":
        daily_path = glob.glob("/data/L0/playback/"+Year+"/"+Month+"/"+Day+"/*c.spi")
        daily_path = sorted(daily_path, key=lambda x: int("".join([i for i in x if i.isdigit()])))
    elif file_type=="dat":
        daily_path = glob.glob("/data/L0/playback/"+Year+"/"+Month+"/"+Day+"/*dat.root")
        daily_path = sorted(daily_path, key=lambda x: int("".join([i for i in x if i.isdigit()])))
    return daily_path



# start = "2018/08/07"
start = "2017/08/22"
end = "2019/02/11"
# end = "2018/12/31"
tracker_row = int(input('enter tracked row number:  '))
# tracker_row = 24382
date1 = dt.datetime.strptime(start,"%Y/%m/%d")
fn = 'no_ped'+str(date1)[0:10]
while os.path.exists(fn):
    date1 = date1 + dt.timedelta(days=1)
    fn = 'no_ped'+str(date1)[0:10]
    continue


date2 = dt.datetime.strptime(end,"%Y/%m/%d")
deltaDays = (date2-date1).days+1
yesterday =  date1 - dt.timedelta(days=1)
today_ped_path = get_path(yesterday)
tomorrow_ped_path = get_path(date1)
for i in range(deltaDays):
    fn = 'no_ped'+str(date1)[0:10]
    try:
        if date1!=dt.datetime.strptime(date_str[tracker_row],"%Y%m%d"):
            if cache == tracker_row:
                tracker_row+=1
            cache = tracker_row
            print('wrong day! tracker_row = {}, day: {} --> {}'.format(tracker_row,str(date1)[:10],date_str[tracker_row]))
            date1 = dt.datetime.strptime(date_str[tracker_row],"%Y%m%d")
    except:
        tracker_row+=1
        continue
    hist = [([],[]) for i in range(1000)]
    histC = [([],[]) for i in range(1000)]
    not_depleted = True
    hist_day, ped_iter, ped_midtime = [], 0, []
    yesterday_ped_path, today_ped_path = today_ped_path, tomorrow_ped_path
    tomorrow_ped_path = get_path(date1 + dt.timedelta(days=1),"spi")
    dat_path = get_path(date1,"dat")
    if len(today_ped_path)==0:
        date1 = date1 + dt.timedelta(days=1)
        print(str(date1)+" has no pedestal files")
        try:
            while date1!=dt.datetime.strptime(date_str[tracker_row],"%Y%m%d"):
                tracker_row+=1
        except:
            tracker_row+=1
        continue
    elif len(tomorrow_ped_path)==0:
        if len(yesterday_ped_path)==0:
            ped_path = today_ped_path
        else:
            ped_path = [yesterday_ped_path[-1]]+ today_ped_path
    else: 
        if len(yesterday_ped_path)==0:
            ped_path = today_ped_path + [tomorrow_ped_path[0]]
        else:
            ped_path = [yesterday_ped_path[-1]] + today_ped_path + [tomorrow_ped_path[0]]
    ped= [[int(j) for j in np.loadtxt(i,usecols = 5)] for i in ped_path]
    # ped= [np.loadtxt(i,usecols = (3,5)) for i in ped_path]
    try:
        data = rpd.read_root(dat_path, columns=['info*','bchCal','badcCal'])
    except:
        print("none of the data files can be opened")
        date1 = date1 + dt.timedelta(days=1)
        try:
            while date1!=dt.datetime.strptime(date_str[tracker_row],"%Y%m%d"):
                tracker_row+=1
        except:
            tracker_row+=1
        continue
    for i in range(len(ped_path)-1):
        ped_time1 = dt.datetime.strptime(ped_path[i][29:44],"%Y%m%d-%H%M%S")
        ped_time2 = dt.datetime.strptime(ped_path[i+1][29:44],"%Y%m%d-%H%M%S")
        diff_time = ped_time1+(ped_time2 - ped_time1)/2
        ped_midtime.append(diff_time)
    PedLen = len(ped_midtime)
    # print('start looping')
    for index, row in data.iterrows():
        # if row['info_trig'][7] == 1 or row['info_trig'][6] == 1:
        # if int(time_us[tracker_row])==row['info_us']:
        if int(time_us[tracker_row])==row['info_us']:
            # print('right microsecond')
            if not (row['info_trig'][7]==1 or row['info_trig'][8]==1):
                tracker_row+=1
                print('not desired event'+str(row['info_trig']))
                try:
                    _ = dt.datetime.strptime(date_str[tracker_row],"%Y%m%d")
                except:
                    tracker_row+=1
                    continue
            # print('right HMS')
            dat_time = int(str(row['info_time'][3])+str(row['info_time'][4])+str(row['info_time'][5]))
            evt_t = dt.datetime.strptime(str(row['info_time']),"[%Y   %m   %d   %H   %M   %S]")
            row_time = dt.datetime.strptime(time_str[tracker_row],"%H%M%S")
            if row_time.hour == evt_t.hour and row_time.minute==evt_t.minute and row_time.second ==evt_t.second:
                noise = False
                ch_list= []
                while not_depleted:
                    if (ped_midtime[ped_iter]-evt_t).total_seconds() < 0 and ped_iter < PedLen:
                        ped_iter+=1
                        if ped_iter < PedLen-1 and (ped_midtime[ped_iter]-evt_t).total_seconds() < 0:
                            ped_iter+=1
                            if ped_iter < PedLen-2 and (ped_midtime[ped_iter]-evt_t).total_seconds() < 0:
                                ped_iter+=1
                        if ped_iter == PedLen:
                            not_depleted = False
                    break
                for i,(ch1, ch2, ch3) in enumerate(all_ribbon): # Check for shortcut
                    if ch1 in row['bchCal'] and ch2 in row['bchCal']:
                        adc1 = row['badcCal'][np.where(row['bchCal']==ch1)[0][0]]
                        adc2 = row['badcCal'][np.where(row['bchCal']==ch2)[0][0]]
                        adc1c = adc1-ped[ped_iter][int(ch1)] 
                        adc2c = adc2-ped[ped_iter][int(ch2)]
                        ch_list.append((i,ch1,ch2,adc1,adc2,adc1c,adc2c))
                # if ch_list==[]:
                #     tracker_row+=1
                #     break
                for l, layer in enumerate(layers):
                    active_ch = 0
                    for j, Ch in enumerate(row['bchCal']):
                        if Ch in layer:
                            chADC= row['badcCal'][j]-ped[ped_iter][Ch] #to be added to ADC sum
                            if chADC>0:
                                active_ch += 1
                    if active_ch>23:
                        noise=True
                        print('Noise!')
                if noise == False:
                    for i,ch1,ch2,adc1,adc2,adc1c,adc2c in ch_list:
                        hist[i][0].append(adc1)
                        hist[i][1].append(adc2)
                        histC[i][0].append(adc1c)
                        histC[i][1].append(adc2c)
            tracker_row+=1
            print(tracker_row)
            try:
                _ = dt.datetime.strptime(date_str[tracker_row],"%Y%m%d")
            except:
                tracker_row+=1
    print(str(tracker_row)+'      '+str(dt.datetime.now())[11:16]+'     ' + str(date1)[0:11])
    with open(fn, 'wb') as fp:
        pickle.dump(hist, fp)
    fn = 'C'+str(date1)[0:10]
    with open('C'+str(date1)[0:10], 'wb') as fp:
        pickle.dump(histC, fp)
    date1 = date1 + dt.timedelta(days=1)



