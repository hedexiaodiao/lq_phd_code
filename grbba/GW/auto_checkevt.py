#!/hxmt/soft/Develop/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 14:31:32 2019

@author: Luoqi
"""
import subprocess
import numpy as np
import astropy.time
from astropy.time import Time
import datetime
from datetime import date
import argparse
import os
#import sys
import time

def run_task(ra,dec,trigtime,obsmode,t,t_bkg_l_b,t_bkg_l_f,t_bkg_r_b,t_bkg_r_f,dirstr):
    print('run_task')
    if os.path.exists('/sharefs/hbkg/user/luoqi/GRB/grbba/GW/'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr))==0:
        os.mkdir('/sharefs/hbkg/user/luoqi/GRB/grbba/GW/'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr))
    os.chdir('/sharefs/hbkg/user/luoqi/GRB/grbba/GW/'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr))
    if os.path.exists('./t'+str(int(t)))==0:
        os.mkdir('t'+str(int(t)))
    
    try:
        #write run_spec.pro
        m1 = 'make_pha_v1, t_trig_utc=\''+str(trigtime.value)+'\', t_bkg_l_b='+str(t_bkg_l_b)+',t_bkg_l_f='+str(t_bkg_l_f)+',t_src_b=0'+\
            ',t_src_f='+str(t)+',t_bkg_r_b='+str(t_bkg_r_b)+',t_bkg_r_f='+str(t_bkg_r_f)+'\n'
        n1 = 'exit\n'
        m1file = open('./t'+str(int(t))+'/run_spec'+str(int(t))+'.pro','a+')
        m1file.write(m1)
        m1file.write(n1)
        m1file.close()
       
        #write submitwork
        l1 = '#!/bin/bash\n'
        q1 = 'cd /sharefs/hbkg/user/luoqi/GRB/grbba/GW/'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr)+'/t'+str(int(t))+'/\n'
        ll1 = 'source /hxmt/work/GRB/ENV/.bashrcBA\n'
        qq1 = 'idl run_spec'+str(int(t))+'.pro\n'
        
        l1file = open('subspec'+str(int(t))+'.sh','a+')
        l1file.write(l1)
        l1file.write(q1)
        l1file.write(ll1)
        l1file.write(qq1)
        l1file.close()
        
        subprocess.Popen('chmod +x subspec'+str(int(t))+'.sh',shell=True)
        subprocess.Popen('sh subspec'+str(int(t))+'.sh',shell=True)
        
        if obsmode == 'GRB':
            os.system('/hxmt/work/GRB/Software/RSPgenerator_v9/HXMTrspg_v9.py --obsmode '+str(obsmode)+' --ra-deg '+str(ra)+' --dec-deg '+str(dec)+' --starttime '+str(trigtime.value)+' --ecfile CsI_GRB_Orbit_v1_EC.txt --fwhmfile CsI_GRB_Orbit_v1_RasE.txt')
        else:
            os.system('/hxmt/work/GRB/Software/RSPgenerator_v9/HXMTrspg_v9.py --obsmode '+str(obsmode)+' --ra-deg '+str(ra)+' --dec-deg '+str(dec)+' --starttime '+str(trigtime.value))
        
        l3 = '#!/bin/bash\n'
        q3 = 'cd /sharefs/hbkg/user/luoqi/GRB/grbba/GW/'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr)+'/\n'
        #ll3 = 'source /hxmt/work/GRB/ENV/.bashrcBA\n'
        qq3 = 'python /sharefs/hbkg/user/luoqi/GRB/grbba/GW/auto_calclim.py --ra-deg '+str(ra)+' --dec-deg '+str(dec)+' --trigtime '+str(trigtime.value)+\
        ' --obsmode '+str(obsmode)+' --t '+str(t)+' --t_bkg_l_b '+str(t_bkg_l_b)+' --t_bkg_l_f '+str(t_bkg_l_f)+' --t_bkg_r_b '+str(t_bkg_r_b)+' --t_bkg_r_f '+str(t_bkg_r_f)+' --dirstr '+str(dirstr)+'\n'

        l3file = open('sub_auto_ca'+str(int(t))+'.sh','a+')
        l3file.write(l3)
        l3file.write(q3)
        l3file.write(qq3)
        l3file.close()
        subprocess.Popen('chmod +x sub_auto_ca'+str(int(t))+'.sh', shell=True)
        subprocess.Popen('sh sub_auto_ca'+str(int(t))+'.sh', shell=True)
    except Exception as e:
        print(str(e))
    return 0
        
def check_file(ra,dec,trigtime,obsmode,t,t_bkg_l_b,t_bkg_l_f,t_bkg_r_b,t_bkg_r_f,dirstr):
    l = os.listdir('/hxmt/work/HXMT-DATA/1K/Y'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'/')
    len_l = len(l)
    if len_l!=0:
        k = -1
        for i in range(0,len_l):
            key1 = l[i]
            p1 = str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)
            result1 = p1 in key1
            if result1:
                k = i
        if k != -1:
            print('change dir /hxmt/work/HXMT-DATA/1K/Y'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'/'+\
                  l[k])
            os.chdir('/hxmt/work/HXMT-DATA/1K/Y'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'/'+\
                  l[k])
            print('try to find HXMT_'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+'T'+\
                              '{:02d}'.format(trigtime.datetime.hour)+'_HE-Evt_FFFFFF_V1_1K.FITS')      
            if os.path.exists('HXMT_'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+'T'+\
                              '{:02d}'.format(trigtime.datetime.hour)+'_HE-Evt_FFFFFF_V1_1K.FITS')==1:
                print('1K_EVT file exist')
                print('try to find HXMT_'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+'T'+\
                              '{:02d}'.format(trigtime.datetime.hour)+'_Att_FFFFFF_V1_1K.FITS')
                ver_att = 1
                while ver_att<6:
                    if os.path.exists('HXMT_'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+'T'+\
                                  '{:02d}'.format(trigtime.datetime.hour)+'_Att_FFFFFF_V'+str(int(ver_att))+'_1K.FITS')==1:
                        print('Att file exist')
                        flag = run_task(ra,dec,trigtime,obsmode,t,t_bkg_l_b,t_bkg_l_f,t_bkg_r_b,t_bkg_r_f,dirstr)
                        break
                    else:
                        ver_att += 1
                        flag = 1
            else:
                print('1K_EVT file not exist')
                flag = 1
        else:
            print('1K_EVT file not exist')
            flag = 1
    else:
        print('1K_EVT file not exist')
        flag = 1
    return flag

def timerun(ra,dec,trigtime,obsmode,t,t_bkg_l_b,t_bkg_l_f,t_bkg_r_b,t_bkg_r_f,dirstr):
    
    j = 0
    flag = 1
    while flag==1:
        now = datetime.datetime.now()
        print(now,j)        
        flag = check_file(ra,dec,trigtime,obsmode,t,t_bkg_l_b,t_bkg_l_f,t_bkg_r_b,t_bkg_r_f,dirstr)
        j = j+1
        time.sleep(180)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--ra-deg', action='store',
                        dest='radeg',
                        help='input Incident Ra(0-360 deg)',
                        type=float,default=None)

    parser.add_argument('--dec-deg', action='store',
                        dest='decdeg',
                        help='input Incident Dec(-90-90 deg)',
                        type=float,default=None)

    parser.add_argument('--trigtime',
                        dest='trigtime',
                        help='input start time (format:YYYY-MM-DDThh:mm:ss)',
                        type=str,
                        default=None)
    parser.add_argument('--obsmode',
                        dest='obsmode',
                        help='input {GRB} or {Nor}',
                        type=str,
                        default='Nor')
    parser.add_argument('--t',
                        dest='t',
                        help='input a timescale(e.g. 1 10)',
                        type=float,
                        default='1')   
    parser.add_argument('--t_bkg_l_b',
                        dest='t_bkg_l_b',
                        help='input a bkg_left_back(e.g. -120),default is -120',
                        type=float,
                        default='-120')
    parser.add_argument('--t_bkg_l_f',
                        dest='t_bkg_l_f',
                        help='input a bkg_left_front(e.g. -20),default is -20',
                        type=float,
                        default='-20')         
    parser.add_argument('--t_bkg_r_b',
                        dest='t_bkg_r_b',
                        help='input a bkg_right_front(e.g. 20),default is 20',
                        type=float,
                        default='20')
    parser.add_argument('--t_bkg_r_f',
                        dest='t_bkg_r_f',
                        help='input a bkg_right_front(e.g. 20),default is 120',
                        type=float,
                        default='120')                                                         
    parser.add_argument('--dirstr',
                        dest='dirstr',
                        help='input a token to distinguish dir',
                        type=str,
                        default='_1')
    args = parser.parse_args()
   
    calsttime = time.time()
    
    
    ra=args.radeg
    dec=args.decdeg
    trigtime=args.trigtime
    obsmode = args.obsmode
    t = args.t
    t_bkg_l_b = args.t_bkg_l_b
    t_bkg_l_f = args.t_bkg_l_f
    t_bkg_r_b = args.t_bkg_r_b
    t_bkg_r_f = args.t_bkg_r_f
    dirstr = args.dirstr
    trigtime=Time(trigtime,scale='utc',format='isot')

    if os.path.exists('/sharefs/hbkg/user/luoqi/GRB/grbba/GW/'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr))==1:
        print('please change the argument \'--dirstr\'!!!\nor if spectum files exists in this dir already,\nyou can use command:\n'+\
              'python /sharefs/hbkg/user/luoqi/GRB/grbba/GW/auto_calclim.py --ra-deg '+str(ra)+' --dec-deg '+str(dec)+' --trigtime '+str(trigtime.value)+\
              ' --obsmode '+str(obsmode)+' --t '+str(t)+' --t_bkg_l_b '+str(t_bkg_l_b)+' --t_bkg_l_f '+str(t_bkg_l_f)+' --t_bkg_r_b '+str(t_bkg_r_b)+' --t_bkg_r_f '+str(t_bkg_r_f)+' --dirstr '+str(dirstr)+'_1'+'\n'+\
              'to calculate the limits directly!!!')
    else:
        timerun(ra,dec,trigtime,obsmode,t,t_bkg_l_b,t_bkg_l_f,t_bkg_r_b,t_bkg_r_f,dirstr)

    caledtime = time.time()
    #print('time cost:{:.5f} seconds'.format(caledtime-calsttime))
    print('this procedure is still running, \n\
          only if you can find results'+str(int(t))+'.txt in dir you set,\n\
            or received the Email can you input other shell command')
    
    
