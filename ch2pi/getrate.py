#!/hxmt/soft/Develop/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 17:23:16 2019

@author: Luoqi
"""




import numpy as np
import argparse
import os
from scipy.interpolate import interp1d
import time
from functools import partial
from astropy.io import fits
import warnings
from astropy.time import Time

serverlibDIR = '/hxmt/work/GRB/Software/RSPgenerator_v9/'

workdir = os.getcwd()
Chan_l_A = np.zeros((18,256))
Chan_h_A = np.zeros((18,256))

dets = [0,1,2,3,4,5,6,7,8,9,10\
        ,11,12,13,14,15,16,17]

theta = 45
phi = 90
ra = 45
dec = 90
nchn_r = 20 #0:nchn_r,(total amount: nchn_r-1) move reverse  ## 255chn后的道数道翻转到0-(ndepE_r-1)道

rate = np.zeros((18,256,256))

def old_EC(serverlibDIR,obsmode=None,sttime=None):
    global Chan_l_A,Chan_h_A
    ecdatestr = np.loadtxt(serverlibDIR+'EC_FWHM/EC_'+obsmode+'_date_list.txt',dtype=str)
    ecdate = list(map(lambda x:Time(x,scale='utc',format='isot'),ecdatestr))
    idd = np.where(Time(sttime,scale='utc',format='isot') > ecdate)[0]

    nchn = 256
    
    if len(idd) != 0:
        matchdate = ecdatestr[idd[-1]]
        ecfile = '{:s}_EC_{:s}.txt'.format(matchdate[0:10],obsmode)
    else:
        matchdate = ecdatestr[0]
        ecfile = '{:s}_EC_{:s}.txt'.format(matchdate[0:10],obsmode)

    #default searching folder: serverlibDIR/EC_FWHM/
            
    if ('/' in ecfile) == False:
        ecfile = '{:s}EC_FWHM/{:s}'.format(serverlibDIR,ecfile)

    if os.path.isfile(ecfile) == False:
        print("cannot find EC file:{}! check available date: EC_FWHM/EC_{}_date_list.txt".format(ecfile,obsmode))
        sys.exit()
    
    ecdata = np.loadtxt(ecfile,usecols=[1,2])
    
    #ecfile = np.loadtxt('/hxmt/work/GRB/Software/RSPgenerator_v9/EC_FWHM/2018-03-10_EC_Nor.txt')

    chnl = np.arange(256)-0.5 #0~255个能道。第n道的下界为n-0.5，对应能量为ecfile[detid,0]*(n-0.5) + ecfile[detid,1]
    chnh = np.arange(256)+0.5  #第n道的上界为n-0.5

    print('deposite energy for each det (detnumber,Edmin,Edmax):')
    for detid in range(0,18): #detid是18个探测器 
        detid = int(detid)
        Chan_l_A[detid,:] = (chnl + nchn_r - 1) * ecdata[detid,0] + ecdata[detid,1]
        Chan_h_A[detid,:] = (chnh + nchn_r - 1) * ecdata[detid,0] + ecdata[detid,1]

        print('det{:d}, min:{:.3f} keV, max:{:.3f} keV'.format(int(detid),Chan_l_A[detid,0],Chan_h_A[detid,-1]))

        
def merge_content(n_chan,e_l,e_h,chan_decimals):
    global dets
    print('dets choosed:',dets)    
    
    #------------------rebin---------------------------
    if (e_l==0 or e_h==0):
        if obsmode == 'Nor':
            e_l = 30.
            e_h = 1100.
        else:
            e_l = 150.
            e_h = 3750.

    dindgen_a = np.arange(0,n_chan+1,1.0)
    chan_t = (dindgen_a/n_chan)*(e_h-e_l) + e_l
    chan_t = np.around(chan_t,decimals = int(chan_decimals))
    chan_high = chan_t[1:n_chan+1]
    chan_low = chan_t[0:n_chan]

    PI_low = chan_low
    PI_high = chan_high

    newbinlist = np.append(chan_low,chan_high[-1])
    #rsp_m_tem = np.zeros((18,128,n_chan))
    print('PI:')
    print(newbinlist)
    
    new_bin = chan_high - chan_low
    

    
    for i in dets:
        Chan_old = (Chan_l_A[i,:] + Chan_h_A[i,:])*0.5
        chan_new = (chan_low + chan_high)*0.5
        old_bin = Chan_h_A[i,:] - Chan_l_A[i,:]


        
        #------------------------handling channel by channel------------------------
        chan_tem = np.rint(np.linspace(0,n_chan-1,n_chan))
        rate_tem = np.zeros((256,256))
        for j in range(0,256):

            
            #----------------------compute every old chan distribute ratio in new chan--------------
            if (Chan_l_A[i,j]<=0 or Chan_h_A[i,j]<=0):
                rate_c = np.zeros(n_chan)
                rate_l = np.zeros(n_chan)
            else:
                #---------------------寻找投影能道起始\结束位置-----------------------
                c_start_dex1 = np.rint(Chan_l_A[i,j]>=chan_low)
                c_start_dex2 = np.zeros(n_chan)
                c_start_dex2[n_chan-1] = np.rint(Chan_l_A[i,j]>=chan_high[n_chan-1])
                c_start_dex2[0:n_chan-1] = c_start_dex1[1:n_chan]
                c_start_dex = c_start_dex1 - c_start_dex2
                c_start_dex = c_start_dex.astype(bool)
                c_start = np.rint(chan_tem[c_start_dex])
                c_start = c_start.astype(int)
               
                c_stop_dex1 = np.rint(chan_high>=Chan_h_A[i,j])
                c_stop_dex2 = np.zeros(n_chan)
                c_stop_dex2[0] = np.rint(chan_low[0]>=Chan_h_A[i,j])
                c_stop_dex2[1:n_chan] = c_stop_dex1[0:n_chan-1]
                c_stop_dex = c_stop_dex1 - c_stop_dex2
                c_stop_dex = c_stop_dex.astype(bool)
                c_stop = np.rint(chan_tem[c_stop_dex])#out_c
                c_stop = c_stop.astype(int)



                rate_c = np.zeros(n_chan)
                rate_l = np.zeros(n_chan)
                #----------------------能道log均匀比例计算-------------------------------
                if (c_start==(n_chan-1) and c_stop.shape==(0,)):         
                    rate_c[c_start] = (np.log(chan_high[c_start]) - np.log(Chan_l_A[i,j]))/(np.log(Chan_h_A[i,j])-np.log(Chan_l_A[i,j]))
                elif (c_stop==0 and c_start.shape==(0,)):
                    rate_c[c_stop] = (np.log(Chan_h_A[i,j])-np.log(chan_low[c_stop]))/(np.log(Chan_h_A[i,j])-np.log(Chan_l_A[i,j]))
                elif (c_start.shape==(1,) and c_stop.shape==(1,)):  
                    for k in range(c_start[0],c_stop[0]+1):
                        if c_stop==c_start:
                            rate_c[k] = 1
                        else:
                            #can be simple by use max{Chan_l_A[i,j],chan_low[k]},min{chan_high[k],Chan_h_A[i,j]},but not must
                            if k==c_start:
                                rate_c[k] = (np.log(chan_high[c_start]) - np.log(Chan_l_A[i,j]))/(np.log(Chan_h_A[i,j])-np.log(Chan_l_A[i,j]))
                            elif k==c_stop:
                                rate_c[k] = (np.log(Chan_h_A[i,j])-np.log(chan_low[c_stop]))/(np.log(Chan_h_A[i,j])-np.log(Chan_l_A[i,j]))
                            else:
                                rate_c[k] = (np.log(chan_high[k])-np.log(chan_low[k]))/(np.log(Chan_h_A[i,j])-np.log(Chan_l_A[i,j]))
                #----------------------linear均匀比例计算--------------------------------
                if (c_start==(n_chan-1) and c_stop.shape==(0,)):         
                    rate_l[c_start] = (np.array(chan_high[c_start]) - np.array(Chan_l_A[i,j]))/(np.array(Chan_h_A[i,j])-np.array(Chan_l_A[i,j]))
                elif (c_stop==0 and c_start.shape==(0,)):
                    rate_l[c_stop] = (np.array(Chan_h_A[i,j])-np.array(chan_low[c_stop]))/(np.array(Chan_h_A[i,j])-np.array(Chan_l_A[i,j]))
                elif (c_start.shape==(1,) and c_stop.shape==(1,)):
                    for k in range(c_start[0],c_stop[0]+1):
                        if c_stop==c_start:
                            rate_l[k] = 1
                        else:
                            #can be simple by use max{Chan_l_A[i,j],chan_low[k]},min{chan_high[k],Chan_h_A[i,j]},but not must
                            if k==c_start:
                                rate_l[k] = (np.array(chan_high[c_start]) - np.array(Chan_l_A[i,j]))/(np.array(Chan_h_A[i,j])-np.array(Chan_l_A[i,j]))
                            elif k==c_stop:
                                rate_l[k] = (np.array(Chan_h_A[i,j])-np.array(chan_low[c_stop]))/(np.array(Chan_h_A[i,j])-np.array(Chan_l_A[i,j]))
                            else:
                                rate_l[k] = (np.array(chan_high[k])-np.array(chan_low[k]))/(np.array(Chan_h_A[i,j])-np.array(Chan_l_A[i,j]))
                
            
            rate_tem[j,:] = rate_c
        rate[i,(nchn_r-1):256,:] = rate_tem[0:(256-nchn_r+1),:]
        rate[i,0:(nchn_r-1),:] = rate_tem[(256-nchn_r+1):256,:]   

        np.savetxt(workdir+'/'+'rate_{:02d}.txt'.format(int(i)),rate[i])#,fmt='%.4f')
    return rate,PI_low,PI_high


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--dets', action='store',
                        dest='dets',
                        help='input the detectors you want merge,e.g. 0 1 2 3,default will use all detectors!',
                        type=int,nargs='+',default=[0,1,2,3,4,5,6,7,8,9,10,11,12,\
                                           13,14,15,16,17])
    parser.add_argument('--n_chan', action='store',
                        dest='n_chan',
                        help='input the total number of PIs you want, default is 256',
                        type=int,default=256)
    
    parser.add_argument('--e_l', action='store',
                        dest='e_l',
                        help='input the lowest energy(keV) of channel left edge you want, if not given,/n\
                        Nor mode default is 30, GRB mode default is 400.',
                        type=float,default=0)
    parser.add_argument('--e_h', action='store',
                        dest='e_h',
                        help='input the highest energy(keV) of channel right edge you want, if not given,/n\
                        Nor mode default is 1100, GRB mode default is 2800.',
                        type=float,default=0)
    parser.add_argument('--chan_decimals', action='store',
                        dest='chan_decimals',
                        help='input a int, it will decide the decimals of channels, e.g. /n\
                        chan_decimals is 2 will make 0.455 be 0.46\
                        default is 3.',
                        type=int,default=3)
    parser.add_argument('--evttime', action='store',
                        dest='evttime',
                        help='input start time (format:YYYY-MM-DDThh:mm:ss)',
                        type=str,default=0)
    parser.add_argument('--obsmode', action='store',
                        dest='obsmode',
                        help='input Observation mode(Nor(default) or LHV)',
                        type=str,default='Nor')

    
    starttime = time.time()
    args = parser.parse_args()
    
    obsmode = args.obsmode
    dets = args.dets

    dets = [int(u) for u in dets]

    evttime = args.evttime
    n_chan = args.n_chan
    e_l = args.e_l
    e_h = args.e_h
    chan_decimals = args.chan_decimals

    #warnings.filterwarnings("ignore")
    
    old_EC(serverlibDIR,obsmode=obsmode,sttime=evttime)
    
    rate,PI_low,PI_high = merge_content(n_chan,e_l,e_h,chan_decimals)

    np.savetxt(workdir+'/'+'PI_low_'+obsmode+'.txt',PI_low)
    np.savetxt(workdir+'/'+'PI_high_'+obsmode+'.txt',PI_high)
    
    '''
    for d in range(18):
        np.savetxt('rate_{:02d}.txt',np.array(rate[d,:,:]))
    '''
    print('done!')
    
    endtime = time.time()
    totaltime = endtime - starttime
    
    print('time cost:{:.5f} seconds'.format(totaltime))
