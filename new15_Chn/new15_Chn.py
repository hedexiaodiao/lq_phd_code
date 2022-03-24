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
        print('det{:d}, 10 Chan (equal 30 Chan before turn):{:.3f} keV, {:.3f} keV'.format(int(detid),Chan_l_A[detid,9],Chan_h_A[detid,9]))


###-----------8 PI setting-------for convenient------
def new_Chn(e_l,e_h):
    n_chan = 8
    # ------------------rebin---------------------------
    if (e_l == 0 or e_h == 0):
        if obsmode == 'Nor':
            e_l = 30.
            e_h = 1100.
        else:
            e_l = 150.
            e_h = 3750.

    dindgen_a = np.arange(0,n_chan+1,1.0)
    chan_t = (dindgen_a/n_chan)*(np.log(e_h)-np.log(e_l)) + np.log(e_l)
    chan_t = np.around(np.exp(chan_t),decimals = int(chan_decimals))
    chan_high = chan_t[1:n_chan+1]
    chan_low = chan_t[0:n_chan]

    PI_low = chan_low
    PI_high = chan_high

    newbinlist = np.append(chan_low, chan_high[-1])
    # rsp_m_tem = np.zeros((18,128,n_chan))
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

        n_chan_old = 256

        a_tem = np.arange(0,256)

        print('the {:d}th Detector'.format(int(i)))

        for j in range(0,n_chan):

            #----------------------compute every old chan distribute ratio in new chan--------------

            #---------------------寻找投影能道起始\结束位置-----------------------
            c_start_dex1 = np.rint(chan_low[j]>=Chan_l_A[i])
            c_start_dex2 = np.zeros(n_chan_old)
            c_start_dex2[n_chan_old-1] = np.rint(chan_high[j]>=Chan_l_A[i,n_chan_old-1])
            c_start_dex2[0:n_chan_old-1] = c_start_dex1[1:n_chan_old]
            c_start_dex = c_start_dex1 - c_start_dex2
            c_start_dex = c_start_dex.astype(bool)
            #c_start = np.rint(chan_tem[c_start_dex])
            #c_start = c_start.astype(int)

            c_stop_dex1 = np.rint(Chan_h_A[i]>=chan_high[j])
            c_stop_dex2 = np.zeros(n_chan_old)
            c_stop_dex2[0] = np.rint(Chan_h_A[i,0]>=chan_low[j])
            c_stop_dex2[1:n_chan_old] = c_stop_dex1[0:n_chan_old-1]
            c_stop_dex = c_stop_dex1 - c_stop_dex2
            c_stop_dex = c_stop_dex.astype(bool)
            #c_stop = np.rint(chan_tem[c_stop_dex])#out_c
            #c_stop = c_stop.astype(int)

            print('the start/stop address of new Channel {:d} in the old electronic Channel:\n'.format(int(j)),a_tem[c_start_dex],a_tem[c_stop_dex])
        print('\n\n')

def get_ECdate():
    if ecfile == None:
        ecdatestr = np.loadtxt(serverlibDIR+'EC_FWHM/EC_'+obsmode+'_date_list.txt',dtype=str)
        ecdate = list(map(lambda x:Time(x,scale='utc',format='isot'),ecdatestr))
        idd = np.where(Time(sttime,scale='utc',format='isot') > ecdate)[0]
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

def merge_content(n_chan,e_l,e_h,chan_decimals):
    global dets
    print('dets choosed:',dets)


    #------------------------256pi------------------------
    if obsmode == 'Nor':
        oldpi_e_l = 30.
        oldpi_e_h = 1100.
    else:
        oldpi_e_l = 150.
        oldpi_e_h = 3750.
    
    dindgen_a = np.arange(0, 256 + 1, 1.0)
    oldpi_chan_t = (dindgen_a / 256) * (oldpi_e_h - oldpi_e_l) + oldpi_e_l
    oldpi_chan_t = np.around(oldpi_chan_t, decimals=int(chan_decimals))
    oldpi_chan_high = oldpi_chan_t[1:256 + 1]
    oldpi_chan_low = oldpi_chan_t[0:256]


    #------------------15pi; rebin---------------------------
    rest_dex = 15
    if (e_l==0 or e_h==0):
        if obsmode == 'Nor':
            e_l = oldpi_chan_low[rest_dex]###前面rest_dex个能道(0、1、2、……rest_dex-1)去除，从第rest_dex的左边起始
            e_h = 1100.
            rest_n_chan = n_chan - 1###后14个能道按照这样分配，0能道按照原PI的0～14合并
        else:
            rest_dex = 0
            e_l = oldpi_chan_low[rest_dex]
            e_h = 3750.
            rest_n_chan = n_chan###15个能道直接分

    #rest_n_chan = n_chan - 1


    if logorlinear =='linear':
        dindgen_a = np.arange(0,rest_n_chan+1,1.0)
        chan_t = (dindgen_a/rest_n_chan)*(e_h-e_l) + e_l
        chan_t = np.around(chan_t,decimals = int(chan_decimals))
        rest_chan_high = chan_t[1:n_chan+1]
        rest_chan_low = chan_t[0:n_chan]
    else:
        dindgen_a = np.arange(0, rest_n_chan + 1, 1.0)
        chan_t = (dindgen_a / rest_n_chan) * (np.log(e_h) - np.log(e_l)) + np.log(e_l)
        chan_t = np.around(np.exp(chan_t), decimals=int(chan_decimals))
        rest_chan_high = chan_t[1:rest_n_chan + 1]
        rest_chan_low = chan_t[0:rest_n_chan]

    if obsmode=='Nor':
        chan_low = np.append(oldpi_chan_low[0],rest_chan_low)
        chan_high = np.append(oldpi_chan_low[rest_dex],rest_chan_high)
    else:
        chan_low = rest_chan_low
        chan_high = rest_chan_high
    
    newpidex = np.array([])
    for i in range(n_chan):
        dex = np.argmin(np.abs(oldpi_chan_low-chan_low[i]))
        print('dex:',dex)
        chan_low[i] = oldpi_chan_low[dex]
        newpidex = np.append(newpidex,dex)

    np.savetxt(workdir+'/'+'{:s}_newpidex.txt'.format(str(obsmode)),newpidex)

    
    '''
    if os.path.exists(workdir+'/'+'LHV_newpidex.txt'):
        ano_dex = np.loadtxt(workdir+'/'+'LHV_newpidex.txt')
        ano_dex = np.rint(ano_dex)
        for i in range(n_chan):
            chan_low[i] = oldpi_chan_low[int(ano_dex[i])]
        np.savetxt(workdir + '/' + '{:s}_newpidex.txt'.format(str(obsmode)), ano_dex)
    else:
        print('please execute LHV mode first, then Nor mode!')
    '''
    
    #chan_high[0:n_chan-1] = chan_low[1:n_chan]

    PI_low = chan_low
    PI_high = chan_high

    newbinlist = np.append(chan_low,chan_high[-1])
    #rsp_m_tem = np.zeros((18,128,n_chan))
    #print('PI_chan_low:')
    #print(chan_low)
    #print('PI_chan_high:')
    #print(chan_high)
    print('PI:')
    print(newbinlist)

    new_bin = chan_high - chan_low



    for i in dets:
        Chan_old = (Chan_l_A[i,:] + Chan_h_A[i,:])*0.5
        chan_new = (chan_low + chan_high)*0.5
        old_bin = Chan_h_A[i,:] - Chan_l_A[i,:]



        #------------------------handling channel by channel------------------------
        chan_tem = np.rint(np.linspace(0,n_chan-1,n_chan))
        rate_tem = np.zeros((256,n_chan))
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
        np.savetxt(workdir+'/'+'{:s}_rate_{:02d}.txt'.format(str(obsmode),int(i)),rate[i])#,fmt='%.4f')
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
                        help='input the total number of PIs you want, default is 15',
                        type=int,default=15)

    parser.add_argument('--e_l', action='store',
                        dest='e_l',
                        help='input the lowest energy(keV) of channel left edge you want, if not given,/n\
                        Nor mode default is 30, GRB mode default is 150.',
                        type=float,default=0)
    parser.add_argument('--e_h', action='store',
                        dest='e_h',
                        help='input the highest energy(keV) of channel right edge you want, if not given,/n\
                        Nor mode default is 1100, GRB mode default is 3750.',
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
                        type=str,default='2019-01-02T13:31:00')
    parser.add_argument('--obsmode', action='store',
                        dest='obsmode',
                        help='input Observation mode(Nor(default) or LHV)',
                        type=str,default='Nor')
    parser.add_argument('--logorlinear', action='store',
                        dest='logorlinear',
                        help='input the bin type of PI, log(default) or linear',
                        type=str, default='log')


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
    logorlinear = args.logorlinear
    #warnings.filterwarnings("ignore")

    
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

    rate = np.zeros((18,256,n_chan))
    

    old_EC(serverlibDIR,obsmode=obsmode,sttime=evttime)

    #new_Chn(e_l,e_h)
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
