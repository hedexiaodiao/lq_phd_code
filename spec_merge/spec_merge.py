#!/hxmt/soft/Develop/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:35:52 2019

@author: Luoqi
"""


import numpy as np
import multiprocessing as mp
import argparse
import os
from scipy.interpolate import interp1d
import time
from functools import partial
from astropy.io import fits
import writersp_fits as wf
import warnings

workdir = os.getcwd()
Ein_l_A = np.zeros((18,128))
Ein_h_A = np.zeros((18,128))
Chan_l_A = np.zeros((18,256))
Chan_h_A = np.zeros((18,256))
rsp_A = np.zeros((18,128,256))
src_A = np.zeros((18,256))
src_err_A = np.zeros((18,256))
bkg_A = np.zeros((18,256))
bkg_err_A = np.zeros((18,256))
SNR_A = np.zeros(18)
dets = [0,1,2,3,4,5,6,7,8,9,10\
        ,11,12,13,14,15,16,17]
obsmode = 'Nor'
weight_m = np.ones(18)
theta = 45
phi = 90
ra = 45
dec = 90


def read_data():
    global Ein_l_A,Ein_h_A,Chan_l_A,Chan_h_A,theta,phi,ra,dec
    global rsp_A,src_A,src_err_A,bkg_A,bkg_err_A,obsmode
    
    for i in range(0,18):
        if (os.path.isfile(workdir+'/D{:02d}.rsp'.format(i)) and\
        os.path.isfile(workdir+'/D{:02d}.bkg'.format(i)) and\
        os.path.isfile(workdir+'/D{:02d}.src'.format(i)))==0:
            print('D{:02d} src and/or bkg and/or rsp not exists!!!'.format(i))
            #print('if D{:02d} was choosed, merged files will be wrong!!!')
            print('try to use D{:1d} src and/or bkg and/or rsp, if also not exists, may be wrong!!!'.format(i))
            #print('loading file D{:02d} ,please wait!'.format(i))
        try:
            if os.path.isfile(workdir+'/D{:02d}.rsp'.format(i)):
                str_tem = 'D{:02d}'.format(i)
            else:
                str_tem = 'D{:1d}'.format(i)
            with fits.open(workdir+'/'+str_tem+'.rsp') as hdl:
                Ein_l_A[i,:] = hdl[2].data['ENERG_LO']
                Ein_h_A[i,:] = hdl[2].data['ENERG_HI']
                Chan_l_A[i,:] = hdl[1].data['E_MIN']
                Chan_h_A[i,:] = hdl[1].data['E_MAX']
                rsp_A[i,:,:] = hdl[2].data['MATRIX']
                obsmode = hdl[0].header['MISSION'][0:3]
                theta = hdl[0].header['THETA']
                phi = hdl[0].header['PHI']
                ra = hdl[0].header['RA']
                dec = hdl[0].header['DEC']
            if os.path.isfile(workdir+'/D{:02d}.src'.format(i)):
                str_tem = 'D{:02d}'.format(i)
            else:
                str_tem = 'D{:1d}'.format(i)
            with fits.open(workdir+'/'+str_tem+'.src') as hdl:
                src_A[i,:] = hdl[1].data['COUNTS']
                src_err_A[i,:] = hdl[1].data['STAT_ERR']

            if os.path.isfile(workdir+'/D{:02d}.bkg'.format(i)):
                str_tem = 'D{:02d}'.format(i)
            else:
                str_tem = 'D{:1d}'.format(i)
            with fits.open(workdir+'/'+str_tem+'.bkg') as hdl:
                bkg_A[i,:] = hdl[1].data['COUNTS']
                bkg_err_A[i,:] = hdl[1].data['STAT_ERR']
        except:
            pass

def decide_dets(SNR_set,SNRrank_detnum):
    global SNR_A,dets,weight_m
    nd_tem = np.rint(np.linspace(0,17,18))
    nd_tem = nd_tem.astype(int)
    SNR_A =(np.sum(src_A,axis=1))/np.sqrt(np.sum(bkg_A,axis=1))
    l2h_index = np.argsort(SNR_A)
    if SNR_set !=0:
        dets = [u for u in nd_tem[np.greater(SNR_A,SNR_set)]]
    if SNRrank_detnum !=0:
        dets = [u for u in nd_tem[l2h_index[-SNRrank_detnum:]]]
    
    
    if snr_weight:
        w_a_tem = 0
        for u in dets:
            w_a_tem += SNR_A[u]
        weight_m = SNR_A/w_a_tem
        weight_m = weight_m*len(dets)
        a = 0
        for u in dets:
            a += weight_m[u]
        print('Total weight is {:.3f}'.format(a))

'''
def myinterp(oldchan_mid,newchan_mid,matrx_in):
    matrx_out = np.zeros((128,len(newchan_mid)))
    for ss in range(0,128):
        funcinterp = interp1d(oldchan_mid,matrx_in[ss,:],bounds_error=False,fill_value=(0,0))
        matrx_out[ss,:] = funcinterp(newchan_mid)
    return matrx_out
'''
    
def rebinChnMatrx(oldbinlist, newbinlist,matrx_in):
    intrgMatrx_old = np.zeros(([len(matrx_in),len(matrx_in[0])+1]))
    intrgMatrx_old[:,1:] = np.cumsum(matrx_in,axis=1)
    func = interp1d(oldbinlist,intrgMatrx_old,bounds_error=False,fill_value=(intrgMatrx_old[:,0],intrgMatrx_old[:,-1]))
    intrgMatrx_new = func(newbinlist)
    #this method, intrgMatrx's shape is (128,256)
    return np.array(intrgMatrx_new[:,1:]-intrgMatrx_new[:,0:-1])

def rebinrsp(filedir,newbinlist,detid):
    energy_ch_det_l = Chan_l_A[detid,:]
    energy_ch_det_h = Chan_h_A[detid,:]
    detbinlist = np.append(energy_ch_det_l,energy_ch_det_h[-1])
    rsp = rebinChnMatrx(detbinlist, newbinlist,rsp_A[detid,:,:]*weight_m[detid])
    return rsp

        
def merge_content(n_chan,e_l,e_h,chan_iflog,chan_decimals):
    global dets
    print('dets choosed:',dets)

    src_m = np.zeros(n_chan)
    src_err_m = np.zeros(n_chan)
    src_err2_m = np.zeros(n_chan)
    bkg_m = np.zeros(n_chan)
    bkg_err_m = np.zeros(n_chan)
    bkg_err2_m = np.zeros(n_chan)
    rsp_m = np.zeros((128,n_chan))

    
    
    #------------------rebin---------------------------
    if (e_l==0 or e_h==0):
        if obsmode == 'Nor':
            e_l = 100.
            e_h = 600.
        else:
            e_l = 400.
            e_h = 2800.
    if chan_iflog==1:
        print('chan_iflog_log:',chan_iflog)
        dindgen_a = np.arange(0,n_chan+1,1.0)
        chan_t = np.exp((dindgen_a/n_chan)*(np.log(e_h)-np.log(e_l))+np.log(e_l))
        chan_t = np.around(chan_t,decimals = int(chan_decimals))
        chan_high = chan_t[1:n_chan+1]
        chan_low = chan_t[0:n_chan]
    else:
        print('chan_iflog_linear:',chan_iflog)
        dindgen_a = np.arange(0,n_chan+1,1.0)
        chan_t = (dindgen_a/n_chan)*(e_h-e_l) + e_l
        chan_t = np.around(chan_t,decimals = int(chan_decimals))
        chan_high = chan_t[1:n_chan+1]
        chan_low = chan_t[0:n_chan]

    newbinlist = np.append(chan_low,chan_high[-1])
    #rsp_m_tem = np.zeros((18,128,n_chan))


    #---------------------compute t of spec,ignoring dettime and bkg flux factor-------------
    t_src_A = np.zeros((18,256))
    t_bkg_A = np.zeros((18,256))
    src_tdex = src_err_A!=0
    t_src_A[src_tdex] = src_A[src_tdex]/(src_err_A[src_tdex]**2)

    bkg_tdex = bkg_err_A!=0
    t_bkg_A[bkg_tdex] = bkg_A[bkg_tdex]/(bkg_err_A[bkg_tdex]**2)


    #--------------------declare a average t, a kind of approximate--------------------------
    t_src_tot = np.zeros(18)
    t_bkg_tot = np.zeros(18)

    
    new_bin = chan_high - chan_low
    

    
    for i in dets:
        src_tem = np.zeros(n_chan)
        src_err2_tem = np.zeros(n_chan)
        bkg_tem = np.zeros(n_chan)
        bkg_err2_tem = np.zeros(n_chan)
    
        t_src_tot[i] = np.sum(t_src_A[i])/((np.where(src_tdex[i]==1)[0]).shape[0])
        t_bkg_tot[i] = np.sum(t_bkg_A[i])/((np.where(bkg_tdex[i]==1)[0]).shape[0])
        #checked

        #------------------------==cnts--------------------------
        src_err_A[i,:] = (src_err_A[i,:]**2)*t_src_tot[i]
        bkg_err_A[i,:] = (bkg_err_A[i,:]**2)*t_bkg_A[i]
        
        '''
        src_A[i,:] = src_A[i,:]*t_src_tot[i]
        src_err_A[i,:] = src_err_A[i,:]*t_src_tot[i]
        bkg_A[i,:] = bkg_A[i,:]*t_bkg_A[i]
        bkg_err_A[i,:] = bkg_err_A[i,:]*t_bkg_A[i]
        '''
        

        print('merging file D{:02d} '.format(i))

        Chan_old = (Chan_l_A[i,:] + Chan_h_A[i,:])*0.5
        chan_new = (chan_low + chan_high)*0.5
        old_bin = Chan_h_A[i,:] - Chan_l_A[i,:]



        #------------------------turn to np.log,avoid down zero--------------------------------
        upzerodex_src = (np.greater(src_A[i,:],np.zeros(len(src_A[i,:])))*np.greater(src_err_A[i,:],np.zeros(len(src_A[i,:]))))
        src_A[i,upzerodex_src] = np.log(src_A[i,upzerodex_src])
        src_err_A[i,upzerodex_src] = np.log(src_err_A[i,upzerodex_src])
        
        upzerodex_bkg = (np.greater(bkg_A[i,:],np.zeros(len(bkg_A[i,:])))*np.greater(bkg_err_A[i,:],np.zeros(len(bkg_A[i,:]))))
        bkg_A[i,upzerodex_bkg] = np.log(bkg_A[i,upzerodex_bkg])
        bkg_err_A[i,upzerodex_bkg] = np.log(bkg_err_A[i,upzerodex_bkg])
        

        #------------------------times oldbin--------------------------------------
        src_A[i,:] = src_A[i,:] * old_bin
        src_err_A[i,:] = src_err_A[i,:] * old_bin
        bkg_A[i,:] = bkg_A[i,:] * old_bin
        bkg_err_A[i,:] = bkg_err_A[i,:] * old_bin        

        
        #------------------------handling channel by channel------------------------
        chan_tem = np.rint(np.linspace(0,n_chan-1,n_chan))
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
                            #即:可以,但是没必要
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
                            #即:可以,但是没必要
                            if k==c_start:
                                rate_l[k] = (np.array(chan_high[c_start]) - np.array(Chan_l_A[i,j]))/(np.array(Chan_h_A[i,j])-np.array(Chan_l_A[i,j]))
                            elif k==c_stop:
                                rate_l[k] = (np.array(Chan_h_A[i,j])-np.array(chan_low[c_stop]))/(np.array(Chan_h_A[i,j])-np.array(Chan_l_A[i,j]))
                            else:
                                rate_l[k] = (np.array(chan_high[k])-np.array(chan_low[k]))/(np.array(Chan_h_A[i,j])-np.array(Chan_l_A[i,j]))
                
                                
            #---------------makes one det dif chans' re-distribute in new chans-------------------------
            src_tem = src_tem + src_A[i,j]*rate_l*weight_m[i]
            src_err2_tem = src_err2_tem + (src_err_A[i,j])*(rate_l*weight_m[i])

            bkg_tem = bkg_tem + bkg_A[i,j]*rate_l*weight_m[i]
            bkg_err2_tem = bkg_err2_tem + (bkg_err_A[i,j])*(rate_l*weight_m[i])

            '''
            for k in range(0,128):
                rsp_m[k,:] = rsp_m[k,:] + rsp_A[i,k,j]*rate_l*weight_m[i]
            '''
            
        #------------------because initial np.log(src_A) is linear to E, integral can be displaced by delt(E)*np.log(src_A)-------------
        '''
        src_tem = src_tem / new_bin
        src_err2_tem = src_err2_tem / new_bin
        bkg_tem = bkg_tem / new_bin
        bkg_err2_tem = bkg_err2_tem / new_bin
        
        
        src_tem = np.exp(src_tem)
        src_err2_tem = np.exp(src_err2_tem)
        bkg_tem = np.exp(bkg_tem)
        bkg_err2_tem = np.exp(bkg_err2_tem)
        '''
        #----------------------will be zeros after i changed--------------------
        src_tem = np.exp(src_tem / new_bin)
        src_err2_tem = np.exp(src_err2_tem / new_bin)

        bkg_tem = np.exp(bkg_tem / new_bin)
        bkg_err2_tem = np.exp(bkg_err2_tem / new_bin)
        

        #---------------------assignment to more glo var-------------------------
        src_m = src_m + src_tem
        src_err2_m = src_err2_m + src_err2_tem
        bkg_m = bkg_m + bkg_tem
        bkg_err2_m = bkg_err2_m + bkg_err2_tem
        

        #===============================one det's re-distributed end==============================
            #if np.sum(rate_c)>1.02 or np.sum(rate_c)<0.98:
                #print(i,j,np.sum(rate_c))
        #rsp_m = rsp_m + myinterp(Chan_old,chan_new,rsp_A[i,:,:])
        #rsp_m_tem[i,:,:] = rebinrsp(workdir,newbinlist,i)
    #rsp_m = np.sum(rsp_m_tem,axis=0)#this method may obviously underestimated rsp_m when n_chan is under 256

    #the dif between dettime factors of dif dets hasn't been considered,average was

    #bkg_revise_dex = np.greater(bkg_m,src_m)
    #bkg_m[bkg_revise_dex] = 0.9*src_m[bkg_revise_dex]
    #print(bkg_revise_dex)
    
    #------------------== cnts/(t**2)-------------
    src_err2_m = src_err2_m/(np.sum(t_src_tot)/18)
    bkg_err2_m = bkg_err2_m/(np.sum(t_bkg_tot)/18)

    #------------------== cnts_err/t--------------
    src_err_m = src_err2_m**0.5
    bkg_err_m = bkg_err2_m**0.5


    rebinrsp1para = partial(rebinrsp,workdir,newbinlist)
    with mp.Pool(2) as p:
        rsp_tem = p.map(rebinrsp1para,dets)
    rsp_m = np.array(rsp_tem)
    rsp_m = np.sum(rsp_m,axis=0)

    '''
    src_m = src_m/(np.sum(t_src_tot)/18)
    src_err_m = src_err_m/(np.sum(t_src_tot)/18)
    bkg_m = bkg_m/(np.sum(t_bkg_tot)/18)
    bkg_err_m = bkg_err_m/(np.sum(t_bkg_tot)/18)
    '''



    return src_m,src_err_m,bkg_m,bkg_err_m,rsp_m,chan_low,chan_high

def write_fits(src_m,src_err_m,bkg_m,bkg_err_m,rsp_m,chan_low,chan_high):
    for i in dets:
        if (os.path.isfile(workdir+'/D{:02d}.rsp'.format(i)) and\
            os.path.isfile(workdir+'/D{:02d}.bkg'.format(i)) and\
            os.path.isfile(workdir+'/D{:02d}.src'.format(i)))==1:
            break
    chan_num = np.rint(np.linspace(1,n_chan,n_chan))
    chan_num = chan_num.astype(int)
    srchdl = fits.open(workdir+'/D{:02d}.src'.format(i))
    srchdl[0].header['FILENAME'] = str(outputname)+'.src'
    srchdl[1].header['BACKFILE'] = str(outputname)+'.bkg'
    srchdl[1].header['RESPFILE'] = str(outputname)+'.rsp'
    srchdl[1].header['DETCHANS'] = int(n_chan)
    srchdl[1].header['TLMAX1'] = int(n_chan)
    srchdl[1].header['NAXIS2'] = int(n_chan)
    srchdl[1].data['CHANNEL'] = chan_num
    srchdl[1].data['COUNTS'] = src_m
    srchdl[1].data['STAT_ERR'] = src_err_m
    if os.path.isfile(workdir+'/'+str(outputname)+'.src'):
        os.remove(workdir+'/'+str(outputname)+'.src')
    srchdl.writeto(workdir+'/'+str(outputname)+'.src')
    
    bkghdl = fits.open(workdir+'/D{:02d}.bkg'.format(i))
    bkghdl[0].header['FILENAME'] = str(outputname)+'.bkg'
    bkghdl[1].header['DETCHANS'] = int(n_chan)
    bkghdl[1].header['TLMAX1'] = int(n_chan)
    bkghdl[1].header['NAXIS2'] = int(n_chan)
    bkghdl[1].data['CHANNEL'] = chan_num
    bkghdl[1].data['COUNTS'] = bkg_m
    bkghdl[1].data['STAT_ERR'] = bkg_err_m
    if os.path.isfile(workdir+'/'+str(outputname)+'.bkg'):
        os.remove(workdir+'/'+str(outputname)+'.bkg')
    bkghdl.writeto(workdir+'/'+str(outputname)+'.bkg')
    
    wf.saveFitsfile(outputname+'.rsp',outputname,n_chan,chan_low,chan_high,Ein_l_A[i,:],Ein_h_A[i,:],rsp_m,theta=theta,phi=phi,ra=ra,dec=dec)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--snrdownlim', action='store',
                        dest='SNR_set',
                        help='input the lowest SNR(total_cnts/sqrt(back_cnts)) you want,/n\
                        if given,arg \'dets\' will not be used!',
                        type=float,default=0)
    
    parser.add_argument('--snrrank_detnum', action='store',
                        dest='SNRrank_detnum',
                        help='input the number of highest SNR_ranked dets you want,/n\
                        if given,arg \'snrdownlim\'\
                        and arg \'dets\' will not be used!',
                        type=int,default=0)

    parser.add_argument('--dets', action='store',
                        dest='dets',
                        help='input the detectors you want merge,e.g. 0 1 2 3,default will use all detectors!',
                        type=int,nargs='+',default=[0,1,2,3,4,5,6,7,8,9,10,11,12,\
                                           13,14,15,16,17])

    parser.add_argument('--outputname', action='store',
                        dest='outputname',
                        help='input the output name,extension name not need,default is DA',
                        type=str,default='DA')
    
    parser.add_argument('--is_snrweight', action='store',
                        dest='snr_weight',
                        help='input 0 or 1, if 1, will merge files weighted by signal to noise ratio!/n\
                        if 0, will merge files with equal weight!/n\
                        default is 0(equal weight).',
                        type=int,default=0)
   
    parser.add_argument('--n_chan', action='store',
                        dest='n_chan',
                        help='input the number of channels you want, default is 256',
                        type=int,default=256)
    
    parser.add_argument('--e_l', action='store',
                        dest='e_l',
                        help='input the lowest energy(keV) of channel left edge you want, if not given,/n\
                        Nor mode default is 100, GRB mode default is 400.',
                        type=float,default=0)
    parser.add_argument('--e_h', action='store',
                        dest='e_h',
                        help='input the highest energy(keV) of channel right edge you want, if not given,/n\
                        Nor mode default is 600, GRB mode default is 2800.',
                        type=float,default=0)
    parser.add_argument('--chan_iflog', action='store',
                        dest='chan_iflog',
                        help='input 0 or 1, if 1, output channels will use log division,/n\
                        if 0, output channels will use linear division,/n\
                        default is 0(linear).',
                        type=bool,default=0)
    parser.add_argument('--chan_decimals', action='store',
                        dest='chan_decimals',
                        help='input a int, it will decide the decimals of channels, e.g. /n\
                        chan_decimals is 2 will make 0.455 be 0.46\
                        default is 2.',
                        type=int,default=2)

    
    starttime = time.time()
    args = parser.parse_args()
    
    SNR_set = args.SNR_set
    SNRrank_detnum = args.SNRrank_detnum
    dets = args.dets
    outputname = args.outputname
    snr_weight = args.snr_weight

    dets = [int(u) for u in dets]
    
    n_chan = args.n_chan
    e_l = args.e_l
    e_h = args.e_h
    chan_iflog = args.chan_iflog
    chan_decimals = args.chan_decimals

    warnings.filterwarnings("ignore")
    print('if log in:',chan_iflog)
    read_data()
    
    decide_dets(SNR_set,SNRrank_detnum)
    
    src_m,src_err_m,bkg_m,bkg_err_m,rsp_m,chan_low,chan_high = merge_content(n_chan,e_l,e_h,chan_iflog,chan_decimals)
    
    write_fits(src_m,src_err_m,bkg_m,bkg_err_m,rsp_m,chan_low,chan_high)
    
    print('done!')
    
    endtime = time.time()
    totaltime = endtime - starttime
    
    print('time cost:{:.5f} seconds'.format(totaltime))
  
