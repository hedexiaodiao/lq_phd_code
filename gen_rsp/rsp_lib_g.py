# -*- coding: utf-8 -*-
import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
from astropy.time import Time
import argparse
from scipy.interpolate import interp1d
import sys
import subprocess
from matplotlib import ticker, cm, colors
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import shutil



angall = np.loadtxt('/hxmt/work/GRB/Software/RSPgenerator_v9/SimInfo/Incident_Angle.txt',usecols=[0,1],dtype=int)

print(angall)
print(angall.shape)

#testidd = np.random.randint(len(angall),size=1)
#print(angall[testidd])
#angall = [[0,0],[60,0],[60,90],[90,0],[90,90],[150,0],[150,90]]
dirhead = '/sharefs/hbkg/user/luoqi/GRB/soft/rspgenerator/rsplib'

path01 = dirhead+'/RSP_CsI/Nor'
path02 = dirhead+'/RSP_CsI/LG'
#plt.ion()
#for i in range(len(angall)):
def mkrsplib(theta,phi):
    if os.path.exists(dirhead+'/tem/{:d}_{:d}'.format(int(theta),int(phi)))==0:
        os.mkdir(dirhead+'/tem/{:d}_{:d}'.format(int(theta),int(phi)))
    if os.path.exists(dirhead+'/tem/{:d}_{:d}/Nor'.format(int(theta),int(phi)))==0:
        os.mkdir(dirhead+'/tem/{:d}_{:d}/Nor'.format(int(theta),int(phi)))
    if os.path.exists(dirhead+'/tem/{:d}_{:d}/LG'.format(int(theta),int(phi)))==0:
        os.mkdir(dirhead+'/tem/{:d}_{:d}/LG'.format(int(theta),int(phi)))
        
    path1 = dirhead+'/tem/{:d}_{:d}/Nor/HXMTrspg_v1_FITS'.format(int(theta),int(phi))
    path2 = dirhead+'/tem/{:d}_{:d}/LG/HXMTrspg_v1_FITS'.format(int(theta),int(phi))

    
    os.chdir(dirhead+'/tem/{:d}_{:d}/Nor'.format(int(theta),int(phi)))
    commandstr = '{:s} --obsmode Nor \
                --theta-deg {:.2f}  --phi-deg {:.2f} \
                '\
                .format('/sharefs/hbkg/user/luoqi/GRB/soft/rspgenerator/RSPgenerator_v0/HXMTrspg_v1.py', float(theta), float(phi))
    print(commandstr)
    subprocess.call(commandstr,shell=True)
    for det in range(0,18):
        shutil.move(path1+'/D{:02d}.rsp'.format(int(det)),path01+'/D{:d}/D{:d}_T{:d}_P{:d}_Nor.fits'.format(int(det),int(det),int(theta),int(phi)))#os.path.join(
        
    os.chdir(dirhead+'/tem/{:d}_{:d}/LG'.format(int(theta),int(phi)))
    commandstr = '{:s} --obsmode LG \
                        --theta-deg {:.2f}  --phi-deg {:.2f} \
                        '\
                .format('/sharefs/hbkg/user/luoqi/GRB/soft/rspgenerator/RSPgenerator_v0/HXMTrspg_v1.py', float(theta), float(phi))
    print(commandstr)
    subprocess.call(commandstr,shell=True)
    for det in range(0,18):
        shutil.move(path2+'/D{:02d}.rsp'.format(int(det)),path02+'/D{:d}/D{:d}_T{:d}_P{:d}_LG.fits'.format(int(det),int(det),int(theta),int(phi)))
    print('{:d}_{:d} done!'.format(int(theta),int(phi)))

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--startnum', action='store',
                        dest='startnum',
                        help='input startnum',
                        type=int,default='0')
    parser.add_argument('--endnum', action='store',
                        dest='endnum',
                        help='input endnum',
                        type=int,default='1')
    args = parser.parse_args()
    startnum = args.startnum
    endnum = args.endnum
    
    for i in range(startnum,endnum):
        mkrsplib(angall[i,0],angall[i,1])
