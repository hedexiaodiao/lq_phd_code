#!/hxmt/soft/Develop/anaconda3/bin/python
# -*- coding: utf-8 -*-

"""
Created on Mon Oct 28 17:23:16 2019

@author: Luoqi
"""
import numpy as np
import argparse
import os,sys
from scipy.interpolate import interp1d
import time
from functools import partial
from astropy.io import fits
import warnings
from astropy.time import Time
import datetime
from datetime import date
from datetime import datetime

def mkdir_try(dirname):
    flag = os.path.exists(dirname)
    if flag == 0:
        try:
            os.makedirs(dirname)
        except OSError:
            print('Wrong with make dir:\n'+dirname)
    return flag

def update_monthly(serverlibDIR,obsmode='Nor',wkdir=os.getcwd()):
    ecdatestr = np.loadtxt(serverlibDIR+'EC_FWHM/EC_'+obsmode+'_date_list.txt',dtype=str)
    ecdate = list(map(lambda x:Time(x,scale='utc',format='isot'),ecdatestr))
    print(ecdate)
    for i,u_ecdate in enumerate(ecdate):
        outdirdate =u_ecdate.datetime
        outdir = wkdir+'/'+obsmode+'/'+str(outdirdate.year)+'_'+str(outdirdate.month)+'_'+str(outdirdate.day)
        evttime = Time(u_ecdate.unix + 3600,scale='utc',format='unix')
        print(evttime.isot)
        mkdir_try(outdir)
        os.chdir(outdir)
        os.system('python {:s} --obsmode {:s} --evttime {:s}'.format(pyfile,obsmode,str(evttime.isot)))
        os.chdir(wkdir)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        wkdir = os.getcwd()
        print('work dir:',wkdir)
    elif len(sys.argv) == 2:
        wkdir = sys.argv[1]
    else:
        print('Wrong number of parameters!')
        os.exit(0)
    serverlibDIR = '/hxmt/work/GRB/Software/RSPgenerator_v9/'
    pyfile = '/sharefs/hbkg/user/luoqi/GRB/soft/new15_Chn/new15_Chn.py'
    for obsmode in ['Nor','LHV']:
        update_monthly(serverlibDIR,obsmode=obsmode,wkdir=wkdir)
