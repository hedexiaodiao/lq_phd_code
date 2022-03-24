#!/hxmt/soft/Develop/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 14:31:32 2019

@author: Luoqi
"""
from xspec import *
import numpy as np
import astropy.time
from astropy.time import Time
import datetime
from datetime import date
import argparse
import os
#import sys
import time
import smtplib
from email.mime.text import MIMEText

mailto_list=["luoqi@ihep.ac.cn","xiongsl@ihep.ac.cn","xiaoshuo@ihep.ac.cn","lick@ihep.ac.cn","caice@ihep.ac.cn","yiqb@ihep.ac.cn"]
mail_host="smtp.139.com" 
mail_user="luoqi_ba"    
mail_pass="luoqi@172"   
mail_postfix="139.com"

alphas = [-1.9,-1.0,-0.0]
betas = [-3.7,-2.3,-1.5]
Etems = [700.,230.,500.]

def send_mail(to_list,sub,content):
    me="<"+mail_user+"@"+mail_postfix+">"#"GW_uplim_REPORT"+
    msg = MIMEText(content,_subtype='plain',_charset='gb2312')
    msg['Subject'] = sub
    msg['From'] = me
    msg['To'] = ";".join(to_list)
    try:
        #server = smtplib.SMTP()
        #server.connect(mail_host)
        #server.login(mail_user,mail_pass)
        #server.sendmail(me, to_list, msg.as_string())
        #server.close()
        return True
    except Exception as e:
        print(str(e))
        return False

def calc_flu(comporband,a,b,Epeak,t,A,Emin,Emax):
    Emin = int(round(Emin))
    Emax = int(round(Emax))
    
    if comporband=='comp':
        Epiv = 100#一般不用改
        Es = np.linspace(Emin,Emax,Emax*5)
    
        ss = []
        for i in range(0,Emax*5-1):
            s = ((Es[i]/Epiv)**a)*np.exp(-(a+2)*Es[i]/Epeak)*t*Es[i]*(Es[i+1]-Es[i])
            ss.append(s)
        f = sum(ss)
        flu = A*f/6.2415/10 #这里的数值来自单位换算，不用改
        
        
        print('time('+str(t)+')_alpha('+str(a)+')_flu:(',flu,')*10^-7 erg/cm^2')
    elif comporband=='band':
        Epiv = 100#一般不用改
        Es = np.linspace(Emin,Emax,Emax*5)

        Eedge = (a-b)*Epeak/(a+2)
        before_dex_tu = np.where(Es<Eedge)
        before_dex = before_dex_tu[0]
        dex_edge = len(before_dex)
        #print(dex_edge)#如果出问题可能是谱太硬、软，而入射Emin、Emax范围过窄，可以print出来检查一下
        
        ss = []
        for i in range(0,dex_edge):
            s = ((Es[i]/Epiv)**a)*np.exp(-(a+2)*Es[i]/Epeak)*t*Es[i]*(Es[i+1]-Es[i])
            ss.append(s)
        for i in range(dex_edge,Emax*5-1):
            s = t*Es[i]*(Es[i+1]-Es[i])*((Es[i]/Epiv)**b)*np.exp(b-a)*((a-b)*Epeak/(100*(a+2)))**(a-b)
            ss.append(s)
        f = sum(ss)
        flu = A*f/6.2415/10 #这里的数值来自单位换算，不用改
        
        
        print('time('+str(t)+')_alpha('+str(a)+')_flu:(',flu,')*10^-7 erg/cm^2')
    return flu


def fitxspec(obsmode,a,b,Etem,t,dirstr,trigtime):
    if os.path.exists('chain.fits'):
        os.system('rm chain.fits')
    AllData.clear()
    AllChains.clear()
    AllData('DA.src')
    s1 = AllData(1)
    #s1.background='DA.bkg'
    #s1.response='DA.rsp'
 
    if obsmode == 'Nor':
        obs = 0
        AllData.ignore("1: **-30,230-**")
    elif obsmode == 'GRB':
        obs = 1
        AllData.ignore("1: **-30,230-**")
    
    s = []
    s.append(AllData(1))#1是起始dex
    
    Model('grbm',setPars={1:a, 2:b, 3:Etem, 4:0.01})
    
    for i in range(1,4):
        par_s = AllModels(1)(i)
        par_s.frozen = 1
    #Fit.statMethod = "cstat"
    Fit.query = "yes"
    Fit.nIterations = 10000
    Fit.perform()
    
    AllChains.defBurn = 1000
    AllChains.defLength = 5000
    AllChains.defProposal = "gaussian fit"
    c = Chain("chain.fits")
    c.run()
    
    Fit.error('5. 4')
    
    par_s = AllModels(1)(4)
    para = par_s.values[0]
    downlim = par_s.error[0]#是参数4下限
    uplim = par_s.error[1]#是参数4上限
    flu = calc_flu('band',a,b,Etem*(a+2),t,uplim,200.,5000.)
    
    
    lstr = str(obs)+' '+str(t)+' '+str(a)+' '+str(b)+' '+str(Etem)+' '+str(downlim)+' '+str(para)+' '+str(uplim)+' '+str(flu)
    lstr += '\n'
    wrtfile = open('/sharefs/hbkg/user/luoqi/GRB/grbba/GW/'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr)+'/log'+str(int(t))+'.txt','a+')
    wrtfile.write(lstr)
    wrtfile.close()
    
    lstr1 = 'time: '+str(trigtime.value)+' timescal:'+str(t)+' a:'+str(a)+' b:'+str(b)+' Epeak:'+str(Etem*(a+2))+'\n'
    lstr1 += 'flulim: {:.4f}'.format(flu)+'*10^-7 erg/cm^2'+'\n'
    lstr1 += 'norm fitting value: {:.7f}'.format(para)+'\n'
    lstr1 += 'norm upper error: {:.7f}'.format(uplim)+'\n'
    wrtfile1 = open('/sharefs/hbkg/user/luoqi/GRB/grbba/GW/'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr)+'/results'+str(int(t))+'.txt','a+')
    wrtfile1.write(lstr1)
    wrtfile1.close()

def run_task_xspec(ra,dec,trigtime,obsmode,t,t_bkg_l_b,t_bkg_l_f,t_bkg_r_b,t_bkg_r_f,dirstr):
    os.chdir('/sharefs/hbkg/user/luoqi/GRB/grbba/GW/'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr)+'/t'+str(int(t))+'/FITS/')
    ##os.system('cp ../../HXMTrspg_v9_FITS/* ./')
    for i in range(0,3):
        a = alphas[i]
        b = betas[i]
        Etem = Etems[i]
        fitxspec(obsmode,a,b,Etem,t,dirstr,trigtime)

    #write email text    
    result_text = np.loadtxt('/sharefs/hbkg/user/luoqi/GRB/grbba/GW/'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr)+'/log'+str(int(t))+'.txt')
    len_result = result_text.shape[0]
    email_str = 'obsmode:'+str(obsmode)+',   ra:'+str(ra)+',   dec:'+str(dec)+',   trigtime:'+str(trigtime)+'\n'
    email_str += 'the used spec time is t_bkg_l_b='+str(t_bkg_l_b)+',t_bkg_l_f='+str(t_bkg_l_f)+',t_src_b=0'+\
        ',t_src_f='+str(t)+',t_bkg_r_b='+str(t_bkg_r_b)+',t_bkg_r_f='+str(t_bkg_r_f)+'\n'
    for j in range(0,len_result):
        l_tem = result_text[j]
        email_str += ' timescal:'+str(l_tem[1])+'   a:'+str(l_tem[2])+'   b:'+str(l_tem[3])+'   Epeak:'+str(l_tem[4]*(2+l_tem[2]))+'   flulim:'+str(l_tem[8])+'*10^-7 erg/cm^2;\n'

    '''
    #send email
    if send_mail(mailto_list,"GW uplim Result "+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr)+'_t'+str(int(t)),email_str):
        print "send email successful"
    else:
        print "send email to failed"
    '''
    
    return 0
    
    
def check_file_xspec(ra,dec,trigtime,obsmode,t,t_bkg_l_b,t_bkg_l_f,t_bkg_r_b,t_bkg_r_f,dirstr):
    print(str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr)+'/t'+str(int(t))+'/FITS/DA.src')
    if os.path.exists('/sharefs/hbkg/user/luoqi/GRB/grbba/GW/'+str(trigtime.datetime.year)+'{:02d}'.format(trigtime.datetime.month)+'{:02d}'.format(trigtime.datetime.day)+str(dirstr)+'/t'+str(int(t))+'/FITS/DA.src')==1:
        print('OK, is running')
        flag = run_task_xspec(ra,dec,trigtime,obsmode,t,t_bkg_l_b,t_bkg_l_f,t_bkg_r_b,t_bkg_r_f,dirstr)
    else:
        flag = 1
    return flag

def timerun_xspec(ra,dec,trigtime,obsmode,t,t_bkg_l_b,t_bkg_l_f,t_bkg_r_b,t_bkg_r_f,dirstr):
    j = 0
    flag = 1
    
    while flag==1:
        now = datetime.datetime.now()  
        print(now,j,'check src')   
        flag = check_file_xspec(ra,dec,trigtime,obsmode,t,t_bkg_l_b,t_bkg_l_f,t_bkg_r_b,t_bkg_r_f,dirstr)
        j = j+1
        time.sleep(1)


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
    print(trigtime)
    obsmode = args.obsmode
    t = args.t
    t_bkg_l_b = args.t_bkg_l_b
    t_bkg_l_f = args.t_bkg_l_f
    t_bkg_r_b = args.t_bkg_r_b
    t_bkg_r_f = args.t_bkg_r_f
    dirstr = args.dirstr
    trigtime=Time(trigtime,scale='utc',format='isot')

    timerun_xspec(ra,dec,trigtime,obsmode,t,t_bkg_l_b,t_bkg_l_f,t_bkg_r_b,t_bkg_r_f,dirstr)

    caledtime = time.time()
    print('time cost:{:.5f} seconds'.format(caledtime-calsttime))
    
    
