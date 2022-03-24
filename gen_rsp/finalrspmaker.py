# /hxmt/soft/Develop/anaconda2/bin/python
# -*- coding: utf-8 -*-
import os
import numpy as np
import copy
import sys
import logging

#用这个处理完后得到的是沉积数目，最终重复提交完作业后，先判断某次是否有问题，如没有，把(128，8000)累加起来
n_in = 128
n_chan = 8000
n_d = 18
num_S=[14, 15, 3, 4, 5, 6, 16, 17, 7, 8, 9,  10, 12, 13, 11, 0,  1,  2 ]
Ein_low = np.loadtxt('/sharefs/hbkg/user/luoqi/GRB/Ein_low.txt')
Ein_high = np.loadtxt('/sharefs/hbkg/user/luoqi/GRB/Ein_high.txt')
#E_depos_low = np.arange(0,4000,0.5)
#E_depos_high = np.arange(0.5,4000.5,0.5)

#文件数据结构，[:,0:20],[:,18]是入射能量，[:,19]是第多少次入射
#Ein_min,Ein_max是使用文件的入射能量范围，没必要用；Einnum是Ein_min的序号/16,范围从0到7(包括)
def exist_hist(theta,phi,Einnum,filenum):
    #1 need change; because the false of rename unzip,it need change when it's rest or first
    #path_file='../OutputData'+'_'+str(theta)+'_'+str(phi)+'/CsI_4_'+str(theta)+'_'+str(phi)+'_'+str(int(Ein_low[16*Einnum]))+'_'+str(int(Ein_high[16*(Einnum+1)-1]))+'_'+str(filenum)+'.dat'
    path_file='../fianl_OutputData/CsI_4_'+str(theta)+'_'+str(phi)+'_'+str(int(Ein_low[16*Einnum]))+'_'+str(int(Ein_high[16*(Einnum+1)-1]))+'_'+str(filenum)+'.dat'
    data0 = np.loadtxt(path_file)
    

    arf = np.ones((n_in,n_d))*0
    rmf = np.ones((n_in,n_chan,n_d))*0
    
    for ein_num in range(0,16):
        CsI= []
        CsI = copy.deepcopy(data0[:,0:19]) #先按照入射能量分
        
        CsI_Ein = CsI[:,18]
        CsI_Ein = np.array(CsI_Ein)
        condi1 = CsI_Ein>=Ein_low[16*Einnum+ein_num]
        condi2 = CsI_Ein<Ein_high[16*Einnum+ein_num]
        condi3 = condi1 & condi2
        dex_tem_tu = np.where(condi3)
        
        dex_tem = dex_tem_tu[0]
        
        CsI_ein_tem = [CsI[dex,:] for dex in dex_tem]
        if CsI_ein_tem ==[]:
            CsI_ein_tem = np.zeros((1,18))
        CsI_ein_tem = np.array(CsI_ein_tem)
        #CsI_ein_tem把指定入射范围的沉积数组获取了
        
        for det_num in range(0,18):
            CsI_tem = CsI_ein_tem[:,num_S[det_num]]
            CsI_tem_tem = [l for l in CsI_tem if l >0]
            CsI_tem = np.array(CsI_tem_tem)
            #这是一维数组，信息仅是沉积能量，把=0的去除了
            
            y_tem,binofy_tem = np.histogram(CsI_tem,bins=8000,range=(0,4000))
            rmf[16*Einnum+ein_num,:,det_num] = y_tem[0:8000]
            arf[16*Einnum+ein_num,det_num] = np.sum(y_tem)
            
    for det_num in range(0,18):
        np.savetxt('depostimes_D'+str(det_num)+'_T'+str(theta)+'_P'+str(phi)+'_'+'Ein_'+str(Einnum)+'file_'+str(filenum)+'.txt', rmf[:,:,det_num])
        np.savetxt('arftem_D'+str(det_num)+'_T'+str(theta)+'_P'+str(phi)+'_'+'Ein_'+str(Einnum)+'file_'+str(filenum)+'.txt',arf[:,det_num])
            

def file_ifexist(theta,phi,Einnum,filenum,noex_num):
    #2 need change; because the false of rename unzip,it need change when it's rest or first
    #path_file='../OutputData'+'_'+str(theta)+'_'+str(phi)+'/CsI_4_'+str(theta)+'_'+str(phi)+'_'+str(int(Ein_low[16*Einnum]))+'_'+str(int(Ein_high[16*(Einnum+1)-1]))+'_'+str(filenum)+'.dat'
    path_file='../fianl_OutputData/CsI_4_'+str(theta)+'_'+str(phi)+'_'+str(int(Ein_low[16*Einnum]))+'_'+str(int(Ein_high[16*(Einnum+1)-1]))+'_'+str(filenum)+'.dat'
    if os.path.exists(path_file):
        logger.info('CsI_4_'+str(theta)+'_'+str(phi)+'_'+str(int(Ein_low[16*Einnum]))+'_'+str(int(Ein_high[16*(Einnum+1)-1]))+'_'+str(filenum)+'.dat exist but may problem')
        exist_hist(theta,phi,Einnum,filenum)
    else:
        logger.info('CsI_4_'+str(theta)+'_'+str(phi)+'_'+str(int(Ein_low[16*Einnum]))+'_'+str(int(Ein_high[16*(Einnum+1)-1]))+'_'+str(filenum)+'.dat not exist')
        shfile = open('./noexist'+str(noex_num)+'/submitgeant4work_'+str(theta)+'_'+str(phi)+'_'+str(Ein_low[16*Einnum])+'_'+str(Ein_high[16*(Einnum+1)-1])+'_'+str(filenum)+'.sh','w')
        shfile.close
        M =str(theta)+' '+str(phi)+' '+str(Einnum)+' '+str(filenum)+' '+str(noex_num)+'\n'
        nomfile = open('./logfile/arraylog6.sh','a+')
        nomfile.write(M)
        nomfile.close()
    #判断文件存在部分完成

if __name__ == "__main__":
    theta,phi,Einnum,filenum = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
    theta = int(theta)
    phi = int(phi)
    Einnum = int(Einnum)
    filenum = int(filenum)
    if os.path.exists('./compressed_matrix')==0:
        os.mkdir('./compressed_matrix')
    os.chdir('./compressed_matrix')
    if os.path.exists('./logfile')==0:
        os.mkdir('./logfile')
    noex_num = 6
    #while (os.path.exists('./logfile/times'+str(noex_num)+'.sh')==1):
    #noex_num = noex_num + 1
    if os.path.exists('./noexist'+str(noex_num))==0:
        os.mkdir('./noexist'+str(noex_num))

    if os.path.exists('./logfile/arraylog6.sh')==0:
        nomfile = open('./logfile/arraylog6.sh','w')
        nomfile.close
    
    # 获取logger对象,取名mylog
    logger = logging.getLogger('mylog')

    # 输出INFO及以上级别的信息
    logger.setLevel(level=logging.INFO)

    # 获取文件处理器并设置日志级别
    handler = logging.FileHandler('./logfile/'+str(theta)+'_'+str(phi)+'_'+str(Ein_low[16*Einnum])+'_'+str(Ein_high[16*(Einnum+1)-1])+'_'+str(filenum)+'log.txt')
    handler.setLevel(logging.INFO)

    # 生成并设置文件处理器的日志格式
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # 获取流处理器并设置日志级别
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    # 为logger对象同时添加文件处理器和流处理器
    logger.addHandler(handler)
    logger.addHandler(console)

    # 打印日志
    logger.info("Start print log")
    logger.debug("Do something")
    logger.warning("Something maybe fail.")
    
    file_ifexist(theta,phi,Einnum,filenum,noex_num)
    
    logger.info("Finish")


