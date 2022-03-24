import numpy as np
from astropy.io import fits
import datetime
import os
import sys
#import global_demo as gl
def saveFitsfile(
                fitsname, detnam, n_chan, 
                energy_ch_l, energy_ch_h,
                energy_in_l, energy_in_h,
                matrix_out, 
                detector='CsI',
                mass_model='v2018', 
                obsmode='Nor', 
                theta=23., phi=12.,
                ra=None, dec=None,
                sttime='2018-09-01T58:22:00', edtime='2018-09-01T58:22:00',
                ecfile = None,
                fwhmfile = None,
                infile = [None],
                w = None,
                thispro=None):
    
    now = datetime.datetime.now()
    nowdate = now.date().strftime("%Y-%m-%d")

    hdu0 = fits.PrimaryHDU()
    hdu0.header['LONGSTRN'] = ('OGIP 1.0','The HEASARC Long String Convention may be used')
    publicheader0(hdu0,nowdate=nowdate,mass_model=mass_model,obs_mode=obsmode,
        detector=detector,detnam=detnam,ecfile=ecfile,fwhmfile=fwhmfile,
        theta=theta, phi=phi,sttime=sttime,edtime=edtime,
        ra=ra,dec=dec,thispro=thispro)

    col1 = fits.Column(name='CHANNEL', format='I', array=np.arange(n_chan)+1)
    col2 = fits.Column(name='E_MIN', format='1E', unit='keV', array=energy_ch_l)
    col3 = fits.Column(name='E_MAX', format='1E', unit='keV', array=energy_ch_h)
    hdu1 = fits.BinTableHDU.from_columns([col1, col2, col3])
    hdu1.header['EXTNAME'] = ('EBOUNDS','Extension name')
    hdu1.header['HDUCLAS2'] = ('EBOUNDS','extension contains a response matrix')
    publicheader0(hdu1,nowdate=nowdate,mass_model=mass_model,obs_mode=obsmode,
        detector=detector,detnam=detnam,ecfile=ecfile,fwhmfile=fwhmfile,
        theta=theta, phi=phi,sttime=sttime,edtime=edtime,
        ra=ra,dec=dec,thispro=thispro)
    publicheader1(hdu1,n_chan=n_chan,infile=infile,weight=w)
    

    col1 = fits.Column(name='ENERG_LO',format='1E',unit='keV',array=energy_in_l)
    col2 = fits.Column(name='ENERG_HI',format='1E',unit='keV',array=energy_in_h)
    col3 = fits.Column(name='N_GRP',format='I',array=np.ones(len(energy_in_l),dtype='int'))
    # col4 = fits.Column(name='F_CHAN',format='PI()',array=np.array(np.ones([len(energy_in_l),1],dtype=int),dtype=np.object))
    # col5 = fits.Column(name='N_CHAN',format='PI()',array=np.array(n_chan*np.ones([len(energy_in_l),1],dtype=int),dtype=np.object))
    # col6 = fits.Column(name='MATRIX',format='PE()',unit='cm2',array=np.array(matrix_out,dtype=np.object))

    col4 = fits.Column(name='F_CHAN',format='1E',array=np.array(np.ones(len(energy_in_l)),dtype=int))
    col5 = fits.Column(name='N_CHAN',format='1E',array=np.array(n_chan*np.ones(len(energy_in_l)),dtype=int))
    col6 = fits.Column(name='MATRIX',format=str(len(matrix_out[0]))+'E',array=matrix_out)

    hdu2 = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
    
    hdu2.header['EXTNAME'] = ('SPECRESP MATRIX','Extension name')
    hdu2.header['HDUCLAS2'] = ('RSP_MATRIX','extension contains a response matrix')
    publicheader0(hdu2,nowdate=nowdate,mass_model=mass_model,obs_mode=obsmode,
        detector=detector,detnam=detnam,ecfile=ecfile,fwhmfile=fwhmfile,
        theta=theta, phi=phi,sttime=sttime,edtime=edtime,
        ra=ra,dec=dec,thispro=thispro)
    publicheader1(hdu2,n_chan=n_chan,infile=infile,weight=w)

    temp = []
    temp.append(hdu0)
    temp.append(hdu1)
    temp.append(hdu2)
    hdulists = fits.HDUList(temp)

    if( os.path.isfile(fitsname)):
        os.remove(fitsname)

    hdulists.writeto(fitsname)

def publicheader0(h,nowdate,mass_model,obs_mode,
        detector,detnam,ecfile,fwhmfile,
        theta, phi,sttime,edtime,ra,dec,thispro):

    #creator = gl.get_value('thispro')

    h.header['DATE'] = (nowdate,'Creation date')
    h.header['TELESCOP'] = ('HXMT','Telescope used')
    h.header['INSTRUME'] = ('HE','Instrument used')
    h.header['DETECTOR'] = (detector,'Detector used')
    h.header['FILTER'] = ('NONE','the instrument filter in use (if any)')
    h.header['MASSMODL'] = (mass_model,'mass model used')
    h.header['DETNAM'] = (detnam, 'Detector used')
    h.header['GRATING'] = ('NONE', 'grating used, if any')
    h.header['ORIGIN'] = ('HXMTrspgenerator','source of FITS file')
    h.header['CREATOR'] = (thispro,'program creating this file')
    h.header['REVISION'] = (1,'Revision')
    h.header['MISSION'] = (obs_mode+' Mode','Mission')
    h.header['ECfile'] = (ecfile,'EC relation used')
    h.header['FWHMfile'] = (fwhmfile,'resolution used')
    h.header['STARTTIM'] = (sttime,'start time')
    h.header['ENDTIME'] = (edtime,'end time')
    h.header['RA'] = (ra,'RA(deg)')
    h.header['DEC'] = (dec,'DEC(deg)')
    h.header['THETA'] = (theta,'THETA(deg)')
    h.header['PHI'] = (phi,'PHI(deg)')

def publicheader1(h,n_chan,infile,weight):
    if infile[0] != None:
        for k in range(len(infile)):
           h.header['INFILE'+str(k)] = infile[k]
        h.header['WEIGHT'] = ('{:.3f} {:.3f} {:.3f} {:.3f}'.format(*weight),'infile weights')
    h.header['CHANTYPE'] = ('PI      ','(PHA or PI) uncorrected or corrected for gain')
    h.header['DETCHANS'] = (n_chan,'raw detector channels')
    h.header['HDUCLASS'] = ('OGIP    ','file format is OGIP standard')
    h.header['HDUCLAS1'] = ('RESPONSE','extension contains response data')
    h.header['HDUVERS'] = ('1.3.0   ','version of the file format')
    h.header['TLMIN4'] = (1,'first channel in the response')
    h.header['RMFVERSN'] = ('1992a   ','obsolete')
    h.header['HDUVERS1'] = ('1.1.0   ','obsolete')
    h.header['HDUVERS2'] = ('1.2.0   ','obsolete')
    h.header['CVSD0001'] = ('2012-01-01','UTC date of beginning applicability')
    h.header['CVST0001'] = ('00:00:00','UTC time of beginning applicability')
    h.header['CDES0001'] = ('RSP     ','brief descriptive summary of this dataset')
