#!/usr/bin/env python
# coding: utf-8

__all__ = ['piszkespipe_from_commandline']

import argparse
import matplotlib
matplotlib.use("Agg")
from pylab import *

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ceres modules
import piszkespipe.piszkesutils as piszkesutils
from piszkespipe.Correlation import correlation
from piszkespipe.GLOBALutils import GLOBALutils
import piszkespipe.Marsh as Marsh

# other useful modules
import ephem
from math import radians as rad
from astropy.io import fits as pyfits
import pickle
import os
import scipy
import ccdproc
import logging

import statsmodels.api as sm
lowess = sm.nonparametric.lowess

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u
from astropy.convolution import convolve, Box1DKernel

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

np.seterr(all='ignore')

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

####### Debug options #####
debug = False
debugfindorders = False
debugfindordersiter = False
debugfitTHARlinesiter = False
debugfitTHARlinesresult = False
debugglobalvelshift = False
debugFLATextraction = False
debugfitflat = False
debugsciencecorrection = False
debugbkgestimation = False
debugscienceorders = False
debugsciencespectrumfit = False
debugfitblaze = False
debugfitblazeiter = False

####### GLOBAL VARIABLES #####
force_pre_process = False
force_flat_extract = False
force_thar_extract = False
force_thar_wavcal = False
force_sci_extract = True
force_stellar_pars = False
force_ofind = True # must be True to always locate orders
#force_bkg = False
#force_tharxc = False
#compute_rv = True
dark_corr = True
back_corr = True

def piszkespipe(dirin,avoid_plot,dirout,DoClass,JustExtract,npools,object2do,
                reffile,onlysimpleextract,remove_cosmic,keepoutliers):

    ########## I/O Path settings ##########
    dirin = os.path.join(dirin,'')
    if not os.path.exists(dirin):
        print(dirin + ' does not exists!')
        exit()

    if dirout == 'default':
        dirout = dirin[:-1]+'_red/'

    dirout = os.path.join(dirout,'')

    if not os.access(dirout,os.F_OK):
        os.system('mkdir '+dirout)
    if os.access(os.path.join(dirout,'proc'),os.F_OK):
        os.system('cp -rp '+ os.path.join(dirout,'proc','')+' '+os.path.join(dirout,'proc_old','') )
        os.system('rm -r '+ os.path.join(dirout,'proc') )
    os.system('mkdir '+ os.path.join(dirout,'proc') )

    f_res = open( os.path.join(dirout,'proc','results.csv') ,'w')
    f_res.write('obname,filename,bjd,RV,RVerr,BiSector,BSerr,RVbary,RVcorr,RVcorrerr,RVerr_systematic,')
    f_res.write('RVtwoG1,RVtwoG1corr,RVtwoG1err,RVtwoG2,RVtwoG2corr,RVtwoG2err,')
    f_res.write('instrument,pipeline,R,Teff,Tefferr,logg,loggerr,Z,Zerr,vsini,XC_min,dispersion,exptime,SNR_5130,ccf_pdf\n')

    if reffile == 'default':
        reffile = os.path.join(dirin,'reffile.txt')


    ####### GLOBAL VARIABLES #####
    Inverse_m = True
    use_cheby = True

    MRMS = 800 # max rms in m/s, global wav solution

    trace_degree  = 4
    Marsh_alg     = 0
    ext_aperture  = 3 # aperture to extract spectra
    NSigma_Marsh  = 5
    NCosmic_Marsh = 5
    S_Marsh       = 0.4
    N_Marsh       = 3 # grado polinomio
    min_extract_col = 20   # cut out edges of frames
    max_extract_col = 2028

    oro0 = 29 # physical m order of 0th aperture

    # coefficients to fit dispersion axis and wavelengths at once
    ncoef_x            = 4
    ncoef_m            = 6
    npar_wsol = (min(ncoef_x,ncoef_m) + 1) * (2*max(ncoef_x,ncoef_m) - min(ncoef_x,ncoef_m) + 2) / 2
    npar_wsol = int(npar_wsol)

    ###### Path to the synthetic models ######
    if DoClass and os.getenv('COELHO') is None:
        print(bcolors.FAIL + 'COELHO model dir is not found, but -do_class is specified!' + bcolors.ENDC)
        print('Please download the models (2.85 GB), e.g. by')
        print('\tmkdir ~/COELHO_MODELS')
        print('\tcd ~/COELHO_MODELS')
        print('\twget https://cloud.konkoly.hu/s/taT2qiSjpBwWWDr/download/coelho_05_red4_R40.tar.gz')
        print('\ttar -xf coelho_05_red4_R40.tar.gz')
        print('And add the location of R_40000b dir to $COELHO path, e.g. by')
        print("\texport COELHO=~/COELHO_MODELS")
        print('Optionally add this line to your .bash_rc file')
        exit()
    elif DoClass:
        models_path = os.path.join( os.getenv('COELHO') , 'R_40000b' )

        # test whether the fits are there
        testfile = os.path.join(models_path,'vsini_0.0','R_0.0_5000_30_p00p00.ms.fits')
        if os.access(testfile,os.F_OK) == False:
            print(bcolors.FAIL + 'COELHO model dir exists, but fits files are missing!' + bcolors.ENDC)
            print('Have you extracted the tar file?')
            exit()

    order_dir = os.path.join(PACKAGEDIR,"wavcals",'')
    n_useful  = 26 # up to which order do we care?

    def gettime(hd):
        date = hd['DATE-OBS']
        ye = float(date[:4])
        mo = float(date[5:7])
        da = float(date[8:10])
        ho = float(date[11:13])
        mi = float(date[14:16])
        se = float(date[17:])
        jd = piszkesutils.jd(ye,mo,da,ho,mi,se)

        return jd - 4./24.

    def get_extension(h):
        if '.fz' in h:
            extension = 1
        else:
            extension = 0
        return extension

    # file containing the log
    logging.basicConfig(filename=os.path.join(dirout,'night.log'), level=logging.INFO)
    log = logging.getLogger(__name__)

    print("\tRAW data is in ",dirin)
    print("\tProducts of reduction will be in",dirout)
    print('\n')

    ###### Download JPL ephemeris if needed ######
    ephempath = os.path.join(PACKAGEDIR,'de421.bsp')
    piszkesutils.check_ephem_file(ephempath)

    ###### Classify all the images according to its imagetype ######
    biases, flats, img_flats, fib_flats, objects, ThAr_ref, darks = piszkesutils.FileClassify(dirin,log,debug=debug)


    ################ Collecting initial informations on calibration data ##############################
    if len(biases) == 0:
        print(bcolors.FAIL + 'Error: No Bias frame found!' + bcolors.ENDC)
        log.warning('Error: No Bias frame found!')
        exit()

    if len(ThAr_ref) == 0:
        print(bcolors.FAIL + 'Error: No ThAr frame found!' + bcolors.ENDC)
        log.warning('Error: No ThAR frame found!')
        exit()

    if len(objects) == 0:
        print(bcolors.FAIL + 'Error: No OBJECT frame found!' + bcolors.ENDC)
        log.warning('Error: No OBJECT frame found!')
        exit()

    if len(flats) > 0:
        have_flat = True
        #JustExtract = False
    else:
        print(bcolors.FAIL + 'Warning: No FLAT frame found!' + bcolors.ENDC)
        print(bcolors.FAIL + 'Simple blaze correction will be performed!' + bcolors.ENDC)
        print(bcolors.FAIL + 'This will affect strong/wide emission lines!' + bcolors.ENDC)
        log.warning('Warning: No FLAT frame found!')
        have_flat = False
        #JustExtract = True

    if len(darks)>0:
        have_darks = True
    else:
        print(bcolors.WARNING + 'Warning: No DARK frame found!' + bcolors.ENDC)
        log.warning('Warning: No DARK frame found!')
        have_darks = False

    if len(fib_flats) > 0:
        have_fib_flats = True
    else:
        have_fib_flats = False

    if ( (os.access(os.path.join(dirout,'trace.pkl'),os.F_OK) == False) or (os.access(os.path.join(dirout,'MasterBias.fits'),os.F_OK) == False) or (force_pre_process) ):
        print("\tNo previous pre-processing files are found")
        pre_process = 1
    else:
        print("\tPre-processing files found, going straight to extraction")
        pre_process = 0

    RON,GAIN = piszkesutils.get_rg(objects[0],biases[0],debug=debug)

    ################ Creating master calibration frames & tracing echelle orders ##############################
    if pre_process == 1 or not os.access( os.path.join(dirout,'MasterBias.fits'), os.F_OK ):
        print("\tGenerating Master calibration frames...")
        # median combine Biases
        MasterBias, roB, gaB = piszkesutils.MedianCombine(biases, False, os.path.join(dirout,'MasterBias.fits'),RON,GAIN)
        hdu = pyfits.PrimaryHDU( MasterBias )
        if (os.access( os.path.join(dirout,'MasterBias.fits') ,os.F_OK)):
            os.remove( os.path.join(dirout,'MasterBias.fits') )
        hdu.writeto( os.path.join(dirout,'MasterBias.fits') )
        print("\t\t-> Masterbias: done!")

        # median combine Darks
        if len(darks)!=0:
            dark_times = []
            each_time = []

            for dark in darks:
                dh = pyfits.getheader(dark,ext=get_extension(dark))
                dtime = dh['EXPTIME']
                each_time.append(dtime)
                esta = False
                for t in dark_times:
                    if dtime == t:
                        esta = True
                        break
                    else:
                        esta = False
                if esta == False:
                    dark_times.append(dtime)

            MasDarl = []

            for t in dark_times:
                sirven = []
                i = 0
                while i < len(darks):
                    if each_time[i] == t:
                        sirven.append(darks[i])
                    i+=1
                Mdark, roD, gaF = piszkesutils.MedianCombine(sirven, True, os.path.join(dirout,'MasterBias.fits'),RON,GAIN)
                nMdark = 'MasterDark_'+str(t)+'s.fits'
                if (os.access( os.path.join(dirout,nMdark) ,os.F_OK)):
                    os.remove( os.path.join(dirout,nMdark) )
                hdu = pyfits.PrimaryHDU( Mdark )
                hdu = GLOBALutils.update_header(hdu,'EXPTIME',t)
                hdu.writeto( os.path.join(dirout,nMdark) )
                MasDarl.append( os.path.join(dirout,nMdark) )

            darks_dict = {'d_names': MasDarl, 'd_times':dark_times}
            pickle.dump( darks_dict, open(  os.path.join(dirout,"darks.pkl"), 'wb' ) )
            print("\t\t-> Masterdarks: done!")

        else:
            dark_times = []
            MasDarl = []
            darks_dict = {'d_names': MasDarl, 'd_times':dark_times}
            pickle.dump( darks_dict, open(  os.path.join(dirout,"darks.pkl"), 'wb' ) )
            print("\t\t-> 0 Masterdarks found!")

        # median combine list of flats
        if len(flats)>0:
            Flat, roF, gaF =piszkesutils.MedianCombine(flats, True, os.path.join(dirout,'MasterBias.fits'),RON,GAIN, dark_bo=have_darks, dlist=MasDarl)
            #if Flat.ravel().min() < 0:
            #    Flat -= Flat.ravel().min() + 1
            hdu = pyfits.PrimaryHDU( Flat )
            if (os.access( os.path.join(dirout,'MasterFlat.fits') ,os.F_OK)):
                os.remove( os.path.join(dirout,'MasterFlat.fits') )
            hdu.writeto( os.path.join(dirout,'MasterFlat.fits') )

            print("\t\t-> Masterflats: done!")
        else:
            print(bcolors.WARNING + 'Warning! No FLAT frames.' + bcolors.ENDC)
            print(bcolors.WARNING + 'Using median of science images to trace orders!' + bcolors.ENDC)
            log.warning('No FLAT frames. Using median of science images to trace orders!')
            Flat, roF, gaF =piszkesutils.MedianCombine(objects, True, os.path.join(dirout,'MasterBias.fits'),RON,GAIN, dark_bo=have_darks, dlist=MasDarl)
            #if Flat.ravel().min() < 0:
            #    Flat -= Flat.ravel().min() + 1
            hdu = pyfits.PrimaryHDU( Flat )
            if (os.access( os.path.join(dirout,'MasterFlat.fits') ,os.F_OK)):
                os.remove( os.path.join(dirout,'MasterFlat.fits') )
            hdu.writeto( os.path.join(dirout,'MasterFlat.fits') )

        print("\tTracing echelle orders...")
        h = pyfits.open( os.path.join(dirout,'MasterFlat.fits') )[0]
        d = h.data
        #d = d.T
        c_all, nord = GLOBALutils.get_them(d,ext_aperture+1,trace_degree,maxords=45,startfrom=20,debug=debugfindorders,debugiter=debugfindordersiter)
        print('\t\t'+str(nord)+' orders found ...')
        trace_dict = {'c_all':c_all, 'nord':nord,'roF':roF,'gaF':gaF}
        pickle.dump( trace_dict, open(  os.path.join(dirout,"trace.pkl"), 'wb' ) )

    else:
        print("\tLoading Master calibration frames...")
        h = pyfits.open( os.path.join(dirout,'MasterBias.fits') )
        MasterBias = h[0].data
        # load orders & tracers
        trace_dict = pickle.load( open(  os.path.join(dirout,"trace.pkl"), 'rb' ) )
        #print trace_dict['c_all'].shape
        c_all = trace_dict['c_all']
        nord = trace_dict['nord']
        roF  = trace_dict['roF']
        gaF  = trace_dict['gaF']
        darks_dict = pickle.load( open(  os.path.join(dirout,"darks.pkl"), 'rb' ) )
        dark_times = darks_dict['d_times']
        MasDarl = darks_dict['d_names']
        h = pyfits.open( os.path.join(dirout,'MasterFlat.fits'))[0]
        d = h.data
        #d = d.T
        Flat = d.copy()


    ################################ ThAr spectra extraction  ##################################################
    thtimes   = []
    nThAr_ref = []
    nthtimes  = []
    print('\n\tExtraction of ThAr calibration frames:')
    #force_thar_extract = True
    for thar in ThAr_ref:
        #print(thar)
        h = pyfits.open(thar)
        mjd,mjd0 = piszkesutils.mjd_fromheader(h[get_extension(thar)])
        hdth     = pyfits.getheader(thar, ext=get_extension(thar))
        thtime   = mjd
        thtimes.append( thtime )

        dthar = pyfits.getdata( thar, ext=get_extension(thar) ) - MasterBias
        #dthar = dthar.T

        #force_thar_extract=False
        thar_fits_simple = os.path.join(dirout,'ThAr_'+os.path.splitext(os.path.basename(thar))[0]+'.spec.simple.fits.S')

        if ( os.access(thar_fits_simple,os.F_OK) == False ) or (force_thar_extract):

            print("\t\tNo previous extraction or extraction forced for ThAr file", os.path.basename(thar), ", extracting...")
            thar_Ss = GLOBALutils.simple_extraction(dthar,c_all,ext_aperture,min_extract_col,max_extract_col,npools)
            thar_Ss = thar_Ss[::-1]

            # save as fits file
            if (os.access(thar_fits_simple,os.F_OK)):
                os.remove( thar_fits_simple )

            hdu = pyfits.PrimaryHDU( thar_Ss )
            hdu.writeto( thar_fits_simple )

        else:
            print("\t\tThAr file", os.path.basename(thar), "already extracted, loading...")
            thar_Ss = pyfits.getdata(thar_fits_simple)

    ################################ ThAr spectra calibration ##################################################
    numt = 0
    print("\n\tWavelength solution of ThAr calibration spectra:")
    for thar in ThAr_ref:
        h = pyfits.open(thar,ext=get_extension(thar))
        mjd,mjd0 = piszkesutils.mjd_fromheader(h[get_extension(thar)])
        hdth     = pyfits.getheader(thar,ext=get_extension(thar))
        thtime   = mjd
        thar_fits_simple =  os.path.join(dirout,'ThAr_'+os.path.splitext(os.path.basename(thar))[0]+'.spec.simple.fits.S')
        thar_fits_simple_wav =  os.path.join(dirout,'ThAr_'+os.path.splitext(os.path.basename(thar))[0]+'.spec.wav.fits.S')
        wavsol_pkl =  os.path.join(dirout,'ThAr_'+os.path.splitext(os.path.basename(thar))[0]+'.wavsolpars.pkl')

        # If forced do ThAr wavelength calibration even if it is already exists
        #force_thar_wavcal = True
        if ( os.access(wavsol_pkl,os.F_OK) == False ) or (force_thar_wavcal):
            thar_Ss = pyfits.getdata(thar_fits_simple)
            print("\t\tWorking on ThAr file", os.path.basename(thar) )

            # Find which aperture is the 14th to make sure we match the input ThAr line list
            lines_thar = thar_Ss.copy()
            if numt == 0:
                if os.access( os.path.join(dirout,'order_find.pkl'),os.F_OK)==False or force_ofind:
                    maxes = 0
                    or08 = 0
                    for order in range(len(lines_thar)):
                        ccf_max, deltar = GLOBALutils.cor_thar(lines_thar[order],span=50,filename=os.path.join(order_dir,'order08.dat'),debug=debugfitTHARlinesiter)
                        if debugfitTHARlinesiter: print('order=',order, 'ccf_max',ccf_max)
                        if ccf_max > maxes:
                            maxes = ccf_max
                            or08  =  order
                            delta = deltar

                    or0 = or08 - 8
                    or25 = or08 + 17
                    #print(or0,or08,or25,nord)
                    if or0 >= 0:
                        orwa = 0
                    else:
                        orwa = - or0
                        or0  = 0

                    if or25 <= nord - 1:
                        orwb = 25
                    else:
                        orwb = 25 - (or25 - nord - 1)
                        or25 = nord - 1
                    #print(or0,or08,or25,nord)
                    if debugfitTHARlinesiter: print("or0,or25, orwa, orwb, delta:",or0,or25, orwa, orwb,delta)

                    pdict = {'orwa':orwa, 'or0':or0, 'orwb':orwb, 'or25':or25, 'delta':delta}
                    pickle.dump( pdict, open( os.path.join(dirout,'order_find.pkl'), 'wb' ) )
                else:
                    pdict = pickle.load(open(os.path.join(dirout,'order_find.pkl'),'rb'))
                    orwa = pdict['orwa']
                    or0 = pdict['or0']
                    orwb = pdict['orwb']
                    or25 = pdict['or25']
                    delta = pdict['delta']

            iv_thar = 1/((lines_thar/GAIN) + (RON**2/GAIN**2))

            All_Pixel_Centers = np.array([])
            All_Wavelengths = np.array([])
            All_Orders = np.array([])
            All_Centroids = np.array([])
            All_Sigmas = np.array([])
            All_Intensities = np.array([])
            All_Residuals = np.array([])
            All_Sigmas = np.array([])

            orre = or0
            order = orwa

            # Fit ThAR lines one by one in each order
            mid_col_wv = []
            OK = []
            OW = []
            nup = or25 - or0 + 1
            trans = np.zeros([nup,4])
            while orre <= or25:
                order_s = str(order)
                if (order < 10):
                    order_s = '0'+str(order)

                thar_order_orig = lines_thar[orre,:]
                IV = iv_thar[orre,:]
                wei = np.sqrt( IV )
                bkg = GLOBALutils.Lines_mBack(thar_order_orig, IV, thres_rel=3000000, line_w=15)
                thar_order = thar_order_orig - bkg
                if debugfitTHARlinesiter:
                    plt.figure(figsize=(20,10))
                    plt.title('Bkg corrected ThAR spectrum')
                    plt.plot(thar_order_orig,label='ThAr order=%d' % orre)
                    plt.plot(thar_order,label='ThAr - bkg')
                    #plt.ylim(-1000,10000)
                    plt.grid()
                    plt.legend()
                    plt.show()
                coeffs_pix2wav, coeffs_pix2sigma, pixel_centers, wavelengths,rms_ms, residuals, centroids, sigmas, intensities =\
                GLOBALutils.Initial_Wav_Calibration(os.path.join(order_dir,'order'+order_s+'.dat'),thar_order,order,wei,rmsmax=1000,
                                                    minlines=6,FixEnds=True,Dump_Argon=False,Dump_AllLines=True, Cheby=use_cheby,
                                                    rough_shift=delta,debug=debugfitTHARlinesiter)

                mid_col_wv.append(GLOBALutils.Cheby_eval(coeffs_pix2wav, int(.5*len(thar_order)), len(thar_order)))

                if (order == 13):
                    if (use_cheby):
                        Global_ZP = GLOBALutils.Cheby_eval( coeffs_pix2wav, int(.5*len(thar_order)), len(thar_order) )
                    else:
                        Global_ZP = np.polyval( coeffs_pix2wav, 0.0 )
                #print residuals
                All_Pixel_Centers = np.append( All_Pixel_Centers, pixel_centers )
                All_Wavelengths = np.append( All_Wavelengths, wavelengths )
                All_Orders = np.append( All_Orders, np.zeros( len(pixel_centers) ) + order )
                All_Centroids = np.append( All_Centroids, centroids)
                All_Sigmas = np.append( All_Sigmas, sigmas)
                All_Intensities = np.append( All_Intensities, intensities )
                All_Residuals = np.append( All_Residuals, residuals)
                All_Sigmas = np.append( All_Sigmas,sigmas)
                trans[order,:] =  coeffs_pix2wav

                order += 1
                orre += 1

            if debugfitTHARlinesresult:
                JJ = np.unique(All_Orders)
                for order in JJ:
                    I = np.where(All_Orders == order)[0]
                    plt.title('ThAr residual after order by order fit')
                    plt.plot(All_Wavelengths[I],All_Residuals[I],'.')
                    plt.xlabel('All Wavelengths')
                    plt.ylabel('All Residuals')
                plt.show()

            # Fit all lines at once assuming physical orders goes by 1/m
            p0 = np.zeros( npar_wsol )
            p0[0] =  (13+oro0) * Global_ZP
            p1, G_pix, G_ord, G_wav, II, rms_ms, G_res = GLOBALutils.Fit_Global_Wav_Solution(All_Pixel_Centers, All_Wavelengths,All_Orders,
                                                                                            np.ones(All_Intensities.shape), p0, Cheby=use_cheby,
                                                                                            maxrms=150,Inv=Inverse_m,minlines=250,order0=oro0,
                                                                                            ntotal=nup,npix=len(thar_order),nx=ncoef_x,nm=ncoef_m,
                                                                                            debug=debugfitTHARlinesresult)

            if debugfitTHARlinesresult:
                plt.figure(figsize=(20,12))
                ejxx = np.arange(thar_order.shape[0])
                for i in np.unique(G_ord):
                    I = np.where(G_ord == i)[0]
                    #print(i,len(I))
                    m = G_ord+oro0
                    m2 = m[I][0]+np.zeros(len(ejxx))
                    chebs = GLOBALutils.Calculate_chebs(ejxx, m2, order0=oro0,
                                                    ntotal=n_useful,
                                                    npix=thar_order.shape[0],nx=ncoef_x,nm=ncoef_m,
                                                    Inverse=Inverse_m)
                    ret = (1.0/m2) * GLOBALutils.Joint_Polynomial_Cheby(p1, chebs,nx=ncoef_x,nm=ncoef_m)
                    plt.plot(ejxx,ret,'k')
                    plt.plot(G_pix,G_wav,'ro')
                    #plt.plot(G_wav[I],G_res[I],'.')
                plt.xlabel('All Pixels')
                plt.ylabel('All Wavelengths')
                plt.show()

            thar_wav_Ss = np.zeros( (2,nup,dthar.shape[1]) )
            equis = np.arange( np.shape(thar_wav_Ss)[2] )
            order = orwa
            orre = 0

            # Save polynomials that are used to locate orders
            while orre < nup:
                m = order + oro0
                chebs = GLOBALutils.Calculate_chebs(equis, m, order0=oro0,ntotal=nup, npix=len(thar_order), Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
                thar_wav_Ss[0,orre,:] = GLOBALutils.ToVacuum( (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(p1,chebs,ncoef_x,ncoef_m) )[::-1]
                thar_wav_Ss[1,orre,:] = thar_Ss[orre][::-1]
                orre += 1
                order+=1

            if (os.access(thar_fits_simple_wav,os.F_OK)):
                os.remove( thar_fits_simple_wav )

            hdu = pyfits.PrimaryHDU( thar_wav_Ss )
            hdu.writeto( thar_fits_simple_wav )

            #if rms_ms/np.sqrt(NL)<50:
            nThAr_ref.append(thar)
            nthtimes.append(thtime)
            pdict = {'p1':p1, 'G_pix':G_pix, 'G_ord':G_ord, 'G_wav':G_wav, 'II':II, 'rms_ms':rms_ms, 'G_res':G_res,\
                    'All_Centroids':All_Centroids, 'All_Sigmas':All_Sigmas, 'trans':trans, 'or0':or0, 'orwa':orwa, 'oro0':oro0}
            pickle.dump( pdict, open( wavsol_pkl, 'wb' ) )
            numt+=1

            debugfitTHARlines = False
        else:
            print("\t\tLoading ThAr file", os.path.basename(thar))

    ######### Collect result from saved files #########
    pdict = pickle.load(open(os.path.join(dirout,'order_find.pkl'),'rb'))
    orwa = pdict['orwa']
    or0 = pdict['or0']
    orwb = pdict['orwb']
    or25 = pdict['or25']
    delta = pdict['delta']
    nup = or25 - or0 + 1
    #ThAr_ref= nThAr_ref
    #print ThAr_ref
    ThAr_ref = np.array(ThAr_ref)
    thtimes = np.array(thtimes)
    I = np.argsort(thtimes)

    thtimes = thtimes[I]
    ThAr_ref= ThAr_ref[I]
    ThAr_ref = list(ThAr_ref)

    ######### Find best ThAr with lowest RMS residual #########
    min_rms = 100000
    thar_min = ThAr_ref[0]
    for thar in ThAr_ref:
        hdth = pyfits.getheader(thar,ext=get_extension(thar))
        wavsol_pkl = os.path.join(dirout,'ThAr_'+os.path.splitext(os.path.basename(thar))[0]+'.wavsolpars.pkl')
        wavsol = pickle.load( open( wavsol_pkl, 'rb' ) )
        if wavsol['rms_ms'] < min_rms:
            min_rms = wavsol['rms_ms']
            thar_min = thar
            if debugglobalvelshift: print('min_rms',min_rms,'for',thar)

    ######### Reference ThAr #########
    hdth = pyfits.getheader(thar_min,ext=get_extension(thar_min))
    thtime = gettime(hdth)
    wavsol_pkl = os.path.join(dirout,'ThAr_'+os.path.splitext(os.path.basename(thar_min))[0]+'.wavsolpars.pkl')
    wavsol = pickle.load( open( wavsol_pkl, 'rb' ) )
    p1_ref = wavsol['p1']
    best_p1 =  wavsol['p1']

    ######### Loop over all ThAr to calculate velocity shift from reference ThAr #########
    f = open( os.path.join(dirout , 'ThAr_shifts.txt') , 'w')
    thshifts = []
    thtimes = []
    for thar in ThAr_ref:
        h = pyfits.open(thar,ext=get_extension(thar))
        h = h[get_extension(thar)]
        mjd,mjd0 = piszkesutils.mjd_fromheader(h)
        hdth = pyfits.getheader(thar,ext=get_extension(thar))
        thtimes.append(mjd)
        wavsol_pkl = os.path.join(dirout,'ThAr_'+os.path.splitext(os.path.basename(thar))[0]+'.wavsolpars.pkl')
        wavsol = pickle.load( open( wavsol_pkl, 'rb' ) )
        G_pix = wavsol['G_pix']
        G_wav = wavsol['G_wav']
        G_ord = wavsol['G_ord']
        oro0 = wavsol['oro0']
        p_shift, pix_centers, orders, wavelengths, I, rms_ms, residuals = GLOBALutils.Global_Wav_Solution_vel_shift(G_pix, G_wav, G_ord,np.ones(G_wav.shape),
                                                                                                                    p1_ref,Cheby=use_cheby,Inv=True,maxrms=150,
                                                                                                                    minlines=250,order0=oro0,ntotal=nup,
                                                                                                                    npix=h.data.shape[1],nx=ncoef_x,nm=ncoef_m,
                                                                                                                    debug=debugglobalvelshift)
        thshifts.append(p_shift)
        f.write(thar+'\t'+str(thtime)+'\t'+str((1e-6*p_shift)*299792.458)+'\n')
        #print p_shift
    thshifts,thtimes = np.array(thshifts),np.array(thtimes)

    ######### Fit polynomial to time vs ThAr velocity shift #########
    Is = np.argsort(thtimes)
    thtimes,thshifts = thtimes[Is],thshifts[Is]
    # If we have >3 ThAr and there time distribution is approx. uniform(ish)
    if len(thtimes) > 3 and not np.any( np.abs(np.diff(thtimes) / np.median(np.diff(thtimes))) > 100 ):
        thtck = scipy.interpolate.splrep(thtimes,thshifts,k=3,s=0)
        ejeje = np.linspace(thtimes[0],thtimes[-1],10000)
        ejeyy = scipy.interpolate.splev(ejeje,thtck,der=0)
        if debugglobalvelshift:
            plt.plot(thtimes,thshifts,'ro')
            plt.plot(ejeje,ejeyy)
            plt.show()
    elif len(thtimes) > 3:
        # If ThAr time distribution is not uniform(ish) use 1st order Bspline
        thtck = scipy.interpolate.splrep(thtimes,thshifts,k=1,s=0)
        ejeje = np.linspace(thtimes[0],thtimes[-1],10000)
        ejeyy = scipy.interpolate.splev(ejeje,thtck,der=0)
        if debugglobalvelshift:
            plt.plot(thtimes,thshifts,'ro')
            plt.plot(ejeje,ejeyy)
            plt.show()
    else:
        thtck = 0
    f.close()

    ################### FLAT frame extraction ##############################
    #print ThAr_ref
    c_all = trace_dict['c_all'][::-1]
    c_all = c_all[or0:or25+1]
    c_all = c_all[::-1]
    print('\n\tExtraction of Flat calibration frames:')
    if not have_flat:
        print(bcolors.WARNING + 'Warning! No FLAT frames.' + bcolors.ENDC)
        print(bcolors.WARNING + 'Using median of science frames instead.' + bcolors.ENDC)
        log.warning('No FLAT frames. Using median of science frames instead.')
    # names of extracted Masterflat (optimal & simple) & of the P flat matrix
    S_flat_fits = os.path.join(dirout,'Masterflat.spec.fits')
    S_flat_fits_simple = os.path.join(dirout,'Masterflat.spec.fits.S')
    flat_P = os.path.join(dirout,'P_Flat.fits')
    sm_flat_fits = os.path.join(dirout,'SmoothMasterflat.spec.fits')
    # name of the extracted background image
    bacfile = os.path.join(dirout,'BkgFlatImage_Flat.fits')
    if ( os.access(S_flat_fits,os.F_OK) == False )  or ( os.access(S_flat_fits_simple,os.F_OK) == False ) or (force_flat_extract) or ( os.access(flat_P,os.F_OK) == False ):

        print("\t\tNo previous extraction or extraction forced for flat file,",  "getting background...")
        #Flat = utils.invert(Flat)

        Flat  = pyfits.getdata(os.path.join(dirout,'MasterFlat.fits'))
        #if Flat.ravel().min() < 0:
        #    Flat -= Flat.ravel().min() + 1
        #Flat  = Flat.T
        if os.access(bacfile,os.F_OK)== False or force_flat_extract:
            Centers = np.zeros((len(c_all),Flat.shape[1]))
            for i in range(c_all.shape[0]):
                Centers[i,:]=np.polyval(c_all[i,:],np.arange(len(Centers[i,:])))
            bac = GLOBALutils.get_scat_fast(Flat,Centers,span=4,debug=debugFLATextraction,allow_neg=True)
            hdbac = pyfits.PrimaryHDU( bac )
            if os.access(bacfile,os.F_OK):
                os.system('rm -r '+bacfile)
            hdbac.writeto(bacfile)

        if debugFLATextraction:
            from astropy.visualization import ZScaleInterval

            vmin,vmax = ZScaleInterval().get_limits(bac)

            plt.figure(figsize=(20,5))
            plt.title('Background of Flat')
            plt.imshow(bac,cmap='gray',origin='lower')
            plt.plot(Centers[0,:],'r--')
            plt.plot(Centers[-1,:],'r--')
            plt.colorbar()
            plt.show()

            vmin,vmax = ZScaleInterval().get_limits(Flat)

            plt.figure(figsize=(20,5))
            plt.title('Original Flat')
            plt.imshow(Flat,vmin=vmin,vmax=vmax,cmap='gray',origin='lower')
            plt.plot(Centers[0,:],'r--')
            plt.plot(Centers[-1,:],'r--')
            plt.colorbar()
            plt.show()

            vmin,vmax = ZScaleInterval().get_limits(Flat-bac)

            plt.figure(figsize=(20,5))
            plt.title('Background Corrected Flat')
            plt.imshow(Flat-bac,vmin=vmin,vmax=vmax,cmap='gray',origin='lower')
            plt.plot(Centers[0,:],'r--')
            plt.plot(Centers[-1,:],'r--')
            plt.colorbar()
            plt.show()

        Flat -= bac

        print("\t\tNo previous extraction or extraction forced for flat file,",  "extracting...")
        # Determination the P matrix
        if os.access(flat_P,os.F_OK):
            P = pyfits.getdata(flat_P)
        else:
            P = GLOBALutils.obtain_P(Flat,c_all,ext_aperture,roF,gaF,NSigma_Marsh, S_Marsh,N_Marsh, Marsh_alg, min_extract_col, max_extract_col, npools)
            hdu = pyfits.PrimaryHDU(P)
            hdu.writeto(flat_P)

        flat_S  = GLOBALutils.optimal_extraction(Flat,P,c_all,ext_aperture,roF,gaF,S_Marsh,NCosmic_Marsh,min_extract_col,max_extract_col,npools)
        sm_flat = flat_S.copy()
        flat_Ss = GLOBALutils.simple_extraction(Flat,c_all,ext_aperture,min_extract_col,max_extract_col,npools)
        if onlysimpleextract:
            sm_flat[:,1,:] = flat_Ss.copy()

        flat_S  = flat_S[::-1]
        flat_Ss = flat_Ss[::-1]
        sm_flat = sm_flat[::-1]
        # save as fits file
        if (os.access(S_flat_fits,os.F_OK)):
            os.remove( S_flat_fits )
        if (os.access(S_flat_fits_simple,os.F_OK)):
            os.remove( S_flat_fits_simple )
        if (os.access(flat_P,os.F_OK)):
            os.remove( flat_P )
        if (os.access(sm_flat_fits,os.F_OK)):
            os.remove( sm_flat_fits )

        hdu = pyfits.PrimaryHDU( flat_S )
        hdu.writeto( S_flat_fits )
        hdu = pyfits.PrimaryHDU( flat_Ss )
        hdu.writeto( S_flat_fits_simple )
        hdu = pyfits.PrimaryHDU( P )
        hdu.writeto( flat_P )
        hdu = pyfits.PrimaryHDU( sm_flat )
        hdu.writeto( sm_flat_fits )

    # recover Flat spectra and P matrix
    else:
        print("\t\tLoading previous extraction for flat file")
        flat_S = pyfits.getdata(S_flat_fits)
        flat_Ss = pyfits.getdata(S_flat_fits_simple)
        P = pyfits.getdata(flat_P)
        sm_flat = pyfits.getdata( sm_flat_fits )

    ##################### Science images extraction ############################################
    sm_flat, norms = GLOBALutils.FlatNormalize_single( sm_flat[:,1,:], mid=int(0.5*sm_flat.shape[2]),span=200)


    ##### Fit Flat with polynomial to smooth out wave-like structures #####
    for ii in range(sm_flat.shape[0]):
        sm_flat[ii,:] = piszkesutils.FlatPolynomial_single(sm_flat[ii,:],nc=7,debug=debugfitflat)
        if not have_flat:
            sm_flat[ii,:] = 1

    ##### Specify which files will be reduced #####
    new_list = []
    new_list_obnames = []

    for i in range(len(objects)):
        fsim = objects[i]
        obname = piszkesutils.search_name(fsim,log=log)
        if (object2do == 'all'):
            new_list.append(fsim)
            new_list_obnames.append( obname )
        else:
            if (os.path.basename(objects[i]) == object2do):
                new_list.append(fsim)
                new_list_obnames.append( obname )


    print('\n\tThe following targets will be processed:')
    for nlisti in range(len(new_list)):
        print('\t\t'+new_list_obnames[nlisti]+' ('+os.path.basename(new_list[nlisti])+')')

    # Does any image have a special requirement for dealing with the moonlight?
    if os.access(dirin + 'moon_corr.txt', os.F_OK):
        fmoon = open(dirin + 'moon_corr.txt','r')
        moon_lns = fmoon.readlines()
        spec_moon = []
        use_moon = []
        for line in moon_lns:
            spec_moon.append(line.split()[0])
            #if line.split()[1] == '0':
            #    use_moon.append(False)
            #else:
            use_moon.append(True)
            log.info('Moonlight considered for %s' % line)
    else:
        spec_moon = []
        use_moon = []

    ################ Calibrate, extract all spectra and do measurements ###############
    fin_spec = []
    have_specph = False
    for obj in sorted(new_list):

        know_moon = False
        if obj.split('/')[-1] in spec_moon:
            I = np.where(np.isin(spec_moon , obj.split('/')[-1]))[0][0]
            know_moon = True
            here_moon = use_moon[I]

        print("\t--> Working on image: ", os.path.basename(obj))
        index1,index2 = -1,-1
        h = pyfits.open(obj)
        h = h[get_extension(obj)]
        mjd,mjd0 = piszkesutils.mjd_fromheader(h)
        hd = pyfits.getheader(obj,ext=get_extension(obj))
        obname = piszkesutils.search_name(obj,log=log)
        exptime = hd['EXPTIME']

        print("\t\tObject name:",obname)

        # Coordinates of RCC in Piszkesteto
        altitude    = 934.6
        latitude    = 47. + 55./60. + 6./3600.
        longitude   = 19. + 53./60. + 41.7/3600.
        epoch = 2000.

        try:
            ra          = h.header['RA']
            dec         = h.header['DEC']
            c = SkyCoord(ra,dec,unit=(u.hourangle, u.deg), frame='icrs')
            ra = c.ra.value
            dec = c.dec.value
            known_coords = True
            log.info('RA DEC loaded from header')
        except:
            log.warning('No RA DEC found in header')
            known_coords = False

        # proper motion corrected coordinates
        if known_coords:
            try:
                ra2,dec2 = piszkesutils.get_coords(ra,dec,mjd,mjd0)
                ra = ra2
                dec = dec2
                print('\t\tUsing coordinates from gaia query.')
                log.info('Using coordinates from gaia query.')
            except:
                print('\t\tUsing coordinates from header.')
                log.info('Using coordinates from header.')

        if known_coords:
            earthloc = EarthLocation.from_geodetic(lat=latitude*u.deg, lon=longitude*u.deg, height=altitude*u.m)
            sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
            barycorr = sc.radial_velocity_correction(obstime=Time(mjd, format='mjd', scale='utc'), location=earthloc)
            bcvel_baryc = barycorr.to(u.km/u.s).value

            print("\t\tBarycentric velocity:", np.round(bcvel_baryc,2))

            times = Time(mjd, format='mjd', scale='utc', location=earthloc)
            ltt_bary = times.light_travel_time(sc)
            mbjd = times.tdb + ltt_bary
            mbjd = mbjd.value
        else:
            bcvel_baryc = 0.
            mbjd = mjd

        # Moon Phase Calculations
        gobs      = ephem.Observer()
        gobs.name = 'piszkes'
        gobs.lat  = rad(latitude)  # lat/long in decimal degrees
        gobs.long = rad(longitude)

        gobs.date = hd['DATE-OBS'].replace('T',' ')

        mephem = ephem.Moon()
        mephem.compute(gobs)
        lunation,moon_state,moonsep,moonvel = piszkesutils.get_lunar_props(ephem,gobs,ra,dec,latitude,longitude,altitude,mjd,mjd0,ephempath=ephempath)
        refvel = bcvel_baryc + moonvel.value
        print('\t\tRadial Velocity of scattered moonlight:',refvel)

        #print utils.search_name(obj)
        #nama = piszkesutils.search_name(obj)+'_'+hd['DATE-OBS'][:10]+'_'+hd['DATE-OBS'][11:13]+'-'+hd['DATE-OBS'][14:16]+'-'+hd['DATE-OBS'][17:]
        nama = os.path.splitext(os.path.basename(obj))[0]
        obj_fits = os.path.join(dirout,nama+'.spec.fits.S')
        bkg_obj_fits = os.path.join(dirout,'Bkg_'+nama+'.fits')
        obj_fits_simple = os.path.join(dirout,nama+'.spec.simple.fits.S')

        # Correct and extract science spectra
        if ( os.access(obj_fits,os.F_OK) == False )  or \
           ( os.access(obj_fits_simple,os.F_OK) == False ) or \
           (force_sci_extract):
            print("\t\tNo previous extraction or extraction forced for science file", os.path.basename(obj), ", extracting...")

            # Raw science image
            dat = pyfits.getdata(obj).astype('float')
            if debugsciencecorrection:
                from astropy.visualization import ZScaleInterval
                vmin,vmax = ZScaleInterval().get_limits(dat)
                plt.figure(figsize=(20,5))
                plt.title('Original science image')
                plt.imshow(dat,cmap='gray',origin='lower',vmin=vmin,vmax=vmax)
                plt.colorbar()
                plt.show()

            # Cosmic ray removal
            if remove_cosmic:
                cosmicsigma = 5
                while True:
                    datclean = ccdproc.cosmicray_lacosmic(dat,sigclip=cosmicsigma,readnoise=RON,
                    gain_apply=False)
                    if np.where(datclean[1])[0].shape[0] < 5000:
                        break
                    else:
                        # if number if cosmics > 5000, increase sigma
                        cosmicsigma += 30
                if debugsciencecorrection:
                    vmin,vmax = ZScaleInterval().get_limits(dat)
                    plt.figure(figsize=(20,5))
                    ax = plt.gca()
                    plt.title('Cosmic rays')
                    plt.imshow(dat,cmap='gray',origin='lower',vmin=vmin,vmax=vmax)
                    for xx,yy in zip(np.where(datclean[1])[1] , np.where(datclean[1])[0]):
                        circle = plt.Circle((xx, yy), 10, color='r', fill=False)
                        ax.add_patch(circle)
                    plt.colorbar()
                    plt.show()

                    vmin,vmax = ZScaleInterval().get_limits(dat)
                    plt.figure(figsize=(20,5))
                    plt.title('Cosmic ray removed science image')
                    plt.imshow(datclean[0],cmap='gray',origin='lower',vmin=vmin,vmax=vmax)
                    plt.colorbar()
                    plt.show()

                dat = datclean[0]
                if isinstance(dat, u.Quantity):
                    dat = dat.value

            # Bias correction
            dat -= MasterBias
            if debugsciencecorrection:
                vmin,vmax = ZScaleInterval().get_limits(dat)
                plt.figure(figsize=(20,5))
                plt.title('Bias corrected science image')
                plt.imshow(dat,cmap='gray',origin='lower',vmin=vmin,vmax=vmax)
                plt.colorbar()
                plt.show()
            # Dark correction
            if len(MasDarl)>0 and dark_corr:
                dark = piszkesutils.get_dark(MasDarl, hd['EXPTIME'])
                dat -= dark
                #dat = dat.T
                if debugsciencecorrection:
                    vmin,vmax = ZScaleInterval().get_limits(dat[50:,])
                    plt.figure(figsize=(20,5))
                    plt.title('Dark corrected science image')
                    plt.imshow(dat,cmap='gray',origin='lower',vmin=vmin,vmax=vmax)
                    #plt.plot(Centers[0,:],'r--')
                    #plt.plot(Centers[-1,:],'r--')
                    plt.colorbar()
                    plt.show()

            # Vertical drift
            trimmeddat = dat.copy() # Cut out lower 50 pixels to remove bad lines
            trimmeddat[:50,] = 0
            drift, c_new = GLOBALutils.get_drift(trimmeddat,P,c_all,pii=dat.shape[1]//2,win=5)
            if debugsciencecorrection: print('\t\ty drift:', drift)

            P_new = GLOBALutils.shift_P(P,drift,c_new,ext_aperture)

            # Background (scattered light) correction
            if back_corr:
                Centers = np.zeros((len(c_new),dat.shape[1]))
                for i in range(c_new.shape[0]):
                    Centers[i,:]=np.polyval(c_new[i,:],np.arange(len(Centers[i,:])))
                bac   = GLOBALutils.get_scat_fast(dat, Centers,span=5,allow_neg=True,debug=debugbkgestimation)
                dat -= bac
                if debugsciencecorrection:
                    vmin,vmax = ZScaleInterval().get_limits(dat[50:,])
                    plt.figure(figsize=(20,5))
                    plt.title('Bkg of science image')
                    plt.imshow(bac,cmap='gray',origin='lower',vmin=vmin,vmax=vmax)
                    plt.colorbar()
                    plt.show()
                    vmin,vmax = ZScaleInterval().get_limits(dat[50:,])
                    plt.figure(figsize=(20,5))
                    plt.title('Bkg corrected science image')
                    plt.imshow(dat,cmap='gray',origin='lower',vmin=vmin,vmax=vmax)
                    plt.plot(Centers[0,:],c='r',ls='--')
                    plt.colorbar()
                    plt.show()

            #dat -= dat.ravel().min() + 1
            # Tracing and extraction of science orders
            if not have_flat:
                P_new = GLOBALutils.obtain_P(dat,c_new,ext_aperture,RON,\
                                        GAIN,NSigma_Marsh, S_Marsh, \
                        N_Marsh, Marsh_alg, min_extract_col, max_extract_col,  npools)
            obj_Ss = GLOBALutils.simple_extraction(dat,c_new,ext_aperture,min_extract_col,max_extract_col,npools)
            obj_S  = GLOBALutils.optimal_extraction(dat,P_new,c_new,ext_aperture,\
                                                        RON,GAIN,S_Marsh,NCosmic_Marsh,min_extract_col,max_extract_col,npools)

            obj_Ss = obj_Ss[::-1]
            obj_S  = obj_S[::-1]
            # 2 point convolution to reduce noise
            boxc = Box1DKernel(2)
            for ii in range(obj_Ss.shape[0]):
                obj_Ss[ii,:]  = convolve(obj_Ss[ii,:],boxc,boundary='wrap')
                obj_S[ii,1,:] = convolve(obj_S[ii,1,:],boxc,boundary='wrap')
            if onlysimpleextract:
                obj_S[:,1,:] = obj_Ss

            # Save as fits file
            if (os.access(obj_fits,os.F_OK)):
                os.remove( obj_fits )
            if (os.access(obj_fits_simple,os.F_OK)):
                os.remove( obj_fits_simple )

            hdu = pyfits.PrimaryHDU( obj_S )
            hdu.writeto( obj_fits )
            hdu = pyfits.PrimaryHDU( obj_Ss )
            hdu.writeto( obj_fits_simple )

        else:
            obj_S = pyfits.getdata(obj_fits)
            obj_Ss = pyfits.getdata(obj_fits_simple)

        # Time delay between scince spectrum and ThAr times
        deltat=1000
        i = 0
        while i < len(thtimes):
            if abs(thtimes[i]-mjd) < deltat:
                index1 = i
                deltat = abs(thtimes[i]-mjd)
            i+=1

        # Collect best ThAr spectra
        if mjd < thtimes[0]:
            print(bcolors.WARNING + "Problem with ThAr and science times:" + bcolors.ENDC)
            print(bcolors.WARNING + "\tNo ThAr taken before science frames" + bcolors.ENDC)
            log.warning("Problem with ThAr and science times: <mjd")
            index1 = 0
            index2 = 0
            indexx = 0
        elif mjd > thtimes[-1]:
            print(bcolors.WARNING + "Problem with ThAr and science times:" + bcolors.ENDC)
            print(bcolors.WARNING + "\tNo ThAr taken after science frames" + bcolors.ENDC)
            log.warning("Problem with ThAr and science times: >mjd")
            index1 = -1
            index2 = -1
            indexx = -1
        else:
            for i in range(len(thtimes)-1):
                if mjd >= thtimes[i] and mjd < thtimes[i+1]:
                    index1 = i
                    index2 = i+1
                    if abs(mjd - thtimes[i]) < abs(mjd - thtimes[i+1]):
                        indexx = i
                    else:
                        indexx = i+1
                    break

        # Load best ThAr spectra
        #print ThAr_ref[index1], obj, ThAr_ref[index2]
        hdth = pyfits.getheader(ThAr_ref[index1],ext=get_extension(ThAr_ref[index1]))
        wavsol_pkl = os.path.join(dirout,'ThAr_'+os.path.splitext(os.path.basename(ThAr_ref[index1]))[0]+'.wavsolpars.pkl')
        pdict = pickle.load(open(wavsol_pkl,'rb'))
        global1 = pdict['p1']
        All_Pixel_Centers = pdict['G_pix']
        All_Orders = pdict['G_ord']
        All_Wavelengths = pdict['G_wav']
        rms_ms = pdict['rms_ms']
        All_Residuals = pdict['G_res']
        All_Centroids = pdict['All_Centroids']
        All_Sigmas = pdict['All_Sigmas']
        trans = pdict['trans']
        or0 = pdict['or0']
        orwa = pdict['orwa']
        oro0 = pdict['oro0']

        hdth = pyfits.getheader(ThAr_ref[index2],ext=get_extension(ThAr_ref[index2]))
        wavsol_pkl2 = os.path.join(dirout,'ThAr_'+os.path.splitext(os.path.basename(ThAr_ref[index2]))[0]+'.wavsolpars.pkl')
        pdict2 = pickle.load(open(wavsol_pkl2,'rb'))
        global2 = pdict2['p1']

        hdth = pyfits.getheader(ThAr_ref[indexx],ext=get_extension(ThAr_ref[indexx]))
        wavsol_pklx = os.path.join(dirout,'ThAr_'+os.path.splitext(os.path.basename(ThAr_ref[indexx]))[0]+'.wavsolpars.pkl')
        pdictx = pickle.load(open(wavsol_pklx,'rb'))
        globalx = pdictx['p1']

        # Arrays for collecting results
        #if JustExtract or not have_flat:
        #    final = np.zeros( [6, nup,np.shape(obj_S)[2]] )
        #else:
        final = np.zeros( [11, nup,np.shape(obj_S)[2]] )
        # Create header
        hdu = pyfits.PrimaryHDU( final )
        hdu = GLOBALutils.update_header(hdu,'OBJECT', obname)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MJD', mjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH MBJD', mbjd)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH BJD', mbjd + 2400000.5)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START DATE', hd['DATE-OBS'][:10] )
        hdu = GLOBALutils.update_header(hdu,'HIERARCH SHUTTER START UT',  hd['DATE-OBS'][:11])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TEXP (S)',hd['EXPTIME'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH BARYCENTRIC CORRECTION (KM/S)', bcvel_baryc)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH TARGET NAME', obname)
        try:
            hdu = GLOBALutils.update_header(hdu,'HIERARCH RA',hd['RA'])
            hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC',hd['DEC'])
            hdu = GLOBALutils.update_header(hdu,'HIERARCH RA BARY',ra)
            hdu = GLOBALutils.update_header(hdu,'HIERARCH DEC BARY',dec)
        except:
            None
        hdu = GLOBALutils.update_header(hdu,'HIERARCH EQUINOX',hd['EPOCH'])
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS LATITUDE',latitude)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH vOBS LONGITUDE',longitude)
        hdu = GLOBALutils.update_header(hdu,'HIERARCH OBS ALTITUDE',altitude)
        try:
            hdu = GLOBALutils.update_header(hdu,'HIERARCH MOON_VEL',refvel.value,'[km/s]')
            hdu = GLOBALutils.update_header(hdu,'HIERARCH MOONST',moon_state)
            hdu = GLOBALutils.update_header(hdu,'HIERARCH LUNATION',lunation)
            hdu = GLOBALutils.update_header(hdu,'HIERARCH MOONSEP',moonsep)
            hdu = GLOBALutils.update_header(hdu,'HIERARCH MOONALT',float(mephem.alt))
            hdu = GLOBALutils.update_header(hdu,'HIERARCH SMOONALT',str(mephem.alt))
        except:
            None

        ######### Loop over all orders in given science image ########
        equis = np.arange( np.shape(obj_S)[2] )
        order = orwa
        orre = 0
        while orre < nup:
            # Get wavelength solutions for ThAr spectra
            m = order + oro0
            chebs = GLOBALutils.Calculate_chebs(equis, m, order0=oro0, npix=len(equis), ntotal=nup, Inverse=Inverse_m,nx=ncoef_x,nm=ncoef_m)
            WavSol1 = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(global1,chebs,ncoef_x,ncoef_m)
            WavSol2 = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(global2,chebs,ncoef_x,ncoef_m)
            WavSolx = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(globalx,chebs,ncoef_x,ncoef_m)
            bestWavSol = (1.0/m) * GLOBALutils.Joint_Polynomial_Cheby(best_p1,chebs,ncoef_x,ncoef_m)

            # Interpolate wavelength solutions
            if thtimes[index2] == thtimes[index1]:
                WavSol = WavSol1
            else:
                pen = (WavSol2 - WavSol1) / (thtimes[index2] - thtimes[index1])
                coe = WavSol2 - pen* thtimes[index2]
                WavSol = pen*mjd + coe
                if thtck == 0.:
                    p_shift = 0.
                else:
                    p_shift = scipy.interpolate.splev(mjd, thtck, der=0)
                    #p_shift = 0
                WavSol = bestWavSol * (1.0 + 1.0e-6*p_shift)
            final[0,orre,:] = GLOBALutils.ToVacuum(WavSol[::-1])

            #final[0,orre,:] = WavSolx[::-1]
            if len(np.where(np.isnan(obj_S[orre,2,:][::-1]))[0]) < 100:
                final[1,orre,:] = obj_S[orre,1,:][::-1]   # ...flux in ADU (optimal extraction)...
                final[2,orre,:] = obj_S[orre,2,:][::-1]   # ...and 1/variance...
                if have_flat:
                    final[3,orre,:] = final[1,orre,:] / sm_flat[orre,:][::-1]   # ...flux (optimal)/flat...
                    final[4,orre,:] = final[2,orre,:] * (sm_flat[orre,:][::-1] ** 2)   # ...and 1/variance * flat**2...
                    cont = GLOBALutils.get_cont_single(final[0,orre,:],final[3,orre,:],final[4,orre,:], nc = 5,lu=3,span=10)
                    ratio = np.polyval(cont,final[0,orre])
                    final[6,orre,:] = final[4,orre,:] * (ratio**2)
                    final[7,orre,:] = ratio
                    final[8,orre,:] = ratio * sm_flat[orre,:][::-1] / np.sqrt( ratio * sm_flat[orre,:][::-1] / GAIN + (RON/GAIN)**2 )
                    final[5,orre,:] = final[3,orre,:] / ratio
                    nJ = np.where(np.isnan(final[5,orre]))[0]
                    nJ2 = np.where(np.isinf(final[5,orre]))[0]
                    final[5,orre,nJ] = 1.0
                    final[5,orre,nJ2] = 1.0
                    #plot(final[8,orre])
                    if not keepoutliers:
                        rI = np.where(final[5,orre] > 1. + 8./final[8,orre])
                        final[5,orre,rI] = 1.
                    spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), final[0,orre,:],k=3)
                    dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
                    NN            = np.average(dlambda_dx)
                    dlambda_dx    /= NN
                    final[9,orre] = final[5,orre] * (dlambda_dx ** 1)
                    final[10,orre] = final[6,orre] / (dlambda_dx ** 2)

                    if debugscienceorders:
                        checkthisorder = orre

                        rw,rf = final[0,orre],obj_Ss[orre,:][::-1]
                        cbl = piszkesutils.fit_blaze(rw,rf,min_extract_col,debug=debugfitblaze,debugiter=debugfitblazeiter)
                        ratio = np.polyval(cbl,rw)
                        finalsimple = rf/ratio

                        print('##### Order ',checkthisorder, '#####')
                        plt.figure(figsize=(20,6))
                        plt.plot(obj_S[checkthisorder,1,:],label='Marsh extracted object')
                        plt.plot(obj_Ss[checkthisorder,:]-500,label='Simple extracted object')
                        plt.legend()
                        plt.show()

                        plt.figure(figsize=(20,6))
                        plt.plot(flat_S[checkthisorder,1,:],label='Marsh flat')
                        plt.plot(flat_Ss[checkthisorder,:]-500,label='Simple flat')
                        plt.plot(sm_flat[orre,:]*np.max(flat_S[checkthisorder,1,:]))
                        plt.legend()
                        plt.show()

                        plt.figure(figsize=(20,5))
                        plt.plot(final[0,checkthisorder,:],final[1,checkthisorder,:],label='Original spectrum')
                        plt.plot(final[0,checkthisorder,:],final[3,checkthisorder,:],label='Blaze corrected spectrum')
                        plt.plot(final[0,checkthisorder,:],final[7,checkthisorder,:],label='Low-order polynomial')
                        plt.ylim(0,final[7,checkthisorder,:].max()+final[1,checkthisorder,:].std())
                        plt.legend()
                        plt.show()
                        plt.figure(figsize=(20,5))
                        plt.plot(final[0,checkthisorder,:],final[5,checkthisorder,:],label='Blaze + poly corrected')
                        plt.plot(final[0,checkthisorder,:],finalsimple-0.5,label='Simple extract+corr')
                        plt.ylim(-0.5,2)
                        plt.grid()
                        plt.legend()
                        plt.show()
                        plt.figure(figsize=(20,5))
                        plt.plot(final[0,checkthisorder,:],final[8,checkthisorder,:])
                        plt.ylabel('S/N')
                        plt.show()
                else:
                    # If there is NO flat
                    rw,rf = final[0,orre],final[1,orre]
                    cbl = piszkesutils.fit_blaze(rw,rf,min_extract_col,debug=debugfitblaze,debugiter=debugfitblazeiter)
                    ratio = np.polyval(cbl,rw)
                    final[3,orre] = rf/ratio
                    final[5,orre] = rf/ratio
                    final[4,orre] = final[2,orre]*(ratio**2)
                    final[6,orre] = final[2,orre]*(ratio**2)
                    final[7,orre] = ratio
                    final[8,orre] = ratio / np.sqrt( ratio / GAIN + (RON/GAIN)**2 )
                    medflx = np.zeros(len(final[3,orre,:]))
                    nI = np.where(np.isfinite(final[3,orre]))[0]
                    medflx = scipy.signal.medfilt(final[3,orre,:][nI],3)
                    res = final[3,orre,:][nI] - medflx
                    dev = np.sqrt(np.var(res))
                    nI = np.where(np.isinf(final[3,orre]))[0]
                    nJ = np.where(np.isnan(final[3,orre]))[0]
                    final[3,orre][nI]=1.
                    final[3,orre][nJ]=1.
                    if not keepoutliers:
                        rI = np.where(final[3,orre,:] > 1. + 5*dev)[0]
                        final[3,orre,rI] = 1.
                        rI = np.where(final[5,orre] > 1. + 8./final[8,orre])
                        final[5,orre,rI] = 1.
                    spl           = scipy.interpolate.splrep(np.arange(WavSol.shape[0]), final[0,orre,:],k=3)
                    dlambda_dx    = scipy.interpolate.splev(np.arange(WavSol.shape[0]), spl, der=1)
                    NN            = np.average(dlambda_dx)
                    dlambda_dx    /= NN
                    final[9,orre] = final[5,orre] * (dlambda_dx ** 1)
                    final[10,orre] = final[6,orre] / (dlambda_dx ** 2)

                    if debugscienceorders:
                        checkthisorder = orre

                        plt.figure(figsize=(20,5))
                        plt.plot(final[0,checkthisorder,:],final[3,checkthisorder,:],label='Simple corr')
                        plt.ylim(0,2)
                        plt.grid()
                        plt.legend()
                        plt.show()
                        plt.figure(figsize=(20,5))
                        plt.plot(final[0,checkthisorder,:],final[5,checkthisorder,:])
                        plt.ylabel('S/N')
                        plt.show()

            order += 1
            orre += 1

        if (not JustExtract):
            spec = final[:,:,::-1].copy()  # reverse wavelength order
            # query SIMBAD by coordinates for spectral type
            print('\t\tQuerying Simbad for spectral type')
            query_success,sp_type_query = piszkesutils.simbad_query_sptype(obj,log)
            if query_success:
                print("\t\t\tSpectral type returned by SIMBAD query:",sp_type_query)
                log.info('Spectral type returned by SIMBAD query')
            else:
                print("\t\t\tSpectral type not found")
                log.info('Spectral type not found')
            hdu = GLOBALutils.update_header(hdu,'HIERARCH SIMBAD SPTYP', sp_type_query)

            #DoClass = False
            if DoClass:
                print('\t\tSpectral Analysis:')
                # spectral analysis

                pars_file = os.path.join(dirout , obj.split('/')[-1][:-4]+'_stellar_pars.txt' )

                if os.access(pars_file,os.F_OK) == False or force_stellar_pars:
                    print("\t\t\tEstimating atmospheric parameters:")
                    T_eff, logg, Z, vsini, vel0, ccf = correlation.CCF(spec,model_path=models_path,npools=npools,
                                                                       velmin=-200.0,velmax=200.0,debug=debugsciencespectrumfit,
                                                                       base=os.path.join(PACKAGEDIR,'Correlation') )
                    line = "%6d %4.1f %4.1f %8.1f %8.1f\n" % (T_eff,logg, Z, vsini, vel0)
                    f = open(pars_file,'w')
                    f.write(line)
                    f.close()

                else:
                    print("\t\t\tAtmospheric parameters loaded from file:")
                    T_eff, logg, Z, vsini, vel0 = np.loadtxt(pars_file,unpack=True)

                print("\t\t\t\tT_eff=",T_eff,"log(g)=",logg,"Z=",Z,"vsin(i)=",vsini,"vel0=",vel0)

            else:
                T_eff, logg, Z, vsini, vel0 = -999,-999,-999,-999,-999

            T_eff_epoch = T_eff
            logg_epoch  = logg
            Z_epoch     = Z
            vsini_epoch = vsini
            vel0_epoch  = vel0
            hdu = GLOBALutils.update_header(hdu,'HIERARCH TEFF', float(T_eff))
            hdu = GLOBALutils.update_header(hdu,'HIERARCH LOGG', float(logg))
            hdu = GLOBALutils.update_header(hdu,'HIERARCH Z', Z)
            hdu = GLOBALutils.update_header(hdu,'HIERARCH VSINI', vsini)
            hdu = GLOBALutils.update_header(hdu,'HIERARCH VEL0', vel0)

            print("\t\tRadial Velocity analysis:")
            # assign mask
            mask_path=os.path.join(PACKAGEDIR,"data/xc_masks")
            sp_type, mask = piszkesutils.get_mask_reffile(obj,reffile=reffile,base=mask_path,sp_type_query=sp_type_query,T_eff=T_eff)
            print("\t\t\tWill use",sp_type,"mask for CCF.")

            ######### Read in mask ########
            ml, mh, weight = np.loadtxt(mask,unpack=True)
            goodlines = (ml>4500) & (mh<6800)
            ml = ml[goodlines]
            mh = mh[goodlines]
            weight = weight[goodlines]
            ml_v = GLOBALutils.ToVacuum( ml )
            mh_v = GLOBALutils.ToVacuum( mh )

            # make mask larger accounting for factor ~2 lower res in CORALIE w/r to HARPS
            av_m = 0.5*( ml_v + mh_v )
            ml_v -= 5*(av_m - ml_v)
            mh_v += 5*(mh_v - av_m)
            mask_hw_kms = (GLOBALutils.Constants.c/1e3) * 0.5*(mh_v - ml_v) / av_m

            #sigma_fout = stellar_pars_dir + obname + '_' +'sigma.txt'

            disp = GLOBALutils.get_disp(obj, reffile=reffile)

            if disp == 0:
                known_sigma = False
                if vsini != -999 and vsini != 0.:
                    disp = vsini
                    known_sigma = True
                else:
                    disp = 30.

                    print(bcolors.WARNING + '\t\tWarning! There is no predefined dispersion of the CCF.' + bcolors.ENDC)
                    print(bcolors.WARNING + '\t\tFor more reliable RVs provide a reffile or use option -do_class.' + bcolors.ENDC)
            else:
                known_sigma = True

            if disp < 5:
                disp = 5
            mask_hw_wide = av_m * disp / (GLOBALutils.Constants.c/1.0e3)
            ml_v = av_m - mask_hw_wide
            mh_v = av_m + mask_hw_wide

            print('\t\t\tComputing the CCF...')
            cond = True
            lbary_ltopo = bcvel_baryc/2.99792458E5 + 1
            while (cond):
                ######### first rough correlation to find the minimum ########
                vels, xc_full, sn, nlines_ccf, W_ccf = \
                    GLOBALutils.XCor(spec, ml_v, mh_v, weight, 0, lbary_ltopo, vel_width=500,vel_step=1.,spec_order=9,iv_order=10,sn_order=8,max_vel_rough=300)

                # minimum S/N for orders to be averaged for CCF
                #if np.sum(sn>15) > 10:
                #    sn_min=15
                #elif np.sum(sn>10) > 10:
                #    sn_min=10
                #elif np.sum(sn>5) > 10:
                #    sn_min=5
                #else:
                #    sn_min=0
                sn_min=5

                #xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=0.0,    Simple=True, W=W_ccf)
                xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=sn_min, Simple=True, W=W_ccf)
                #print W_ccf
                # Normalize the continuum of the CCF
                yy = scipy.signal.medfilt(xc_av,11)
                pred = piszkesutils.get_CCF_bkg(vels,yy,order=1)
                tck1 = scipy.interpolate.splrep(vels,pred,k=1)
                xc_av_orig = xc_av.copy()
                xc_av /= pred
                vel0_xc = vels[ np.argmin( xc_av ) ]
                rvels, rxc_av, rpred, rxc_av_orig, rvel0_xc = vels.copy(), \
                xc_av.copy(), pred.copy(), xc_av_orig.copy(), vel0_xc

                ######## Zoom in on minimum of wide range correlation then fitting CCF ########
                xc_av_rough = xc_av
                vels_rough  = vels
                if disp > 50:
                    disp = 50.
                vel_width = np.maximum( 40.0, 6*disp )

                vels, xc_full, sn, nlines_ccf, W_ccf =\
                GLOBALutils.XCor(spec, ml_v, mh_v, weight, vel0_xc, lbary_ltopo, vel_width=vel_width,vel_step=1.,  spec_order=9,iv_order=10,sn_order=8,max_vel_rough=300)
                #print W_ccf

                xc_av = GLOBALutils.Average_CCF(xc_full, sn, sn_min=sn_min, Simple=True, W=W_ccf)
                pred = scipy.interpolate.splev(vels,tck1)
                xc_av /= pred

                if sp_type == 'M5':
                    moon_sig = 5.0
                elif sp_type == 'K5':
                    moon_sig = 6.6
                else:
                    moon_sig = 9.0

                ######## Fitting CCF with a Gaussian ########
                p1,XCmodel,p1gau,XCmodelgau,Ls2 = GLOBALutils.XC_Final_Fit( vels, xc_av ,\
                                      sigma_res = 4, horder=8, moonv = refvel, moons = moon_sig, moon = False)

                twoGmodel = piszkesutils.XCfit_two_gaussians(rvels, rxc_av, disp, debug=debugsciencespectrumfit)

                #ldc = CoralieUtils.get_ldc(T_eff, logg, Z, 1.0, ldfile = 'lin_coe_sloan2.dat')
                #p1R, ROTmodel = CoralieUtils.XC_Final_Fit_Rot( vels, xc_av, ldc = ldc, vsini = vsini )
                moonmatters = False
                if (know_moon and here_moon):
                    moonmatters = True
                    ismoon = True
                    confused = False
                    p1_m,XCmodel_m,p1gau_m,XCmodelgau_m,Ls2_m = GLOBALutils.XC_Final_Fit( vels, xc_av , sigma_res = 4, horder=8, moonv = refvel, moons = moon_sig, moon = True)
                    moon_flag = 1
                else:
                    confused = False
                    ismoon = False
                    p1_m,XCmodel_m,p1gau_m,XCmodelgau_m,Ls2_m = p1,XCmodel,p1gau,XCmodelgau,Ls2
                    moon_flag = 0

                ######### BiSpector Span ########
                bspan = GLOBALutils.calc_bss(vels,xc_av,debug=debugsciencespectrumfit,log=log)
                if np.isnan(bspan[0]):
                    SP = -999.0
                else:
                    SP = bspan[0]
                if debugsciencespectrumfit: print('Bisector span:', SP)

                # If dispersion was not known, update it using Gaussian fit to CCF
                if (not known_sigma):
                    disp = np.floor(p1gau[2])
                    if (disp < 3.0):
                        disp = 3.0
                    mask_hw_wide = av_m * disp / (GLOBALutils.Constants.c/1.0e3)
                    ml_v = av_m - mask_hw_wide
                    mh_v = av_m + mask_hw_wide
                    known_sigma = True
                else:
                    cond = False

            xc_dict = {'vels':vels,'xc_av':xc_av,'XCmodelgau':XCmodelgau,'Ls2':Ls2,'refvel':refvel,\
               'rvels':rvels,'rxc_av':rxc_av,'rpred':rpred,'rxc_av_orig':rxc_av_orig,\
               'rvel0_xc':rvel0_xc,'xc_full':xc_full, 'p1':p1, 'sn':sn, 'p1gau':p1gau,\
               'p1_m':p1_m,'XCmodel_m':XCmodel_m,'p1gau_m':p1gau_m,'Ls2_m':Ls2_m,\
               'XCmodelgau_m':XCmodelgau_m,'twoGmodel':twoGmodel}

            moon_dict = {'moonmatters':moonmatters,'moon_state':moon_state,'moonsep':moonsep,\
                     'lunation':lunation,'mephem':mephem,'texp':exptime}

            pkl_xc = os.path.join( dirout , obj.split('/')[-1][:-4]+obname+'_XC_'+sp_type+'.pkl' )
            pickle.dump( xc_dict, open( pkl_xc, 'wb' ) )

            SNR_5130 = np.median(spec[8,15,700:901])
            airmass  = -999
            seeing   = -999

            if sp_type == 'G2':
                if T_eff < 6000:
                    A = 0.06544
                    B = 0.11066 #0.00146 (from betaCVn scatter)
                    D = 0.77210 #0.24416 (original*sqrt(n points))
                    C = 0.00181
                else:
                    A = 0.09821
                    B = 0.11066 #0.00014
                    D = 1.05908 #0.33491
                    C = 0.00113
            elif  sp_type == 'K5':
                A = 0.05348
                B = 0.11066 #0.00147
                D = 0.65443 #0.20695
                C = 0.00321
            else:
                A = 0.05348
                B = 0.11066 #0.00147
                D = 0.65443 #0.20695
                C = 0.00321

            RVerr =  B + ( 1.6 + 0.2 * p1gau[2] ) * A / np.round(SNR_5130)
            RVtwoG1err =  B + ( 1.6 + 0.2 * twoGmodel.stddev_0.value ) * A / np.round(SNR_5130)
            RVtwoG2err =  B + ( 1.6 + 0.2 * twoGmodel.stddev_1.value ) * A / np.round(SNR_5130)
            depth_fact = 1. + p1gau[0]/(p1gau[2]*np.sqrt(2*np.pi))
            depth_fact_twoG1 = 1 + (-1-twoGmodel.amplitude_0.value)/(twoGmodel.stddev_0.value*np.sqrt(2*np.pi))
            depth_fact_twoG2 = 1 + (-1-twoGmodel.amplitude_1.value)/(twoGmodel.stddev_1.value*np.sqrt(2*np.pi))

            # single Gaussian error
            if depth_fact >= 1.:
                RVerr2 = -999.000
            else:
                if sp_type == 'G2':
                    depth_fact = (1 - 0.62) / (1 - depth_fact)
                else:
                    depth_fact = (1 - 0.59) / (1 - depth_fact)
                RVerr2 = RVerr * depth_fact
                if (RVerr2 <= 0.009):
                    RVerr2 = 0.009
            # two Gaussians error 1
            if depth_fact_twoG1 >= 1.:
                RVtwoG1err2 = -999.000
            else:
                if sp_type == 'G2':
                    depth_fact_twoG1 = (1 - 0.62) / (1 - depth_fact_twoG1)
                else:
                    depth_fact_twoG1 = (1 - 0.59) / (1 - depth_fact_twoG1)
                RVtwoG1err2 = RVtwoG1err * depth_fact_twoG1
                if (RVtwoG1err2 <= 0.009):
                    RVtwoG1err2 = 0.009
            # two Gaussians error 2
            if depth_fact_twoG2 >= 1.:
                RVtwoG2err2 = -999.000
            else:
                if sp_type == 'G2':
                    depth_fact_twoG2 = (1 - 0.62) / (1 - depth_fact_twoG2)
                else:
                    depth_fact_twoG2 = (1 - 0.59) / (1 - depth_fact_twoG2)
                RVtwoG2err2 = RVtwoG2err * depth_fact_twoG2
                if (RVtwoG2err2 <= 0.009):
                    RVtwoG2err2 = 0.009

            if RVerr2 < 0.3:
                RVerr2 = 0.3
            if RVtwoG1err2 < 0.3:
                RVtwoG1err2 = 0.3
            if RVtwoG2err2 < 0.3:
                RVtwoG2err2 = 0.3
            try:
                BSerr = D / float(np.round(SNR_5130)) + C
            except ZeroDivisionError as e:
                BSerr = -999.000
                log.error('BS error calculation failed')
                log.error(e)

            RV     = np.around(p1gau_m[1],4)
            RVtwoG1     = np.around(twoGmodel.mean_0.value,4)
            RVtwoG2     = np.around(twoGmodel.mean_1.value,4)
            RVerr2 = np.around(RVerr2,4)
            RVtwoG1err2 = np.around(RVtwoG1err2,4)
            RVtwoG2err2 = np.around(RVtwoG2err2,4)
            BS     = np.around(SP,4)
            BSerr = np.around(BSerr,4)


            print('\t\t\tRV = '+str(RV)+' +- '+str(RVerr2))
            print('\t\t\tBS = '+str(BS)+' +- '+str(BSerr))
            print('\t\t\t2 gauss RV 1 = '+str(RVtwoG1)+' +- '+str(RVtwoG1err2))
            print('\t\t\t2 gauss RV 2 = '+str(RVtwoG2)+' +- '+str(RVtwoG2err2))

            bjd_out = 2400000.5 + mbjd
            T_eff_err = 300
            logg_err = 0.5
            Z_err = 0.5
            vsini_err = 2
            XC_min = np.abs(np.around(np.min(XCmodel),2))

            SNR_5130 = np.around(SNR_5130)
            SNR_5130_R = np.around(SNR_5130*np.sqrt(2.9))

            disp_epoch = np.around(p1gau_m[2],1)
            hdu = GLOBALutils.update_header(hdu,'RV', RV)
            try:
                hdu = GLOBALutils.update_header(hdu,'RV_E', RVerr2)
            except ValueError:
                hdu = GLOBALutils.update_header(hdu,'RV_E', -999.000)
            RVsyserror = 1.89114
            hdu = GLOBALutils.update_header(hdu,'RV_E_SYS', RVsyserror)
            hdu = GLOBALutils.update_header(hdu,'BS', BS)
            hdu = GLOBALutils.update_header(hdu,'BS_E', BSerr)
            hdu = GLOBALutils.update_header(hdu,'DISP', disp_epoch)
            hdu = GLOBALutils.update_header(hdu,'XC_MIN', XC_min)
            hdu = GLOBALutils.update_header(hdu,'GaussRV1', RVtwoG1)
            hdu = GLOBALutils.update_header(hdu,'GaussRV2', RVtwoG1)
            try:
                hdu = GLOBALutils.update_header(hdu,'RV1_E', RVtwoG1err2)
            except ValueError:
                hdu = GLOBALutils.update_header(hdu,'RV1_E', -999.000)
            try:
                hdu = GLOBALutils.update_header(hdu,'RV2_E', RVtwoG1err2)
            except ValueError:
                hdu = GLOBALutils.update_header(hdu,'RV2_E', -999.000)
            hdu = GLOBALutils.update_header(hdu,'SNR', SNR_5130)
            hdu = GLOBALutils.update_header(hdu,'SNR_R', SNR_5130_R)
            hdu = GLOBALutils.update_header(hdu,'INST', 'RCC1m')
            hdu = GLOBALutils.update_header(hdu,'RESOL', '20000')
            hdu = GLOBALutils.update_header(hdu,'PIPELINE', 'CERES')

            # Save CCF plot
            ccf_pdf = os.path.join(dirout , 'proc' , nama + '_XCs_' + sp_type + '.pdf' )
            if not avoid_plot:
                GLOBALutils.plot_CCF(xc_dict,moon_dict,path=ccf_pdf)

            # Save CCF data
            ccf_fits = os.path.join(dirout , 'proc' , nama + '_XCs_' + sp_type + '.fits')
            piszkesutils.save_CCF(xc_dict,ccf_fits,obname,mjd,mbjd,bjd_out,bcvel_baryc,hd,ra,dec,refvel,moon_state,lunation,moonsep,mephem,
                                  sp_type,RV,RVerr2,BS,BSerr,disp_epoch,SNR_5130,SNR_5130_R,XC_min)

            line_out = "%s,%s,%.8f,%.4f,%.4f,%.3f,%.3f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,rcc1m,piszkespipe,20000,%d,%d,%.2f,%.2f,%.2f,%.2f,%.1f,%.2f,%.2f,%.1f,%d,%s\n"%\
                      (obname, nama, bjd_out, RV, RVerr2, BS, BSerr, bcvel_baryc, RV+bcvel_baryc, RVerr2, RVsyserror,\
                      RVtwoG1, RVtwoG1+bcvel_baryc, RVtwoG1err2, RVtwoG2, RVtwoG2+bcvel_baryc, RVtwoG2err2,\
                      T_eff_epoch, T_eff_err, logg_epoch, logg_err, Z_epoch, Z_err, vsini_epoch, XC_min, disp_epoch,exptime, SNR_5130_R, ccf_pdf)
            f_res.write(line_out)

        # Save combined spectrum
        if not JustExtract:
            wcombined,fcombined,errors = piszkesutils.combine_orders(final,lbary_ltopo)
            combined_fits = os.path.join(dirout , 'proc' , nama + '_combined_barycorrected.fits' )
            piszkesutils.save_combined_spectrum(combined_fits,wcombined,fcombined,errors,obname,mjd,mbjd,bjd_out,bcvel_baryc,hd,ra,dec,SNR_5130,SNR_5130_R)
        else:
            wcombined,fcombined,errors = piszkesutils.combine_orders(final,None)
            combined_fits = os.path.join(dirout , 'proc' , nama + '_combined.fits' )
            piszkesutils.save_combined_spectrum(combined_fits,wcombined,fcombined,errors,obname,mjd,mbjd,2400000.5+mbjd,None,hd,ra,dec,None,None)

        # Save final spectrum orders
        fout = 'proc/'+nama+'_final.fits'
        if (os.access( os.path.join(dirout , fout) ,os.F_OK)):
            os.remove( os.path.join(dirout , fout))
        hdu.writeto( os.path.join(dirout , fout) )


    f_res.close()


    print('Done')

    return 0


#########################
# Command-line interfaces
#########################

def piszkespipe_from_commandline(args=None):
    """Use piszkespipe from command-line."""
    # Recieve input parameters
    parser = argparse.ArgumentParser(
                description="piszkespipe reduces echelle spectrum obtained in Piszkesteto.",
                formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('directorio',
                        help="Path to the raw data.")
    parser.add_argument('-dirout',default='default',
                        help="Path to the directory where the pipeline products will "
                            "be placed. The default path will be a new directory with "
                            "the same name that the input directory but followed by a '_red' suffix.")
    parser.add_argument('-do_class', action="store_true", default=False,
                        help="This option will enable the estimation of atmospheric parameters.")
    parser.add_argument('-just_extract', action="store_true", default=False,
                        help="If activated, the code will not compute the CCF and atmospheric parameters.")
    parser.add_argument('-o2do',default='all',
                        help="If you want to process just one particular science object you have to "
                            "enter this option followed by the filename of the object.")
    parser.add_argument('-reffile',default='default',
                        help="Name of the auxiliary file. The default is './reffile.txt', "
                        "a file located in the directory where the raw data is.")
    parser.add_argument('-npools', default=-1,
                        help="Number of CPU cores to be used by the code. Default is all.")
    parser.add_argument('-marsch', action="store_true", default=False,
                        help="If enabled, Marsch optimized raw extraction algorithm will be saved.")
    parser.add_argument('-nocosmic', action="store_true", default=False,
                        help="If activated, no cosmic ray identification will be performed.")
    parser.add_argument('-keepoutliers', action="store_true", default=False,
                        help="If activated, strong upper (emission line-like) outliers won't be removed.")
    parser.add_argument('-avoid_plot', action="store_true", default=False,
                        help="If activated, the code will not generate a pdf file "
                            "with the plot of the computed CCF.")

    args   = parser.parse_args()

    # Input parameters
    dirin             = args.directorio
    avoid_plot        = args.avoid_plot
    dirout            = args.dirout
    DoClass           = args.do_class
    JustExtract       = args.just_extract
    npools            = int(args.npools)
    object2do         = args.o2do
    reffile           = args.reffile
    keepoutliers      = args.keepoutliers
    if args.marsch:
        onlysimpleextract = False
    else:
        onlysimpleextract = True
    if args.nocosmic:
        remove_cosmic     = False
    else:
        remove_cosmic     = True

    print("\n\n\tPISZKESTETO RCC1.0m  PIPELINE\n")

    _ = piszkespipe(dirin,
                    avoid_plot=avoid_plot,
                    dirout=dirout,
                    DoClass=DoClass,
                    JustExtract=JustExtract,
                    npools=npools,
                    object2do=object2do,
                    reffile=reffile,
                    onlysimpleextract=onlysimpleextract,
                    remove_cosmic=remove_cosmic,
                    keepoutliers=keepoutliers)
