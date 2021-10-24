import matplotlib
#matplotlib.use("Agg")

base = '../'

import piszkespipe.GLOBALutils as GLOBALutils

import numpy as np
import scipy
from astropy.io import fits as pyfits
import os
import glob
import matplotlib.pyplot as plt
import requests

from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
from astropy.time import Time

np.seterr(all='ignore')
import warnings
warnings.filterwarnings("ignore",category=RuntimeWarning)

def is_there(string, word):
    l=len(word)
    i=0
    ist = False
    while i < len(string)-l:
        if string[i:i+l] == word:
            ist = True
        i+=1
    return ist

def search_name(file,log):
    from astroquery.simbad import Simbad
    from astroquery.exceptions import TableParseError
    from pathlib import Path

    h = pyfits.open(file)[get_extension(file)].header
    ra  = h['RA']
    dec = h['DEC']

    try:
        result_table = Simbad.query_region(SkyCoord(ra,dec,unit=(u.hourangle, u.deg), frame='icrs'), radius=40*u.arcsec)
        name = result_table['MAIN_ID'][0].strip('* ').strip('V*').replace(' ','')
    except (TypeError,TableParseError,IndexError):
        name = Path(file).stem
        log.warning('%s %s not found in Simbad! Using %s instead.' % (ra,dec,name))

    return name

def FileClassify(path,log,debug=True):
    biases   = []
    flats = []
    img_flats = []
    fib_flats = []
    objects  = []
    darks    = []
    thars    = []


    lines = []
    dates = []
    archs = glob.glob(path+'*.fit*')

    bad_files = []
    if os.access(path+'bad_files.txt',os.F_OK):
        bf = open(path+'bad_files.txt')
        linesbf = bf.readlines()
        for line in linesbf:
            bad_files.append(path+line[:-1])

        bf.close()
    for arch in archs:
        dump = False
        for bf in bad_files:
            if arch == bf:
                dump = True
                break
        if not dump:

            h = pyfits.open(arch)

            if '.fz' in arch:
                extension = 1
            else:
                extension = 0

            try:
                #if h[0].header['CCDXBIN'] == 1 and h[0].header['CCDXBIN'] == 1:
                if h[extension].header['IMAGETYP'].upper() == 'OBJECT':
                    if 'flat' in arch.lower():
                        flats.append(arch)
                    else:
                        ra = h[extension].header['RA']
                        dec = h[extension].header['DEC']
                        expt = h[extension].header['EXPTIME']
                        date = h[extension].header['DATE-OBS']
                        JD = h[extension].header['JD']
                        line = "%10s %10s %8.2f %8s %s\n" % (ra, dec, expt, date, arch)

                        ye = float(date[:4])
                        mo = float(date[5:7])
                        da = float(date[8:10])
                        ho = float(date[11:13])-4.0
                        mi = float(date[14:15])
                        se = float(date[17:])

                        lines.append(line)
                        dates.append( JD )
                        #f.write(line)

                        if is_there(arch.lower(),'thar')  or is_there(arch.lower(),'th_ar'):
                            thars.append(arch)
                        else:
                            objects.append(arch)

                elif  (h[extension].header['IMAGETYP'] == 'Bias Frame' or h[extension].header['IMAGETYP'].upper() == 'BIAS') and 'MasterBias' not in arch:
                    biases.append(arch)
                elif  (h[extension].header['IMAGETYP'] == 'Flat Frame' or h[extension].header['IMAGETYP'].upper() == 'FLAT') and 'MasterFlat' not in arch:
                    # Now check which kind of flat it is.
                    # Maybe a "surface" flat...
                    if(is_there(arch.lower(),'imgflat')):
                        img_flats.append(arch)
                    # ...a fibre flat...
                    elif(is_there(arch.lower(),'fibre')):
                        fib_flats.append(arch)
                        # (use them for traces, blaze and col-to-col)
                        flats.append(arch)
                    # ...else, it is a screen flat (w/difussor):
                    else:
                        flats.append(arch)
                elif  (h[extension].header['IMAGETYP'] == 'Dark Frame' or h[extension].header['IMAGETYP'].upper() == 'DARK') and 'MasterDark' not in arch:
                    if h[extension].header['EXPTIME']!=0.0:
                        darks.append(arch)
            except KeyError:
                print('\tWarning! There is no IMAGETYP keyword. Classifying file based on filename!')
                log.warning('There is no IMAGETYP keyword. Classifying file based on filename!')

                if 'bias' in arch.lower() and 'MasterBias' not in arch:
                    biases.append(arch)
                elif 'flat' in arch.lower() and 'MasterFlat' not in arch:
                    # Now check which kind of flat it is.
                    # Maybe a "surface" flat...
                    if(is_there(arch.lower(),'imgflat')):
                        img_flats.append(arch)
                    # ...a fibre flat...
                    elif(is_there(arch.lower(),'fibre')):
                        fib_flats.append(arch)
                        # (use them for traces, blaze and col-to-col)
                        flats.append(arch)
                    # ...else, it is a screen flat (w/difussor):
                    else:
                        flats.append(arch)
                elif 'dark' in arch.lower() and 'MasterDark' not in arch:
                    darks.append(arch)
                else:
                    ra = h[extension].header['RA']
                    dec = h[extension].header['DEC']
                    expt = h[extension].header['EXPTIME']
                    date = h[extension].header['DATE-OBS']
                    JD = h[extension].header['JD']
                    line = "%10s %10s %8.2f %8s %s\n" % (ra, dec, expt, date, arch)

                    ye = float(date[:4])
                    mo = float(date[5:7])
                    da = float(date[8:10])
                    ho = float(date[11:13])-4.0
                    mi = float(date[14:15])
                    se = float(date[17:])

                    lines.append(line)
                    dates.append( JD )
                    #f.write(line)

                    objects.append(arch)

            h.close()
    lines = np.array(lines)
    dates = np.array(dates)
    I = np.argsort(dates)
    lines = lines[I]
    for line in lines:
        log.info('%s',line)

    if debug:
        print('BIAS:')
        for b in biases:
            print('\t'+b)
        print('DARK:')
        for b in darks:
            print('\t'+b)
        print('FLAT:')
        for b in flats+img_flats+fib_flats:
            print('\t'+b)
        print('OBJECT:')
        for b in objects:
            print('\t'+b)
        print('ThAr:')
        for b in thars:
            print('\t'+b)

    return biases,flats,img_flats,fib_flats,objects,thars,darks

def get_extension(h):
    if '.fz' in h:
        extension = 1
    else:
        extension = 0
    return extension

def get_rg(h,b, debug=False):

    gain = float( pyfits.open(h)[get_extension(h)].header['GAIN'] )
    ron  = float( np.sqrt(pyfits.open(b)[get_extension(b)].data.mean()) )

    if debug:
        print('RON=',ron,'GAIN=',gain)

    return ron,gain

def MedianCombine(ImgList,zero_bo,zero,ronoise,gain,dark_bo=False, dlist = []):
    """
    Median combine a list of images
    """

    hf = pyfits.getheader(ImgList[0],ext=get_extension(ImgList[0]))

    if zero_bo:
        Master = pyfits.getdata(zero,ext=get_extension(zero))
    if dark_bo:
        Dark = get_dark(dlist,hf['EXPTIME'])

    n = len(ImgList)

    if n==0:
        raise ValueError("empty list provided!")

    d = pyfits.getdata(ImgList[0],ext=get_extension(ImgList[0]))

    if zero_bo:
        d = d - Master
    if dark_bo:
        d = d - Dark

    if (n == 1):
        return d, ronoise, gain

    else:
        for i in range(n-1):
            h = pyfits.getdata(ImgList[i+1],get_extension(ImgList[i+1]))
            if zero_bo:
                h = h-Master
            if dark_bo:
                h = h-Dark
            d = np.dstack((d,h))

        return np.median(d,axis=2), ronoise/np.sqrt(n), gain

def get_dark(darks,t):
    exact = 0
    dts = []
    for dark in darks:
        hd = pyfits.getheader(dark)
        dt = hd['EXPTIME']
        dts.append(dt)
        if dt == t:
            #print 'dark:',dark
            DARK = pyfits.getdata(dark)
            exact = 1

    dts = np.array(dts)
    if exact == 0:
        if t < dts.min():
            I = np.where( dts == dts.min() )[0]
            DARK = pyfits.getdata(darks[I[0]])*t/dts[I[0]]
        elif t > dts.max():
            I = np.where( dts == dts.max() )[0]
            DARK = pyfits.getdata(darks[I[0]])*t/dts[I[0]]
            #print darks[I[0]]
        else:
            tmin = dts.min()
            tmax = dts.max()
            I = np.where( dts == dts.min() )[0]
            Dmin = pyfits.getdata(darks[I[0]])
            Dminname=darks[I[0]]
            I = np.where( dts == dts.max() )[0]
            Dmax = pyfits.getdata(darks[I[0]])
            Dmaxname = darks[I[0]]

            i = 0
            while i < len(dts):
                if dts[i] < t and dts[i] > tmin:
                    tmin = dts[i]
                    Dminname = darks[i]
                    Dmin = pyfits.getdata(darks[i])
                elif dts[i] > t and dts[i] < tmax:
                    tmax = dts[i]
                    Dmaxname = darks[i]
                    Dmax = pyfits.getdata(darks[i])
                i+=1

            num = Dmax - Dmin
            den = tmax-tmin
            m = num/den
            n = Dmax - m*tmax
            DARK = m*t+n

    return DARK


def jd(y,m,d,h,mins,s):
       "Julian day is calculated here if it's needed"
       MY = (m-14)/12
       y = MY+y
       return ( 1461 * ( y + 4800 ) ) / 4  + (  367 * ( m - 2 - 12*MY ) ) / 12 - ( 3 * ( ( y + 4900 ) / 100 ) ) / 4 + d -32077.5 + htosec(h,mins,s)/86400.0

def htosec(h,m,s):
       "transform from hour,minute and seconds, to seconds"
       return s+60.0*(m+60.0*h)

def fit_blaze(w,f,min_extract_col,n=5,debug=False,debugiter=False):
    w = w[min_extract_col:-min_extract_col]
    f = f[min_extract_col:-min_extract_col]
    worig = w.copy()
    forig = f.copy()
    li = len(w)

    # Pre-cleaning of emission lines
    co = np.polyfit(w,f,2)
    res = f - np.polyval(co,w)
    dev = np.sqrt(np.var(res))
    I = np.where(res < 2.5*dev)[0]
    cond = True
    if len(I) < .3*li:
        cond = False
    while cond:
        w,f = w[I],f[I]
        ntps = len(f)
        co = np.polyfit(w,f,2)
        if debugiter:
            plt.title('Elliminating emission lines')
            plt.plot(worig,forig,'k')
            plt.plot(w,f,'r')
            plt.plot(worig,np.polyval(co,worig),'g')
            plt.show()
        res = f - np.polyval(co,w)
        dev = np.sqrt(np.var(res))
        I = np.where(res < 2.5*dev)[0]
        if ntps==len(I) or len(I) < .3*li:
            cond = False

    # Iterative fitting
    co = np.polyfit(w,f,n)
    if debug:
        plt.title('Fitting blaze initial')
        plt.plot(worig,forig)
        plt.plot(w,f)
        plt.plot(worig,np.polyval(co,worig),'r')
        plt.show()
    res = f - np.polyval(co,w)
    dev = np.sqrt(np.var(res))
    J1 = np.where(res < -1.5*dev)[0]
    J2 = np.where(res > 3*dev)[0]
    J = np.hstack((J1,J2))
    J = np.sort(J)
    I = np.where( (res >= -1.5*dev) & (res <= 3*dev) )[0]
    cond = True
    if len(J)==0 or len(I) < .3*li:
        cond = False
    while cond:
        w,f = w[I],f[I]
        co = np.polyfit(w,f,n)
        if debugiter:
            plt.plot(worig,forig)
            plt.plot(w,f)
            plt.plot(worig,np.polyval(co,worig),'r')
            plt.show()
        res = f - np.polyval(co,w)
        dev = np.sqrt(np.var(res))
        J1 = np.where(res < -1.5*dev)[0]
        J2 = np.where(res > 3*dev)[0]
        J = np.hstack((J1,J2))
        J = np.sort(J)
        I = np.where( (res >= -1.5*dev) & (res <= 3*dev) )[0]
        cond = True
        if len(J)==0 or len(I) < .3*li:
            cond = False

    if debug:
        plt.title('Fitting blaze final')
        plt.plot(worig,forig)
        plt.plot(w,f)
        plt.plot(worig,np.polyval(co,worig),'r')
        plt.show()

    # Fill edges if missing
    # IMPORT: wavelengths are in reversed order!
    refit = False
    missing_wav = 10 # AA
    if w.max() < worig.max()-missing_wav:
        if debug: print('Filling end')
        # If last 20 AA has been removed extrapolate it back
        fillfrom = np.argmin(np.abs(worig - (worig.max()-missing_wav) )) + 1
        filledflux = np.poly1d(np.polyfit(w,f,2))(worig[:fillfrom])
        w = np.append(worig[:fillfrom], w)
        f = np.append(filledflux , f)
        refit = True

    if w.min() > worig.min()+missing_wav:
        if debug: print('Filling beginning')
        # If first 20 AA has been removed extrapolate it back
        fillto = np.argmin(np.abs(worig - (worig.min()+missing_wav) ))
        filledflux = np.poly1d(np.polyfit(w,f,2))(worig[fillto:])
        w = np.append(w , worig[fillto:])
        f = np.append(f , filledflux)
        refit = True

    if refit:
        co = np.polyfit(w,f,n)
        if debug:
            plt.title('After edge filling')
            plt.plot(worig,forig,'k')
            plt.plot(w,f,'C1')
            plt.plot(worig,np.polyval(co,worig),'r')
            plt.show()

    return co

def mjd_fromheader(h):
    """
    return modified Julian date from header
    """

    datetu = h.header['TIME-BEG']
    t = Time(datetu, format='isot', scale='utc')
    jd  = t.jd
    mjd = t.mjd
    mjd0 = jd-mjd

    secinday = 24*3600.0
    fraction = 0.5
    texp     = h.header['EXPTIME'] #sec

    mjd = mjd + (fraction * texp) / secinday

    return mjd, mjd0

def get_coords(ra,dec,mjd,mjd0):
    jd = mjd + mjd0
    # Get the positions of the source
    c1 = SkyCoord(ra, dec, frame='icrs', unit='deg')
    rad = 5*30.0  # arcseconds
    from astroquery.vizier import Vizier
    Vizier.ROW_LIMIT = -1
    result = Vizier.query_region(c1, catalog=["I/345/gaia2"],
                                radius=Angle(rad, "arcsec"))
    no_targets_found_message = ValueError('Either no sources were found in the query region '
                                          'or Vizier is unavailable')
    if result is None or len(result) == 0:
        raise no_targets_found_message
    dist = []
    for r in result:
        c2 = SkyCoord(r['RA_ICRS'],r['DE_ICRS'], frame='icrs')
        sep = c1.separation(c2)
        dist.append( sep.value )
    result = result["I/345/gaia2"].to_pandas()
    result = result.assign( separation = np.array(dist).ravel())
    result.sort_values('separation',inplace=True)
    result = result.iloc[0]
    radecs = np.vstack([result['RA_ICRS'], result['DE_ICRS']]).T
    year = ((jd - 2457206.375) * u.day).to(u.year)
    pmra = ((np.nan_to_num(np.asarray(result.pmRA)) * u.milliarcsecond/u.year) * year).to(u.deg).value
    pmdec = ((np.nan_to_num(np.asarray(result.pmDE)) * u.milliarcsecond/u.year) * year).to(u.deg).value
    result.RA_ICRS += pmra
    result.DE_ICRS += pmdec

    return result['RA_ICRS'],result['DE_ICRS']

def get_lunar_props(ephem,gobs,RA,DEC,latitude,longitude,altitude,mjd,mjd0,ephempath='./'):
    mephem = ephem.Moon()
    mephem.compute(gobs)
    pfm = ephem.previous_full_moon(gobs.date)
    nfm = ephem.next_full_moon(gobs.date)
    pnm = ephem.previous_new_moon(gobs.date)
    nnm = ephem.next_new_moon(gobs.date)
    if gobs.date - pnm < gobs.date - pfm:
        moon_state = 'crescent'
        lunation   = (gobs.date-pnm)/(nfm-pnm)
    else:
        moon_state = 'waning'
        lunation   = 1. - (gobs.date-pfm)/(nnm-pfm)

    # now compute the radial velocity sun-moon
    ram  = mephem.ra
    decm = mephem.dec

    from skyfield.api import N, W, load, load_file, wgs84

    ts = load.timescale()

    planets = load_file(ephempath)
    earth = planets['earth']
    moon = planets['moon']
    sun = planets['sun']

    # Altitude and azimuth in the sky of a
    # specific geographic location

    piszkes = earth + wgs84.latlon(latitude * N, longitude * W, elevation_m=altitude)
    t = ts.ut1_jd(mjd+mjd0)

    astro = piszkes.at(t).observe(moon)
    app = astro.apparent()

    Mvx,Mvy,Mvz = app.velocity.km_per_s
    Mx,My,Mz = app.position.km
    moon_obs = (Mx*Mvx + My*Mvy +Mz*Mvz)/np.sqrt(Mx**2 + My**2 +Mz**2)
    moon_obs = moon_obs * u.km/u.s
    #print('moon_obs',moon_obs)

    Mx,My,Mz = moon.at(t).position.km
    Mvx,Mvy,Mvz = moon.at(t).velocity.km_per_s
    Sx,Sy,Sz = sun.at(t).position.km
    Svx,Svy,Svz = sun.at(t).velocity.km_per_s

    SMvx  = Svx - Mvx
    SMvy  = Svy - Mvy
    SMvz  = Svz - Mvz
    SMx   = Sx - Mx
    SMy   = Sy - My
    SMz   = Sz - Mz
    sun_moon = (SMx*SMvx + SMy*SMvy +SMz*SMvz)/np.sqrt(SMx**2 + SMy**2 +SMz**2)        # This is the radial velocity between the sun and the moon
    sun_moon = sun_moon * u.km/u.s
    #print('sun_moon',sun_moon)

    moonvel =  sun_moon + moon_obs
    #print('moonvel',moonvel)

    rasep    = np.sqrt((RA - ram*360./24.)**2)
    decsep   = np.sqrt((DEC - decm)**2)
    if rasep > 180:
        rasep = 360 - rasep
    if decsep > 180:
        decsep = 360 -decsep
    moonsep2 = np.sqrt( (rasep)**2 + (decsep)**2 )
    moonvel =  sun_moon + moon_obs
    return lunation,moon_state,moonsep2,moonvel

def FlatPolynomial_single(flat,nc=7,debug=False):
    flatcorr = scipy.ndimage.median_filter(flat,151)
    if debug:
        plt.figure(figsize=(20,5))
        plt.plot(flat,label='Flat')
        plt.plot(flatcorr,label='Medfilt',ls="--")

    flatx = np.arange(flat.shape[0])
    flatcorr = np.poly1d( np.polyfit(flatx,flatcorr,nc) )(np.arange(flat.shape[0]))

    if debug:
        plt.plot(flatcorr,'r',lw=2,label='Final Fit')
        plt.legend()
        plt.show()

    return flatcorr

def simbad_query_sptype(file,log):
    from astroquery.simbad import Simbad
    Simbad.add_votable_fields('sptype')
    from astroquery.exceptions import TableParseError

    h = pyfits.open(file)[get_extension(file)].header
    ra  = h['RA']
    dec = h['DEC']

    try:
        try:
            result_table = Simbad.query_region(SkyCoord(ra,dec,unit=(u.hourangle, u.deg), frame='icrs'), radius=40*u.arcsec)
        except (TableParseError,IndexError):
            result_table = None

        query_success = False
        if result_table is None or len(result_table['SP_TYPE'][0]) == 0:
            query_success = False
            sp_type_query = None
        else:
            query_success = True
            sp_type_query = result_table['SP_TYPE'][0]
    except requests.exceptions.ConnectTimeout as e:
        print('Warning! Simbad is not available.')
        log.warning('Simbad is not available.')
        log.warning(e)
        query_success = False
        sp_type_query = None

    return query_success, sp_type_query

def get_mask_reffile(obname,reffile='reffile.txt',base='../../xc_masks/',sp_type_query=None,T_eff=-999.):
    import re

    xc_masks = [os.path.join(base,'G2.mas'),\
                os.path.join(base,'K5.mas'),\
                os.path.join(base,'M2.mas') ]
    sp_type = 'G2'

    '''
    try:
        f = open(reffile,'r')
        lines = f.readlines()
        f.close()

        found = False
        for line in lines:
            cos = line.split(',')
            msk = cos[6].strip()
            if cos[0] == obname and (msk=='G2' or msk=='K5' or msk=='M2'):
                sp_type = cos[6][:2]
                found = True
                break
        if not found:
            print( '\t\tWarning! Target not found in reference mask file.')
    except:
        try:
            sp_type = re.match("[O,B,A,F,G,K,M][0-9]",sp_type_query)[0]
            if sp_type[0] in ['O','B','A','F','G']:
                sp_type = 'G2'
            elif sp_type[0] in ['K']:
                sp_type = 'K5'
            elif sp_type[0] in ['M']:
                sp_type = 'M2'
            else:
                sp_type = 'G2'
        except TypeError:
            if T_eff > 0:
                if T_eff < 4000:
                    sp_type = 'M2'
                elif T_eff < 5400:
                    sp_type = 'K5'
                else:
                    sp_type = 'G2'
            else:
                print( '\t\tWarning! Problem with reference mask file and Simbad spectral type. Forcing to G2 mask')
    '''

    found = False
    try:
        f = open(reffile,'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            cos = line.split(',')
            msk = cos[6].strip()
            if cos[0] == os.path.basename(obname) and (msk=='G2' or msk=='K5' or msk=='M2'):
                sp_type = cos[6][:2]
                found = True
                break
    except:
        pass
    if not found:
        print( '\t\tWarning! Target not found in reference mask file.')

    if sp_type == 'G2':
        xc_mask = xc_masks[0]
    elif sp_type == 'K5':
        xc_mask = xc_masks[1]
    else:
        xc_mask = xc_masks[2]

    return sp_type, xc_mask

def save_CCF(xc_dict,path,obname,mjd,mbjd,bjd,bcvel_baryc,hd,ra,dec,refvel,moon_state,lunation,moonsep,mephem,sp_type,RV,RVerr2,
                BS,BSerr,disp_epoch,SNR_5130,SNR_5130_R,XC_min):

    rvels        = xc_dict['rvels']
    rxc_av       = xc_dict['rxc_av']
    rxc_av_orig  = xc_dict['rxc_av_orig']

    data = np.array([rvels,rxc_av_orig,rxc_av])
    # Create header
    hdu = pyfits.PrimaryHDU( data )
    hdu = GLOBALutils.update_header(hdu,'TTYPE1', 'VEL')
    hdu = GLOBALutils.update_header(hdu,'TTYPE2', 'XC ORIG')
    hdu = GLOBALutils.update_header(hdu,'TTYPE3', 'XC CORR')
    hdu = GLOBALutils.update_header(hdu,'TUNIT1', 'KM/S')
    hdu = GLOBALutils.update_header(hdu,'TUNIT2', '1-CCF')
    hdu = GLOBALutils.update_header(hdu,'TUNIT3', '1-CCF')
    hdu = GLOBALutils.update_header(hdu,'TARGET', obname)
    hdu = GLOBALutils.update_header(hdu,'MJD', mjd)
    hdu = GLOBALutils.update_header(hdu,'MBJD', mbjd)
    hdu = GLOBALutils.update_header(hdu,'BJD', bjd)
    try:
        hdu = GLOBALutils.update_header(hdu,'RA',hd['RA'])
        hdu = GLOBALutils.update_header(hdu,'DEC',hd['DEC'])
        hdu = GLOBALutils.update_header(hdu,'RA BARY',ra)
        hdu = GLOBALutils.update_header(hdu,'DEC BARY',dec)
    except:
        None
    hdu = GLOBALutils.update_header(hdu,'EQUINOX',hd['EPOCH'])
    try:
        hdu = GLOBALutils.update_header(hdu,'MOON_VEL',refvel.value,'[km/s]')
        hdu = GLOBALutils.update_header(hdu,'MOONST',moon_state)
        hdu = GLOBALutils.update_header(hdu,'LUNATION',lunation)
        hdu = GLOBALutils.update_header(hdu,'MOONSEP',moonsep)
        hdu = GLOBALutils.update_header(hdu,'MOONALT',float(mephem.alt))
        hdu = GLOBALutils.update_header(hdu,'SMOONALT',str(mephem.alt))
    except:
        None
    hdu = GLOBALutils.update_header(hdu,'HIERARCH BARYCENTRIC CORRECTION (KM/S)', bcvel_baryc)
    hdu = GLOBALutils.update_header(hdu,'MASK',sp_type)
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
    hdu = GLOBALutils.update_header(hdu,'SNR', SNR_5130)
    hdu = GLOBALutils.update_header(hdu,'SNR_R', SNR_5130_R)
    hdu = GLOBALutils.update_header(hdu,'INST', 'RCC1m')
    hdu = GLOBALutils.update_header(hdu,'RESOL', '20000')
    hdu = GLOBALutils.update_header(hdu,'PIPELINE', 'CERES')

    if os.access(path,os.F_OK):
        os.remove( path )
    hdu.writeto( path )

    pass

def weights(err):
    """ generate observation weights from uncertainties """
    w = np.power(err, -2)
    w[~np.isfinite(w)] = 0
    return w/np.sum(w)

def combine_orders(spec,lbary_ltopo=None):
    # get wavelength range
    w1 = spec[0,:,:].ravel().min()
    w2 = spec[0,:,:].ravel().max()
    dw = np.abs(np.diff(spec[0,:,:].ravel())).mean()

    # define combined wavelegth grid
    wgrid = np.arange(w1,w2,dw)

    # interpolate orders to new wavelength grid
    # and cut out orders' edges
    newflux    = np.empty((spec.shape[1],wgrid.shape[0]))
    newfluxerr = np.empty((spec.shape[1],wgrid.shape[0]))
    for ii in range(spec.shape[1]):
        try:
            overlap1 = (spec[0,ii+1,:].max() - spec[0,ii,:].min())/2
        except IndexError:
            overlap1 = 0
        try:
            if ii<=0: raise IndexError
            overlap2 = (spec[0,ii,:].max()-spec[0,ii-1,:].min())/2
        except IndexError:
            overlap2 = 0
        if overlap1 < 0: overlap1 = 0
        if overlap2 < 0: overlap2 = 0

        um = np.where( (spec[0,ii,:] >= spec[0,ii,:].min()+overlap1) & (spec[0,ii,:] <= spec[0,ii,:].max()-overlap2) )[0]
        f  = scipy.interpolate.interp1d(spec[0,ii,um],spec[5,ii,um],bounds_error=False)
        fs = scipy.interpolate.interp1d(spec[0,ii,um],weights(spec[2,ii,um]),bounds_error=False)

        newflux[ii,:]    = f(wgrid)
        newfluxerr[ii,:] = fs(wgrid)

    # weighted average orders
    combined = np.nansum( newflux*newfluxerr,axis=0 )/np.nansum(newfluxerr,axis=0)
    errors = np.sqrt(1. / np.nansum(1./newfluxerr**2., axis=0) )

    # convert wavelengths to air
    wgrid = GLOBALutils.ToAir( wgrid )
    # apply barycentric correction
    if lbary_ltopo is not None:
        wgrid = wgrid*lbary_ltopo

    return wgrid,combined,errors

def get_CCF_bkg(t,y,order):
    z = np.polyfit( t, y, order )

    residualtest = y- np.poly1d(z)(t)
    npoints = 0
    while True:
        umtest = residualtest>-0.5*np.std(residualtest)
        try:
            z = np.polyfit( t[umtest],  y[umtest], order )
        except TypeError:
            # Probably too low number of points
            z = np.polyfit( t,   residualtest, order )
            break
        residualtest = y-np.poly1d(z)(t)
        if npoints==len(umtest):
            # Break if no more pts below 0.5 std
            break
        npoints=len(umtest)
    del residualtest
    del npoints

    return np.poly1d(z)(t)

def save_combined_spectrum(path,wcombined,fcombined,errors,obname,mjd,mbjd,bjd,bcvel_baryc,hd,ra,dec,SNR_5130,SNR_5130_R):
    goodpts = np.isfinite(fcombined)
    data = np.array([wcombined[goodpts],fcombined[goodpts],errors[goodpts]])
    # Create header
    hdu = pyfits.PrimaryHDU( data )
    hdu = GLOBALutils.update_header(hdu,'TTYPE1', 'WAV')
    hdu = GLOBALutils.update_header(hdu,'TTYPE2', 'FLUX')
    hdu = GLOBALutils.update_header(hdu,'TTYPE3', 'FLUX_E')
    hdu = GLOBALutils.update_header(hdu,'TUNIT1', 'AA')
    hdu = GLOBALutils.update_header(hdu,'TUNIT2', '')
    hdu = GLOBALutils.update_header(hdu,'TUNIT3', '')
    hdu = GLOBALutils.update_header(hdu,'TARGET', obname)
    hdu = GLOBALutils.update_header(hdu,'MJD', mjd)
    hdu = GLOBALutils.update_header(hdu,'MBJD', mbjd)
    hdu = GLOBALutils.update_header(hdu,'BJD', bjd)
    try:
        hdu = GLOBALutils.update_header(hdu,'RA',hd['RA'])
        hdu = GLOBALutils.update_header(hdu,'DEC',hd['DEC'])
        hdu = GLOBALutils.update_header(hdu,'RA BARY',ra)
        hdu = GLOBALutils.update_header(hdu,'DEC BARY',dec)
    except:
        None
    hdu = GLOBALutils.update_header(hdu,'EQUINOX',hd['EPOCH'])
    if bcvel_baryc is not None:
        hdu = GLOBALutils.update_header(hdu,'HIERARCH BARYCENTRIC CORRECTION (KM/S)', bcvel_baryc)
        hdu = GLOBALutils.update_header(hdu,'SNR', SNR_5130)
        hdu = GLOBALutils.update_header(hdu,'SNR_R', SNR_5130_R)
    hdu = GLOBALutils.update_header(hdu,'INST', 'RCC1m')
    hdu = GLOBALutils.update_header(hdu,'RESOL', '20000')
    hdu = GLOBALutils.update_header(hdu,'PIPELINE', 'CERES')

    if os.access(path,os.F_OK):
        os.remove( path )
    hdu.writeto( path )

    pass

def XCfit_two_gaussians(ccfxvals,ccfyvals,disp,debug=False):
    from astropy.modeling import models, fitting

    disp *= 1.5

    mean_init = ccfxvals[ np.argmax(1-ccfyvals) ]
    amp_init = np.ptp(ccfyvals)
    g_init = models.Gaussian1D(amplitude=amp_init, mean=mean_init, stddev=disp)

    fit_g = fitting.LevMarLSQFitter()
    g_init.amplitude.fixed = True
    g_init.mean.fixed = True
    gm = fit_g(g_init, ccfxvals,1-ccfyvals)

    g_init = models.Gaussian1D(amplitude=amp_init, mean=mean_init, stddev=gm.stddev.value)
    residual = (1-ccfyvals-g_init(ccfxvals))

    if debug:
        plt.figure(figsize=(10,5))
        plt.title('First Gaussian')
        plt.plot(ccfxvals,ccfyvals,'k',lw=2)
        plt.plot(ccfxvals,1-g_init(ccfxvals),'r--',lw=3)
        plt.xlabel("Velocity (km/s)")
        plt.ylabel('XC = 1 - CCF')
        plt.show()

    mean_init2 = ccfxvals[ np.argmax(residual) ]
    amp_init2 = np.ptp(residual)
    g_init =  models.Gaussian1D(amplitude=amp_init, mean=mean_init, stddev=disp) + models.Gaussian1D(amplitude=amp_init2, mean=mean_init2, stddev=disp) + models.Linear1D(slope=-1e-5)

    if debug:
        plt.figure(figsize=(10,5))
        plt.title('Residual and second Gaussian')
        plt.plot(ccfxvals,residual,'k',lw=2)
        plt.plot(ccfxvals,models.Gaussian1D(amplitude=amp_init2, mean=mean_init2, stddev=disp)(ccfxvals),'r--',lw=3)
        plt.xlabel("Velocity (km/s)")
        plt.ylabel('Residual XC')
        plt.show()

    gm = fit_g(g_init, ccfxvals,1-ccfyvals)

    if debug:
        plt.figure(figsize=(10,5))
        plt.title('Fitted two Gaussians + linear')
        plt.plot(ccfxvals,ccfyvals,'k',lw=2)
        plt.plot(ccfxvals,1-gm(ccfxvals),'r',lw=1)
        plt.xlabel("Velocity (km/s)")
        plt.ylabel('XC = 1 - CCF')
        plt.show()

    return gm

def check_ephem_file(ephempath):
    if not os.access(ephempath,os.F_OK):
        from skyfield.iokit import download

        print('Downloading JPL ephemeris...')

        url = 'ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/de421.bsp'
        download(url,ephempath,verbose=True)
