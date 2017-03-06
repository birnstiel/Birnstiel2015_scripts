import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

from IPython.display import HTML

pi           = np.pi;                   # PI
c_light      = sc.c*1e2;                # speed of light in cm/s
arcsec       = (pi/180./3600.)	        # 1 arcsecond in radian
arcsec_sq    = (pi/180./3600.)**2	# 1 square arcsecond in sterad
AU           = sc.au*1e2;               # astronomical unit in cm
k_b          = sc.k*1e7;                # Boltzmann constant in erg/K
mu           = 2.3e0;                   # mean molecular mass in proton masses
m_p          = sc.proton_mass*1e3;      # proton mass in g
Grav         = sc.G*1e3;                # gravitational constant in cm^3 g^-1 s^-2
year         = sc.Julian_year;          # year in s
sig_h2       = 2e-15;                   # cross section of H2 [cm^2]
PC           = sc.parsec*1e2;           # parsec in cm
M_sun        = 1.9891e+33;              # mass of the sun in g

def dlydlx(x,R):
    """
    calculates the log-derivative
    
     dlog(y)
    -------- = dlydlx(x,y)
     dlog(x)
    """
    from numpy import zeros,shape,interp,log10 
    #
    # define the interpolation function (one for each row)
    #
    r = zeros(shape(R))
    if len(shape(R))>1:
        for i,row in enumerate(R):
            R_int = lambda x_int: 10**(interp(log10(x_int),log10(x),log10(row)))
            h = x/100.
            r[i] = x/row*(R_int(x+h)-R_int(x-h))/(2.*h)
    else:
        R_int = lambda x_int: 10**(interp(log10(x_int),log10(x),log10(R)))
        h = x/100.
        r = x/R*(R_int(x+h)-R_int(x-h))/(2.*h)
    return r

def show_pdf(filename,width=0.5,aspect=None):
    """
    Displays the specified pdf file in a ipython/jupiter notebook.
    The width is given as screen width, the height is given via the
    aspect ratio.
    
    Arguments:
    ----------
    
    filename : string
        The path of the pdf to be shown
       
    Keywords:
    ---------
    
    width : float
        The width where 1 = full width
        
    aspect : float
        The aspect ratio width/height. Defaults to last figure's
        aspect ratio or to 4./3. if no figure present.
        
    Returns:
    --------
    A HTML object
    """
    if aspect is None:
        if plt.get_fignums()==[]:
            aspect = aspect or 4./3.
        else:
            aspect = plt.gcf().get_size_inches()
            aspect = aspect[0]/aspect[1]
    return HTML('<div style="position:relative;width:{:g}%;height:0;padding-bottom:{:g}%">'.format(width*100,width*100/aspect+2)+\
         '<iframe src="'+filename+'" style="width:100%;height:100%"></iframe></div>')

def get_image_data(filename):
    """
    Loads the specified fits image and returns x, y in arcsec
    and the intensity in Jy*arcsec^-2, as well as wavelength and
    the fits header.
    
    Arguments:
    ----------
    
    filename : string
    :   the path to the fits image
    
    Output:
    -------
    
    x,y,img,img_lam,h
    
    x : array
    :   x coordinate [arcsec]
    
    x : array
    :   y coordinate [arcsec]
    
    img : array
    :   the 2D image data [Jy/arcsec^2]
    
    img_lam : float
    :   wavelengh [cm]
    
    h : fits header
    :   the header information of the fits file
    """
    #
    # open fits file and define x and y in arcsec
    #
    from astropy.io import fits
    f = fits.open(filename)
    h = f[0].header
    if h['CUNIT1']!='deg' or h['CUNIT2']!='deg' \
    or h['NAXIS1']!=h['NAXIS2'] or h['BUNIT']!='JY/PIXEL':                                                            
        raise NameError('Something wrong with the image, check units & shape!')
    x = (np.arange(h['NAXIS1'])-h['CRPIX1'])*h['CDELT1']*pi/180./arcsec
    y = (np.arange(h['NAXIS2'])-h['CRPIX2'])*h['CDELT2']*pi/180./arcsec
    #
    # get image data (in Jy/pix) and convert to Jy/arcsec^2
    #
    img     = f[0].data.copy()/(h['CDELT1']*h['CDELT2'])*(180./pi)**2*arcsec_sq
    img_lam = c_light/h['RESTFREQ']
    #
    # close fits file
    #
    f.close()
    return x,y,img,img_lam,h

def makeimage(img_lam,dirname='./',incl=0.,PA=0.,npix=512,sizeau=1000.,\
              img_name=None,dpc=10.0,**kwargs):
    """
    Make an image at the specified wavelength
    
    Arguments:
    ----------
    
    img_lam : float
        Wavelength of the image [cm]
        
    Keywords:
    ---------
    
    dirname : string
        Directory where to call RADMC3D
    
    incl : float
        Inclination of the source
        
    PA : float
        Position angle of the source
        
    npix : float
        Pixel size of the image
        
    sizeau : float
        Size of the image in AU
    
    img_name : str
        Name under which to save the image (gets deleted otherwise)
        
    dpc : float
        Distance of the source in PC
        
    Other keywords should be passed as nphot='100000' for example and will be
    appended to the radmc call
        
    Output:
    -------
    returns result from readimage: image is in erg/(s cm^2 Hz ster)
    """
    import os,subprocess,shutil
    delete_image = False
    if img_name is None:
        delete_image = True
        img_name = 'temp_image'
    if os.path.exists(dirname+os.sep+img_name):
        os.unlink(dirname+os.sep+img_name)
    
    incl_str = '{:.0f}'.format(round(incl))
    PA_str   = '{:.0f}'.format(round(PA-90.))
    npix_str = '{:.0f}'.format(round(npix))
    szau_str = '{:.0f}'.format(round(sizeau))
    wl_str   = '{:.4f}'.format(img_lam*1e4)
    dpc_str  = '{:.4f}'.format(dpc)
    #
    # create image at that wavelength
    #
    params = [item for kv in kwargs.iteritems() for item in kv]
    subprocess.call(['nice','radmc3d','image','lambda',wl_str,\
                     'incl',incl_str,'posang',PA_str,'npix',npix_str,\
                     'sizeau',szau_str,'dpc',dpc_str]+params,cwd=dirname)
    if os.path.exists(dirname+os.sep+'image.fits'):
        os.unlink(dirname+os.sep+'image.fits')
    radmcimage_to_fits(dirname+os.sep+'image.out',\
                       dirname+os.sep+'image.fits',dpc)
    #
    # read in image
    #
    im=readimage(filename=dirname+os.sep+'image.out')
    #
    # delete if necessary
    #
    os.unlink(dirname+os.sep+'image.out')
    if delete_image:
        os.unlink(dirname+os.sep+'image.fits')
    else:
        shutil.move(dirname+os.sep+'image.fits',\
                    dirname+os.sep+img_name+'.fits')
    #
    # return result
    #
    return im
    
def radmcimage_to_fits(imagename,fitsname,dpc,arcsec=None,mas=None):
    """
    CONVERT RARMC-3D IMAGE TO FITS
    
    RADMC-3D Images are text files (the standard is image.out). This routine
    converts such a file to the FITS standard.
    
    Arguments:
    ----------
    
    imagename : string
        Name of the image file that RADMC-3D produces. The standard file name
        RADMC-3D produces is: 'image.out'.
                
    fitsname : string
        Name of the output fits file you want to produce.
        
    dpc : float
        Distance of observer to object in units of parsec.
    
    Keywords:
    ---------
    
    arcsec : float
        Passed as keyword with the same name to `wirtefitsimage`

    mas : float
        Passed as keyword with the same name to `wirtefitsimage`
    """
    im  = readimage(filename=imagename)
    pixdeg_x = 180.0/pi * (im['sizepix_x']/(dpc*PC))
    pixdeg_y = 180.0/pi * (im['sizepix_y']/(dpc*PC))
    #
    # XXX NOTE XXX
    # there seems to be some difference between idl and python which causes
    # the idl image to be the transpose of the python image. Should find out
    # where this is happening 
    #
    writefitsimage(im['image'],fitsname,1e4*c_light/im['lamb'],\
                   pixdeg_x,pixdeg_y,arcsec=arcsec,mas=mas)

def writefitsimage(image,filename,freqhz,pixdeg_x,pixdeg_y,\
                   radeg=None,decdeg=None,arcsec=False,mas=False):
    """
    ROUTINE FOR WRITING FITS IMAGE
    
    Arguments:
    ----------
    
    image : array 
        Image data
        
    filename : string
        Name of the output file
    
    freqhz : float
        The frequency of at which the image was taken/calculated
        
    pixdeg_x,pixdeg_y : float
        The extend of the image on the sky in degree
        
    radeg : float
        RA of the image on the sky
    
    decdeg : float
        DEC of the image on the sky
        
    arcsec : bool
        use units of arcsecs for the image header instead of the default degree
    
    mas : bool
        use units of milli arc seconds for the image header
    
    Output:
    -------
    
    Writes the data to the fits file specified by `filename`, converted to
    Jansky / pixel
    """
    from numpy import pi,squeeze,shape
    from astropy.io import fits
    #
    # check
    #
    if mas and arcsec:
        raise NameError('Seti EITHER mas OR arcsec to true, not both!')
    #
    # get image size
    #
    nx,ny = shape(squeeze(image))
    #
    # Compute the conversion factor from erg/cm^2/s/Hz/ster to Jy/pixel
    #
    pixsurf_ster = pixdeg_x*pixdeg_y/((180/pi)**2)
    factor = 1e+23 * pixsurf_ster
    #
    # Make FITS header information, reverse order
    #
    header    = fits.Header()
    #
    # ...Rest frequency
    #
    header['RESTFREQ'] = freqhz
    #
    # ...Zero point of coordinate system
    #
    header['CRPIX2']  =  (ny+1)/2
    #
    # ...Pixel scale
    #
    if arcsec  is not None:
        header['CUNIT2'] = 'arcsec'
        header['CDELT2'] = 3.6e3*pixdeg_y
    else:
        if mas  is not None:
            header['CUNIT2'] = 'mas'
            header['CDELT2'] = 3.6e6*pixdeg_y
        else:
            header['CUNIT2'] = 'deg'
            header['CDELT2'] = pixdeg_y
    #
    # ...CRVAL2: value of y-axis at critical pixel
    #
    if decdeg  is not None:
        header['CRVAL2'] = decdeg
        header['CTYPE2'] = 'DEC--SIN'
    #
    # ...Zero point of coordinate system
    #
    header['CRPIX1'] = (nx+1)/2
    #
    # ...Pixel scale
    #
    if arcsec  is not None:
        header['CUNIT1'] = 'arcsec  '
        header['CDELT1'] = 3.6e3*pixdeg_x
    else:
        if mas  is not None:
            header['CUNIT1'] = 'mas     '
            header['CDELT1'] = 3.6e6*pixdeg_x
        else:
            header['CUNIT1'] = 'deg     '
            header['CDELT1'] = pixdeg_x
    #
    # ...CRVAL1: value of x-axis at critical pixel
    #
    if radeg  is not None:
        header['CRVAL1'] = radeg
        header['CTYPE1'] = 'RA---SIN'
        header['LONPOLE'] = 1.8e+02
        header['EPOCH'] = 2e+03
    #
    # ...Unit of intensity
    #
    header['BUNIT'] = 'JY/PIXEL'
    #
    # ...BZERO
    #
    header['BZERO'] = 0e0
    #
    # ...BSCALE
    #
    header['BSCALE'] = 1e0
    #
    # ...Type of data
    #
    header['BTYPE'] = 'Intensity'
    #
    # Make a FITS file
    #
    imjypix  = factor * image
    hdu      = fits.PrimaryHDU(imjypix,header=header)
    thdulist = fits.HDUList([hdu])
    thdulist.writeto(filename)

def readimage(ext=None,filename=None):
    """
    Reads the rectangular telescope image produced by RADMC3D. The file name of
    the image is assumed to be image.out if no keyword is given. If keyword
    `ext` is given, the filename  'image_'+ext+'.out' is used. If keyword
    `filename` is given, it is used as the file name.
    
    Keywords:
    ---------
    
    ext : string
    :   Filename extension of the image file, see above
        
    filename : string
    :   file name of the image file 
        
    Output:
    -------
    
    Returns a dictionary containing the image data with the following entries:
    nx,ny,nrfr,sizepix_x,sizepix_y,image,flux,x,y,lamb,radian,stokes
    
    The image units are erg/(s cm^2 Hz ster)

    """
    from numpy import fromfile,product,arange
    import glob
    #
    # Read from normal file, so make filename
    #
    if filename is None:
        if ext is None:
            filename = 'image.out'
        else:
            filename = 'image_'+str(ext)+'.out'
    fstr = glob.glob(filename)
    if len(fstr) == 0:
        print('Sorry, cannot find '+filename)
        print('Presumably radmc3d exited without succes.')
        print('See above for possible error messages of radmc3d!')
        raise NameError('File not found')
    funit = open(filename)
    #
    # Read the image
    #
    iformat = fromfile(funit,dtype='int',count=1,sep=' ')[0]
    if iformat < 1 or iformat > 4:
        raise NameError('ERROR: File format of '+filename+' not recognized.')
    if iformat == 1 or iformat == 3:
        radian = False
    else:
        radian = True
    if iformat == 1 or iformat == 2:
        stokes = False
    else:
        stokes = True
        
    nx,ny               = fromfile(funit,dtype=int,count=2,sep=' ')
    nf                  = fromfile(funit,dtype=int,count=1,sep=' ')[0]
    sizepix_x,sizepix_y = fromfile(funit,dtype=float,count=2,sep=' ')
    lamb                = fromfile(funit,dtype=float,count=nf,sep=' ')
    if nf==1:
        lamb = lamb[0]
    if stokes:
        image_shape = [4,nx,ny,nf]
    else:
        image_shape = [nx,ny,nf]
    image = fromfile(funit,dtype=float,count=product(image_shape),\
                     sep=' ').reshape(image_shape,order='F')
    funit.close()
    #
    # If the image contains all four Stokes vector components,
    # then it is useful to transpose the image array such that
    # the Stokes index is the third index, so that the first
    # two indices remain x and y
    #
    if stokes:
        if nf > 1:
            image = image[[1,2,0,3]]
        else:
            image = image[[1,2,0]]
    #
    # Compute the flux in this image as seen at 1 pc
    #
    flux=0.0
    if stokes:
        for ix in arange(nx):
            for iy in arange(ny):
                flux=flux+image[ix,iy,0,:]
    else:
        for ix in arange(nx):
            for iy in arange(ny):
                flux=flux+image[ix,iy,:]
    flux=flux*sizepix_x*sizepix_y
    if not radian: flux=flux/PC**2
    #
    # ADDED 13.12.06:
    # Compute the x- and y- coordinates
    #
    x=((arange(nx)+0.5)/(nx*1.)-0.5)*sizepix_x*nx
    y=((arange(ny)+0.5)/(ny*1.)-0.5)*sizepix_y*ny
    #
    # Return all
    #
    return {'nx':nx,'ny':ny,'nrfr':nf,'sizepix_x':sizepix_x,\
            'sizepix_y':sizepix_y,'image':image.squeeze(),'flux':flux,\
            'x':x,'y':y,'lamb':lamb,'radian':radian,'stokes':stokes}            
    
def better_plots(back='w',front='k',fs=15,cmap='spectral',lw=1,sans=False,\
                 usetex=False,brewer=True):
    """
    Changes matplotlib default parameters to get an improved look.
    
    Keywords:
    ---------
    
    back : [*'k'* | color]
        the background color of the axes and the figure
    
    front : [*'w'* | color]
        the foreground color of the axes, lines, font, ...
    
    fs : [*12* | int]
        the default font size 

    lw : float
        line width modifier (1 is already somewhat thicker than default) 

    brewer : bool
        wether to use normal colors or brewer colors 
    
    cmap : [*'spectral'* | colormap]
        the default colormap to be used

    sans : [*False* | True]
        whether to use sans-serif font family or not
    
    usetex : [*False* | True]
        whether tex-fonts are used by default. This also sets the fonts to have
        a consistent look, but makes plotting quite a lot slower
    
    Example:
    
        for usetex in [False,True]:
            for sans in [False,True]:
                better_plots(usetex=usetex,sans=sans)
                figure()
                title('usetex = {:}, sans = {:}'.format(usetex,sans))
                plot(sin(linspace(0,2*pi,100)))
                xlabel(r'$x$ [m]')
                ylabel(r'$\alpha$ [$\mu$m]')
                draw()
                pause(5)
    """
    from matplotlib.pyplot import rcParams, rcParamsDefault, rc
    import brewer2mpl
    dark2_colors = brewer2mpl.get_map('Dark2', 'Qualitative', 8).mpl_colors
    
    rcParams['figure.facecolor']            = back
    rcParams['axes.edgecolor']              = front
    rcParams['axes.facecolor']              = back
    rcParams['axes.linewidth']              = 1.5*lw
    rcParams['axes.labelcolor']             = front
    rcParams['axes.color_cycle']            = [front, 'g', 'r', 'c', 'm', 'y']*(not brewer)+brewer*dark2_colors
    rcParams['axes.formatter.limits']       = [-10000,10000]
    rcParams['axes.formatter.use_mathtext'] = True
    rcParams['xtick.color']                 = front
    rcParams['ytick.color']                 = front
    rcParams['xtick.major.size']            = 6*lw
    rcParams['ytick.major.size']            = 6*lw
    rcParams['ytick.major.width']           = 1*lw
    rcParams['xtick.major.width']           = 1*lw
    rcParams['xtick.minor.size']            = 3*lw
    rcParams['ytick.minor.size']            = 3*lw
    rcParams['ytick.minor.width']           = 0.75*lw
    rcParams['xtick.minor.width']           = 0.75*lw
    rcParams['lines.linewidth']             = 1.5*lw
    rcParams['image.cmap']                  = cmap
    rcParams['font.size']                   = fs
    rcParams['text.color']                  = front
    rcParams['savefig.facecolor']           = back
    #
    # avoid tick labels overlapping with the axes for large fonts
    #
    if fs > 16:
        rcParams['xtick.major.pad']='6'
        rcParams['ytick.major.pad']='6'
    #
    # make tex and non-tex text look similar
    #
    if sans > 0:
        if sans==1:
            rcParams['mathtext.fontset']='stixsans'
        elif sans==2:
            rcParams['mathtext.fontset'] = 'custom'
            rcParams['font.family'] = 'sans-serif'
            rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
            rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
            rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
    else:
        rcParams['mathtext.fontset'] = 'stix'
        rcParams['font.family']      = 'STIXGeneral'
        rcParams['mathtext.rm']      = rcParamsDefault['mathtext.rm']
        rcParams['mathtext.it']      = rcParamsDefault['mathtext.it'] 
        rcParams['mathtext.bf']      = rcParamsDefault['mathtext.bf'] 
    #
    # the tex settings
    #
    rc('text', usetex=usetex)
    if usetex:
        rcParams['text.latex.preamble']         = [r"\usepackage{amsmath}"]
        if sans:
            rcParams['text.latex.preamble'] += [r"\usepackage{cmbright}"]
            rc('font',**{'family':'sans-serif',\
                         'sans-serif':['Bitstream Vera Sans']})
        else:
            rc('font',**{'family':'serif',\
                         'serif':'Computer Modern Roman',\
                         'sans-serif':'Computer Modern Sans serif',\
                         'monospace':'Computer Modern Typewriter'})
    else:
        rcParams['text.latex.preamble']         = ['']
        
def trace_line_though_grid(xi,yi,f,x=None,y=None):
    """
    Returns the cell indices through which the curve moves

    """
    if x is None: x=0.5*(xi[1:]+xi[:-1])
    if y is None: y=0.5*(yi[1:]+yi[:-1])
    
    def fill_cells(x,y,yi,ixstart,iystart,yend):
        """
        Takes an index, returns, all cells from (including) this index
        up until the last cell that includes the value yend (towards yend)
        
        """
        iyend     = np.searchsorted(yi,yend)-1
        direction = int(np.sign(yend-y[iystart]))
        return [(ixstart,iy) for iy in range(iystart,min(max(0,iyend),len(y)-1)+direction,direction)]
    #
    # begin function
    #
    fx = f(x)
    result = set()
    #
    # find first cell center where the function value is on the grid
    #
    mask = np.where((fx<=yi[-1]) & (fx>=yi[0]))[0]
    if len(mask) == 0: return result
    ix0 = mask[0]
    y_interface = f(xi[ix0])
    
    if y_interface>yi[-1]:
        iy0    = len(y)-1
        dum    = fill_cells(x,y,yi,ix0,iy0,y_interface)
        result = result.union(dum)
    elif y_interface>yi[-2]:
        ix0    = max(0,ix0-1)
        iy0    = len(y)-1
        dum    = fill_cells(x,y,yi,ix0,iy0,y_interface)
        result = result.union(dum)
        ix0   += 1
    elif y_interface>yi[0]:
        ix0 = max(0,ix0-1)
        iy0 = np.searchsorted(yi,y_interface)-1
        dum    = fill_cells(x,y,yi,ix0,iy0,y_interface)
        result = result.union(dum)
        ix0   += 1
    else:
        iy0 = 0
        dum    = fill_cells(x,y,yi,ix0,iy0,y_interface)
        result = result.union(dum)
        
    iy0    = dum[-1][-1]
    dum    = fill_cells(x,y,yi,ix0,iy0,fx[ix0])
    result = result.union(dum)
    iy0    = dum[-1][-1]
    
    while ix0<=mask[-1]:
        #
        # evaluate function at next interface
        #
        y_interface = f(xi[ix0+1])
        #
        # fill until interface value is reached
        #
        dum    = fill_cells(x,y,yi,ix0,iy0,y_interface)
        result = result.union(set(dum))
        #
        # go right
        #
        ix0 += 1  # update 1
        iy0  = dum[-1][-1]
        if ix0>mask[-1]: break
        #
        # fill until final value is reached
        #
        dum    = fill_cells(x,y,yi,ix0,iy0,fx[ix0])
        result = result.union(set(dum))
        iy0    = dum[-1][-1] 
    #
    # transform into sorted list
    #
    result=[list(i) for i in list(result)]
    result.sort()
    return result