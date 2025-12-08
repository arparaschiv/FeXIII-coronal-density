import numpy as np                       ## in DKIST_DENS, COMP_DENS
import multiprocessing                   ## in DKIST_DENS, COMP_DENS
from tqdm import tqdm                    ## in DKIST_DENS, COMP_DENS
from scipy.io import readsav             ## in DKIST_DENS, COMP_DENS ; Only if using the deprecated sswidl table
from datetime import datetime,timedelta  ## in COMP_DENS
from scipy import interpolate            ## in work_1pix 

## Automatic parallel computing of raster slit scan timeseries densities of DKIST Cryo-NIRSP or DL-NIRSP, or uCOMP observation.
## Input data are matching pairs of Fe XIII 1074.6nm and 1079.8nm reduced observations of >= Level-1.

## The data can be single slit positions or raster scans of observations (single-dimension to three-dimensional arrays).

# The **DKIST_DENS** and *COMP_DENS* functions need a CHIANTI derived calculation help file. It is included in this repository.
# The function can process peak line emission or integrated line observations, provided the required parameters for each case are provided.

# The **work_1pix_dkist** and *work_1pix_ucomp* worker function is required for the parallel implementation.

# a multitemp implementation for DKIST_dens is prototyped at the end of this script. Temperature effects were found to exist, but  too small to be significant.
# These scripts are provided "as-is", are not free of bugs, and nor safe from improper usage. Please mind inherent assumptions.

# More information on the method and software codes can be found in [Paraschiv & Judge, 2022](https://ui.adsabs.harvard.edu/abs/2022SoPh..297...63P/abstract) and [Schad & Dima, 2020](https://ui.adsabs.harvard.edu/abs/2020SoPh..295...98S/abstract).

# Contact: Alin Paraschiv, NSO ---- arparaschiv at nso edu


##########################################################################
####################### DKIST_DENS implementation ########################
##########################################################################


def DKIST_DENS(i1074,i1079,chianti_link,header_1074,header_1079,hpxy=None,integrated=False):
## compute the density from two stokes I observations (serialized parallel runs for 1 pix at a time) of Fe XIII assuming:
##    - a line ratio look-up table created using chianti atomic data. 
##    - accurate header pointing information. shape [number_of_steps]
##    - array of pointing information (hpxy) computed using the Cryo-NIRSP public community notebooks https://bitbucket.org/dkist-community-code/cryonirsp-notebooks/src/main/first_release_specFitting/ (optional). shape [2,number_of_steps,slit_size]
##    - that the arrays are aligned properly. shape [number_of_steps,slit_size]
##    - others ??

## function should process single-dimensional, or bi-dimensional input arrays of 1 OR 3 emission line fit parameters. Input array dimension should be between 1 and 3.
##       The input arrays contain EITHER peak intensity, OR peak intensity, wavelength center position, and distribution sigma if integrated=True is passed; e.g. 
##       [slitsize]                               --single slit position, only peak intensity
##       [slitsize,fit_params]                    --single slit position, peak intensity, center, and sigma.
##       [number_of_steps,slit_size]              --raster of slit positions, only peak intensity
##       [number_of_steps,slit_size,fit_params]   --raster of slit positions,  peak intensity, center, and sigma.

## - Two chianti look-up tables are computed via PyCELP or SSWIDL(deprecated) and imported via the "chianti_link" path given to this function.
## - Requires full observation headers (1 header per slit step) from both Fe XIII line datasets.
## - Requires the numpy, multiprocessing, tqdm, scipy's interpolate and io.readsav imports 
    
## input ingestion and return array. Sanity checks will return 0 in case something is incompatible.
    
    ## single slit position with only peak intensity.
    if np.ndim(i1074) <= 1:                              ## allows for 1 value calculations
        density=np.zeros(i1074.shape,dtype=np.float32)   ## single-dimensional spatial array of peak intensities
        if integrated:
            print("Not enough parameter inputs for integrated line calculations. Aborting!")
            return density    
        if (i1074.shape != i1079.shape):   ## sanity check 
            print("Arrays are not of the same physical shape. Aborting!")
            print("Shapes are: Array1074=", str(i1074.shape)," Array1079=", str(i1079.shape))
            return density
        print('Single slit position of peak intensity ingested. Continuing!')        
    
    ## two sub-combinations of inputs.
    elif np.ndim(i1074) == 2:
        if integrated:                                                                ## understand a single slit with 3 fit parameters
            density=np.zeros(i1074.shape[0],dtype=np.float32)   
            if (i1074.shape[0] != i1079.shape[0]):   ## sanity check 
                print("Arrays are not of the same physical shape. Aborting!")
                print("Shapes are: Array1074=", str(i1074.shape)," Array1079=", str(i1079.shape))
                return density
            print('Single slit position with peak, center and standard deviation parameters are ingested. Continuing!')
        else:                                                                         ## understand a raster of slits with only peak intensity.
            density=np.zeros(i1074.shape,dtype=np.float32)   
            if (i1074.shape != i1079.shape):   ## sanity check 
                print("Arrays are not of the same physical shape. Aborting!")
                print("Shapes are: Array1074=", str(i1074.shape)," Array1079=", str(i1079.shape))
                return density
            print('Raster of slit positions of only peak intensity ingested. Continuing!')
    
    ## bi-dimensional spatial array with added 3rd dimension for peak, center an standard deviation of coronal line.
    ## not capturing if more than three fit parameters are ingested. 
    elif np.ndim(i1074) == 3:      
        density=np.zeros(i1074.shape[:2],dtype=np.float32)
        if (i1074.shape[:2] != i1079.shape[:2]):   ## sanity check 
            print("Arrays are not of the same physical shape. Aborting!")
            print("Shapes are: Array1074=", str(i1074.shape)," Array1079=", str(i1079.shape))
            return density
        if i1074.shape[2] < 3:
            print("Not enough parameter inputs for integrated line calculations. Aborting!")
            return density
        print('Raster of slit positions with peak, center and standard deviation parameters are ingested. Continuing!')
        
    ## Input not compatible with this function
    else:
        print("input array dimensions not understood. Aborting!")
        return np.zeros((2,2))


    ## theoretical chianti ratio calculations created via PyCELP or SSWIDL        
    ## read the chianti table         
    if (chianti_link[-4:]    == ".npz"):                     ## default
        chianti_table        =  dict(np.load(chianti_link))  ## variables (h,den,rat) directly readable by work_1pix_dkist .files required for loading the data directly.
    elif (chianti_link[-4:]  == ".sav"):                     ## deprecated
        print("The SSWIDL look-up table is deprecated. The PyCELP look-up table is better supported and recommended at this stage.")
        chianti_table        =  readsav(chianti_link)        ## variables  (h,den,rat) directly readable by work_1pix_dkist
        chianti_table["den"] =  chianti_table["den"][0,:]    ## density sample values are the same at all heights. Use just one instance.
    else:
        print("CHIANTI look-up table not found. Is the path correct?")
        return density                                ## a zero array at this point
    #print(chianti_table.keys())               ## Debug - check the arrays
    
    ##  in the save files:
    ## 'h' is an array of heights in solar radii ranged 1.01 - 2.00(pycelp) or 1.01 - 1.50(sswidl);  
    ## 'den' is an array of density (ranged 6.00 to 12.00 (pycelp) or 5.0 to 13.0 (sswidl); the broader interval is not really recommended here. 
    ## 'rat' is a 2D array containing line ratios corresponding to the density range values at each distinct height.
    
    ##  To query the look-up table:
    ##  print(chia['h'],chia['h'].shape,chia['rat'].shape,chia['den'].shape)

    ## look-up table resolutions
    ## pycelp: h is of shape [99]; den and rat are arrays of shape [99, 120] corresponding to the 99 h height values and 120 density values
    ## sswidl: h is of shape [50]; den and rat are arrays of shape [50,  33] corresponding to the 50 h height values and  33 density values

## set up the cpu worker and argument arrays
    p          = multiprocessing.Pool(processes=multiprocessing.cpu_count()-2,maxtasksperchild = 10000)     ##ncpus = 54  ## dynamically defined from system query as total CPU core number - 2
    ## argument index keeper for splitting tasks to cpu cores
    ## two branches to separate slingle slits vs rasters
    arg_array  = []

    ## check if hpxy exists:
    if hpxy is None: hpxy=np.zeros((2,i1074.shape[0],i1074.shape[1]))   # dummy to skip hpxy-based calculations

    
    if np.ndim(i1074) <= 1 or (np.ndim(i1074) == 2 and integrated):                                        ## Single slit/value branch
        for ypos in range(i1074.shape[0]): 
            arg_array.append((0,ypos,i1074[ypos],i1079[ypos],header_1074[0],header_1079[0],chianti_table,integrated))        
        rs     = p.starmap(work_1pix_dkist,arg_array)
        p.close()
        pbar   = tqdm(total=len(arg_array))  ## progress bar implmentation via tqdm
        for i,res in enumerate(rs):
            step,ypos,density[ypos] = res    ## step is just a dummy in this case; Density defined as a single-dimension vector in this case.
            pbar.update()
        pbar.close()
   
    else:                                                                                                  ## Raster slit branch
        for step in range(i1074.shape[0]): 
            for ypos in range(i1074.shape[1]): 
                arg_array.append((step,ypos,i1074[step,ypos],i1079[step,ypos],header_1074[step],header_1079[step],hpxy[0,step,ypos],hpxy[1,step,ypos],chianti_table,integrated))

        rs     = p.starmap(work_1pix_dkist,arg_array)
        p.close()
        pbar   = tqdm(total=len(arg_array))  ## progress bar implmentation via tqdm
        for i,res in enumerate(rs):
            step,ypos,density[step,ypos] = res 
            pbar.update()
        pbar.close()
    
    density[density <= 0] = np.nan

    return density

def work_1pix_dkist(xx,yy,a_obs,b_obs,hdr_a_1pix,hdr_b_1pix,hpxy_a,hpxy_b,chianti_table,integrated):
    ## compute the ratio of the observation ## 1074/1079 fraction, same as rat component of the chianti_table 
    
    ## Sanity checks at pixel level
    if integrated and b_obs[0] == 0:                    ## don't divide by 0
        return xx,yy,0                                  ## return 0 value
    elif integrated == False and b_obs == 0:
        return xx,yy,0                                  ## return 0 value
   
    ## line ratio calculations for peak and integrated quantities
    ## Ratios of integrated line counts, requires amplitude, distribution center and sigma parameters
    ## There are two ways of doing this; setting a wavelength linearspace and then sum the gaussian over that inverval, OR take the analytical gaussian integral (seems faster).
    if integrated:
        #rg_a = np.linspace(hdr_a_1pix["CRVAL1"]-(hdr_a_1pix["NAXIS1"]/2)*hdr_a_1pix["CDELT1"],hdr_a_1pix["CRVAL1"]+(hdr_a_1pix["NAXIS1"]/2)*hdr_a_1pix["CDELT1"], hdr_a_1pix["NAXIS1"]) ## Tom cnPipeline version, makes the assumption of the wavelength reference being in the center of the wavelength axis space
        #rg_b = np.linspace(hdr_b_1pix["CRVAL1"]-(hdr_b_1pix["NAXIS1"]/2)*hdr_b_1pix["CDELT1"],hdr_b_1pix["CRVAL1"]+(hdr_b_1pix["NAXIS1"]/2)*hdr_b_1pix["CDELT1"], hdr_b_1pix["NAXIS1"])
        # rg_a = np.linspace(hdr_a_1pix["CRVAL1"]-(hdr_a_1pix["CRPIX1"])*hdr_a_1pix["CDELT1"],hdr_a_1pix["CRVAL1"]+(hdr_a_1pix["NAXIS1"]-hdr_a_1pix["CRPIX1"])*hdr_a_1pix["CDELT1"],hdr_a_1pix["NAXIS1"]) ## FITS conforming version (Datacenter?), uses the CRVAL1 keyword for reference. Unclear implementation. 
        # rg_b = np.linspace(hdr_b_1pix["CRVAL1"]-(hdr_b_1pix["CRPIX1"])*hdr_b_1pix["CDELT1"],hdr_b_1pix["CRVAL1"]+(hdr_b_1pix["NAXIS1"]-hdr_b_1pix["CRPIX1"])*hdr_b_1pix["CDELT1"],hdr_b_1pix["NAXIS1"])
        # rat_obs = np.sum(a_obs[0]*np.exp(-(rg_a-a_obs[1])**2/(2*a_obs[2]**2)))/np.sum(b_obs[0]*np.exp(-(rg_b-b_obs[1])**2/(2*b_obs[2]**2))) ##ratios of 3 point gausian integral of input fitting parameters.

        rat_obs = (a_obs[0]*np.sqrt(2*np.pi*(a_obs[2])**2)) / (b_obs[0]*np.sqrt(2*np.pi*(b_obs[2])**2))  ## Analytical gaussian integral formula.


    
    ##ratios of peak intensity only. Ensure that the ratio is scalar regardless of vector input
    elif np.ndim(a_obs) != 0 or np.ndim(b_obs) != 0:   
        rat_obs = a_obs[0]/b_obs[0]  ## requires only one input number for each line 
    else:
        rat_obs = a_obs/b_obs        ## requires only one input number for each line 
    #rat_obs_noise = rat_obs*np.sqrt((np.sqrt(a_obs)/a_obs)**2+(np.sqrt(b_obs)/b_obs)**2)    ## error propagation
   
    
    if (np.isnan(rat_obs) or np.isinf(rat_obs)):        ## discard nans and infs
        return xx,yy,0                                  ## return 0 value
    
    else:                                               ## main loop for valid "rat_obs" value
        ## find the corresponding height (in solar radii) for each pixel
        ## Changed to the manual hpxy helioprojective coordinate calculation
        if (hpxy_a, hpxy_b) == (0,0):
            rpos = np.sqrt( (hdr_a_1pix['CRVAL2']+hdr_a_1pix['CDELT2']*(-hdr_a_1pix['CRPIX2'] + xx))**2 + (hdr_a_1pix['CRVAL3']+hdr_a_1pix['CDELT3']*(-hdr_a_1pix['CRPIX3'] + yy))**2 ) / hdr_a_1pix['SOLARRAD']
        else:
            rpos= np.sqrt(hpxy_a**2+hpxy_b**2)/hdr_a_1pix["SOLARRAD"]
        subh = np.argwhere(chianti_table['h'] > rpos)                                                                                  
        
        ## if height is greater than maximum h (2.0R_sun as in the currently implemented table) just use the 2.0R_sun corresponding ratios.
        if len(subh) == 0: 
            subh = [-1]       
        
        ## make the interpolation function; Quadratic as radial density drop is usually not linear
        ifunc = interpolate.interp1d(chianti_table['rat'][subh[0],:].flatten(),chianti_table['den'], kind="quadratic",fill_value="extrapolate")  
        
        ## apply the interpolation to the data 
        dens_1pix = ifunc(rat_obs) 
        #dens_1pix_noise = ifunc(rat_obs+rat_obs_noise) - dens_1pix
        
        ## debug prints
        #print("Radius from limb: ",rpos," at pixel positions (",xx,yy,")")         ## debug
        #print(len(subh),rpos/rsun)                                                 ## debug          
        return xx,yy,dens_1pix

##########################################################################
####################### COMP_DENS implementation ########################
##########################################################################

def COMP_DENS(i1074,i1079,chianti_link,header_1074,header_1079,integrated=False):
## compute the density from two stokes I observations (serialized parallel runs for 1 pixel at a time) of Fe XIII assuming:
##    - a line ratio look-up table created using chianti atomic data.
##    - accurate header pointing information
##    - that the data arrays are aligned properly
##    - others ??
    
## - One of two CHIANTI look-up tables, computed via PyCELP or SSWIDL(deprecated) aare imported via the "chianti_link" path given when calling this function.
## - Requires full observation headers (1 header per map) from both Fe XIII line datasets.
## - Requires the numpy, multiprocessing, tqdm, scipy's interpolate and io.readsav imports for loading the older IDL look-up table.
    
## Ingesting the inputs and shaping the return array. Sanity checks will return 0 in case something is incompatible.
    if np.ndim(i1074) == 2: ## Function understands a set of two 2D maps containing Stokes I measurements.
        density=np.zeros(i1074.shape,dtype=np.float32)   
        if (i1074.shape != i1079.shape):           ## sanity check 
            print("Input arrays are not of the same physical shape. Aborting!")
            print("Shapes are: Array1074= ", str(i1074.shape)," Array1079= ", str(i1079.shape))
            return density                         ## a zero array at this point
        print('Intensity only maps ingested. Continuing!')
    elif np.ndim(i1074) == 3: ## Function understands a set of two 2D maps containing Stokes I measurements.
        density=np.zeros(i1074.shape[:2],dtype=np.float32)   
        if (i1074.shape != i1079.shape):           ## sanity check 
            print("Input arrays are not of the same physical shape. Aborting!")
            print("Shapes are: Array1074= ", str(i1074.shape)," Array1079= ", str(i1079.shape))
            return density                        ## a zero array at this point
        print('Intensity, Doppler, and Linewidth maps ingested. Continuing!')

    ## Input not compatible with this function
    else:
        print("Input arrays dimensions not understood. Aborting!")
        return np.zeros((2,2))

## theoretical chianti ratio calculations created via PyCELP or SSWIDL        
    ## read the chianti table         
    if (chianti_link[-4:]    == ".npz"):                     ## default
        chianti_table        =  dict(np.load(chianti_link))  ## variables (h,den,rat) directly readable by work_1pix_comp .files required for loading the data directly.
    elif (chianti_link[-4:]  == ".sav"):                     ## deprecated
        print("The sswidl look-up table is deprecated. The pycelp look-up table is better supported and recommended at this stage")
        chianti_table        =  readsav(chianti_link)        ## variables  (h,den,rat) directly readable by work_1pix_comp
        chianti_table["den"] =  chianti_table["den"][0,:]    ## density sample values are the same at all heights. Use just one instance.
    else:
        print("Chianti look-up table not found. Is the path correct?")
        return density                                ## a zero array at this point  
    #print(chianti_table.keys())                      ## Debug - check the arrays
    
    ##  in the save files:
    ## 'h' is an array of heights in solar radii ranged 1.01 - 2.00(pycelp) or 1.01 - 1.50(sswidl);  
    ## 'den' is an array of density (ranged 6.00 to 12.00 (pycelp) or 5.0 to 13.0 (sswidl); the broader interval is not really recommended here. 
    ## 'rat' is a 2D array containing line ratios corresponding to the density range values at each distinct height.
    
    ##  To query the look-up table:
    ##  print(chia['h'],chia['h'].shape,chia['rat'].shape,chia['den'].shape)

    ## look-up table resolutions
    ## pycelp: h is of shape [99]; den and rat are arrays of shape [99, 120] corresponding to the 99 h height values and 120 density values
    ## sswidl: h is of shape [50]; den and rat are arrays of shape [50,  33] corresponding to the 50 h height values and  33 density values

## Print the difference in time between the two frames if you want
    FMT =  '%Y-%m-%dT%H:%M:%S.%f'                 ## format for the dates that are read
    ## Compute the seconds difference between each frame from the shorter dataseries and the entire longer dataseries, 
    if header_1074['INSTRUME']   =="UCoMP": 
        occmask = 0
        aa0 = datetime.strptime(header_1074['DATE-OBS'], FMT) #datetime.strptime(header_1074['DATE-OBS']+'T'+header_1074['TIME-OBS']+'.00', FMT) 
        bb0 = datetime.strptime(header_1079['DATE-OBS'], FMT) #datetime.strptime(header_1079['DATE-OBS']+'T'+header_1079['TIME-OBS']+'.00', FMT)
        print("Time between selected frames: " + str((np.abs( aa0 - bb0 )).total_seconds()) + ' [s]' ) 
    
    
    elif header_1074['INSTRUME'] =="COMP":        ## This includes leap seconds that are not compatible with datetime module. Solve "60"s by changing to "59"s    
        if integrated == True and header_1074["POL_LIST"] == "IQU":
            print("CoMP polarization data does not have the information for an integrated search. Aborting!")
            return density
        occmask =  header_1074['CDELT1']* header_1074['ORADIUS']/ header_1074['RSUN']
        print(occmask)
        aa0 = header_1074['DATE-OBS']+'T'+header_1074['TIME-OBS']+'.00'
        bb0 = header_1079['DATE-OBS']+'T'+header_1079['TIME-OBS']+'.00'
        if aa0[-5:-3] == "60":   aa0 = header_1074['DATE-OBS']+'T'+header_1074['TIME-OBS'][:-2]+'59.00'
        if bb0[-5:-3] == "60":   bb0 = header_1079['DATE-OBS']+'T'+header_1079['TIME-OBS'][:-2]+'59.00'
        aa0 = datetime.strptime(aa0, FMT) 
        bb0 = datetime.strptime(bb0, FMT)
        print("Time between selected frames: " + str((np.abs( aa0 - bb0 )).total_seconds()) + ' [s]' )        
  
## set up the cpu worker and argument arrays
    p         = multiprocessing.Pool(processes=multiprocessing.cpu_count()-2,maxtasksperchild = 10000)     ## dynamically defined from system query as total CPU core number - 2
    ## argument index keeper for splitting tasks to cpu cores
    ## two branches to separate slingle slits vs rasters
    arg_array = []
                                                                                              ## Raster slit branch
    for xx in range(i1074.shape[0]): 
        for yy in range(i1074.shape[1]): 
            arg_array.append((xx,yy,i1074[xx,yy],i1079[xx,yy],header_1074,header_1079,occmask,chianti_table,integrated)) ## Only one header instance goes in for oa set of maps. pointing should be the same in both.

    rs        = p.starmap(work_1pix_comp,arg_array)
    p.close()
    pbar      = tqdm(total=len(arg_array))  ## progress bar implmentation via tqdm
    for i,res in enumerate(rs):
        xx,yy,density[xx,yy] = res 
        pbar.update()
    pbar.close()
    
    density[(density <= 0) | (density > 10.**12)] = np.nan ## negative density is unphysical, Maximum value in lookup table is log_Ne=12. Anything above is ratio errors.

    return density


def work_1pix_comp(xx,yy,a_obs,b_obs,a_hdr,b_hdr,occmask,chianti_table,integrated):
    ## compute the ratio of the observation ## 1074/1079 fraction, same as rat component of the chianti_table 
    
    ## Sanity checks at pixel level
    if integrated and b_obs[0] == 0:                    ## don't divide by 0
        return xx,yy,0                                  ## return 0 value        
    elif  integrated == False and b_obs == 0:           ## don't divide by 0
        return xx,yy,0                                  ## return 0 value

    ## line ratio calculations for peak and integrated quantities
    ## line ratio calculations of integrated line counts, requires amplitude, distribution center and sigma parameters.
    ## There are two ways of doing this; setting a wavelength linearspace and then sum the gaussian over that inverval, OR take the analytical gaussian integral (seems faster).
    ## to not have to load header extensions, we hardcode the central wavelngth sampling to 1074.7 and 1079.8 respectively.
    ## Distribution center and sigma need to be changed from km/s units to nm units.
    if integrated: 
        # rg_a = np.linspace(1074.1,1075.3,50) 
        # rg_b = np.linspace(1079.2,1080.4,50)
        # rat_obs  = np.sum(a_obs[0]*np.exp(-(rg_a-(1074.7+a_obs[1]*1074.7/3e5))**2/(2*(a_obs[2]*1074.7/3e5/2.355)**2)))/np.sum(b_obs[0]*np.exp(-(rg_b-(1079.8+b_obs[2]*1079.8/3e5))**2/(2*(b_obs[2]*1079.8/3e5/2.355)**2))) ##ratios of 3 point gausian integral of input fitting parameters.
        rat_obs = (a_obs[0]*np.sqrt(2*np.pi*(a_obs[2]*1074.7/3e5/2.355)**2)) / (b_obs[0]*np.sqrt(2*np.pi*(b_obs[2]*1079.8/3e5/2.355)**2))  ## Analytical gaussian integral formula.

    ## line ratio calculations for peak quantities only
    else:  
        rat_obs       = a_obs/b_obs        ## requires only one input number for each line 
        #rat_obs_noise = rat_obs*np.sqrt((np.sqrt(a_obs)/a_obs)**2+(np.sqrt(b_obs)/b_obs)**2)    ## error propagation
    
    ## another sanity check for numerical issues
    if (np.isnan(rat_obs) or np.isinf(rat_obs)):        ## discard nans and infs ratios
        return xx,yy,0                                  ## return 0 value
    
    else:                                               ## main loop for valid "rat_obs" value
        ## find the corresponding height (in solar radii) for each pixel
        if a_hdr['instrume']   == "UCoMP":
            rpos = np.sqrt( (a_hdr['CRVAL1'] + a_hdr['CDELT1']*(yy-a_hdr['CRPIX1']-1))**2 + (a_hdr['CRVAL2'] + a_hdr['CDELT2']*(xx-a_hdr['CRPIX2']-1))**2 ) / a_hdr['RSUN_OBS']
            
        elif a_hdr['instrume'] == "COMP":                
            rpos = np.sqrt( (a_hdr['CRVAL1'] + a_hdr['CDELT1']*(yy-a_hdr['CRPIX1']-1))**2 + (a_hdr['CRVAL2'] + a_hdr['CDELT2']*(xx-a_hdr['CRPIX2']-1))**2 ) / a_hdr['RSUN']
        else:                                           ## here you can modify add other image-like datasources
            return xx,yy,0                              ## return 0 value
        
        if rpos <=occmask:    # CEDLT1*ORADIUS/RSUN(COMP) or (uCOMP)
            return xx,yy,0  
            
        subh = np.argwhere(chianti_table['h'] > rpos)                                                                                  
        
        ## if height is greater than maximum h (2.0R_sun as in the currently implemented table) just use the 2.0R_sun corresponding ratios.
        if len(subh) == 0: 
            subh = [-1]       
        
        ## make the interpolation function; Quadratic as radial density drop is usually not linear
        ifunc = interpolate.interp1d(chianti_table['rat'][subh[0],:].flatten(),chianti_table['den'], kind="quadratic",fill_value="extrapolate")  
        
        ## apply the interpolation to the data 
        dens_1pix       = ifunc(rat_obs)                                                                                      
        #dens_1pix_noise = ifunc(rat_obs+rat_obs_noise) - dens_1pix
       
        ## debug prints
        #print("Radius from limb: ",rpos," at pixel positions (",xx,yy,")")         ## debug
        #print(len(subh),rpos/rsun)                                                 ## debug          
        return xx,yy,dens_1pix

##########################################################################
####################### Multithermal implementation of the DKIST_DENS#####
##########################################################################

# ## Usually not required! Effective temperature effects are insignificant.
# ## A ucomp version can also be adapted from this model

# def DKIST_DENS_MULTITEMP(i1074,i1079,tmap,chianti_link,header_1074,header_1079,hpxy=None,integrated=False):
# ## compute the density from two stokes I observations (serialized parallel runs for 1 pix at a time) of Fe XIII assuming:
# ##    - a line ratio look-up table created using chianti atomic data. 
# ##    - accurate header pointing information. shape [number_of_steps]
# ##    - array of pointing information (hpxy) computed using the Cryo-NIRSP public community notebooks https://bitbucket.org/dkist-community-code/cryonirsp-notebooks/src/main/first_release_specFitting/ (optional). shape [2,number_of_steps,slit_size]
# ##    - that the arrays are aligned properly. shape [number_of_steps,slit_size]
# ##    - others ??

# ## function should process single-dimensional, or bi-dimensional input arrays of 1 OR 3 emission line fit parameters. Input array dimension should be between 1 and 3.
# ##       The input arrays contain EITHER peak intensity, OR peak intensity, wavelength center position, and distribution sigma if integrated=True is passed; e.g. 
# ##       [slitsize]                               --single slit position, only peak intensity
# ##       [slitsize,fit_params]                    --single slit position, peak intensity, center, and sigma.
# ##       [number_of_steps,slit_size]              --raster of slit positions, only peak intensity
# ##       [number_of_steps,slit_size,fit_params]   --raster of slit positions,  peak intensity, center, and sigma.

# ## - Two chianti look-up tables are computed via PyCELP or SSWIDL(deprecated) and imported via the "chianti_link" path given to this function.
# ## - Requires full observation headers (1 header per slit step) from both Fe XIII line datasets.
# ## - Requires the numpy, multiprocessing, tqdm, scipy's interpolate and io.readsav imports 
    
# ## input ingestion and return array. Sanity checks will return 0 in case something is incompatible.
    
#     ## single slit position with only peak intensity.
#     if np.ndim(i1074) <= 1:                              ## allows for 1 value calculations
#         density=np.zeros(i1074.shape,dtype=np.float32)   ## single-dimensional spatial array of peak intensities
#         if integrated:
#             print("Not enough parameter inputs for integrated line calculations. Aborting!")
#             return density    
#         if (i1074.shape != i1079.shape):   ## sanity check 
#             print("Arrays are not of the same physical shape. Aborting!")
#             print("Shapes are: Array1074=", str(i1074.shape)," Array1079=", str(i1079.shape))
#             return density
#         print('Single slit position of peak intensity ingested. Continuing!')        
    
#     ## two sub-combinations of inputs.
#     elif np.ndim(i1074) == 2:
#         if integrated:                                                                ## understand a single slit with 3 fit parameters
#             density=np.zeros(i1074.shape[0],dtype=np.float32)   
#             if (i1074.shape[0] != i1079.shape[0]):   ## sanity check 
#                 print("Arrays are not of the same physical shape. Aborting!")
#                 print("Shapes are: Array1074=", str(i1074.shape)," Array1079=", str(i1079.shape))
#                 return density
#             print('Single slit position with peak, center and standard deviation parameters are ingested. Continuing!')
#         else:                                                                         ## understand a raster of slits with only peak intensity.
#             density=np.zeros(i1074.shape,dtype=np.float32)   
#             if (i1074.shape != i1079.shape):   ## sanity check 
#                 print("Arrays are not of the same physical shape. Aborting!")
#                 print("Shapes are: Array1074=", str(i1074.shape)," Array1079=", str(i1079.shape))
#                 return density
#             print('Raster of slit positions of only peak intensity ingested. Continuing!')
    
#     ## bi-dimensional spatial array with added 3rd dimension for peak, center an standard deviation of coronal line.
#     ## not capturing if more than three fit parameters are ingested. 
#     elif np.ndim(i1074) == 3:      
#         density=np.zeros(i1074.shape[:2],dtype=np.float32)
#         if (i1074.shape[:2] != i1079.shape[:2]):   ## sanity check 
#             print("Arrays are not of the same physical shape. Aborting!")
#             print("Shapes are: Array1074=", str(i1074.shape)," Array1079=", str(i1079.shape))
#             return density
#         if i1074.shape[2] < 3:
#             print("Not enough parameter inputs for integrated line calculations. Aborting!")
#             return density
#         print('Raster of slit positions with peak, center and standard deviation parameters are ingested. Continuing!')
        
#     ## Input not compatible with this function
#     else:
#         print("input array dimensions not understood. Aborting!")
#         return np.zeros((2,2))


#     ## theoretical chianti ratio calculations created via PyCELP or SSWIDL        
#     ## read the chianti table         
#     if (chianti_link[-4:]    == ".npz"):                     ## default
#         chianti_table        =  dict(np.load(chianti_link))  ## variables (h,den,rat) directly readable by work_1pix_dkist .files required for loading the data directly.
#     elif (chianti_link[-4:]  == ".sav"):                     ## deprecated
#         print("The SSWIDL look-up table is deprecated. The PyCELP look-up table is better supported and recommended at this stage.")
#         chianti_table        =  readsav(chianti_link)        ## variables  (h,den,rat) directly readable by work_1pix_dkist
#         chianti_table["den"] =  chianti_table["den"][0,:]    ## density sample values are the same at all heights. Use just one instance.
#     else:
#         print("CHIANTI look-up table not found. Is the path correct?")
#         return density                                ## a zero array at this point
#     #print(chianti_table.keys())               ## Debug - check the arrays
    
#     ##  in the save files:
#     ## 'h' is an array of heights in solar radii ranged 1.01 - 2.00(pycelp) or 1.01 - 1.50(sswidl);  
#     ## 'den' is an array of density (ranged 6.00 to 12.00 (pycelp) or 5.0 to 13.0 (sswidl); the broader interval is not really recommended here. 
#     ## 'rat' is a 2D array containing line ratios corresponding to the density range values at each distinct height.
    
#     ##  To query the look-up table:
#     ##  print(chia['h'],chia['h'].shape,chia['rat'].shape,chia['den'].shape)

#     ## look-up table resolutions
#     ## pycelp: h is of shape [99]; den and rat are arrays of shape [99, 120] corresponding to the 99 h height values and 120 density values
#     ## sswidl: h is of shape [50]; den and rat are arrays of shape [50,  33] corresponding to the 50 h height values and  33 density values

# ## set up the cpu worker and argument arrays
#     p          = multiprocessing.Pool(processes=multiprocessing.cpu_count()-2,maxtasksperchild = 10000)     ##ncpus = 54  ## dynamically defined from system query as total CPU core number - 2
#     ## argument index keeper for splitting tasks to cpu cores
#     ## two branches to separate slingle slits vs rasters
#     arg_array  = []

#     ## check if hpxy exists:
#     if hpxy is None: hpxy=np.zeros((2,i1074.shape[0],i1074.shape[1]))   # dummy to skip hpxy-based calculations

    
#     if np.ndim(i1074) <= 1 or (np.ndim(i1074) == 2 and integrated):                                        ## Single slit/value branch
#         for ypos in range(i1074.shape[0]): 
#             arg_array.append((0,ypos,i1074[ypos],i1079[ypos],tmap[ypos],header_1074[0],header_1079[0],chianti_table,integrated))        
#         rs     = p.starmap(work_1pix_dkist_multitemp,arg_array)
#         p.close()
#         pbar   = tqdm(total=len(arg_array))  ## progress bar implmentation via tqdm
#         for i,res in enumerate(rs):
#             step,ypos,density[ypos] = res    ## step is just a dummy in this case; Density defined as a single-dimension vector in this case.
#             pbar.update()
#         pbar.close()
   
#     else:                                                                                                  ## Raster slit branch
#         for step in range(i1074.shape[0]): 
#             for ypos in range(i1074.shape[1]): 
#                 arg_array.append((step,ypos,i1074[step,ypos],i1079[step,ypos],tmap[step,ypos],header_1074[step],header_1079[step],hpxy[0,step,ypos],hpxy[1,step,ypos],chianti_table,integrated))

#         rs     = p.starmap(work_1pix_dkist_multitemp,arg_array)
#         p.close()
#         pbar   = tqdm(total=len(arg_array))  ## progress bar implmentation via tqdm
#         for i,res in enumerate(rs):
#             step,ypos,density[step,ypos] = res 
#             pbar.update()
#         pbar.close()
    
#     density[density <= 0] = np.nan

#     return density

# def work_1pix_dkist_multitemp(xx,yy,a_obs,b_obs,eft_obs,hdr_a_1pix,hdr_b_1pix,hpxy_a,hpxy_b,chianti_table,integrated):
#     ## compute the ratio of the observation ## 1074/1079 fraction, same as rat component of the chianti_table 
    
#     ## Sanity checks at pixel level
#     if integrated and b_obs[0] == 0:                    ## don't divide by 0
#         return xx,yy,0                                  ## return 0 value
#     elif integrated == False and b_obs == 0:
#         return xx,yy,0                                  ## return 0 value
   
#     ## line ratio calculations for peak and integrated quantities
#     ## Ratios of integrated line counts, requires amplitude, distribution center and sigma parameters
#     ## There are two ways of doing this; setting a wavelength linearspace and then sum the gaussian over that inverval, OR take the analytical gaussian integral (seems faster).
#     if integrated:
#         #rg_a = np.linspace(hdr_a_1pix["CRVAL1"]-(hdr_a_1pix["NAXIS1"]/2)*hdr_a_1pix["CDELT1"],hdr_a_1pix["CRVAL1"]+(hdr_a_1pix["NAXIS1"]/2)*hdr_a_1pix["CDELT1"], hdr_a_1pix["NAXIS1"]) ## Tom cnPipeline version, makes the assumption of the wavelength reference being in the center of the wavelength axis space
#         #rg_b = np.linspace(hdr_b_1pix["CRVAL1"]-(hdr_b_1pix["NAXIS1"]/2)*hdr_b_1pix["CDELT1"],hdr_b_1pix["CRVAL1"]+(hdr_b_1pix["NAXIS1"]/2)*hdr_b_1pix["CDELT1"], hdr_b_1pix["NAXIS1"])
#         # rg_a = np.linspace(hdr_a_1pix["CRVAL1"]-(hdr_a_1pix["CRPIX1"])*hdr_a_1pix["CDELT1"],hdr_a_1pix["CRVAL1"]+(hdr_a_1pix["NAXIS1"]-hdr_a_1pix["CRPIX1"])*hdr_a_1pix["CDELT1"],hdr_a_1pix["NAXIS1"]) ## FITS conforming version (Datacenter?), uses the CRVAL1 keyword for reference. Unclear implementation. 
#         # rg_b = np.linspace(hdr_b_1pix["CRVAL1"]-(hdr_b_1pix["CRPIX1"])*hdr_b_1pix["CDELT1"],hdr_b_1pix["CRVAL1"]+(hdr_b_1pix["NAXIS1"]-hdr_b_1pix["CRPIX1"])*hdr_b_1pix["CDELT1"],hdr_b_1pix["NAXIS1"])
#         # rat_obs = np.sum(a_obs[0]*np.exp(-(rg_a-a_obs[1])**2/(2*a_obs[2]**2)))/np.sum(b_obs[0]*np.exp(-(rg_b-b_obs[1])**2/(2*b_obs[2]**2))) ##ratios of 3 point gausian integral of input fitting parameters.

#         rat_obs = (a_obs[0]*np.sqrt(2*np.pi*(a_obs[2])**2)) / (b_obs[0]*np.sqrt(2*np.pi*(b_obs[2])**2))  ## Analytical gaussian integral formula.


    
#     ##ratios of peak intensity only. Ensure that the ratio is scalar regardless of vector input
#     elif np.ndim(a_obs) != 0 or np.ndim(b_obs) != 0:   
#         rat_obs = a_obs[0]/b_obs[0]  ## requires only one input number for each line 
#     else:
#         rat_obs = a_obs/b_obs        ## requires only one input number for each line 
#     #rat_obs_noise = rat_obs*np.sqrt((np.sqrt(a_obs)/a_obs)**2+(np.sqrt(b_obs)/b_obs)**2)    ## error propagation
   
    
#     if (np.isnan(rat_obs) or np.isinf(rat_obs)):        ## discard nans and infs
#         return xx,yy,0                                  ## return 0 value
    
#     else:                                               ## main loop for valid "rat_obs" value
#         ## find the corresponding height (in solar radii) for each pixel
#         ## Changed to the manual hpxy helioprojective coordinate calculation
#         if (hpxy_a, hpxy_b) == (0,0):
#             rpos = np.sqrt( (hdr_a_1pix['CRVAL2']+hdr_a_1pix['CDELT2']*(-hdr_a_1pix['CRPIX2'] + xx))**2 + (hdr_a_1pix['CRVAL3']+hdr_a_1pix['CDELT3']*(-hdr_a_1pix['CRPIX3'] + yy))**2 ) / hdr_a_1pix['SOLARRAD']
#         else:
#             rpos= np.sqrt(hpxy_a**2+hpxy_b**2)/hdr_a_1pix["SOLARRAD"]
#         subh = np.argwhere(chianti_table['h'] > rpos)                                                                                  
        
#         ## if height is greater than maximum h (2.0R_sun as in the currently implemented table) just use the 2.0R_sun corresponding ratios.
#         if len(subh) == 0: 
#             subh = [-1]       
        
#         ## make the interpolation function; Quadratic as radial density drop is usually not linear
#         k = np.argmin(np.abs( np.array([6.10,6.15,6.20,6.23,6.26,6.31,6.36]) - eft_obs))              ## Temperatures capcured in the multitemp lookup
#         ifunc = interpolate.interp1d(chianti_table['rat'][k,subh[0],:].flatten(),chianti_table['den'], kind="quadratic",fill_value="extrapolate")  ## Searck for k in multitemp lookup
        
#         ## apply the interpolation to the data 
#         dens_1pix = ifunc(rat_obs) 
#         #dens_1pix_noise = ifunc(rat_obs+rat_obs_noise) - dens_1pix
        
#         ## debug prints
#         #print("Radius from limb: ",rpos," at pixel positions (",xx,yy,")")         ## debug
#         #print(len(subh),rpos/rsun)                                                 ## debug          
#         return xx,yy,dens_1pix
