import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import multivariate_normal
import scipy
# import spextra as sp

import scipy as sp

from scopesim_templates.extragalactic import galaxy
from scopesim_templates.extragalactic import galaxy3d
from scopesim_templates.extragalactic import spiral_two_component
from scopesim_templates.extragalactic import elliptical
from scopesim_templates.misc import point_source
from scopesim_templates.misc import source_from_file
from scopesim_templates.misc import source_from_array
from scopesim_templates.stellar import star,clusters
from scopesim_templates.stellar import stars
from scopesim_templates.calibration import empty_sky
import os
import astropy.units as u 
from astropy.io.fits import Header
import yaml
from astropy.io import fits

import scopesim as sim
from scopesim_templates.calibration import empty_sky
import os
import yaml

from __config__ import *

from command import Command

#always ensure this has the correct relative path to inst_pkgs
# sim.rc.__config__["!SIM.file.local_packages_path"] = "..\inst_pkgs"
# yaml_dicts = []
# with open("..\inst_pkgs\MICADO\default.yaml") as f:
#     yaml_dicts += [dic for dic in yaml.full_load_all(f)]
# mspec = yaml_dicts

#plt.rcParams.update({"font.family": "serif", "font.serif": "Computer Modern", "text.usetex": True})

def makeFits2(customCommand : Command):
    source = point_source(customCommand.sed, customCommand.amp, customCommand.filterCurve, customCommand.coordinates[0], customCommand.coordinates[1])

    file_name = f"point{customCommand.amp}_{customCommand.filterCurve}_{customCommand.fieldOfView}_{customCommand.adaptiveOptics}_espt{customCommand.exposureTime}"

    defaults = Command()
    for i in filter(lambda prop: "__" not in prop and prop != "exposureTime", customCommand.__dir__()):
        #walrus operator ":=" since we never need the 'value' variable again
        if (value := customCommand.__getattribute__(i)) != defaults.__getattribute__(i):
            #modify how the non-default args are added to filestring here:
            file_name += f"_{i}-{value}"
           
    emptysky = empty_sky()
    sky_file_name = file_name + f"_sky" 
    
    cmds = customCommand.asCMD()
    micado = sim.OpticalTrain(cmds)
    micado.observe(source)
    hdus = micado.readout(filename=f"../artifacts/fitsOutput/{file_name}.fits")

    micado.observe(emptysky)
    hdussky = micado.readout(filename=f"../artifacts/fitsOutput/{sky_file_name}.fits")
    
    #plt.imshow(hdus[0][1].data)
    #plt.show()
    for i in range(len(source.spectra)):
        try:
            os.remove(f"../artifacts/fitsOutput/spec{i}_{file_name}.fits")
        except:
            pass  
        source.spectra[i].to_fits(f"../artifacts/fitsOutput/spec{i}_{file_name}.fits", wavelengths=None, flux_unit=None, area=None, vegaspec=None)

    return file_name

def makeFits(sed="sne/sn1a", amp=18, filt_cur="Ks", xcoord=0, ycoord=0, exposure_time=1000, fov="IMG_4mas", AOsystem="SCAO", airmass=1.2, fw1="open", fw2="Ks", humidity=0.1, temp=7, pressure=0.755):
    argnames = ["sed", "amp", "filt_cur", "xcoord", "ycoord", "exposure_time", "fov", "AOsystem", "airmass", "fw1", "fw2", "humidity", "temp", "pressure"]
    defaults = ["sne/sn1a", 18, "Ks", 0, 0, 1000, "IMG_4mas", "SCAO", 1.2, "open", "Ks", 0.1, 7, 0.755]
    vals = [sed, amp, filt_cur, xcoord, ycoord, exposure_time, fov, AOsystem, airmass, fw1, fw2, humidity, temp, pressure]
    #checkvalsindex = [0, 3, 4, 8, 9, 10, 11, 12, 13] # index values of args that would be appended to file name if they differ from defaults 
    checkvalsindex = [0, 3, 4, 8, 11, 12, 13]

    source = point_source(sed, amp, filt_cur, xcoord, ycoord)
    file_name = f"point{amp}_{filt_cur}_{fov}_{AOsystem}_espt{exposure_time}"

    for x in checkvalsindex:
        if vals[x] != defaults[x]:
            file_name += f"_{argnames[x]}_{vals[x]}"
           
    emptysky = empty_sky()
    sky_file_name = file_name + f"_sky" 
    
    cmds = sim.UserCommands(use_instrument="MICADO",
                    set_modes=[AOsystem,fov],  
                    properties={"!OBS.dit": exposure_time,
                                "!OBS.airmass": airmass,
                                "!OBS.filter_name_fw1": fw1,
                                "!OBS.filter_name_fw2": fw2,
                                "!ATMO.humidity" : humidity,
                                "!ATMO.temperature" : temp,
                                "!ATMO.pressure" : pressure
                                })
    micado = sim.OpticalTrain(cmds)
    micado.observe(source)
    hdus = micado.readout(filename=f"../artifacts/fitsOutput/{file_name}.fits")

    micado.observe(emptysky)
    hdussky = micado.readout(filename=f"../artifacts/fitsOutput/{sky_file_name}.fits")
    
    for i in range(len(source.spectra)):
        try:
            os.remove(f"../artifacts/fitsOutput/spec{i}_{file_name}.fits")
        except:
            pass  
        source.spectra[i].to_fits(f"../artifacts/fitsOutput/spec{i}_{file_name}.fits", wavelengths=None, flux_unit=None, area=None, vegaspec=None)

    return file_name
    
def plot_img(file_discriminator, image_cmap = "magma"):
    # image = None
    with fits.open(f'../artifacts/fitsOutput/{file_discriminator}.fits') as file:

        image = file[1].data.copy()
        file.close()

    fig, ax = plt.subplots(1,1, figsize=(20,10))
    im = ax.imshow(image,cmap=image_cmap,vmin=0,vmax=150000)
    
    plt.plot()
    plt.colorbar(im, ax=ax)
    plt.show() 
 


def intensity_1D(file_discriminator, y0, x0, R, unit): #unit variable is a string, "pixels", "4mas", or "1.5mas", coordinates and R are given as pixels in the frame, or mas with 0,0 at the center
    #with fits.open(f'fitsOutput/{file_discriminator}.fits') as file: #use this if function is given a fits file name and not an array 

    #    image = file[1].data.copy()
    image = file_discriminator # use this line if the function is given an array instead of a fits file name 
    
    calculate = False
    if unit == "pixels":
        calculate = True
    elif unit == "4mas":
        x0 = int(512 + 250 * x0)
        y0 = int(512 + 250 * y0)
        R = int(R*250)
        calculate = True
    elif unit == "1.5mas":
        x0 = int(511 + 665 * x0)
        y0 = int(511 + 665 * y0)
        R = int(R*665)
        calculate  = True
    else:
        print("invalid input, need unit type")

    if calculate == True:
        limitlist = [R,y0,x0,len(image)-x0,len(image)-y0]
        R2 = min(limitlist)

        ymax = y0+R2
        ymin = y0-R2
        xmax = x0+R2
        xmin = x0-R2
        xlist = image[y0,xmin:xmax]
        ylist = image[ymin:ymax,x0]
    #print("x0 = ", x0, " y0 = ", y0, " R = ", R2)
    return xlist, ylist 

def FWHM(vals):
    peak = max(vals)
    end = min(vals)
    mid = (peak + end)/2
    midlist = []
    count = 0
    while len(midlist) == 0:
        if vals[count] >= mid:
            midlist.append(count)
        count += 1

    count = len(vals) -1
    while len(midlist) == 1:
        if vals[count] >= mid:
            midlist.append(count)
        count -= 1

    diff = midlist[1] - midlist[0] +1

    #print("peak is ", peak, ", min is ", end, ", peak range is ", midlist)
    if diff > len(vals)/2:
        return "none"
    else:
        return diff

def PSF_FWHM_set(nameset, pointcoordy, pointcoordx, radius, coordtype):
    hPSFset = []
    vPSFset = []
    hFWHMset = []
    vFWHMset = []
    skyname = None
    if type(pointcoordy) is int and type(pointcoordx) is int:
        for x in range(len(nameset)):
            skyname = nameset[x] + f"_sky"
            reducedarry = reduce(nameset[x], skyname) 
            hPSF, vPSF = intensity_1D(reducedarry, pointcoordy, pointcoordx, radius, coordtype)
            hPSFset.append(hPSF)
            vPSFset.append(vPSF)
            hFWHMset.append(FWHM(hPSF))
            vFWHMset.append(FWHM(vPSF))
    namecount = 0
    if type(pointcoordy) is list and type(pointcoordx) is list:
        for x in pointcoordx:
            for y in pointcoordy:
                skyname = nameset[namecount] +f"_sky"
                reducedarry = reduce(nameset[namecount], skyname)
                hPSF,vPSF = intensity_1D(reducedarry, y, x, radius, coordtype)
                hPSFset.append(hPSF)
                vPSFset.append(vPSF)
                hFWHMset.append(FWHM(hPSF))
                vFWHMset.append(FWHM(vPSF))
                namecount +=1

    return hPSFset, vPSFset, hFWHMset, vFWHMset 

def plot_PSFs(PSFset, paramname, paramvals, HorV, photoband, AOsystem, fov, filename):
    plt.figure(figsize=[11,6])
    FWHMval = 0
    for i in range(len(PSFset)):
        FWHMval = FWHM(PSFset[i])
        #plt.scatter(range(0,len(PSFset[i])), PSFset[i], label=f"{paramname} = {round(paramvals[i],1)}",s=3)
        plt.plot(range(0,len(PSFset[i])), PSFset[i], alpha = 0.5, label=f"{paramname} = {paramvals[i]}, FWHM = {FWHMval}")
    plt.xlabel("pixels")
    if HorV == "H" or HorV == "h":
        plt.title(f"Horizontal cross-section of PSF, Photometric Filter = {photoband}, AO system = {AOsystem}, {fov} Pixels")
    if HorV =="V" or HorV =="v":
        plt.title(f"Vertical cross-section of PSF, Photometric Filter = {photoband}, AO system = {AOsystem}, {fov} Pixels")
    plt.legend()
    plt.savefig(f"../thesis_images/{filename}")
    plt.show()

def plot_FWHM_vs_param(FWHMlist, paramname, paramvals, filename):
    plt.figure(figsize=[12,10])
    plt.plot(paramvals,FWHMlist)
    plt.xlabel(paramname)
    plt.ylabel("FWHM")
    plt.savefig(f"../thesis_images/{filename}")
    plt.show()

def reduce(unreduced,sky):
    """
    Subtract a sky HDU from a full image HDU
    """
    with fits.open(f'../artifacts/fitsOutput/{unreduced}.fits') as file:
        image1 = file[1].data.copy()
        file.close()
        
    with fits.open(f'../artifacts/fitsOutput/{sky}.fits') as file:
        sky = file[1].data.copy()
        file.close()

    reducedarr = np.zeros((len(sky),len(sky)))
    for i in range(len(image1)):
        for j in range(len(image1[0])):
            reducedarr[i][j] = image1[i][j] - sky[i][j]

    return reducedarr 


def reduce_array(image_array, emptySky_array):
    reducedarr = np.zeros((len(emptySky_array),len(emptySky_array)))
    for i in range(len(image_array)):
        for j in range(len(image_array[0])):
            reducedarr[i][j] = image_array[i][j] - emptySky_array[i][j]

    return reducedarr 

def strehl_ratio(amp,AOsystem,pixelscale,photofilter,a = 19.5,espt=1000):
    if photofilter=="Ks":
        endname = ""
    if photofilter=="J":
        endname= "" #"_fw1_J_fw2_open"
    if photofilter=="Y":
        endname= "" #"_fw1_Y_fw2_open"

    file_discriminator = f"point{amp}_{photofilter}_IMG_{pixelscale}mas_{AOsystem}_espt{espt}{endname}"
    Ireal = reduce(file_discriminator,f'{file_discriminator}_sky')
    center = int(len(Ireal)/2)
    radius4 = 140
    radius1p5 = 350
    if pixelscale == 4:
        cropIreal = Ireal[center - radius4:center + radius4,center - radius4:center + radius4]
        pix_scale = 4 * 4.84814e-9 # convert pixel scale to radians per pixel

    if pixelscale == 1.5:
        cropIreal = Ireal[center - radius1p5:center + radius1p5,center - radius1p5:center + radius1p5]
        pix_scale = 1.5 * 4.84814e-9 # convert pixel scale to radians per pixel

    width = len(cropIreal)
    realtot = sum(sum(cropIreal))
    realmax = cropIreal.max() #cropIreal[int(width/2),int(width/2)] 
    

    if photofilter == "Ks":
        k = 2*np.pi/(2.16e-6) #wavenumber 2pi/wavelength(m) Ks band

    if photofilter == "J":
        k = 2*np.pi/(1.240e-6) #j band

    if photofilter == "Y": #Y band
        k = 2*np.pi/(1.03e-6)

    ka = k*a
    airydisk = np.zeros([width,width])
    for i in range(width):
        for j in range(width):
            R0 = ((i-width/2)**2+(j-width/2)**2)**(.5)
            theta = R0*pix_scale
            u = ka*np.sin(theta)
            Jval = sp.special.jv(1,u)
            if u == 0:
                airydisk[i][j] = 1
            else:
                airydisk[i][j] = (2*Jval/u)**2

    idealtot = sum(sum(airydisk))
    idealmax = airydisk.max()

    realheight = realmax/realtot
    idealheight = idealmax/idealtot

    strehl = realheight/idealheight
    #print("realtot = ", realtot)
    #print("realmax = ", realmax)
    #print("idealtot = ", idealtot)
    #print("idealmax = ", idealmax)

    return strehl

def strehl_ratio2(file_discriminator,pixelscale,photofilter,a = 19.5):

    Ireal = reduce(file_discriminator,f'{file_discriminator}_sky')
    center = int(len(Ireal)/2)
    radius4 = 140
    radius1p5 = 350
    if pixelscale == 4:
        cropIreal = Ireal[center - radius4:center + radius4,center - radius4:center + radius4]
        pix_scale = 4 * 4.84814e-9 # convert pixel scale to radians per pixel

    if pixelscale == 1.5:
        cropIreal = Ireal[center - radius1p5:center + radius1p5,center - radius1p5:center + radius1p5]
        pix_scale = 1.5 * 4.84814e-9 # convert pixel scale to radians per pixel

    width = len(cropIreal)
    realtot = sum(sum(cropIreal))
    realmax = cropIreal.max() #cropIreal[int(width/2),int(width/2)] 
    

    if photofilter == "Ks":
        k = 2*np.pi/(2.16e-6) #wavenumber 2pi/wavelength(m) Ks band

    if photofilter == "J":
        k = 2*np.pi/(1.240e-6) #j band

    if photofilter == "Y": #Y band
        k = 2*np.pi/(1.03e-6)

    ka = k*a
    airydisk = np.zeros([width,width])
    for i in range(width):
        for j in range(width):
            R0 = ((i-width/2)**2+(j-width/2)**2)**(.5)
            theta = R0*pix_scale
            u = ka*np.sin(theta)
            Jval = sp.special.jv(1,u)
            if u == 0:
                airydisk[i][j] = 1
            else:
                airydisk[i][j] = (2*Jval/u)**2

    idealtot = sum(sum(airydisk))
    idealmax = airydisk.max()

    realheight = realmax/realtot
    idealheight = idealmax/idealtot

    strehl = realheight/idealheight

    return strehl

def strehl_ratio_and_arrays(amp,AOsystem,pixelscale,photofilter,a = 19.5,espt=1000):
    if photofilter=="Ks":
        fw1 = "open"
        fw2 = "Ks"

    if photofilter=="J":
        fw1 = "J"
        fw2 = "open"

    if photofilter=="Y":
        fw1 = "Y"
        fw2 = "open"

    file_discriminator = makeFits(amp=amp,AOsystem=AOsystem,fov=f'IMG_{pixelscale}mas',filt_cur=photofilter,fw1=fw1,fw2=fw2,exposure_time=espt)
    Ireal = reduce(file_discriminator,f'{file_discriminator}_sky')
    center = int(len(Ireal)/2)
    radius4 = 140
    radius1p5 = 350
    if pixelscale == 4:
        cropIreal = Ireal[center - radius4:center + radius4,center - radius4:center + radius4]
        pix_scale = 4 * 4.84814e-9 # convert pixel scale to radians per pixel

    if pixelscale == 1.5:
        cropIreal = Ireal[center - radius1p5:center + radius1p5,center - radius1p5:center + radius1p5]
        pix_scale = 1.5 * 4.84814e-9 # convert pixel scale to radians per pixel

    width = len(cropIreal)
    realtot = sum(sum(cropIreal))
    realmax = cropIreal.max() #cropIreal[int(width/2),int(width/2)] 
    

    if photofilter == "Ks":
        k = 2*np.pi/(2.16e-6) #wavenumber 2pi/wavelength(m) Ks band

    if photofilter == "J":
        k = 2*np.pi/(1.240e-6) #j band

    if photofilter == "Y": #Y band
        k = 2*np.pi/(1.03e-6)

    ka = k*a
    airydisk = np.zeros([width,width])
    for i in range(width):
        for j in range(width):
            R0 = ((i-width/2)**2+(j-width/2)**2)**(.5)
            theta = R0*pix_scale
            u = ka*np.sin(theta)
            Jval = sp.special.jv(1,u)
            if u == 0:
                airydisk[i][j] = 1
            else:
                airydisk[i][j] = (2*Jval/u)**2

    idealtot = sum(sum(airydisk))
    idealmax = airydisk.max()

    realheight = realmax/realtot
    idealheight = idealmax/idealtot

    strehl = realheight/idealheight
    #print("realtot = ", realtot)
    #print("realmax = ", realmax)
    #print("idealtot = ", idealtot)
    #print("idealmax = ", idealmax)

    return strehl, cropIreal, airydisk*(realtot/idealtot)