from typing import List, Tuple
import matplotlib.pyplot as plt
import scopesim_templates as st
import scopesim as sim
import GalaxySN as gsn
import numpy as np
from command import Command
from fitsFunctions import makeFits, reduce
from time import process_time
from enum import Enum
import scipy as sp

class PhotoFilter(Enum):
    Ks = 1
    J = 2
    Y = 3

def calculate_strehl_from_fits(file_discriminator, pixelscale, photofilter, telescope_radius=19.5):
    Ireal = reduce(file_discriminator,f'{file_discriminator}_sky')
#    center = int(len(Ireal)/2)
#    radius4 = 140
#    radius1p5 = 350
#    if pixelscale == 4:
#        cropIreal = Ireal[center - radius4:center + radius4,center - radius4:center + radius4]
#        pix_scale = 4 * 4.84814e-9 # convert pixel scale to radians per pixel
#
#    if pixelscale == 1.5:
#        cropIreal = Ireal[center - radius1p5:center + radius1p5,center - radius1p5:center + radius1p5]
#        pix_scale = 1.5 * 4.84814e-9 # convert pixel scale to radians per pixel
#
#    width = len(cropIreal)
#    realtot = sum(sum(cropIreal))
#    realmax = cropIreal.max() #cropIreal[int(width/2),int(width/2)] 
#    
#
#    if photofilter == "Ks":
#        k = 2*np.pi/(2.16e-6) #wavenumber 2pi/wavelength(m) Ks band
#
#    if photofilter == "J":
#        k = 2*np.pi/(1.240e-6) #j band
#
#    if photofilter == "Y": #Y band
#        k = 2*np.pi/(1.03e-6)
#
#    ka = k*telescope_radius
#    airydisk = np.zeros([width,width])
#    for i in range(width):
#        for j in range(width):
#            R0 = ((i-width/2)**2+(j-width/2)**2)**(.5)
#            theta = R0*pix_scale
#            u = ka*np.sin(theta)
#            Jval = sp.special.jv(1,u)
#            if u == 0:
#                airydisk[i][j] = 1
#            else:
#                airydisk[i][j] = (2*Jval/u)**2
#
#    idealtot = sum(sum(airydisk))
#    idealmax = airydisk.max()
#
#    realheight = realmax/realtot
#    idealheight = idealmax/idealtot

    return calculate_strehl_from_array(array=Ireal,pixelscale=pixelscale,photofilter=photofilter,telescope_radius=telescope_radius)

def calculate_strehl_from_array(array, pixelscale, photofilter, telescope_radius=19.5):
    Ireal = array
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

    ka = k*telescope_radius
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

    return realheight/idealheight, cropIreal, airydisk*(realtot/idealtot)

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

def make_source_return_strehl(amp,AOsystem,pixelscale,photofilter,telescope_radius = 19.5,espt=1000):
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

    return calculate_strehl_from_fits(file_discriminator, pixelscale, photofilter, telescope_radius)

def generate_plots(magnitudes, photofilter):
    
    strehls = []
    # for pf in ["Ks", "J", "Y"]
    pixelScales = [4, 1.5]
    aoSystems = ["SCAO", "MCAO"]
    for ps in pixelScales:
        for ao in aoSystems:
            #ps == 4, ao == "SCAO" | ps == 4, ao == "MCAO" 
            temp = [] #a set of strehl ratios for a given (ps, ao) over the range of magnitudes
            for i in magnitudes:
                s, _, _ = strehl_ratio_and_arrays(i, ao, ps, photofilter)
                temp.append(s)
            strehls.append(zip(magnitudes, temp))

    fig, ax = plt.subplots(1,1, figsize=(14,8))

    # magnitudes.reverse()
    # strehls.reverse()

    fig.suptitle("Strehl Ratio vs Magnitude")
    ax.set_ylabel("SR")
    ax.set_xlabel("Magnitude")
    ax.set_title(f"{photofilter} Photo Band")
    i = 0
    for ps in pixelScales:
        for ao in aoSystems:
            #filter out any value that isn't in (0,1)
            filteredData = list(filter(lambda coord : 0 < coord[1] < 1, strehls[i]))

            #transpose [(x0,y0), (x1,y1),...] into [[x0,x1,x2,x3,...], [y0,y1,y2,y3,...]] because mpl is dumb
            x,y = np.array(filteredData).T

            ax.plot(x,y, linestyle="dashed", alpha=0.4)
            ax.scatter(x,y, label= f"{ao} | {ps}")
            i+=1
    
    plt.legend(title="AO System | Pixel Scale (mas)")
    ax.set_xlim(max(magnitudes)+.5, min(magnitudes)-.5)
    ax.set_ylim(0, 1)
    plt.savefig(f"../thesis_images/Strehl_study_{photofilter}.png")


def generate_plots_by_photoband(magnitude, photofilters : List[PhotoFilter]):
    
    photoFilterNums = [i.value for i in photofilters]
    photoFilterNames = [i.name for i in photofilters]
    strehls = []
    # for pf in ["Ks", "J", "Y"]
    pixelScales = [4, 1.5]
    aoSystems = ["SCAO", "MCAO"]
    for ps in pixelScales:
        for ao in aoSystems:
            #ps == 4, ao == "SCAO" | ps == 4, ao == "MCAO" 
            temp = [] #a set of strehl ratios for a given (ps, ao) over the range of magnitudes
            for i in photofilters:
                s, _, _ = strehl_ratio_and_arrays(magnitude, ao, ps, i.name)
                temp.append(s)
            strehls.append(zip(photoFilterNums, temp))

    fig, ax = plt.subplots(1,1, figsize=(14,8))

    # magnitudes.reverse()
    # strehls.reverse()

    fig.suptitle("Strehl Ratio vs Photo Band")
    ax.set_ylabel("SR")
    ax.set_xlabel("Photo Band")
    ax.set_title(f"Magnitude {magnitude}")
    
    i = 0
    for ps in pixelScales:
        for ao in aoSystems:
            #filter out any value that isn't in (0,1)
            filteredData = list(strehls[i]) #list(filter(lambda coord : 0 < coord[1] < 1, strehls[i]))

            #transpose [(x0,y0), (x1,y1),...] into [[x0,x1,x2,x3,...], [y0,y1,y2,y3,...]] because mpl is dumb
            x,y = np.array(filteredData).T
            #print(filteredData)
            # ax.plot(x,y, linestyle="dashed", alpha=0.4)
            ax.scatter(x,y, label= f"{ao} | {ps}")
            print(f"{ao}_{ps}_{x}_{y}")
            i+=1
    
    ax.grid(True, linestyle='--', color='gray', alpha=0.3)
    plt.legend(title="AO System | Pixel Scale (mas)")
    ax.set_ylim(0, .8)
    ax.set_xticks(photoFilterNums, photoFilterNames)
    plt.savefig(f"../thesis_images/Strehl_study_allBands_{magnitude}.png")


if __name__ == "__main__":
    # print(PhotoFilter.Ks == PhotoFilter.Ks)
    # print(PhotoFilter.Ks.name)
    magnitude = 22
    generate_plots_by_photoband(magnitude, [PhotoFilter.Ks, PhotoFilter.J, PhotoFilter.Y]) #, , 
    print(f"Run time: {process_time()}")
    plt.show()
