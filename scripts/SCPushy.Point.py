
import csv
import os
import sys
from matplotlib.colors import LogNorm
from time import process_time

from result import Ok, Err
import GalaxySN as gsn
import fitsFunctions as ff
import SNR as SNR
import Strehl as strl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scopesim_templates as st
import scopesim as sim
from astropy.visualization import simple_norm, MinMaxInterval, SqrtStretch, ImageNormalize

import numpy as np

from __memos__ import memoize
from __config__ import *

CLUSTER_SEED = 68770

ARCSEC_PER_DEGREE = 3600

MY_SEED = 0x5EED

DATA_DIR = f"{ROOT_DATA_DIR}\scpush"
IMG_DIR = f"{ROOT_IMG_DIR}\scpush"

def move_tableSource(source, x, y):
    for j in range(0,len(source.fields[0])):
        source.fields[0]["x"][j] += x
        source.fields[0]["y"][j] += y

def clusterSN_distanceShift_analysis(SNmag, clusterDistance,diffx,dit,ndit,filter,fw1,fw2,AOsystem,pix_scale,ex_count,mass):
    
    possibilities = [i for i in np.arange(-2,2,0.1)]
    getem = lambda poss : np.random.Generator(np.random.PCG64(MY_SEED)).choice(poss, (4,ex_count), False)
    xlist,ylist,xlistgal,ylistgal = getem(possibilities)
    print(xlist, ylist, xlistgal, ylistgal)
    R=500
    finalarr = np.zeros((2*R,2*R))
    finalgalarr = np.zeros((2*R,2*R))
    unreducedfinalgalarr = np.zeros((2*R,2*R))
    unreducedfinalarr = np.zeros((2*R,2*R))
    emptysky = st.empty_sky()

    #gal_filter is either a normal filter band like "J" or "Ks" (which ScopeSim can handle) or is a local path/svo filter id that specifies what scopesim should use for the galaxy spectrum (scopesim is dumb)
    gal_filter = filter
    if filter == "Y":
        instrument_path = sim.rc.__config__["!SIM.file.local_packages_path"]
        gal_filter = rf"{instrument_path}\MICADO\filters\TC_filter_Y.dat"

    maximum_displacement = max(abs(np.array(list(xlist) + list(ylist) + list(xlistgal) + list(ylistgal)))) * 1000 #mas
    detwidth = np.ceil(2*(maximum_displacement/pix_scale + R))
    print(f"Detector Width: {detwidth}")
    #fov_"size" is half the width/length of the image, converted to angle
    fov_size_as = (R) * (pix_scale / 1000) #arcseconds
    fov_size_radians = (fov_size_as/ARCSEC_PER_DEGREE)*np.pi/180
    #mass = 100_000 #Msol
    distance = clusterDistance #pc
    core_radius = 7 #distance * (fov_size_radians) / 2 #pc
    # hdussky = None
    print(f'core radius = {core_radius}')
    #tan fov_radius_radians = core_radius / distance
    cmds = sim.UserCommands(use_instrument="MICADO",
                    set_modes=[AOsystem,f'IMG_{pix_scale}mas'],  
                    properties={"!OBS.dit": dit,
                                "!OBS.ndit": ndit,
                                "!OBS.filter_name_fw1": fw1,
                                "!OBS.filter_name_fw2": fw2,
                                "!DET.width":detwidth,
                                "!DET.height":detwidth
                                })
    
    micado = sim.OpticalTrain(cmds)
    micado.observe(emptysky)
    hdussky = micado.readout()#filename=f"../artifacts/fitsOutput/{sky_file_name}.fits")
    skydata = np.array(hdussky[0][1].data)
    for i in range(ex_count):
        #scluster = st.cluster(mass = mass, distance = distance, core_radius=core_radius, seed= CLUSTER_SEED)
        #move_tableSource(scluster, xlist[i]+diffx, ylist[i])
        #scluster2 = st.cluster(mass = mass, distance = distance, core_radius=core_radius, seed= CLUSTER_SEED)
        #move_tableSource(scluster2, xlistgal[i]+diffx, ylistgal[i])
        
        SN = st.misc.point_source(sed="sne/sn1a",amplitude=SNmag,filter_curve=filter,x=xlist[i],y=ylist[i])
        source = SN
        cmds = sim.UserCommands(use_instrument="MICADO",
                    set_modes=[AOsystem,f'IMG_{pix_scale}mas'],  
                    properties={"!OBS.dit": dit,
                                "!OBS.ndit": ndit,
                                "!OBS.filter_name_fw1": fw1,
                                "!OBS.filter_name_fw2": fw2,
                                "!DET.width":detwidth,
                                "!DET.height":detwidth
                                })
        micado = sim.OpticalTrain(cmds)
        micado.observe(source)
        hdus = micado.readout()#filename=f"../artifacts/fitsOutput/{file_name}.fits")
        
        cmds = sim.UserCommands(use_instrument="MICADO",
                    set_modes=[AOsystem,f'IMG_{pix_scale}mas'],  
                    properties={"!OBS.dit": dit,
                                "!OBS.ndit": ndit,
                                "!OBS.filter_name_fw1": fw1,
                                "!OBS.filter_name_fw2": fw2,
                                "!DET.width":detwidth,
                                "!DET.height":detwidth
                                })
        #micado = sim.OpticalTrain(cmds)
        #micado.observe(scluster2)
        #hdus2 = micado.readout()#filename=f"../artifacts/fitsOutput/{gal_file_name}.fits")

        # if i==0:
        #     micado3.observe(emptysky)
        #     hdussky = micado3.readout()#filename=f
        # "../artifacts/fitsOutput/{sky_file_name}.fits")
        #subtract empty sky from the two images
        # M_1 + M_2 = Sum_ij M_1i^j + M_2i^j  
        reducedarr = np.array(hdus[0][1].data) - skydata#ff.reduce_array(hdus[0][1].data, hdussky[0][1].data)#reduce(file_name,sky_file_name)
        #reducedgal = np.array(hdus2[0][1].data) - skydata#ff.reduce_array(hdus2[0][1].data, hdussky[0][1].data)#reduce(gal_file_name,sky_file_name)
        widthtot = len(reducedarr)
        SNindexx = int((xlist[i]/(pix_scale*.001))+widthtot/2) # -3.7 * 250 + 500
        SNindexy = int((ylist[i]/(pix_scale*.001))+widthtot/2)

        SNindexxgal = int((xlistgal[i]/(pix_scale*.001))+widthtot/2)
        SNindexygal = int((ylistgal[i]/(pix_scale*.001))+widthtot/2)
        
        croppedarr = reducedarr[SNindexy-R:SNindexy+R,SNindexx-R:SNindexx+R]
        #croppedgal = reducedgal[SNindexygal-R:SNindexygal+R,SNindexxgal-R:SNindexxgal+R]

        finalarr += croppedarr
        #finalgalarr += croppedgal

        #unreducedarr = np.array(hdus[0][1].data)
        #unreducedgal = np.array(hdus2[0][1].data)

        #unreducedcroppedarr = unreducedarr[SNindexy-R:SNindexy+R,SNindexx-R:SNindexx+R]
        #unreducedcroppedgal = unreducedgal[SNindexygal-R:SNindexygal+R,SNindexxgal-R:SNindexxgal+R]

        #unreducedfinalarr += unreducedcroppedarr
        #unreducedfinalgalarr += unreducedcroppedgal

    #galsubtractedarr = finalarr - finalgalarr
    #unreducedpoint = unreducedfinalarr - unreducedfinalgalarr
    print(finalarr)
    return finalarr,

def memoizer(func, filenames):
    output_arrays = []
    if ((not os.path.isfile(filenames[0])) or sys.argv.count("--no-cache") > 0):
        output_arrays = func()
        for i,arr in enumerate(output_arrays):
            with open(filenames[i], "w+") as f :
                w = csv.writer(f, delimiter="|",lineterminator="\n")
                for row in arr:
                    w.writerow(row)

    else:
        for file in filenames:
            with open(file, "r") as f :
                r = csv.reader(f, delimiter="|",lineterminator="\n")
                temp = []
                for row in r:
                    temp.append(row)
                output_arrays.append(temp)
    
    return np.array(output_arrays).astype(np.float32)

def plot_results(data, filePrim):
    fig, ax = plt.subplots(1,1, figsize = (8,8))
    #SN2008BK: 8.5 Msol, 2.5 Mpc, ApMag 22 (progenitor), ApMag 12.5 (SN), AbMag -15 (SN)
    # fig.suptitle("Supernova, SN with Galaxy, Galaxy")
    norm = ImageNormalize(data[0], interval=MinMaxInterval(),stretch=SqrtStretch(), vmin=0, vmax=15000)
    # norm = simple_norm(arr[0:450][0:450], 'sqrt')
    
    im = ax.imshow(data[0], origin='lower', norm=norm, cmap="magma")
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im, cax=cax, orientation="vertical")
    ax.title.set_text("Point Source")
    
    plt.savefig(f"./imgOutputs/{IMG_DIR}/{filePrim}_ps.png")

if __name__ == "__main__":
    #real measured distance to 2008bk
    REF_DISTANCE = 3.5E6
    PIX_SCALES = [1.5, 4]
    #real measured apparent magnitudes in various bands
    snmags = {
                'Y' : 20.05,
                'J' : 19.45,
                'Ks' : 18.39
            }
    
    #measure apparent mag -> Absolute magnitude (ref) + distance -> modified apparent magnitude
    filterBandlist = ["Y", "J", "Ks"]

    #variable distance for simulation
    dlist = [(REF_DISTANCE*3) * i for i in range(0,10)] #definitely want to start higher than 2E6
    dlist[0] = REF_DISTANCE
    totmasslist = [1E5, 5E5]
    
    NDIT_BASE = 50
    EXPOSURE_TIME_CONST = 4
    for i,d in enumerate(dlist):
        for totmass in totmasslist:
            for pix_scale in PIX_SCALES:
                for filterBand in filterBandlist:
                    ex_count = 4
                    ndit = 900
                    # ndit = 900
                    snmag = snmags[filterBand] + 5 * np.log10(d/REF_DISTANCE)
                    print(snmag)
                    filePrim = f"scpushy_{filterBand}_ps{pix_scale}_d{d:.2E}_mass{totmass:.1E}_et{ex_count*ndit}"
                    fileFunc = lambda i : rf"{DATA_DIR}\{filePrim}_{i}.csv"
                    filenames = fileFunc(4)
                    fw1 = filterBand
                    fw2 = "open"
                    if filterBand == "Ks":
                        fw1 = "open"
                        fw2 = filterBand

                    output_arrays = memoize(lambda : clusterSN_distanceShift_analysis(SNmag=snmag,clusterDistance = d,diffx=0,dit=1,ndit=ndit,filter=filterBand,fw1=fw1,fw2=fw2,AOsystem="SCAO",pix_scale=pix_scale,ex_count=ex_count,mass=totmass), filenames)
                    
                    match output_arrays:
                        case Ok(value):
                            plot_results(value.astype(np.float32), filePrim)
                        case Err(e):
                            print(e)
                    
                    print(process_time())
                    
                    
        