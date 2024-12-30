import csv
import os
import sys
import matplotlib.pyplot as plt
import scopesim_templates as st
import scopesim as sim
import numpy as np
from astropy.visualization import simple_norm, MinMaxInterval, SqrtStretch, ImageNormalize
from dataclasses import dataclass
from fitsFunctions import reduce, reduce_array

from __config__ import *


@dataclass
class GalaxySupernova:
    SNmag=20
    galmag=13
    galr=2
    diffx=.5
    dit=1
    ndit=900
    filter="Ks"
    fw1="open"
    fw2="Ks"
    AOsystem="SCAO"
    pix_scale=4
    ex_count=4

DATA_DIR = "./dataOutputs/galsubtract"

def gal_SN_galsubtract(SNmag,galmag,galr,diffx,dit,ndit,filter,fw1,fw2,AOsystem,pix_scale,ex_count=4):
    
    exposure_time = dit*ndit*ex_count
    file_part = f'SN{SNmag}_gal{galmag}_galR{galr}_dx{diffx}_{filter}_{AOsystem}_{pix_scale}mas_et{exposure_time}_ndit{ndit}'
    file_func = lambda disc: f"{DATA_DIR}/{file_part}_{disc}.csv"
    file_discs = ["reduced", "raw", "galaxy"]
    results = []
    if ((not os.path.isfile(file_func(file_discs[0]))) or sys.argv.count("--no-cache") > 0): 
        xlist0 = np.arange(-1,1,2/ex_count) #[-1, -.5, 0, .5]
        ylist0 = xlist0
        xlist = [round(x,1) for x in xlist0] 
        ylist = [round(x,1) for x in ylist0] 
        xlist0gal = np.arange(-1,1,2/ex_count) 
        ylist0gal = xlist0gal[::-1] 
        xlistgal = [round(x,1) for x in xlist0gal] 
        ylistgal = [round(x,1) for x in ylist0gal]
        R=500
        finalarr = np.zeros((2*R,2*R))
        finalgalarr = np.zeros((2*R,2*R))
        emptysky = st.empty_sky()
        #gal_filter is either a normal filter band like "J" or "Ks" (which ScopeSim can handle) or is a local path/svo filter id that specifies what scopesim should use for the galaxy spectrum (why...)
        gal_filter = filter
        if filter == "Y":
            instrument_path = sim.rc.__config__["!SIM.file.local_packages_path"]
            gal_filter = rf"{instrument_path}\MICADO\filters\TC_filter_Y.dat"
        #detwidth >= 2*(max(abs(all_list))*1000)/ps + R
        maximum_displacement = max(abs(np.array(xlist + ylist + xlistgal + ylistgal))) * 1000 #mas
        detwidth = np.ceil(2*(maximum_displacement/pix_scale + R))
        # if pix_scale == 4:
        #     detwidth = 2000
        # if pix_scale == 1.5:
        #     detwidth = 3668

        for i in range(ex_count):
            galellip = st.elliptical(r_eff=galr,pixel_scale=(pix_scale*.001),amplitude=galmag,filter_name=gal_filter,x_offset=(xlist[i]+diffx),y_offset=ylist[i])
            galellip2 = st.elliptical(r_eff=galr,pixel_scale=(pix_scale*.001),amplitude=galmag,filter_name=gal_filter,x_offset=(xlistgal[i]+diffx),y_offset=ylistgal[i])
            SN = st.misc.point_source(sed="sne/sn1a",amplitude=SNmag,filter_curve=filter,x=xlist[i],y=ylist[i])
            source = galellip + SN
            
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
            micado.observe(galellip2)
            hdus2 = micado.readout()#filename=f"../artifacts/fitsOutput/{gal_file_name}.fits")

            if i==0:
                micado.observe(emptysky)
                hdussky = micado.readout()#filename=f"../artifacts/fitsOutput/{sky_file_name}.fits")
            #subtract empty sky from the two images
            reducedarr = reduce_array(hdus[0][1].data, hdussky[0][1].data)#reduce(file_name,sky_file_name)
            reducedgal = reduce_array(hdus2[0][1].data, hdussky[0][1].data)#reduce(gal_file_name,sky_file_name)
            widthtot = len(reducedarr)
            SNindexx = int((xlist[i]/(pix_scale*.001))+widthtot/2)
            SNindexy = int((ylist[i]/(pix_scale*.001))+widthtot/2)

            SNindexxgal = int((xlistgal[i]/(pix_scale*.001))+widthtot/2)
            SNindexygal = int((ylistgal[i]/(pix_scale*.001))+widthtot/2)
            
            croppedarr = reducedarr[SNindexy-R:SNindexy+R,SNindexx-R:SNindexx+R]
            croppedgal = reducedgal[SNindexygal-R:SNindexygal+R,SNindexxgal-R:SNindexxgal+R]

            for j in range(len(croppedarr)):
                for k in range(len(croppedarr[0])):
                    finalarr[j][k] += croppedarr[j][k]
                    finalgalarr[j][k] += croppedgal[j][k]

        galsubtractedarr = np.zeros((len(croppedarr),len(croppedarr)))

        print(galsubtractedarr.shape, croppedarr.shape, croppedgal.shape)

        for j in range(len(galsubtractedarr)):
            for k in range(len(galsubtractedarr[0])):
                galsubtractedarr[j][k]=finalarr[j][k]-finalgalarr[j][k]

        results = [galsubtractedarr, finalarr, finalgalarr]

        for j, res in enumerate(results):
            with open(file_func(file_discs[j]), "w+") as f :
                w = csv.writer(f, delimiter="|",lineterminator="\n")
                for row in res:
                    w.writerow(row)
    else:
        for disc in file_discs:
            with open(file_func(disc), "r") as f :
                reader = csv.reader(f, delimiter="|",lineterminator="\n")
                result = []
                for row in reader:
                    result.append(row)
                results.append(result)
    
    return tuple(np.array(results).astype(np.float32))


if __name__ == "__main__":
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    output_arrays = gal_SN_galsubtract(SNmag=20,galmag=13,galr=2,diffx=.5,dit=1,ndit=900,filter="Ks",fw1="open",fw2="Ks",AOsystem="SCAO",pix_scale=4,ex_count=4)

    fig, ax = plt.subplots(1,3, figsize = (14,8))
    
    # fig.suptitle("Supernova, SN with Galaxy, Galaxy")
    for i,arr in enumerate(output_arrays):
        norm = ImageNormalize(arr, interval=MinMaxInterval(),stretch=SqrtStretch())
        
        im = ax[i].imshow(arr, origin='lower', norm=norm, cmap="magma")
        
        divider = make_axes_locatable(ax[i])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax, orientation="vertical")
    ax[0].title.set_text("Supernova")
    ax[1].title.set_text("Supernova and Galaxy")
    ax[2].title.set_text("Galaxy")
    plt.savefig("galaxysn_test.png")
    # plt.colorbar()
    plt.show()