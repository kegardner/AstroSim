import sys 
sys.path.append("..")

import csv
import os
import sys
from matplotlib.colors import LogNorm
from time import process_time
import GalaxySN as gsn
import fitsFunctions as ff
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scopesim_templates as st
import scopesim as sim
import Cluster as c
from astropy.visualization import simple_norm, MinMaxInterval, SqrtStretch, ImageNormalize
import numpy as np

CLUSTER_SEED = 68770

ARCSEC_PER_DEGREE = 3600
sim.rc.__config__["!SIM.file.local_packages_path"] = "..\..\inst_pkgs"

def move_tableSource(source, x, y):
    for j in range(0,len(source.fields[0])):
        source.fields[0]["x"][j] += x
        source.fields[0]["y"][j] += y

def ndit_analysis(SNmag, dit,ndit,filter,fw1,fw2,AOsystem,pix_scale):
    R=500
    emptysky = st.empty_sky()

    #gal_filter is either a normal filter band like "J" or "Ks" (which ScopeSim can handle) or is a local path/svo filter id that specifies what scopesim should use for the galaxy spectrum (spextra is dumb)
    gal_filter = filter
    if filter == "Y":
        instrument_path = sim.rc.__config__["!SIM.file.local_packages_path"]
        gal_filter = rf"{instrument_path}\MICADO\filters\TC_filter_Y.dat"

    SN = st.misc.point_source(sed="sne/sn1a",amplitude=SNmag,filter_curve=filter)
    source = SN
    
    cmds = sim.UserCommands(use_instrument="MICADO",
                set_modes=[AOsystem,f'IMG_{pix_scale}mas'],  
                properties={"!OBS.dit": dit,
                            "!OBS.ndit": ndit,
                            "!OBS.filter_name_fw1": fw1,
                            "!OBS.filter_name_fw2": fw2,
                            "!DET.width":1000,
                            "!DET.height":1000
                            })
    micado = sim.OpticalTrain(cmds)
    micado.observe(source)
    hdus = micado.readout()

    micado.observe(emptysky)
    hdussky = micado.readout()
    #subtract empty sky from the two images
    reducedarr = ff.reduce_array(hdus[0][1].data, hdussky[0][1].data)
    
    croppedarr = reducedarr#reducedarr[SNindexy-R:SNindexy+R,SNindexx-R:SNindexx+R]

    
    return croppedarr, hdussky[0][1].data

if __name__ == "__main__":
    d = 1E8
    snMag = -3.2
    snmag = 5* np.log10(d/10) + snMag
    print(snmag)

    output_arrays = ndit_analysis(SNmag=snmag,dit=1,ndit=6*3600,filter="Ks",fw1="open",fw2="Ks",AOsystem="SCAO",pix_scale=4)
    # arr = output_arrays[0] 
    fig, ax = plt.subplots(1,2, figsize = (14,8))
    #SN2008BK: 8.5 Msol, 2.5 Mpc, ApMag 22 (progenitor), ApMag 12.5 (SN), AbMag -15 (SN)
    # fig.suptitle("Supernova, SN with Galaxy, Galaxy")
    for i,arr in enumerate(output_arrays):
        norm = ImageNormalize(arr, interval=MinMaxInterval(),stretch=SqrtStretch())
    # norm = simple_norm(arr[0:450][0:450], 'sqrt')
    
        im = ax[i].imshow(arr, origin='lower', norm=norm, cmap="magma")
    
        divider = make_axes_locatable(ax[i])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax, orientation="vertical")
    ax[0].title.set_text("Supernova")
    ax[1].title.set_text("Empty Sky")
    # ax[2].title.set_text("Star Cluster")
    plt.savefig("dangoltest.png")
    # plt.colorbar()
    print(process_time())
    plt.show()