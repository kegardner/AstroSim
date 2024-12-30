
import csv
import os
import sys
from matplotlib.colors import LogNorm
import GalaxySN as gsn
import fitsFunctions as ff
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.visualization import simple_norm, MinMaxInterval, LogStretch, ImageNormalize

import numpy as np

LOG_STRETCH = 100_000

def plot_airy_disk(data):
    fig, main_ax = plt.subplots(figsize=(10, 10))
    divider = make_axes_locatable(main_ax)
    top_ax = divider.append_axes("top", 2.5, pad=0.2, sharex=main_ax)
    right_ax = divider.append_axes("right", 2.5, pad=0.2, sharey=main_ax)

    # make some labels invisible
    top_ax.xaxis.set_tick_params(labelbottom=False)
    right_ax.yaxis.set_tick_params(labelleft=False)

    main_ax.set_xlabel('Airy Disk')
    #main_ax.set_ylabel('dim 2')
    top_ax.set_ylabel('Horizontal Cross-Section')
    right_ax.set_xlabel('Vertical Cross-Section')

    center = int(len(data[0])/2)
    
    z = data #[center-150:center+150,center-150:center+150]
    z_max = max(data[center])
    norm = ImageNormalize(data, interval=MinMaxInterval(),stretch=LogStretch(LOG_STRETCH), vmin=0, vmax=z_max)

    cur_x = 140
    cur_y = 140

    main_ax.imshow(z,cmap="magma",norm=norm)
    main_ax.autoscale(enable=False)
    #plt.colorbar(main_ax.imshow(z,cmap="magma",vmin=0), ax=main_ax)
    right_ax.autoscale(enable=False)
    top_ax.autoscale(enable=False)

    #set subplot limits to be log10
    right_ax.set_xlim(left = 1, right=z_max*1.015)
    top_ax.set_ylim(bottom = 1, top=z_max*1.015)

    v_line = main_ax.axvline(cur_x, color='orangered',linewidth = 1)
    h_line = main_ax.axhline(cur_y, color='forestgreen',linewidth = 1)

    #plot the log10 of the data
    v = np.array(z[:,int(cur_x)])#list(filter(lambda x : x > 1,z[:,int(cur_x)])))
    h = np.array(z[int(cur_y),:])#list(filter(lambda x : x > 1,z[int(cur_y),:])))
    right_ax.set_xscale("log")
    top_ax.set_yscale("log")
    v_prof, = right_ax.plot(v,np.arange(len(v)), 'r-',linewidth = 1,)
    h_prof, = top_ax.plot(np.arange(len(h)),h, 'g-',linewidth = 1)
    right_ax.grid(True, linestyle="--")
    top_ax.grid(True, linestyle="--")

    plt.savefig(r"..\thesis_images\AiryDisk_study_18.png")

    return z_max


def plot_psf(data,z_max):
    fig, main_ax = plt.subplots(figsize=(10, 10))
    divider = make_axes_locatable(main_ax)
    top_ax = divider.append_axes("top", 2.5, pad=0.2, sharex=main_ax)
    right_ax = divider.append_axes("right", 2.5, pad=0.2, sharey=main_ax)
    
    norm = ImageNormalize(data, interval=MinMaxInterval(),stretch=LogStretch(LOG_STRETCH), vmin=0, vmax=z_max)
    # make some labels invisible
    top_ax.xaxis.set_tick_params(labelbottom=False)
    right_ax.yaxis.set_tick_params(labelleft=False)

    main_ax.set_xlabel('PSF')
    #main_ax.set_ylabel('dim 2')
    top_ax.set_ylabel('Horizontal Cross-Section')
    right_ax.set_xlabel('Vertical Cross-Section')

    center = int(len(data[0])/2)

    z = data #[center-150:center+150,center-150:center+150]
    #z_max = max(data[center])

    cur_x = 140
    cur_y = 140

    main_ax.imshow(z,cmap="magma",norm=norm)#,norm=LogNorm(vmin=.01,vmax=z_max, clip=True)
    main_ax.autoscale(enable=False)
    #plt.colorbar(main_ax.imshow(z,cmap="magma",vmin=0), ax=main_ax)
    right_ax.autoscale(enable=False)
    top_ax.autoscale(enable=False)

    #set subplot limits to be log10
    right_ax.set_xlim(left = 1, right=z_max*1.015)
    top_ax.set_ylim(bottom = 1, top=z_max*1.015)

    v_line = main_ax.axvline(cur_x, color='orangered',linewidth = 1)
    h_line = main_ax.axhline(cur_y, color='forestgreen',linewidth = 1)

    #plot the log10 of the data
    v = np.array(z[:,int(cur_x)])#list(filter(lambda x : x > 1,z[:,int(cur_x)])))
    h = np.array(z[int(cur_y),:])#list(filter(lambda x : x > 1,z[int(cur_y),:])))
    right_ax.set_xscale("log")
    top_ax.set_yscale("log")
    v_prof, = right_ax.plot(v,np.arange(len(v)), 'r-',linewidth = 1,)
    h_prof, = top_ax.plot(np.arange(len(h)),h, 'g-',linewidth = 1)
    right_ax.grid(True, linestyle="--")
    top_ax.grid(True, linestyle="--")

    plt.savefig(r"..\thesis_images\PSF_study_18.png")

if __name__ == "__main__":
    dir = r".\dataOutputs"

    a_file = rf"{dir}\airy18.csv"
    e_file = rf"{dir}\ex18.csv"
    airy_result = []
    ex18_result = []
    if ((not os.path.isfile(a_file)) or sys.argv.count("--no-cache") > 0):

        # sntest, sngaltest, galtest = gsn.gal_SN_galsubtract(SNmag=24,galmag=13,galr=2,diffx=.5,dit=1,ndit=900,filter="Ks",fw1="open",fw2="Ks",AOsystem="SCAO",pix_scale=4,ex_count=4)

        strehl18, ex18, Airy18 = ff.strehl_ratio_and_arrays(amp=18,AOsystem="SCAO",pixelscale=4,photofilter="Ks",espt=1000)

        with open(a_file, "w+") as f :
            w = csv.writer(f, delimiter="|",lineterminator="\n")
            for row in Airy18:
                w.writerow(row)

        with open(e_file, "w+") as f :
            w = csv.writer(f, delimiter="|",lineterminator="\n")
            for row in ex18:
                w.writerow(row)
    else:
        with open(a_file, "r") as f:
            reader = csv.reader(f, delimiter="|",lineterminator="\n")
            result = []
            for row in reader:
                result.append(row)
            Airy18 = np.array(result).astype(np.float64)

        with open(e_file, "r") as f:
            reader = csv.reader(f, delimiter="|",lineterminator="\n")
            result = []
            for row in reader:
                result.append(row)
            ex18 = np.array(result).astype(np.float64)
    z_max = plot_airy_disk(Airy18)
    print(f"PSF min:{np.min(ex18)}, max:{np.max(ex18)}")
    plot_psf(ex18,z_max)
    plt.show()