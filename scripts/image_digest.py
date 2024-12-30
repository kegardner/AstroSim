from os import listdir
from os.path import isfile, join
import csv
import numpy as np
import re
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.visualization import simple_norm, MinMaxInterval, SqrtStretch, ImageNormalize

from __config__ import *

DATA_DIR = f"{ROOT_DATA_DIR}/scpush"
IMG_DIR = f"{ROOT_IMG_DIR}/scpush" 
PATTERN = re.compile("^scpushy_(Y|Ks|J)_ps([0-9.]+)_d([0-9.E+]+)_mass([0-9.E+]+)_et(3600)_([0124]).csv$")

def generate_grouped_comparisons():
    data = dict(dict())
    for f in filter(lambda n: isfile(join(DATA_DIR, n)) and PATTERN.match(n) is not None, listdir(DATA_DIR)):
        print(f)
        with open(f"{DATA_DIR}/{f}", "r") as d:
            r = csv.reader(d, delimiter="|",lineterminator="\n")
            temp = []
            for row in r:
                temp.append(row)
            g = PATTERN.match(f).groups()
            # print(type(g[0]))
            k = g[0:4]
            k2 = g[4]
            if data.get(k, None) is None:
                data[k] = {}
            data[k][k2] = np.array(temp).astype(np.float32)

    for key, value in data.items():
        print(key)
        print(value.keys())
        print(len(value.keys()))
        fig, ax = plt.subplots(len(value),1, figsize = (8, 8*len(value)))
        file_name = "_".join(key)
        fig.suptitle(f"{key}")
        if type(ax) is plt.Axes:
            ax = [ax]
        for i,(key2,value2) in enumerate(value.items()):
            norm = ImageNormalize(value2, interval=MinMaxInterval(),stretch=SqrtStretch())
            # norm = simple_norm(arr[0:450][0:450], 'sqrt')
            
            im = ax[i].imshow(value2, origin='lower', norm=norm, cmap="magma")
            ax[i].set_title(f"{key2}")
            divider = make_axes_locatable(ax[i])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(im, cax=cax, orientation="vertical")
        plt.savefig(f"{IMG_DIR}/scpushy_study/{file_name}.png")
        plt.close()

def make_group_grid_imgs():
    data = dict(dict())
    for f in filter(lambda n: isfile(join(DATA_DIR, n)) and PATTERN.match(n) is not None, listdir(DATA_DIR)):
        print(f)
        with open(f"{DATA_DIR}/{f}", "r") as d:
            r = csv.reader(d, delimiter="|",lineterminator="\n")
            temp = []
            for row in r:
                temp.append(row)
            g = PATTERN.match(f).groups()
            print(g)
            k = g[0:5]
            k2 = g[5]
            if data.get(k, None) is None:
                data[k] = {}
            data[k][k2] = np.array(temp).astype(np.float32)

    for key, value in data.items():
        # print(value)
        try:
            subdata = [value["1"], value["2"], value["0"], value["4"]] 
            fig, axs = plt.subplots(2,2, figsize = (14, 14))
            file_name = "_".join(key)
            normmax = np.max((subdata[1]))
            for i,value2 in enumerate(subdata):
                norm = ImageNormalize(value2, interval=MinMaxInterval(),stretch=SqrtStretch(), vmin=0, vmax = normmax)
                
                im = axs[(i >> 1) & 1][i%2].imshow(value2, origin='lower', norm=norm, cmap="magma")
                
                divider = make_axes_locatable(axs[(i >> 1) & 1][i%2])
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(im, cax=cax, orientation="vertical")
            axs[0][0].set_title("Progenitor and Cluster")
            axs[0][1].set_title("Cluster")
            axs[1][0].set_title("Difference Image")
            axs[1][1].set_title("Isolated Progenitor")
            
            plt.savefig(f"{IMG_DIR}/scpushy_final/scpushy_{file_name}_final2.png")
            plt.close()
        except:
            continue

if __name__ == "__main__":
    make_group_grid_imgs()