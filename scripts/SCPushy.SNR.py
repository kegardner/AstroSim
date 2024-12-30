from os import listdir
from os.path import isfile, join
import csv
import re
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterExponent
from time import process_time
import GalaxySN as gsn
import fitsFunctions as ff
import SNR as SNR
import Strehl as strl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scopesim_templates as st
import scopesim as sim
from astropy.visualization import simple_norm, MinMaxInterval, SqrtStretch, ImageNormalize
import json
import numpy as np

import SNR as snr
from tqdm import tqdm

from __config__ import ROOT_DATA_DIR, ROOT_IMG_DIR

DATA_DIR = f"{ROOT_DATA_DIR}/scpush"
IMG_DIR = f"{ROOT_IMG_DIR}/scpushy_snr" 
#PATTERN = re.compile("^scpushy_(J|Y|Ks)_ps([0-9.]+)_d([0-9.E+]+)_mass([0-9.E+]+)_et([0-9]+)_([04]).csv$")
PATTERN = re.compile("^scpushy_(J|Y|Ks)_ps([0-9.]+)_d([0-9.E+]+)_mass([0-9.E+]+)_et(3600)_([04]).csv$")
data = {}

for f in (pbar := tqdm(list(filter(lambda n: isfile(join(DATA_DIR, n)) and PATTERN.match(n) is not None, listdir(DATA_DIR))))):
    pbar.set_description(f)
    file = f
    g = PATTERN.match(f).groups()
    fb, ps, d, m, et = g[0:5]

    #for d=3.5E6, use the 400 et file, for fb=Ks and ps=1.5, also use it for distance 1.05E7 
    # one = lambda d, fb, ps: fb == "Ks" and d == "1.05E+07" and ps == "1.5" #use 400, do not use 3600
    two = lambda d : d=="3.50E+06"
    # if ((et == "3600" and not two(d) )#and not one(d, fb, ps)) 
    # or  (et == "200" and (two(d)))): #one(d, fb, ps) or 
    tqdm.write(f"reading {fb} {ps} {d} {m} {et}")
    with open(f"{DATA_DIR}/{f}", "r") as file:
        r = csv.reader(file, delimiter="|",lineterminator="\n")
        temp = []
        for row in r:
            temp.append(row)
        k2 = g[5]
        if data.get(fb, None) is None:
            data[fb] = {}
        if data.get(fb, None).get(ps, None) is None:
            data[fb][ps] = {}
        if data.get(fb, None).get(ps, None).get(d, None) is None:
            data[fb][ps][d] = {}
        if data.get(fb, None).get(ps, None).get(d, None).get(m, None) is None:
            data[fb][ps][d][m] = {}
        data[fb][ps][d][m][k2] = np.array(temp).astype(np.float32)
    # else:
    #     tqdm.write(f"skipping {fb} {ps} {d} {m} {et}")

for fb,sets1 in data.items():
    for ps, sets2 in sets1.items():
        xy = [] #like [(d, snr_m1, snr_m2, snr_i)]
        for d, sets3 in sets2.items():
            try:
                isolate = snr.get_SNR3(list(sets3.items())[0][1]["4"], float(ps), fb)
            except KeyError:
                print(f" failed on: {fb} {ps} {d}")
                exit(1)
            snrs = []
            print(f"{fb} {ps} {d} - {list(sets3.keys())}")
            for m, values in sets3.items():
                snrs.append(snr.get_SNR3(values["0"], float(ps), fb))

            xy.append([float(d), snrs[0], snrs[1], isolate])
            # except IndexError:
            #     print(f" failed on: {fb} {ps} {d}")
            #     exit(1)
        fig, ax = plt.subplots(1,1, figsize = (5,5))
        newd = np.array(xy)
        newd = (newd[newd[:,0].argsort()]).T
        #ax.scatter(newd[0], np.log10(newd[1]), label="Reduced 10^5")
        #ax.plot(newd[0], np.log10(newd[1]), "--")
        #ax.scatter(newd[0], np.log10(newd[2]), label="Reduced 5*10^5")
        #ax.plot(newd[0], np.log10(newd[2]), "--")
        #ax.scatter(newd[0], np.log10(newd[3]), label="Isolated")
        #ax.plot(newd[0], np.log10(newd[3]), "--")
        x = newd[0]

        ax.scatter(x, newd[1], label="10^5 M☉")
        ax.plot(x, newd[1], "--")
        ax.scatter(x, newd[2], label="5*10^5 M☉")
        ax.plot(x, newd[2], "--")
        ax.scatter(x, newd[3], label="0 M☉")
        ax.plot(x, newd[3], "--")
        #ax.set_yscale('log')
        #ax.yaxis.set_major_formatter(LogFormatterExponent())
        #ax.set_ylabel("log(SNR)")
        ax.set_ylabel("SNR")
        ax.set_xlabel("distance (pc)")
        ax.legend(title="Host Star Cluster Mass")
        fig.suptitle(f"{fb} Band Supernova Progenitor Observations")
        ax.set_title("Using Difference Image and Lone-observation")
        plt.savefig(f"{IMG_DIR}/{fb}_{ps}snr_final.png", dpi=300)