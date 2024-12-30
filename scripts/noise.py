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

import numpy as np
from __config__ import *

CLUSTER_SEED = 68770

ARCSEC_PER_DEGREE = 3600

MY_SEED = 0x5EED

def move_tableSource(source, x, y):
    for j in range(0,len(source.fields[0])):
        source.fields[0]["x"][j] += x
        source.fields[0]["y"][j] += y
mass = 100_000 #Msol
distance = 1E7 #pc
core_radius = 7
scluster = st.cluster(mass = mass, distance = distance, core_radius=core_radius, seed= CLUSTER_SEED)
# move_tableSource(scluster, xlistgal[i]+diffx, ylistgal[i])
AOSystem = "SCAO"
filterBand = "Ks"
pix_scale="4"
fw1 = "open"
fw2 = filterBand

dit=1
ndit=int(900)
fw1=fw1
fw2=fw2

# SN = st.misc.point_source(sed="sne/sn1a",amplitude=SNmag,filter_curve=filter,x=xlist[i],y=ylist[i])
source = scluster# + SN
sky = st.empty_sky()
cmds = sim.UserCommands(use_instrument="MICADO",
            set_modes=[AOSystem,f'IMG_{pix_scale}mas'],  
            properties={"!OBS.dit": dit,
                        "!OBS.ndit": ndit,
                        "!OBS.filter_name_fw1": fw1,
                        "!OBS.filter_name_fw2": fw2,
                        "!DET.width":3400,
                        "!DET.height":3400
                        })

micado = sim.OpticalTrain(cmds)
micado.observe(source)
clust = np.array(micado.readout()[0][1].data)
micado.observe(sky)
skydata = np.array(micado.readout()[0][1].data)



fig, axs = plt.subplots(1,3,figsize=(14,8))

norm = ImageNormalize(clust, interval=MinMaxInterval(),stretch=SqrtStretch(), vmin=0, vmax=15000)
axs[0].imshow(clust, cmap="magma")
norm = ImageNormalize(skydata, interval=MinMaxInterval(),stretch=SqrtStretch(), vmin=0, vmax=15000)
axs[1].imshow(skydata, cmap="magma")
norm = ImageNormalize(clust - skydata, interval=MinMaxInterval(),stretch=SqrtStretch(), vmin=0, vmax=15000)
im = axs[2].imshow(clust - skydata, cmap="magma")
divider = make_axes_locatable(axs[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
fig.colorbar(im, cax=cax, orientation="vertical")
plt.show()