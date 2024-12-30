from scopesim_templates.extragalactic import galaxy
from scopesim_templates.misc import point_source
from scopesim_templates.stellar import star,clusters
import matplotlib.pyplot as plt
import os
import numpy as np 
from scopesim_templates.stellar import stars
import astropy.units as u 

from matplotlib.colors import LogNorm

import scopesim as sim

starset = stars(filter_name="Ks", amplitudes = np.linspace(18,23,6)*u.mag, spec_types= ["A0V", "G2V", "K0III"] *2, x=np.linspace(-1,1,6), y=[0]*6)

#cluster = clusters.cluster(mass=1000,         # Msun
#                             distance=50000,    # parsec
#                             core_radius=0.3,   # parsec
#                             seed=68770)


#gal = galaxy(sed="brown/Arp118", amplitude=15, filter_curve="I", theta=45)  # This will create a galaxy with an S0 SED from the Kinney-Calzetti library (see speXtra)
#agn = point_source(sed="sne/sn2n", amplitude=15, filter_curve="I") # and this an AGN
#agn2 = point_source(sed="sne/sn1a", amplitude=15, filter_curve="I", x=.5, y=.5)
#agn3 = point_source(sed="sne/hyper", amplitude=15, filter_curve="I", x=1, y=1)

#nova1 = point_source(sed="sne/sn1a", amplitude=15, filter_curve="H",   x = -1.8, y=-1.8)
#nova2 = point_source(sed="sne/sn1a", amplitude=15.5, filter_curve="H", x = -0.9, y=-1.8)
#nova3 = point_source(sed="sne/sn1a", amplitude=16, filter_curve="H",   x = 0, y=-1.8)
#nova4 = point_source(sed="sne/sn1a", amplitude=16.5, filter_curve="H", x = 0.9, y=-1.8)
#nova5 = point_source(sed="sne/sn1a", amplitude=17, filter_curve="H",   x = 1.8, y=-1.8)
#nova6 = point_source(sed="sne/sn1a", amplitude=17.5, filter_curve="H", x = -1.8, y=-0.9)
#nova7 = point_source(sed="sne/sn1a", amplitude=18, filter_curve="H",   x = -0.9, y=-0.9)
#nova8 = point_source(sed="sne/sn1a", amplitude=18.5, filter_curve="H", x = 0, y=-0.9)
#nova9 = point_source(sed="sne/sn1a", amplitude=19, filter_curve="H",   x = 0.9, y=-0.9)
#nova10 = point_source(sed="sne/sn1a", amplitude=19.5, filter_curve="H", x = 1.8, y=-0.9)
#nova11 = point_source(sed="sne/sn1a", amplitude=20, filter_curve="H",   x = -1.8, y=0)
#nova12 = point_source(sed="sne/sn1a", amplitude=20.5, filter_curve="H", x = -0.9, y=0)
#nova13 = point_source(sed="sne/sn1a", amplitude=21, filter_curve="H",   x = 0, y=0)
#nova14 = point_source(sed="sne/sn1a", amplitude=21.5, filter_curve="H", x = 0.9, y=0)
#nova15 = point_source(sed="sne/sn1a", amplitude=22, filter_curve="H",   x = 1.8, y=0)
#nova16 = point_source(sed="sne/sn1a", amplitude=22.5, filter_curve="H", x = -1.8, y=0.9)
#nova17 = point_source(sed="sne/sn1a", amplitude=23, filter_curve="H",   x = -0.9, y=0.9)
#nova18 = point_source(sed="sne/sn1a", amplitude=23.5, filter_curve="H", x = 0, y=0.9)
#nova19 = point_source(sed="sne/sn1a", amplitude=24, filter_curve="H",   x = 0.9, y=0.9)
#nova20 = point_source(sed="sne/sn1a", amplitude=24.5, filter_curve="H", x = 1.8, y=0.9)
#nova21 = point_source(sed="sne/sn1a", amplitude=25, filter_curve="H",   x = -1.8, y=1.8)
#nova22 = point_source(sed="sne/sn1a", amplitude=25.5, filter_curve="H", x = -0.9, y=1.8)
#nova23 = point_source(sed="sne/sn1a", amplitude=26, filter_curve="H",   x = 0, y=1.8)
#nova24 = point_source(sed="sne/sn1a", amplitude=26.5, filter_curve="H", x = 0.9, y=1.8)
#nova25 = point_source(sed="sne/sn1a", amplitude=27, filter_curve="H",   x = 1.8, y=1.8)

#novaA = point_source(sed="sne/sn1a", amplitude=21, filter_curve="H",   x = 0, y=-1.3)
#novaB = point_source(sed="sne/sn1a", amplitude=21, filter_curve="H",   x = .01, y=-1.3)
#novaC = point_source(sed="sne/sn1a", amplitude=21, filter_curve="H",   x = 0, y=-0.3)
#novaD = point_source(sed="sne/sn1a", amplitude=21, filter_curve="H",   x = .02, y=-0.3)
#novaE = point_source(sed="sne/sn1a", amplitude=21, filter_curve="H",   x = 0, y=.7)
#novaF = point_source(sed="sne/sn1a", amplitude=21, filter_curve="H",   x = .03, y=.7)
#novaG = point_source(sed="sne/sn1a", amplitude=21, filter_curve="H",   x = 0, y=1.7)
#novaH = point_source(sed="sne/sn1a", amplitude=21, filter_curve="H",   x = .04, y=1.7)

# print(vars(agn))
# print(dir(agn))
#print(agn.spectra[0])
#print(vars(agn.spectra[0]))
#print(dir(agn.spectra[0]))
# print(agn.distance)

#source = agn + agn2 + agn3 + gal # + cluster + 

#source = nova1 + nova2 + nova3 + nova4 + nova5 + nova6 + nova7 + nova8 + nova9 + nova10 + nova11 + nova12 + nova13 + nova14 + nova15 + nova16 + nova17 + nova18 + nova19 + nova20 + nova21 + nova22 + nova23 + nova24 + nova25 
#source = novaA + novaB + novaC + novaD + novaE + novaF +novaG + novaH
source = starset

exposure_time = 1000
file_name = f"starstest_espt{exposure_time}"

micado = sim.OpticalTrain("MICADO")
micado.cmds["!OBS.dit"] = exposure_time

micado.observe(source)
hdus = micado.readout(filename=f"{file_name}.fits")

try:
    os.remove(f"fitsOutput/spec_{file_name}.fits")
except:
    pass  
source.spectra[0].to_fits(f"fitsOutput/spec_{file_name}.fits", wavelengths=None, flux_unit=None, area=None,
                vegaspec=None)


#fig, axes = plt.subplots(2,2,figsize=(10,8))
#axes[0][0].imshow(source.fields[0].data)
# axes[0][1] 
#source.spectra[0].plot()
#agn.spectra[0].plot()
#im = axes[1][0].imshow(hdus[0][1].data, norm=LogNorm(vmax=3E4, vmin=3E3), cmap="hot")
#axes[1][0].set_title("Title")
#plt.colorbar(mappable=im, ax=axes[1][0])
#plt.show()


#mylist = []

#mylist.remove