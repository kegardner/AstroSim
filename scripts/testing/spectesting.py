from scopesim_templates.extragalactic import galaxy
from scopesim_templates.misc import point_source
from scopesim_templates.misc import source_from_file
from scopesim_templates.misc import source_from_array
from scopesim_templates.stellar import star,clusters
import matplotlib.pyplot as plt
import os
import numpy as np 
from scopesim_templates.stellar import stars
import astropy.units as u 
from astropy.io.fits import Header
import yaml
from astropy.io import fits

from matplotlib.colors import LogNorm

import scopesim as sim

#starset = stars(filter_name="Ks", amplitudes = np.linspace(18,23,6)*u.mag, spec_types= ["A0V", "G2V", "K0III"] *2, x=np.linspace(-1,1,6), y=[0]*6)

#cluster = clusters.cluster(mass=1000,         # Msun
#                             distance=50000,    # parsec
#                             core_radius=0.3,   # parsec
#                             seed=68770)


#fitsource does not work yet

gal = galaxy(sed="brown/Arp118", amplitude=15, filter_curve="I", theta=45)
nova1 = point_source(sed="sne/sn1a", amplitude=12, filter_curve="Ks",   x = 0, y=0)
# print(nova1.cube_fields[0].wave)
#nova2 = point_source(sed="sne/sn1a", amplitude=15.5, filter_curve="Ks", x = -0.9, y=-1.8)
#nova3 = point_source(sed="sne/sn1a", amplitude=16, filter_curve="Ks",   x = 0, y=-1.8)
#nova4 = point_source(sed="sne/sn1a", amplitude=16.5, filter_curve="Ks", x = 0.9, y=-1.8)
#nova5 = point_source(sed="sne/sn1a", amplitude=17, filter_curve="Ks",   x = 1.8, y=-1.8)
#nova6 = point_source(sed="sne/sn1a", amplitude=17.5, filter_curve="Ks", x = -1.8, y=-0.9)
#nova7 = point_source(sed="sne/sn1a", amplitude=18, filter_curve="Ks",   x = -0.9, y=-0.9)
#nova8 = point_source(sed="sne/sn1a", amplitude=18.5, filter_curve="Ks", x = 0, y=-0.9)
#nova9 = point_source(sed="sne/sn1a", amplitude=19, filter_curve="Ks",   x = 0.9, y=-0.9)
#nova10 = point_source(sed="sne/sn1a", amplitude=19.5, filter_curve="Ks", x = 1.8, y=-0.9)
#nova11 = point_source(sed="sne/sn1a", amplitude=20, filter_curve="Ks",   x = -1.8, y=0)
#nova12 = point_source(sed="sne/sn1a", amplitude=20.5, filter_curve="Ks", x = -0.9, y=0)
#nova13 = point_source(sed="sne/sn1a", amplitude=21, filter_curve="Ks",   x = 0, y=0)
#nova14 = point_source(sed="sne/sn1a", amplitude=21.5, filter_curve="Ks", x = 0.9, y=0)
#nova15 = point_source(sed="sne/sn1a", amplitude=22, filter_curve="Ks",   x = 1.8, y=0)
#nova16 = point_source(sed="sne/sn1a", amplitude=22.5, filter_curve="Ks", x = -1.8, y=0.9)
#nova17 = point_source(sed="sne/sn1a", amplitude=23, filter_curve="Ks",   x = -0.9, y=0.9)
#nova18 = point_source(sed="sne/sn1a", amplitude=23.5, filter_curve="Ks", x = 0, y=0.9)
#nova19 = point_source(sed="sne/sn1a", amplitude=24, filter_curve="Ks",   x = 0.9, y=0.9)
#nova20 = point_source(sed="sne/sn1a", amplitude=24.5, filter_curve="Ks", x = 1.8, y=0.9)
#nova21 = point_source(sed="sne/sn1a", amplitude=25, filter_curve="Ks",   x = -1.8, y=1.8)
#nova22 = point_source(sed="sne/sn1a", amplitude=25.5, filter_curve="Ks", x = -0.9, y=1.8)
#nova23 = point_source(sed="sne/sn1a", amplitude=26, filter_curve="Ks",   x = 0, y=1.8)
#nova24 = point_source(sed="sne/sn1a", amplitude=26.5, filter_curve="Ks", x = 0.9, y=1.8)
#nova25 = point_source(sed="sne/sn1a", amplitude=27, filter_curve="Ks",   x = 1.8, y=1.8)

#source = cluster 
#source = nova1 + nova2 + nova3 + nova4 + nova5 + nova6 + nova7 + nova8 + nova9 + nova10 + nova11 + nova12 + nova13 + nova14 + nova15 + nova16 + nova17 + nova18 + nova19 + nova20 + nova21 + nova22 + nova23 + nova24 + nova25 
source = nova1


# print(sim.rc.__config__)

exposure_time = 1000

file_name = f"arraysource_test_espt{exposure_time}"

sim.rc.__config__["!SIM.file.local_packages_path"] = "..\inst_pkgs\\"

# print(mspec)
cmds = sim.UserCommands(use_instrument="MICADO",
                    # mode_yamls=mspec,
                    set_modes=["SCAO", "SPEC_3000x48"], #, 
                    properties={"!OBS.dit": exposure_time,
                                "!OBS.airmass": 2.25,
                                #"!OBS.filter_name_fw1": "J",
                                #"!OBS.filter_name_fw2": "open",
                                "!SIM.spectral.wave_mid": 1.578,
                                "!SIM.spectral.spectral_resolution": 1/.00001,
                                #"!TEL.area" : 5000,
                                #"!TEL.pixel_scale" : 1,
                                #"!TEL.temperature" : 20,
                                #"!ATMO.humidity" : 0.8,
                                #"!ATMO.temperature" : 10,
                                #"!ATMO.pressure" : 1,
                                # "!DET.height": 2048,
                                # "!DET.width": 800
                                })
cmds["!DET.height"] = 1000
cmds["!OBS.dit"] = 3600
cmds["!INST.psf.strehl"]= 0.2
cmds["!INST.aperture.x"] = 0             # [arcsec]
cmds["!INST.aperture.y"] = 0
cmds["!INST.aperture.width"] = 3         # [arcsec]
cmds["!INST.aperture.height"] = 0.048  
print(cmds)
micado = sim.OpticalTrain(cmds)
# for i in micado.effects:
#     print(i)
# print(micado["slit_wheel"])
# print(micado["slit_wheel"].include)
# micado["slit_wheel"].include = True
micado.observe(source)
hdus = micado.readout()
# print(micado.effects)

plt.imshow(hdus[0][1].data)
plt.show()
# try:
#     os.remove(f"fitsOutput/spec_{file_name}.fits")
# except:
#     pass  
# source.spectra[0].to_fits(f"fitsOutput/spec_{file_name}.fits", wavelengths=None, flux_unit=None, area=None,
#                 vegaspec=None)