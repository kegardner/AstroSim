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
import astropy.table as table
import yaml
from astropy.io import fits

from matplotlib.colors import LogNorm

import scopesim as sim

# how many stars
n = 200

# random coordinated in a 100 arcsec box
x, y  = 100 * np.random.random(size=(2, n)) - 50

# All stars reference the Vega spectrum
ref = np.zeros(n)
# Each star is given a different weighting, i.e. magnitude
weight = np.ones(n)

# Note: The Pyckles and SpeXtra libraries contain many more stellar and galactic spectra
vega = sim.source.source_templates.vega_spectrum(mag=15)

tbl = table.Table(names=["x", "y", "ref", "weight"],
                  data= [x,    y,   ref,   weight],
                  units=[u.arcsec, u.arcsec, None, None])

table_source = sim.Source(table=tbl, spectra=[vega])
nova1 = point_source(sed="sne/sn1a", amplitude=20, filter_curve="Ks",   x = -1.8, y=-1.8)
source = table_source #+ nova1


exposure_time = 100

file_name = f"arraysource_test_espt{exposure_time}"

sim.rc.__config__["!SIM.file.local_packages_path"] = ".\inst_pkgs\\"
yaml_dicts = []
with open("inst_pkgs\MICADO\default.yaml") as f:
    yaml_dicts += [dic for dic in yaml.full_load_all(f)]
mspec = yaml_dicts

# print(mspec)
cmds = sim.UserCommands(use_instrument="MICADO",
                    # mode_yamls=mspec,
                    # set_modes=["SCAO","SPEC_3000x16"], #"SCAO", 
                    properties={"!OBS.dit": exposure_time,
                                "!OBS.airmass": 2.25,
                                #"!OBS.filter_name_fw1": "J",
                                #"!OBS.filter_name_fw2": "open",
                                #"!SIM.spectral.wave_mid": 1.578,
                                #"!SIM.spectral.spectral_resolution": 1/.00001,
                                #"!TEL.area" : 5000,
                                #"!TEL.pixel_scale" : 1,
                                #"!TEL.temperature" : 20,
                                #"!ATMO.humidity" : 0.8,
                                #"!ATMO.temperature" : 10,
                                #"!ATMO.pressure" : 1,
                                # "!DET.height": 2048,
                                # "!DET.width": 800
                                })
print(cmds.cmds)
micado = sim.OpticalTrain(cmds)
micado.observe(source)
hdus = micado.readout(filename=f"fitsOutput/{file_name}.fits")
print(micado.effects)

plt.imshow(hdus[0][1].data)
plt.show()

