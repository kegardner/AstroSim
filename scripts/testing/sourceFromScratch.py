import numpy as np
import matplotlib.pyplot as plt
import scopesim_templates as st
from astropy.io import fits
from astropy import units 
from scopesim import Source
from spextra import Spextrum
from synphot import SourceSpectrum
import scopesim as sim
import yaml

arr = np.zeros((100,100))
arr[50:60,50:60] = 5

def _make_spectrum(sed):
    if isinstance(sed, (Spextrum, SourceSpectrum)):
        return Spextrum(modelclass=sed)
    try:
        return Spextrum(sed)
    except TypeError as err:
        raise TypeError("Parameter 'sed' must be Spextrum, SourceSpectrum or "
                        f"str, found {type(sed)=}.") from err


def _make_scaled_spectrum(sed, amplitude, filter_curve):
    spec = _make_spectrum(sed)
    if amplitude is None:
        return spec

    # This might be redundant, check spextra method!
    if not isinstance(amplitude, units.Quantity):
        amplitude = amplitude * units.ABmag

    return spec.scale_to_magnitude(amplitude, filter_curve)

def source_from_array(arr, sed, pixel_scale, amplitude, filter_curve,
                      ra=-10, dec=10) -> Source:
    pixel_scale <<= units.arcsec

    header = st.utils.general_utils.make_img_wcs_header(ra=ra, dec=dec, pixel_scale=pixel_scale,
                                 image_size=arr.shape)
    header.update({"SPEC_REF": 0})

    data = arr / np.sum(arr)
    hdu = fits.ImageHDU(data=data, header=header)

    spec = _make_scaled_spectrum(sed, amplitude, filter_curve)
    src = Source(image_hdu=hdu, spectra=spec)
    return src


arraysource = source_from_array(arr=arr,sed="sne/sn1a", amplitude = 15, pixel_scale=10, filter_curve="Ks")
# pointsource = st.misc.point_source(sed="sne/sn1a", amplitude=15, filter_curve="Ks")
# print(dir(pointsource))
# print(pointsource.image)
print(arraysource.fields)
print(arraysource.fields[0].header)
# plt.imshow(pointsource)
# arraysource.plot()
# pointsource.plot()
# plt.show()

source = arraysource



exposure_time = 10


sim.rc.__config__["!SIM.file.local_packages_path"] = "AstroSim/inst_pkgs/"
yaml_dicts = []
with open("AstroSim/inst_pkgs/MICADO/default.yaml") as f:
    yaml_dicts += [dic for dic in yaml.full_load_all(f)]
mspec = yaml_dicts

# print(mspec)
cmds = sim.UserCommands(use_instrument="MICADO",
                    # mode_yamls=mspec,
                    # set_modes=["SCAO","SPEC_3000x16"], #"SCAO", 
                    properties={"!OBS.dit": exposure_time,
                                # "!OBS.airmass": 2.25,
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

hdus = micado.readout()


plt.imshow(hdus[0][1].data)
plt.show()