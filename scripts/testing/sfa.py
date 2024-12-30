from scopesim_templates.misc import source_from_array
import matplotlib.pyplot as plt
import numpy as np 
import astropy.units as u 
import yaml
import scopesim as sim


testarray = np.zeros((1000,1000)) 
#testarray[0:999,0:999] = 5000
testarray[497:502,497:502] = 50000
testarray[499,499]=0
testarray[497,497]=0
testarray[501,501]=0
arraysource = source_from_array(arr=testarray,sed="sne/sn1a", amplitude = 15, pixel_scale=1 * u.arcsec, filter_curve="Ks")

# testarray2 = np.ones((9,9)) *100
# testarray2[4,4] = 5000
# arraysource2 = source_from_array(arr=testarray2,sed="sne/sn1a", amplitude = 15, pixel_scale=1 * u.arcsec, filter_curve="Ks")

# testarray3 = np.ones((10,10)) *100
# testarray3[4,4] = 5000
# arraysource3 = source_from_array(arr=testarray3,sed="sne/sn1a", amplitude = 15, pixel_scale=1 * u.arcsec, filter_curve="Ks")


source = arraysource #+ arraysource2 + arraysource3

exposure_time = 1000

file_name = f"arraysource_test_espt{exposure_time}"

sim.rc.__config__["!SIM.file.local_packages_path"] = ".\inst_pkgs\\"

cmds = sim.UserCommands(use_instrument="MICADO",
                    properties={"!OBS.dit": exposure_time,
                                "!OBS.airmass": 2.25,
                                })

micado = sim.OpticalTrain(cmds)
micado.observe(source)
hdus = micado.readout(filename=f"fitsOutput/{file_name}.fits")


plt.imshow(hdus[0][1].data)
plt.show()
