from scopesim.effects import SeeingPSF
import scopesim as sim
import matplotlib.pyplot as plt
import scopesim_templates as st
import scopesim as sim

from __config__ import *
# arraysource = source_from_array(arr=arr,sed="sne/sn1a", amplitude = 15, pixel_scale=10, filter_curve="Ks")
pointsource = st.misc.point_source(sed="sne/sn1a", amplitude=15, filter_curve="Ks")

cmds = sim.UserCommands( use_instrument = "MICADO"
                       , set_modes = ["SCAO", "IMG_4mas"]
                       , properties = {
                           "!OBS.dit" : 1000
                       } 
                        )
cmds["!OBS.modes"] = ["IMG_4mas"]
ot = sim.OpticalTrain(cmds)
ot.optics_manager.add_effect(SeeingPSF())
ot.observe(pointsource)
hdu = ot.readout()
print(hdu[0][1].data)
plt.imshow(hdu[0][1].data)
plt.show()