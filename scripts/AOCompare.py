import matplotlib.pyplot as plt
import scopesim_templates as st
import scopesim as sim
from scopesim.effects import SeeingPSF



pointsource = st.misc.point_source(sed="sne/sn1a", amplitude=15, filter_curve="Ks")
source = pointsource

exposure_time = 1000

cmds1 = sim.UserCommands(use_instrument="MICADO",
                    # mode_yamls=mspec,
                    set_modes=["SCAO","IMG_4mas"],
                    properties={"!OBS.dit": exposure_time,
                                })
cmds2 = sim.UserCommands(use_instrument="MICADO",
                    # mode_yamls=mspec,
                    set_modes=["MCAO","IMG_4mas"],
                    properties={"!OBS.dit": exposure_time,
                                })
cmds3 = sim.UserCommands(use_instrument="MICADO",
                    # mode_yamls=mspec,
                    set_modes=["IMG_4mas"],
                    properties={"!OBS.dit": exposure_time,
                                })
micado = sim.OpticalTrain(cmds1)
micado.observe(source)

hdus = micado.readout()


micado2 = sim.OpticalTrain(cmds2)
micado2.observe(source)

hdus2 = micado2.readout()

cmds3["!OBS.modes"] = ["IMG_4mas"]
print(cmds3)

micado3 = sim.OpticalTrain(cmds3)
micado3.observe(source)

hdus3 = micado3.readout()

micado4 = sim.OpticalTrain(cmds1)
micado4.optics_manager.add_effect(SeeingPSF())
micado4.observe(source)

hdus4 = micado4.readout()


micado5 = sim.OpticalTrain(cmds2)
micado5.optics_manager.add_effect(SeeingPSF())
micado5.observe(source)

hdus5 = micado5.readout()

micado6 = sim.OpticalTrain(cmds3)
micado6.optics_manager.add_effect(SeeingPSF())
micado6.observe(source)

hdus6 = micado6.readout()

fig, axes = plt.subplots(2,3,figsize=(14,8))#

fig.suptitle("Comparison between MICADO Adaptive Optics Systems")
axes[0][0].imshow(hdus[0][1].data)
axes[0][0].title.set_text("SCAO")
axes[0][1].imshow(hdus2[0][1].data)
axes[0][1].title.set_text("MCAO")
axes[0][2].imshow(hdus3[0][1].data)
axes[0][2].title.set_text("No AO")
axes[1][0].imshow(hdus4[0][1].data)
axes[1][0].title.set_text("SCAO - SeeingPSF")
axes[1][1].imshow(hdus5[0][1].data)
axes[1][1].title.set_text("MCAO - SeeingPSF")
axes[1][2].imshow(hdus6[0][1].data)
axes[1][2].title.set_text("No AO - SeeingPSF")
plt.show()
plt.savefig("../thesis_images/AOcompareison.png")