import matplotlib.pyplot as plt
import scopesim_templates as st
import scopesim as sim

pointsource = st.misc.point_source(sed="sne/sn1a", amplitude=15, filter_curve="Ks")

source = pointsource

exposure_time = 10


sim.rc.__config__["!SIM.file.local_packages_path"] = "../inst_pkgs/"

# print(mspec)
cmds1 = sim.UserCommands(use_instrument="MICADO",
                    # mode_yamls=mspec,
                    set_modes=["SCAO","IMG_4mas"],
                    properties={"!OBS.dit": exposure_time,
                                "!OBS.ndit": 100
                                })
cmds2 = sim.UserCommands(use_instrument="MICADO",
                    # mode_yamls=mspec,
                    set_modes=["SCAO","IMG_4mas"],
                    properties={"!OBS.dit": exposure_time*10,
                                "!OBS.ndit": 10
                                })
cmds3 = sim.UserCommands(use_instrument="MICADO",
                    # mode_yamls=mspec,
                    set_modes=["SCAO","IMG_4mas"],
                    properties={"!OBS.dit": exposure_time*100,
                                "!OBS.ndit": 1
                                })
micado = sim.OpticalTrain(cmds1)
micado.observe(source)

hdus = micado.readout()


micado2 = sim.OpticalTrain(cmds2)
print(micado2.optics_manager.all_effects)
micado2.observe(source)

hdus2 = micado2.readout()

micado3 = sim.OpticalTrain(cmds3)
print(micado3.optics_manager.all_effects)
micado3.observe(source)

hdus3 = micado3.readout()

print(hdus3)
fig, axes = plt.subplots(1,3,figsize=(14,4))#

fig.suptitle("SN Point Source with varying DIT/NDIT")
axes[0].imshow(hdus[0][1].data)
axes[0].title.set_text("(10/100)")
axes[1].imshow(hdus2[0][1].data)
axes[1].title.set_text("(100/10)")
axes[2].imshow(hdus3[0][1].data)
axes[2].title.set_text("(1000/1)")
plt.show()