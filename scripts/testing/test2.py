import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import scopesim as sim
import scopesim_templates as sim_tp

cluster = sim_tp.stellar.clusters.cluster(mass=1000,         # Msun
                                          distance=50000,    # parsec
                                          core_radius=0.3,   # parsec
                                          seed=9002)


micado = sim.OpticalTrain("MICADO")
micado.observe(cluster)
hdus = micado.readout()
# micado.readout(filename="TEST.fits")

fig, ax = plt.subplots(1,1,figsize=(10,8))
im = ax.imshow(hdus[0][1].data, norm=LogNorm(vmax=3E4, vmin=3E3), cmap="hot")
ax.set_title("Big Dank")
plt.colorbar(mappable=im, ax=ax)
plt.show()