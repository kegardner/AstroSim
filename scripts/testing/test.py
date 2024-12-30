from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

import scopesim as sim
from scopesim.source import source_templates as st

src = st.star_field(n=100,
                    mmax=15,      # [mag]
                    mmin=20,
                    width=200)    # [arcsec]

opt = sim.load_example_optical_train()
opt.cmds["!OBS.dit"] = 60         # [s]
opt.cmds["!OBS.ndit"] = 10

opt.observe(src)
hdulist = opt.readout()[0]



fig,ax = plt.subplots(1,1, figsize=(10,8))
im = ax.imshow(hdulist[1].data, norm=LogNorm(vmin=1))
# plt.colorbar()
plt.colorbar(mappable=im,ax = ax)

# fig.savefig("test.jpg")
plt.show()

