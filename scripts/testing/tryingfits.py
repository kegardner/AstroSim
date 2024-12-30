import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

file = fits.open('fitsOutput/TEST.fits')
#getting information about fits file
#file.info()
#read data
data = file[1].data

# print(vars(file[0]))
print(np.shape(data)) #shape of data

plt.imshow(data,cmap='hot')
plt.colorbar()
plt.show()