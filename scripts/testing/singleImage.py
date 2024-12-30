import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm



#utils

def plot_img_and_spec(file_discriminator, image_cmap = "magma"):
    # image = None
    with fits.open(f'fitsOutput/{file_discriminator}.fits') as file:

        image = file[1].data.copy()
        file.close()


    with fits.open(f'fitsOutput/spec_{file_discriminator}.fits') as file:

        spectra = file[1].data.copy()
        file.close()

    wavelength = np.array(list(spectra))[0:,0]
    amplitude = np.array(list(spectra))[0:,1]

    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8,15))
    im = ax1.imshow(image,cmap=image_cmap)
    plt.colorbar(im, ax=ax1)
    ax2.plot(wavelength, amplitude)
    # plt.xlim(1000,10000)
    plt.show() 
    
def plot_img(file_discriminator, image_cmap = "magma"):
    # image = None
    with fits.open(f'fitsOutput/{file_discriminator}.fits') as file:

        image = file[1].data.copy()
        file.close()

    fig, ax = plt.subplots(1,1, figsize=(20,10))
    im = ax.imshow(image,cmap=image_cmap)
    plt.colorbar(im, ax=ax)
    # plt.xlim(1000,10000)
    plt.show() 

plot_img("arraysource_test_espt1000")