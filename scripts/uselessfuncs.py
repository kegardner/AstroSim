#graveyard of functions no longer used (archive for possible later use)

import numpy as np
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import multivariate_normal
import scipy
# import spextra as sp

from scopesim_templates.extragalactic import galaxy
from scopesim_templates.extragalactic import galaxy3d
from scopesim_templates.extragalactic import spiral_two_component
from scopesim_templates.extragalactic import elliptical
from scopesim_templates.misc import point_source
from scopesim_templates.misc import source_from_file
from scopesim_templates.misc import source_from_array
from scopesim_templates.stellar import star,clusters
from scopesim_templates.stellar import stars
from scopesim_templates.calibration import empty_sky
import os
import astropy.units as u 
from astropy.io.fits import Header
import yaml
from astropy.io import fits

import scopesim as sim
from scopesim_templates.calibration import empty_sky
import os
import yaml

def plot_img_and_spec(file_discriminator, image_cmap = "magma"):
    # image = None
    with fits.open(f'../fitsOutput/{file_discriminator}.fits') as file:

        image = file[1].data.copy()
        file.close()


    with fits.open(f'../fitsOutput/spec_{file_discriminator}.fits') as file:

        spectra = file[1].data.copy()
        file.close()

    wavelength = np.array(list(spectra))[0:,0]
    amplitude = np.array(list(spectra))[0:,1]

    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8,15))
    im = ax1.imshow(image,cmap=image_cmap,vmin=0)
    plt.colorbar(im, ax=ax1)
    ax2.plot(wavelength, amplitude)
    plt.xlim(1000,10000)
    plt.show() 

def plot_spec(file_discriminator):

    with fits.open(f'../fitsOutput/spec_{file_discriminator}.fits') as file:

        spectra = file[1].data.copy()
        file.close()

    wavelength = np.array(list(spectra))[0:,0]
    amplitude = np.array(list(spectra))[0:,1]

    fig, ax = plt.subplots(1,1, figsize=(10,10))
    ax.plot(wavelength, amplitude)
    plt.xlabel("Wavelength (Angstrom)")
    plt.xlim(1000,10000)
    plt.show() 

def plot_spec_set(specfiles, labels):

    fig, ax = plt.subplots(1,1, figsize=(10,10))
    plt.xlabel("Wavelength (Angstrom)")
    plt.xlim(1000,10000)
    wavelength = []
    amplitude = []
    wavelist = []
    amplist = []
    for i in range(len(specfiles)):
        with fits.open(f'../fitsOutput/spec_{specfiles[i]}.fits') as file:
            spectra = file[1].data.copy()
            file.close()

        wavelength = list(np.array(list(spectra))[0:,0])
        amplitude = list(np.array(list(spectra))[0:,1])
        wavelist.append(wavelength)
        amplist.append(amplitude)

        ax.plot(wavelist[i], amplist[i], label=labels[i], linestyle = (i, (1, len(specfiles))), linewidth = 2)
    
    plt.legend()
    plt.show()

def plot_spec_subtraction(spec1,spec2):
    with fits.open(f'../fitsOutput/spec_{spec1}.fits') as file:
        spectra1 = file[1].data.copy()
        file.close()

    with fits.open(f'../fitsOutput/spec_{spec2}.fits') as file2:
        spectra2 = file2[1].data.copy()
        file2.close()

    wavelength1 = list(np.array(list(spectra1))[0:,0])
    amplitude1 = list(np.array(list(spectra1))[0:,1])
    wavelength2 = list(np.array(list(spectra2))[0:,0])
    amplitude2 = list(np.array(list(spectra2))[0:,1])
    ampdiff = []
    d = 0

    if wavelength1[0]==wavelength2[0] and wavelength1[-1]==wavelength2[-1]:
        for i in range(len(amplitude1)):
            d = amplitude1[i] - amplitude2[i]
            ampdiff.append(d)
    else:
        print("wavelength lists don't match")

    fig, ax = plt.subplots(1,1, figsize=(10,10))
    ax.plot(wavelength1, ampdiff, label="difference")
    ax.plot(wavelength1, amplitude1, label = "spec1")
    ax.plot(wavelength1, amplitude2, label = "spec2")
    plt.xlabel("Wavelength (Angstrom)")
    plt.xlim(1000,10000)
    plt.legend()
    plt.show() 