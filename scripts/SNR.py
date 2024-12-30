from typing import List, Tuple
from warnings import warn
import matplotlib.pyplot as plt
import scopesim_templates as st
import scopesim as sim
import GalaxySN as gsn
import numpy as np
from command import Command
from fitsFunctions import makeFits, reduce
from time import process_time

def get_SNR(sn, sws=10, nws=10, debug = False): 
    """
    Takes a two dimensional array with a source of interest in the center and calculates the signal to noise ratio.

    sn = the two dimensional array
    sws = the half length of the signal window
    nws = the half length of the noise window

    TODO (maybe): add logic to be able to specify where the source is, as opposed to assuming that it is in the center.
    """
    wx, wy = np.shape(sn)
    wx, wy = (int(wx/2), int(wy/2)) #coords of the point source, center of signal window
    wxx, wyy = (wx+150+nws, wy+150+nws) #center of the noise window
    
    signalarea = sn[wx-sws:wx+sws,wy-sws:wy+sws]
    noisearea = sn[wxx-nws:wxx+nws,wyy-nws:wyy+nws]
    S = np.mean(signalarea)
    N = np.std(noisearea)
    SNR = S/N
    if debug:
        print(SNR)
        print(signalarea)
        print(noisearea)
        fig, ax = plt.subplots(1,3, figsize=(14,8))

        ax[0].imshow(sn)
        ax[1].imshow(signalarea)
        ax[2].imshow(noisearea)
        plt.savefig("../thesis_images/SNR_testing.png")
        plt.show()

    return SNR

def get_SNR2(sn, nws=10, ps = 4, filter = "Ks", debug = False): 
    """
    Takes a two dimensional array with a source of interest in the center and calculates the signal to noise ratio.

    sn = the two dimensional array
    sws = the half length of the signal window
    nws = the half length of the noise window

    TODO (maybe): add logic to be able to specify where the source is, as opposed to assuming that it is in the center.
    """
    wx, wy = np.shape(sn)
    wx, wy = (int(wx/2), int(wy/2)) #coords of the point source, center of signal window
    wxx, wyy = (wx+400+nws, wy+nws) #center of the noise window
    signal = 0
    signalsize =0

    #choose signal radius based on filter and pixel scale
    if filter == "Ks":
        if ps == 4:
            R = 2

        if ps == 1.5:
            R = 5

    if filter == "J" or filter == "Y":
        if ps == 4:
            R = 1

        if ps == 1.5:
            R = 2

    

    for i in range(len(sn)):
        for j in range(len(sn[0])):
            if (i-wx)**2+(j-wy)**2 <= R*R:
                signal += sn[i][j]
                signalsize +=1

    
    noisearea = sn[wxx-nws:wxx+nws,wyy-nws:wyy+nws]
    print("noise average = ", np.mean(noisearea))
    N = np.std(noisearea)
    S = signal/signalsize
    SNR = S/N
    if debug:
        print(SNR)
        #print(signalarea)
        print(noisearea)
        fig, ax = plt.subplots(1,3, figsize=(14,8))

        ax[0].imshow(sn)
        #ax[1].imshow(signalarea)
        ax[2].imshow(noisearea)
        plt.savefig("../thesis_images/SNR_testing.png")
        plt.show()

    return SNR

def get_SNR3(sn, ps = 4, filter = "Ks", debug = False): 
    """
    Takes a two dimensional array with a source of interest in the center and calculates the signal to noise ratio.

    sn = the two dimensional array
    ps = pixel scale of image
    filter = photo filter of image

    TODO (maybe): add logic to be able to specify where the source is, as opposed to assuming that it is in the center.
    """
    wx, wy = np.shape(sn)
    wx, wy = (int(wx/2), int(wy/2)) #coords of the point source, center of signal window
    signal = 0
    signalsize =0
    noiselist = []

    #choose signal radius based on filter and pixel scale
    if filter == "Ks":
        if ps == 4:
            R = 5
            R0 = 3

        if ps == 1.5:
            R = 10
            R0 = 9

    if filter == "J" or filter == "Y":
        if ps == 4:
            R = 4
            R0 = 3

        if ps == 1.5:
            R = 6
            R0 = 4

    for i in range(len(sn)):
        for j in range(len(sn[0])):
            if (i-wx)**2+(j-wy)**2 <= R0*R0:
                signal += sn[i][j]
                signalsize +=1

            elif (7*R)**2 <= (i-wx)**2+(j-wy)**2 <= (7*R+10)**2:
                noiselist.append(sn[i][j])

    
    
    print("noise average = ", np.mean(noiselist), " ps = ", ps, " filter = ", filter)
    N = np.std(noiselist)
    S = signal/signalsize
    SNR = S/N
    if debug:
        print(SNR)
        #print(signalarea)
        fig, ax = plt.subplots(1,3, figsize=(14,8))

        ax[0].imshow(sn)
        #ax[1].imshow(signalarea)
        #ax[2].imshow(noisearea)
        plt.savefig("../thesis_images/SNR_testing.png")
        plt.show()

    return SNR

def SNR_unreduced(sn,R=6, nws=10):
    wx, wy = np.shape(sn)
    wx, wy = (int(wx/2), int(wy/2)) #coords of the point source, center of signal window
    wxx, wyy = (wx+400+nws, wy+400+nws) #center of the noise window
    
    noisearea = sn[wxx-nws:wxx+nws,wyy-nws:wyy+nws]
    N = np.std(noisearea)
    Navg = np.mean(noisearea)
    signal = 0
    signalsize =0
    for i in range(len(sn)):
        for j in range(len(sn[0])):
            if (i-wx)**2+(j-wy)**2 <= R*R:
                signal += sn[i][j]
                signalsize +=1
    signalreduce = Navg*signalsize
    signal = signal - signalreduce
    S = signal/signalsize
    ratio = S/N 
    return ratio

def snr_overMagnitudes(dataFunc, magnitudes : List[int], windows : Tuple[int]) -> List[List[float]]:    
    sn_array = []
    for j in magnitudes:
        data = dataFunc(j)
        sn_array.append(get_SNR(data, windows[0], windows[1]))
    return sn_array

def  _snr_manyWindows():
    y_arrays = []
    windows = [i for i in range(5,16)]
    magnitudes = [i for i in range(20,31)]
    def dataFunc(mag):
        sn_fn = makeFits(amp=mag)
        return reduce(sn_fn, f"{sn_fn}_sky")
    for i in windows:
        y_arrays.append(
            snr_overMagnitudes(dataFunc, magnitudes, (i, 20-i))
        )
    fig, ax = plt.subplots(1,2, figsize=(14,8))

    fig.suptitle("SNR vs Magnitude")
    ax[0].set_ylabel("S/N")
    ax[1].set_ylabel("S/N (dB)")
    ax[0].set_xlabel("Magnitude")
    ax[1].set_xlabel("Magnitude")
    for window, snr in zip(windows,y_arrays):
        ax[0].plot(magnitudes, snr, linestyle="dashed", alpha=0.4, zorder=0)
        ax[0].scatter(magnitudes, snr, zorder=1, label=f"{window}|{20-window}")
        
        ax[1].plot(magnitudes, [20*np.log10(i) for i in snr], linestyle="dashed", alpha=0.4, zorder=0)
        ax[1].scatter(magnitudes, [20*np.log10(i) for i in snr], zorder=1, label=f"{window}|{20-window}")

    plt.legend(title="Window Sizes")
    plt.savefig("../thesis_images/SNR_study.png")
    print(f"Run time: {process_time()}")
    plt.show()

if __name__ == "__main__":
    _snr_manyWindows()