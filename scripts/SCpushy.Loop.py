from typing import Tuple
from time import process_time
import SNR as SNR
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scopesim_templates as st
import scopesim as sim
from astropy.visualization import MinMaxInterval, SqrtStretch, ImageNormalize

from result import Ok, Err
import numpy as np

from __memos__ import memoize
from __config__ import *

CLUSTER_SEED = 68770

ARCSEC_PER_DEGREE = 3600

MY_SEED = 0x5EED

DATA_DIR = f"{ROOT_DATA_DIR}/scpush"
IMG_DIR = f"{ROOT_IMG_DIR}/scpush"
CACHING_ARG = "scpush-loop"

def move_tableSource(source, x, y):
    for j in range(0,len(source.fields[0])):
        source.fields[0]["x"][j] += x
        source.fields[0]["y"][j] += y

def crop(data : np.ndarray, center : Tuple[int, int], half_width : int) -> np.ndarray:
    """
    Crop a given 2D image around a given center, to given half_width

    Parameters
    ----------
    data : np.ndarray
        The image to crop
    center : tuple of int
        The indices specifying the part of the image to crop around
    half_width : int
        The half_width of the output image
    """
    return data[ center[1]-half_width : center[1]+half_width
               , center[0]-half_width : center[0]+half_width ]

def calculate_shifted_center_indices(shift: Tuple[float, float], pixel_scale: float, width: int):
    """
    With the center of an image being (0,0), generate the 2D index given the shift from center, the width of the image, and the pixel scale

    Parameters
    ----------
    shift : tuple of float
        The shift from the center of an image (in as).
    pixel_scale : float
        The size of each pixel in mas
    width : int
        The width of the image (image _must_ be square)
    """
    return ( int((shift[0]/(pixel_scale*.001))+width/2)
           , int((shift[1]/(pixel_scale*.001))+width/2)
           )

def image_differencing_simulation_SNCluster(SNmag, clusterDistance,diffx,dit,ndit,filter,fw1,fw2,AOsystem,pix_scale,ex_count,mass):
    """
    Performs an image differencing simulation for a supernova nested at the center of a star cluster

    Parameters
    ----------
    SNmag : float
        The apparent magnitude of the supernova (should be calculated to ensure that it is representative of that supernova at the desired distance relative to the cluster)
    clusterDistance : float (parsecs)
        The distance to set the star cluster at
    diffx : float (mas)
        The horizontal bias between the center of the SN and the star cluster
    dit : float (s)
        The detector integration time to use for each exposure
    ndit : int
        The number of detector integrations to use for each exposure
    filter : string
        The filter band to use in simulation (must be supported by ScopeSIM, spextra, etc); generally one of "Ks", "J", or "Y"
    fw1 : string
        The primary filter wheel setting
    fw2 : string
        The secondary filter wheel setting
    AOsystem : string
        The optics system to use for MICADO (either "SCAO" or "MCAO")
    pix_scale : float (mas)
        The pixel scale to use with MICADO (must be 1.5 or 4)
    ex_count : int
        The number of exposures to perform
    mass : float
        The total mass of the star cluster
    """
    possibilities = [i for i in np.arange(-2,2,0.1)]
    getem = lambda poss : np.random.Generator(np.random.PCG64(MY_SEED)).choice(poss, (4,ex_count), False)
    xlist,ylist,xlistgal,ylistgal = getem(possibilities)
    print(xlist, ylist, xlistgal, ylistgal)
    R=500
    subject = np.zeros((2*R,2*R))
    isolate = np.zeros((2*R,2*R))
    isolate_cropped = np.zeros((2*R,2*R))
    subject_cropped = np.zeros((2*R,2*R))
    emptysky = st.empty_sky()

    #gal_filter is either a normal filter band like "J" or "Ks" (which ScopeSim can handle) or is a local path/svo filter id that specifies what scopesim should use for the galaxy spectrum (scopesim is dumb)
    gal_filter = filter
    if filter == "Y":
        instrument_path = sim.rc.__config__["!SIM.file.local_packages_path"]
        gal_filter = rf"{instrument_path}\MICADO\filters\TC_filter_Y.dat"

    maximum_displacement = max(abs(np.array(list(xlist) + list(ylist) + list(xlistgal) + list(ylistgal)))) * 1000 #mas
    detwidth = np.ceil(2*(maximum_displacement/pix_scale + R))
    print(f"Detector Width: {detwidth}")
    #fov_"size" is half the width/length of the image, converted to angle
    fov_size_as = (R) * (pix_scale / 1000) #arcseconds
    fov_size_radians = (fov_size_as/ARCSEC_PER_DEGREE)*np.pi/180
    #mass = 100_000 #Msol
    distance = clusterDistance #pc
    core_radius = 7 #distance * (fov_size_radians) / 2 #pc
    # hdussky = None
    print(f'core radius = {core_radius}')
    #tan fov_radius_radians = core_radius / distance
    cmds = lambda : sim.UserCommands(
                        use_instrument="MICADO", 
                        set_modes=[AOsystem,f'IMG_{pix_scale}mas'],  
                        properties={"!OBS.dit": dit,
                                    "!OBS.ndit": ndit,
                                    "!OBS.filter_name_fw1": fw1,
                                    "!OBS.filter_name_fw2": fw2,
                                    "!DET.width":detwidth,
                                    "!DET.height":detwidth
                                    }
                    )
    
    micado = sim.OpticalTrain(cmds())
    micado.observe(emptysky)
    hdussky = micado.readout()#filename=f"../artifacts/fitsOutput/{sky_file_name}.fits")
    skydata = np.array(hdussky[0][1].data)


    def _observe(source):
        #build the optical train, observe the given source, return the hdus
        ot = sim.OpticalTrain(cmds())
        ot.observe(source)
        return ot.readout()
    for i in range(ex_count):
        scluster = st.cluster(mass = mass, distance = distance, core_radius=core_radius, seed= CLUSTER_SEED)
        move_tableSource(scluster, xlist[i]+diffx, ylist[i])
        scluster2 = st.cluster(mass = mass, distance = distance, core_radius=core_radius, seed= CLUSTER_SEED)
        move_tableSource(scluster2, xlistgal[i]+diffx, ylistgal[i])

        SN = st.misc.point_source(sed="sne/sn1a",amplitude=SNmag,filter_curve=filter,x=xlist[i],y=ylist[i])
        source = scluster + SN

        hdus = _observe(source)
        hdus2 = _observe(scluster2)


        # "_" indicates temporary, i.e. one exposure (e.g. _subject is the subject for one exposure, "subject" is the sum of all "_subject")
        _subject = np.array(hdus[0][1].data)
        _isolate = np.array(hdus2[0][1].data)
        print(f"Raw: {_subject} {_subject.shape}")
        _subject_reduced = _subject - skydata
        _isolate_reduced = _isolate - skydata
        print(f"Reduced: {_subject_reduced} {_subject_reduced.shape}")
        widthtot = len(_subject_reduced)
        #where in the image data that the center of the subject/isolate is
        # _subject_center_x = int((xlist[i]/(pix_scale*.001))+widthtot/2) # -3.7 * 250 + 500
        # _subject_center_y = int((ylist[i]/(pix_scale*.001))+widthtot/2)
        _subject_center = calculate_shifted_center_indices((xlist[i], ylist[i]), pix_scale, widthtot)
        # _isolate_center_x = int((xlistgal[i]/(pix_scale*.001))+widthtot/2)
        # _isolate_center_y = int((ylistgal[i]/(pix_scale*.001))+widthtot/2)
        _isolate_center = calculate_shifted_center_indices((xlistgal[i], ylistgal[i]), pix_scale, widthtot)
        print(f"Center: {_subject_center}")
        _subject_cropped_reduced = crop(_subject_reduced, _subject_center, R)
        print(f"Cropped Reduced: {_subject_cropped_reduced}")
        
        #_subject_reduced[_subject_center_y-R:_subject_center_y+R,_subject_center_x-R:_subject_center_x+R]
        _isolate_cropped_reduced = crop(_isolate_reduced, _isolate_center, R)
        #_isolate_reduced[_isolate_center_y-R:_isolate_center_y+R,_isolate_center_x-R:_isolate_center_x+R]
        #compose observations together
        subject += _subject_cropped_reduced
        isolate += _isolate_cropped_reduced

        # #cropping
        # _subject_cropped = crop(_subject, _subject_center, R)
        # #_subject[_subject_center_y-R:_subject_center_y+R,_subject_center_x-R:_subject_center_x+R]
        # _isolate_cropped = crop(_isolate, _isolate_center, R)
        # #_isolate[_isolate_center_y-R:_isolate_center_y+R,_isolate_center_x-R:_isolate_center_x+R]

        # subject_cropped += _subject_cropped
        # isolate_cropped += _isolate_cropped

    difference = subject - isolate
    print(difference.shape, subject.shape, isolate.shape)
    return difference, subject, isolate

def plot_results(data, filePrim):
    fig, ax = plt.subplots(1,4, figsize = (16,8))
    #SN2008BK: 8.5 Msol, 2.5 Mpc, ApMag 22 (progenitor), ApMag 12.5 (SN), AbMag -15 (SN)
    # fig.suptitle("Supernova, SN with Galaxy, Galaxy")
    for j,arr in enumerate(data):
        norm = ImageNormalize(arr, interval=MinMaxInterval(),stretch=SqrtStretch(), vmin=0, vmax=15000)
        # norm = simple_norm(arr[0:450][0:450], 'sqrt')
        
        im = ax[j].imshow(arr, origin='lower', norm=norm, cmap="magma")
        
        divider = make_axes_locatable(ax[j])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax, orientation="vertical")
    ax[0].title.set_text("Supernova Progenitor")
    ax[1].title.set_text("Progenitor and Star Cluster")
    ax[2].title.set_text("Star Cluster")
    # ax[3].title.set_text("Progenitor Without Empty Sky Subtracted")
    plt.savefig(f"{IMG_DIR}/{filePrim}.png")

if __name__ == "__main__":
    #real measured distance to 2008bk
    REF_DISTANCE = 3.5E6
    PIX_SCALES = [1.5, 4]

    #real measured apparent magnitudes in various bands, for 2008bk
    snmags = {
                'Y' : 20.05,
                'J' : 19.45,
                'Ks' : 18.39
            }
    
    #measure apparent mag -> Absolute magnitude (ref) + distance -> modified apparent magnitude
    filterBandlist = ["Ks"]#, "Y",  "J"] # , 

    #variable distance for simulation
    dlist = [(REF_DISTANCE*3) * i for i in range(0,10)] #definitely want to start higher than 2E6
    dlist[0] = REF_DISTANCE
    totmasslist = [1E5]#, 5E5]
    
    EXPOSURE_TIME_CONST = 4
    NDIT_BASE = 50

    for i,d in enumerate(dlist):
        for totmass in totmasslist:
            for pix_scale in PIX_SCALES:
                for filterBand in filterBandlist:
                    ex_count = 4
                    ndit = 900
                    snmag = snmags[filterBand] + 5 * np.log10(d/REF_DISTANCE)
                    print(snmag)

                    #file naming stuff
                    filePrim = f"scpushy_{filterBand}_ps{pix_scale}_d{d:.2E}_mass{totmass:.1E}_et{ex_count*ndit}"
                    fileFunc = lambda file_index : rf"{DATA_DIR}/{filePrim}_{file_index}.csv"
                    filenames = [fileFunc(file_index) for file_index in range(0,4)]

                    #filter wheel handling
                    fw1 = filterBand
                    fw2 = "open"
                    output_arrays = []
                    if filterBand == "Ks":
                        fw1 = "open"
                        fw2 = filterBand

                    output_arrays = memoize(lambda : image_differencing_simulation_SNCluster(SNmag=snmag,clusterDistance = d,diffx=0,dit=1,ndit=ndit,filter=filterBand,fw1=fw1,fw2=fw2,AOsystem="SCAO",pix_scale=pix_scale,ex_count=ex_count,mass=totmass), filenames, CACHING_ARG)
                                    
                    match output_arrays:
                        case Ok(value):
                            plot_results(value.astype(np.float64), filePrim)
                        case Err(e):
                            print(e)

                    print(process_time())
                    
                    
        