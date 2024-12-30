import csv
from time import process_time
from typing import List
import numpy as np
import GalaxySN as gsn
import fitsFunctions as ff
import SNR as snr
import Strehl as st
import matplotlib.pyplot as plt
from scopesim_templates.misc import point_source
from scopesim_templates.extragalactic import elliptical
from matplotlib.ticker import LogFormatterExponent
import scopesim as sim
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import os

import sys
from dataclasses import dataclass, field

import DbDriver as sql

H0 = 69.6*u.km/u.s/u.Mpc

O_vac = 0.714

O_mat = 0.286

COSMO = FlatLambdaCDM(H0=69.6, Om0=O_mat)
DATA_DIR = "./dataOutputs/absmag"

@dataclass
class PhotometricBandRedshift:
    from_band: str
    to_band: str
    absolute_magnitude: float
    redshift: float

    @staticmethod
    def from_db_result(t):
        return PhotometricBandRedshift(t[0], t[1], float(t[2]), float(t[3]))
@dataclass
class FilterBandCosmology:
    M : float #absolute magnitude
    z : float #the redshift to 
    pbr: PhotometricBandRedshift
    d : float = field(init=False) #distance associated with z
    ad : float = field(init=False) #pc per mas for -band
    m : float = field(init=False) #apparent mag 
    

    def __post_init__(self):
        self.d = COSMO.luminosity_distance(self.z).to(u.pc).value
        self.ad = 1./COSMO.arcsec_per_kpc_proper(self.z).value
        self.m = 5*np.log10(self.d/10)+self.M 

    @staticmethod
    def from_PBR(pbr : list[PhotometricBandRedshift]):
        return FilterBandCosmology(pbr.absolute_magnitude, pbr.redshift, pbr)

@dataclass
class GalacticReductionResult:
    """
    container for the results of a galaxy-point-source reduction simulation (the result of a single gal_SN_galsubtract call)
    """
    point: list[list[float]]
    raw : list[list[float]]
    galaxy : list[list[float]]
    fbc : FilterBandCosmology

@dataclass
class GalacticReductionAnalysisResult:
    """
    container for the results of a galaxy-point-source reduction simulation and analysis (over a series of FilterBandCosmologies)
    """
    analysis: list[list[float]]
    data: list[GalacticReductionResult]

# KS = [('I',-18,1.7 ), ('R',-19,2.28 ), ('V',-19,2.97 ), ('B',-19,3.94), ('U',-19.2 ,5.24), ('uvw1',-18,7.31), ('uvw2',-16,10.1)]
# J  = [('I',-18,0.55), ('R',-19, .884), ('V',-19,1.27 ), ('B',-19,1.84), ('U',-19.2 ,2.58), ('uvw1',-18,3.77), ('uvw2',-16, 5.39)]
# Y  = [('I',-18, .29), ('R',-19, .568), ('V',-19, .897), ('B',-19,1.36), ('U',-19.2 ,1.98), ('uvw1',-18,2.96), ('uvw2',-16, 4.31)]
@dataclass
class absMag:
    """ 
    
    """
    M : float #absolute magnitude
    
    zKs : float #the redshift to Ks
    zJ : float #the redshift to J-band    
    zY : float #the redshift to Y-band 

    
    #luminosity distances
    dKs : float = field(init=False) #distance associated with zKs
    dJ : float = field(init=False) #distance associated with zJ
    dY : float = field(init=False) #distance associated with zY


    #angular scales
    adKs : float = field(init=False) #pc per mas for Ks-band
    adJ : float = field(init=False) #pc per mas for J-band
    adY : float = field(init=False) #pc per mas for Y-band

    #apparent magnitudes
    mKs :float = field(init=False) #apparent mag Ks
    mJ :float = field(init=False) #apparent mag J
    mY :float = field(init=False) #apparent mag Y

    def __post_init__(self):
        
        self.dKs = COSMO.luminosity_distance(self.zKs) * 1e6
        self.dJ = COSMO.luminosity_distance(self.zJ) * 1e6
        self.dY = COSMO.luminosity_distance(self.zY) * 1e6

        self.adKs = 1./COSMO.arcsec_per_kpc_proper(self.zKs)
        self.adJ = 1./COSMO.arcsec_per_kpc_proper(self.zJ)
        self.adY = 1./COSMO.arcsec_per_kpc_proper(self.zY)

        self.mKs = 5*np.log10(self.dKs/10)+self.M 
        self.mJ = 5*np.log10(self.dJ/10)+self.M
        self.mY = 5*np.log10(self.dY/10)+self.M

class absMagLong:

    def __init__():
        pass


# @dataclass
# class PhotometricBandRedshifts:
#     I : tuple[float] # M, z
#     R : tuple[float]
#     V : tuple[float]
#     B : tuple[float]
#     U : tuple[float]
#     uvw1 : tuple[float]
#     uvw2 : tuple[float]

#     def as_list(self):
#         return [self.I, self.R, self.V, self.B, self.U, self.uvw1, self.uvw2]




def plot_3(results, filter_band):
    fig, ax = plt.subplots(3,1, figsize=(8,30), sharex = True)
    fig.subplots_adjust(hspace=0)
    fig.suptitle("Strehl Ratio, SNR, FWHM vs Redshift")
    ax[0].grid(True, linestyle='--', color='gray', alpha=0.3)
    ax[1].grid(True, linestyle='--', color='gray', alpha=0.3)
    ax[2].grid(True, linestyle='--', color='gray', alpha=0.3)
    ax[2].set_xlabel("Redshift")
    tres = np.array(results).astype(np.float32).T
    z, strehl = np.array(list(filter(lambda pair : 0 < pair[1] < 1, tres[(0,1),...].T))).T
    ax[0].scatter(z, strehl)
    ax[0].set_ylabel("SR")
    z, snr = tres[(0,2),...]
    ax[1].scatter(z, snr)
    ax[1].set_ylabel("S/N")
    z, fwhm = np.array(list(filter(lambda pair : pair[1] != -69, tres[(0,3),...].T))).T
    ax[2].scatter(z, fwhm)
    ax[2].set_ylabel("FWHM (px)")

    plt.savefig("../thesis_images/AbsMag_study.png")
    print(f"Run time: {process_time()}")

def plot_snr(results, filter_band, pixel_scales, galactic_displacements, file_disc = ""):
    """
    filter_band - like "Y", "Ks", etc.
    pixel_scale
    """
    fig, ax = plt.subplots(1,len(galactic_displacements), figsize=(15,6), sharey=True)
    fig.suptitle(f"SNR vs Redshift - {filter_band}-Band")
    fig.subplots_adjust(wspace=0)
    ax[0].set_ylabel("S/N")
    results = np.array(results).astype(np.float32)
    z = results[0][(0),...].T
    
    for i, gd in enumerate(galactic_displacements):
        ax[i].grid(True, linestyle='--', color='gray', alpha=0.3)
        ax[i].set_title(f"Galactic Displacement (pc): {gd}")
        
        for j, ps in enumerate(pixel_scales): 
            tres = results[i*len(pixel_scales) + j].T
            z, snr = tres[(0,2),...]
            ax[i].scatter(z, snr, label=f"{ps}", alpha = 0.6)
            
            #ax[i].set_yscale('log')
            #ax[i].yaxis.set_major_formatter(LogFormatterExponent())
            #TODO!: what works for one doesn't work for another...
            # ax[i].set_xlim(0, max(z)+1)
            try:
                strehl, fwhm = np.array(list(filter(lambda pair : pair[1] != -69, tres[(1,3),...].T))).T 
                print(f"Full Width at Half Maximum: {fwhm}, ps={ps},gd={gd},filter={filter_band}")
                print(f"Strehl: {strehl}, mu {np.mean(strehl)} | sigma {np.std(strehl)}, ps={ps},gd={gd},filter={filter_band}")
            except ValueError:
                print(f"No valid FWHM values: {tres[(1,3),...].T}")
        # plt.figtext(0.9, 0.9, f"Galactic Displacement (pc): {galactic_displacement: <20}\nPixel Scale (mas): {pixel_scale: <4}", wrap=True, horizontalalignment='right', fontsize=12)
    ax[int((len(galactic_displacements)-1)/2)].set_xlabel("Redshift")
    plt.legend(title="Pixel Scale (mas)")
    plt.savefig(f"../thesis_images/AbsMag_snr_study{file_disc}_{filter_band}.png")
    print(f"Run time: {process_time()}")

def band_snr_redshift_data(filter_band, fbc : list[FilterBandCosmology], pixel_scale, galactical_displacement) -> list[GalacticReductionResult]:
    fw1 = filter_band
    fw2 = "open"

    if filter_band == "Ks":
        fw1 = "open"
        fw2 = filter_band

    data_results = []
    for X in fbc:
        print(f"Processing for: {X.pbr.from_band} - {X.pbr.to_band}")
        point,pointgal,gal = gsn.gal_SN_galsubtract(SNmag=X.m,galmag=X.m,galr=600/(1000*X.ad),diffx=galactical_displacement/(1000*X.ad),
        dit=1,ndit=900,filter=filter_band,fw1=fw1,fw2=fw2,AOsystem="SCAO",pix_scale=pixel_scale)

        data_results.append(GalacticReductionResult(point,pointgal,gal, X))
    return data_results

def band_snr_redshift_analysis(data: list[GalacticReductionResult], pixel_scale, filter_band) -> GalacticReductionAnalysisResult:
    
    analysis_results = []
    for res in data:
        strehl,_,_ = st.calculate_strehl_from_array(array=res.point, pixelscale=pixel_scale, photofilter=filter_band)
        SNratio = snr.get_SNR3(res.point,filter = filter_band,ps = pixel_scale)
        fwhmval = ff.FWHM(res.point[int(len(res.point)/2)])
        analysis_results.append((res.fbc.z,strehl,SNratio,fwhmval if fwhmval != 'none' else -69))

    return GalacticReductionAnalysisResult(analysis_results, data)

if __name__ == "__main__":
    print(sim.rc.__search_path__)
    # FlatLambdaCDM().
    
    #filter_band = "Y"
    galactical_displacements = [500, 1000, 2000]
    fbcs = {
        'Ks' : None, #[fbc_from_PBR(i) for i in zKs],
        'J'  : None, #[fbc_from_PBR(i) for i in zJs],
        'Y'  : None, #[fbc_from_PBR(i) for i in zYs],
    }

    with sql.connect() as conn:
        cur = conn.cursor()
        query = "select * from fb_redshift where to_band = :to_band order by redshift"
        for i in fbcs.keys():
            cur.execute(query,{"to_band" : i})
            brs = map(PhotometricBandRedshift.from_db_result, cur.fetchall())
            fbc = map(FilterBandCosmology.from_PBR, brs)
            fbcs[i] = list(fbc)

    pixel_scales = [1.5, 4]
    
    print(fbcs)
    total_runs = len(galactical_displacements) * len(pixel_scales)
    for key in fbcs.keys():
        results : list[GalacticReductionResult] = []
        file = rf"data_{key}"
        for i,gd in enumerate(galactical_displacements):
            for j,ps in enumerate(pixel_scales):
                print(f"-----{i*len(pixel_scales) + j + 1}/{total_runs}----- ({key}, {gd}, {ps})")
                file_part = f"abs_mag_{file}_{ps}_{gd}" 
                file_discs = ["analysis2", "reduced", "raw", "galaxy"] 
                file_func = lambda disc: f"{DATA_DIR}/{file_part}_{disc}.csv"
                #if the analysis file is present, we will assume the rest are too
                print(f"Determining Cache State for {file_func(file_discs[0])}")

                #get the galsubtract data (cached or otherwise)
                data = band_snr_redshift_data(key, fbcs[key], ps, gd)
                if ((not os.path.isfile(file_func(file_discs[0]))) or sys.argv.count("--no-cache-analysis") > 0):
                    
                    #analysis not present in cache, calculate from data:
                    result : GalacticReductionResult = band_snr_redshift_analysis(data, ps, key)
                    

                    with open(file_func(file_discs[0]), "w+") as f :
                        w = csv.writer(f, delimiter="|",lineterminator="\n")
                        for row in result.analysis:
                            w.writerow(row)
                    results.append(result)
                else:
                    #the analysis is cached, read it into memory:
                    results_temp = []
                    with open(file_func(file_discs[0]), "r") as f :
                        reader = csv.reader(f, delimiter="|",lineterminator="\n")
                        result = []
                        for row in reader:
                            result.append(row)
                        results.append(GalacticReductionAnalysisResult(result, data))
        plot_snr(list(map(lambda r: r.analysis, results)), key, pixel_scales, galactical_displacements, file_disc="_snrNoiseCalcChange2tallplot")
    print(f"Run time: {process_time()}")
    # plot_3(results, filter_band)
    plt.show()