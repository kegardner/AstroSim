



import scopesim as sim
from scopesim.effects import AnisocadoConstPSF


psf = AnisocadoConstPSF(filename="test_AnisoCADO_rms_map.fits",
                        strehl=0.5,
                        wavelength=2.15,
                        convolve_mode="same",
                        psf_side_length=512)