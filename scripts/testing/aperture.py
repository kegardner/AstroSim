import scopesim as sc
aperture_dict = {"x":["arcsec",1.,2.],"y":["arcsec",1.,2.],"id":1,"conserve_image":True}
mask = sc.effects.apertures.ApertureMask(array_dict=aperture_dict)
# data = []
# colnames = []
# for key, val in aperture_dict.items():
#             data += [val]
#             colnames += [key]
print(mask)