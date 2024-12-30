import scopesim as sc
wcs = {
    "CDELT1" : 1
    , "CRVAL1" : 1
    , "CRPIX1" : 1
    , "NAXIS1" : 1
    , "NAXIS2" : 2
    , "CDELT1D" : 1
    , "CRVAL1D" : 1
    , "CRPIX1D" : 1
    , "CUNIT1" : "arcsec"
    ,"CDELT2" : 2
    , "CRVAL2" : 2
    , "CRPIX2" : 2
    , "CUNIT2" : "arcsec"
}
field = sc.optics.fov.FieldOfView(header=wcs, waverange=[1,2])

print(field.hdu)

from astropy.table import Table
import astropy.units as u
from scopesim.optics import image_plane as imp

my_point_source_table = Table(names=["x", "y", "flux"],
                              data=[(0,  1,  2)*u.mm,
                                    (0, -5, 10)*u.mm,
                                    (100,50,25)*u.ph/u.s])
# imp.
hdr = imp.make_image_plane_header([my_point_source_table],
                                   pixel_size=0.015*u.mm)
img_plane = imp.ImagePlane(hdr)
img_plane.add(my_point_source_table)

print(img_plane.image)