
from astropy.io import fits

image_file = "/pebble_tmp/models/quiescent-galaxy-colour-diagram/HorizonAGN_LAIGLE-DAVIDZON+2019_SPECTRA_z2-4.fits"
info = fits.info(image_file)
image_data = fits.getdata(image_file, ext=0)

def return_info():
    
    shape = image_data.shape

    return info, shape
