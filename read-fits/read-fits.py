
from astropy.io import fits

image_file = "/pebble_tmp/models/quiescent-galaxy-colour-diagram/HorizonAGN.fits"
# image_file = "../data/HorizonAGN.fits"
info = fits.info(image_file)
image_data = fits.getdata(image_file, ext=0)

def return_info():
    
    shape = image_data.shape

    return info, shape

if __name__ == "__main__":
    info, shape = return_info()
    print(shape)