from astropy.io import fits
import numpy as np

image_file = "/pebble_tmp/models/quiescent-galaxy-colour-diagram/HorizonAGN.fits"
# image_file = "../data/HorizonAGN.fits"

image_data = np.array(fits.getdata(image_file))

def return_info():
    shape = image_data.shape
    print(shape)

    return shape

def return_data(layer=0, idx=-1):
    data = image_data[layer][idx]

    return data

if __name__ == "__main__":
    data = return_data()
    shape = return_info()

