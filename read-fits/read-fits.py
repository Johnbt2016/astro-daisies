from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
from io import BytesIO
from summary import *
from PIL import Image
image_data = ... # byte values of the image
image = Image.open(io.BytesIO(image_data))

image_file = "/pebble_tmp/models/quiescent-galaxy-colour-diagram/HorizonAGN.fits"
# image_file = "../data/HorizonAGN.fits"

image_data = np.array(fits.getdata(image_file))

def spectra_hist(spectra):
    fig= plt.figure()
    if len(spectra) > 1:
        title = "Comparing spectra steps distribution for galaxies "
    else:
        title = "Spectra steps distribution for galaxy "
    nb = 0.
    for key, item in spectra.items():
        if item is not None:
            plt.hist(item,bins = 50, alpha= 1. - nb/3, rwidth=1)
            title += str(key) + ", "
            nb += 1
    title = title.strip().strip(",")
    fig.suptitle(title, fontsize = 14 , fontweight = 'bold')
    return fig

def return_spectra(galaxy_idx=0):
    '''
    Return the spectra for a given galaxy

    Parameters :
    - galaxy_idx (int) : index of the galaxy to query

    Returns :
    - an array with the 1132 wavelength steps from 506 to 1e+05 Angstroms
    '''
    if galaxy_idx < 359075:
        data = image_data[galaxy_idx][-1]
        return data
    else:
        return None

def return_hist(indices):
    spectra = dict()
    for i in indices:
        spectra[i] = return_spectra(galaxy_idx=i)
    fig = spectra_hist(spectra)

    return fig

def st_ui():
    st.title("Querying galaxies spectra")
    st.subheader("From the HorizonAGN_LAIGLE-DAVIDZON+2019_SPECTRA_z2-4.fits catalog")
    g1 = int(st.sidebar.text_input("Galaxy 1 - Index (max = 359074)", 456))
    g2 = int(st.sidebar.text_input("Galaxy 2 - Index (max = 359074)", 67845))

    indices = [g1, g2]
    for i in indices:
        if i > 359074:
            st.error(f"Galaxy {i} is not in the catalog")
    st.markdown(get_summary())

    fig = return_hist(indices)

    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches='tight', transparent = True)
    st.image(buf, use_column_width=False)

if __name__ == "__main__":
    st_ui()

