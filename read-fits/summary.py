def get_summary():
    text = '''
    This Daisi query the **HorizonAGN_LAIGLE-DAVIDZON+2019_SPECTRA_z2-4.fits** catalog 
    which contains spectra at 1132 wavelength steps from 506 to 1e+05 Angstroms,
    for all the galaxies in **HorizonAGN_LAIGLE-DAVIDZON+2019_COSMOS_v1.5.fits**.

    The computation of spectra is described in Laigle, Davidzon, Ilbert+2019
    and Davidzon+2019. They are extracted from the Horizon-AGN simulation
    described in Dubois et al. 2014. Please cite these 3 papers if you use these spectra.
    Would you have any question on these spectra, please contact:
    clotilde.laigle@iap.fr and iary.davidzon@nbi.ku.dk.

    You can query the spectra for any galaxy index or obtain this histogram programmatically in Python:

    ```python
    import pydaisi as pyd

    query_spectra = pyd.Daisi("Query Galaxies Spectra")

     # Query the spectra for a galaxy:
    spectra = query_spectra.return_spectra(galaxy_idx=0).value

    # Get the histogram of the distribution:
    indices = [1223, 65434]
    hist = query_spectra.return_hist(indices).value
    ```
    
    '''

    return text