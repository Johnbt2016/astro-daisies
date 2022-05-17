import streamlit as st
from io import BytesIO

import astropy.coordinates as coord
from astropy.table import QTable
import astropy.units as u
from astroquery.gaia import Gaia

# Third-party imports
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import RendererAgg

import numpy as np

# gala imports
import gala.coordinates as gc
import gala.dynamics as gd
import gala.potential as gp
from gala.units import galactic

import pickle

_lock = RendererAgg.lock


def summary(nb = 4096, parallax_over_error = 8, parallax = 10, BP_RP_range = [0.5, 0.7], M_G_range = [2.0, 3.75], time_step = 1., time = 500.):
    text = '''
    Return orbits of stars with data from the Gaia mission.
    The query gets sky positions, distances (parallaxes), proper motions, and radial velocities 
    for a set of stars that are close to the Sun. 
    The observed, heliocentric kinematic measurements are transforme dto Galactocentric Cartesian coordinates 
    and we use the positions and velocities as initial conditions to compute the orbits of these stars in the galaxy using the gala Python package.
    
    You can call it programatically (The code snippet below is up to date with the current app values, 
    uses the High Mass stars definition range !):  
    ```python
    import pydaisi as pyd

    orbits_query = pyd.Daisi("Gaia 2 Orbits Query")
    result = orbits_query.get_orbits(nb = '''
    text += str(nb) + ", "
    text += "parallax_over_error = " + str(parallax_over_error) + ", "
    text += "parallax = " + str(parallax) + ", " 
    text += "BP_RP_range = " + str(BP_RP_range) + ", "
    text += "M_G_range = " + str(M_G_range) + ", "
    text += "time_step = " + str(time_step) + ", "
    text += "time = " + str(time) + ").value"	
    '''```

    '''
    return text

def query(nb = 4096, parallax_over_error = 8, parallax = 10):

    print("Querying Data")

    query_text = '''SELECT TOP '''
    query_text += str(int(nb))
    query_text += ''' ra, dec, parallax, pmra, pmdec, radial_velocity,
    phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag
    FROM gaiadr2.gaia_source
    WHERE parallax_over_error > '''
    query_text += str(int(parallax_over_error))
    query_text += ''' AND
        parallax > '''
    query_text += str(int(parallax))
    query_text += ''' AND
        radial_velocity IS NOT null
    ORDER BY random_index
    '''

    job = Gaia.launch_job(query_text)
    gaia_data = job.get_results()
    
    return gaia_data


def define_mask(galcen, BP_RP, M_G, BP_RP_high, M_G_high, BP_RP_low, M_G_low):
    np.seterr(invalid="ignore")
    hi_mass_mask = ((BP_RP*u.mag > float(BP_RP_high[0])*u.mag) & (BP_RP*u.mag < float(BP_RP_high[1])*u.mag) & 
                    (M_G*u.mag > float(M_G_high[0])*u.mag) & (M_G*u.mag < float(M_G_high[1])*u.mag) & 
                    (np.abs(galcen.v_y - 220*u.km/u.s) < 50*u.km/u.s))

    lo_mass_mask = ((BP_RP*u.mag > BP_RP_low[0]*u.mag) & (BP_RP*u.mag < BP_RP_low[1]*u.mag) & 
                    (M_G*u.mag > M_G_low[0]*u.mag) & (M_G*u.mag < M_G_low[1]*u.mag) &
                    (np.abs(galcen.v_y - 220*u.km/u.s) < 50*u.km/u.s))
    
    return hi_mass_mask, lo_mass_mask

def data_preparation(gaia_data, distance_sun_gc = 8.1):
    dist = coord.Distance(parallax=u.Quantity(gaia_data['parallax']))
    c = coord.SkyCoord(ra=gaia_data['ra'], 
                   dec=gaia_data['dec'],
                   distance=dist,
                   pm_ra_cosdec=gaia_data['pmra'], 
                   pm_dec=gaia_data['pmdec'],
                   radial_velocity=gaia_data['radial_velocity'])
    galcen = c.transform_to(coord.Galactocentric(z_sun=0*u.pc, galcen_distance=distance_sun_gc*u.kpc))

    M_G = gaia_data['phot_g_mean_mag'] - dist.distmod
    BP_RP = gaia_data['phot_bp_mean_mag'] - gaia_data['phot_rp_mean_mag']

    return galcen, BP_RP, M_G
    
def integrate_orbits(galcen, hi_mass_mask, lo_mass_mask, time_step, time):
    milky_way = gp.MilkyWayPotential()
    different_disk_potential = gp.MilkyWayPotential(disk=dict(m=8e10*u.Msun))
    H = gp.Hamiltonian(milky_way)
    w0_hi = gd.PhaseSpacePosition(galcen[hi_mass_mask].cartesian)
    w0_lo = gd.PhaseSpacePosition(galcen[lo_mass_mask].cartesian)

    orbits_hi = H.integrate_orbit(w0_hi, dt=time_step*u.Myr, t1=0*u.Myr, t2=time*u.Myr)

    orbits_lo = H.integrate_orbit(w0_lo, dt=time_step*u.Myr, t1=0*u.Myr, t2=time*u.Myr)

    return orbits_hi, orbits_lo

def get_orbits(nb = 4096, parallax_over_error = 8, parallax = 10, BP_RP_range = [0.5, 0.7], M_G_range = [2.0, 3.75], time_step = 1., time = 500.):
    '''
    Return orbits of stars with data from the Gaia mission.
    The query gets sky positions, distances (parallaxes), proper motions, and radial velocities 
    for a set of stars that are close to the Sun. 
    The observed, heliocentric kinematic measurements are transforme dto Galactocentric Cartesian coordinates 
    and we use the positions and velocities as initial conditions to compute the orbits of these stars in the galaxy using the gala Python package.

    Parameters:
    - nb (int): sample size to be queried. Default = 4096
    - parallax_over_error (int) : minimum threshold for query. Default = 8
    - parallax (int) : minimum threshold for query. Default = 10
    - BP_RP_range (list of two floats) : Blue - Red diff bandpass magnitude filter. Default = [0.5, 0.7]
    - M_G_range (list of two floats) : Broadband bandpass magnitude filter. Default = [2.0, 3.75]
    - time_step (float) : time step for orbits integration. Default = 1 Myr
    - time (float) : duration for orbits integration. Default = 500 Myr

    Returns:
    - a tuple (orbit.pos, orbit.vel)
    '''
    gaia_data = query(nb, parallax_over_error, parallax)
    galcen, BP_RP, M_G = data_preparation(gaia_data)

    mask = ((BP_RP*u.mag > float(BP_RP_range[0])*u.mag) & (BP_RP*u.mag < float(BP_RP_range[1])*u.mag) & 
                    (M_G*u.mag > float(M_G_range[0])*u.mag) & (M_G*u.mag < float(M_G_range[1])*u.mag) & 
                    (np.abs(galcen.v_y - 220*u.km/u.s) < 50*u.km/u.s))
    
    milky_way = gp.MilkyWayPotential()

    H = gp.Hamiltonian(milky_way)
    w0_hi = gd.PhaseSpacePosition(galcen[mask].cartesian)

    orbits = H.integrate_orbit(w0_hi, dt=time_step*u.Myr, t1=0*u.Myr, t2=time*u.Myr)

    return pickle.dumps(orbits.pos), pickle.dumps(orbits.vel)

def st_ui():
    st.set_page_config(layout = "wide")

    st.title("Computing Galactic Orbits of Stars with Gala")
    st.sidebar.header("Data Collection")
    st.sidebar.write("Allow a few seconds to update when you change a parameter in this section")

    if 'query' not in st.session_state:
        st.session_state.query = [0,0,0]
    
    if 'gaia_data' not in st.session_state:
        st.session_state.gaia_data = 0
    
    nb = st.sidebar.slider("Nb of samples", 1000, 10000, 4096)
    parallax_over_error = st.sidebar.slider("Parallax over error threshold", 5, 20, 8)
    parallax = st.sidebar.slider("Parallax", 5, 20, 10)

    query_data = [nb, parallax_over_error, parallax]

    if query_data != st.session_state.query:
        # st.info("Querying new data, please wait")
        st.session_state.gaia_data = query(query_data[0], query_data[1], query_data[2])
        st.session_state.query = query_data

    gaia_data = st.session_state.gaia_data
    st.sidebar.header("Data Preparation and Filtering")
    distance_sun_gc = st.sidebar.slider('Distance Sun to Galaxy center (kpc)', 7.0, 9.0, 8.1)
    
    
    BP_RP_high = st.sidebar.slider('High Mass stars - Blue / Red bandpass diff BP_RP filter', 0.0, 2.0, (0.5, 0.7))
    M_G_high = st.sidebar.slider('High Mass stars - Broadband bandpass magnitude M_G filter', 0.0, 10.0, (2.0, 3.75))

    BP_RP_low = st.sidebar.slider('Low Mass stars - Blue / Red bandpass diff BP_RP filter', 0.0, 3.0, (2.0, 2.4))
    M_G_low = st.sidebar.slider('Low Mass stars - Broadband bandpass magnitude M_G filter', 0.0, 10.0, (8.2, 9.7))

    st.sidebar.header("Orbits Integration")
    time_step = st.sidebar.slider('Time step for orbits integration (Myr)', 0.1, 3.0, 1.0)
    time = st.sidebar.slider('Time duration for orbits integration (Myr)', 100, 1000, 500)

    
    galcen, BP_RP, M_G = data_preparation(gaia_data, distance_sun_gc)
    hi_mass_mask, lo_mass_mask = define_mask(galcen, BP_RP, M_G, BP_RP_high, M_G_high, BP_RP_low, M_G_low)

    with st.expander("Summary"):
        st.markdown(summary(nb = nb, 
                            parallax_over_error = parallax_over_error, 
                            parallax = parallax, 
                            BP_RP_range = BP_RP_high, 
                            M_G_range = M_G_high, 
                            time_step = time_step, 
                            time = time))

    hi_mass_color = 'tab:red'
    lo_mass_color = 'tab:purple'

    with _lock:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))

        ax.plot(BP_RP, M_G, 
                marker='.', linestyle='none', alpha=0.1)

        for mask, color in zip([lo_mass_mask, hi_mass_mask],
                            [lo_mass_color, hi_mass_color]):
            ax.plot(BP_RP[mask], M_G[mask], 
                    marker='.', linestyle='none', 
                    alpha=0.5, color=color)

        ax.set_xlim(0, 3)
        ax.set_ylim(11, 1)

        ax.set_xlabel('$G_{BP}-G_{RP}$')
        ax.set_ylabel('$M_G$')

        buf = BytesIO()
        fig.savefig(buf, format="png", bbox_inches='tight', transparent = True)
        st.image(buf, use_column_width=False, caption=' two CMD selections')
    
    orbits_hi, orbits_lo = integrate_orbits(galcen, hi_mass_mask, lo_mass_mask, time_step, time)

    with _lock:
        fig = orbits_hi[:, 0].plot(color=hi_mass_color)
        _ = orbits_lo[:, 0].plot(axes=fig.axes, color=lo_mass_color)

        buf = BytesIO()
        fig.savefig(buf, format="png", bbox_inches='tight', transparent = True)
        st.image(buf, use_column_width=False, caption=' orbits')

    with _lock:
        fig = orbits_hi[:, 0].plot(['x', 'v_x'], 
                           auto_aspect=False, 
                           color=hi_mass_color)
        buf = BytesIO()
        fig.savefig(buf, format="png", bbox_inches='tight', transparent = True)
        st.image(buf, use_column_width=False, caption=' orbits')
    
    with _lock:
        fig = orbits_hi[:, 0].cylindrical.plot(['rho', 'z'], 
                                       color=hi_mass_color,
                                       label='high mass')
        _ = orbits_lo[:, 0].cylindrical.plot(['rho', 'z'], color=lo_mass_color,
                                     axes=fig.axes,
                                     label='low mass')

        fig.axes[0].legend(loc='upper left')
        fig.axes[0].set_ylim(-0.3, 0.3)

        buf = BytesIO()
        fig.savefig(buf, format="png", bbox_inches='tight', transparent = True)
        st.image(buf, use_column_width=False, caption=' orbits')





if __name__ == "__main__":
    st_ui()