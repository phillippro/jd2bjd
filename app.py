"""
JD → BJD Converter Service
High-precision astrometry calculations using Astropy

Consider:
1. Earth orbital motion (barycenter correction)
2. Stellar parallax and proper motion (optional)
3. Light-travel time effect
"""

from flask import Flask, render_template_string, request, jsonify
from astroquery.simbad import Simbad
from astropy.constants import c as const
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, get_body
from astropy import units as u
import math



def lookup_object(query):
    """
    Look up object by name or parse coordinates.
    
    Args:
        query: Object name (e.g., "Tau Ceti") or coordinates "RA, Dec"
        
    Returns:
        dict with 'ra' (degrees), 'dec' (degrees), 'name' (if found)
    """
    # Try to parse as coordinates first
    if ',' in query:
        try:
            parts = query.split(',')
            ra = float(parts[0].strip())
            dec = float(parts[1].strip())
            return {'ra': ra, 'dec': dec, 'name': None}
        except:
            pass
    
    # Try SIMBAD lookup
    try:
        result = Simbad.query_object(query)
        if result is not None and len(result) > 0:
            ra = float(result['ra'][0])
            dec = float(result['dec'][0])
            name = str(result['main_id'][0])
            return {'ra': ra, 'dec': dec, 'name': name}
    except Exception as e:
        raise ValueError(f"Could not find object '{query}': {e}")
    
    raise ValueError(f"Could not find object '{query}'")


app = Flask(__name__)

# Default Observatory (Shanghai Observatory)
DEFAULT_OBSERVATORY = {
    'name': 'Shanghai',
    'lon': 121.54,  # Longitude
    'lat': 31.22,   # Latitude
    'height': 7     # Elevation (meters)
}

# Famous Observatories
OBSERVATORIES = {
    'Shanghai': {'lon': 121.54, 'lat': 31.22, 'height': 7},
    'Mauna Kea': {'lon': -155.48, 'lat': 19.82, 'height': 4207},
    'Paranal': {'lon': -70.40, 'lat': -24.63, 'height': 2635},
    'La Silla': {'lon': -70.73, 'lat': -29.26, 'height': 2400},
    'McDonald': {'lon': -104.02, 'lat': 30.67, 'height': 2070},
    'Siding Spring': {'lon': 149.04, 'lat': -31.27, 'height': 1164},
    'Cerro Tololo': {'lon': -70.81, 'lat': -30.17, 'height': 2215},
}


def jd_to_bjd(jd, ra, dec, obs_name='Shanghai', parallax=0, pmra=0, pmdec=0, rv_star=0):
    """
    Convert JD to BJD
    
    Parameters:
    - jd: Julian Date (can be a single value or list)
    - ra: Right Ascension (degrees or "HH:MM:SS")
    - dec: Declination (degrees or "DD:MM:SS")
    - obs_name: Observatory name
    - parallax: Parallax (mas)
    - pmra: Proper motion RA (mas/yr)
    - pmdec: Proper motion Dec (mas/yr)
    - rv_star: Radial Velocity (km/s)
    
    Returns:
    - BJD_TDB (Barycentric Julian Date，barycentric coordinate system)
    """
    
    # Get Observatory location
    obs = OBSERVATORIES.get(obs_name, DEFAULT_OBSERVATORY)
    location = EarthLocation(lon=obs['lon']*u.deg, 
                            lat=obs['lat']*u.deg, 
                            height=obs['height']*u.m)
    
    # Create target coordinates
    # Numeric RA is interpreted as HOURS (like OSU), string can be HMS or degrees
    # Get RA/dec values first
    if isinstance(ra, str):
        if ':' in ra:
            ra_val = ra
            unit_ra = u.hourangle
        else:
            try:
                ra_val = float(ra)
                unit_ra = u.hourangle if ra_val <= 24 else u.deg
            except:
                ra_val = ra
                unit_ra = u.hourangle
    else:
        ra_val = ra * u.hourangle
        unit_ra = u.hourangle
    
    # Create base coordinate
    coord = SkyCoord(ra=ra_val, dec=dec, unit=(unit_ra, u.deg))
    
    # Add parallax (distance) and proper motion - must create new coord
    if parallax > 0 or pmra != 0 or pmdec != 0:
        distance = (1000/parallax) * u.pc if parallax > 0 else None
        pm_ra = pmra * u.mas/u.yr if pmra != 0 else 0 * u.mas/u.yr
        pm_dec = pmdec * u.mas/u.yr if pmdec != 0 else 0 * u.mas/u.yr
        coord = SkyCoord(ra=coord.ra, dec=coord.dec, 
                        distance=distance,
                        pm_ra_cosdec=pm_ra, pm_dec=pm_dec,
                        frame='icrs')
    
    # Convert to time object
    if isinstance(jd, (list, tuple)):
        t = Time(jd, format='jd', scale='utc', location=location)
    else:
        t = Time(jd, format='jd', scale='utc', location=location)
    
    # Convert to TDB (barycentric dynamical time) - includes Earth orbital correction
    t_tdb = t.tdb
    
    # Calculate light-travel time (LTT)
    # This is due to Earth-barycenter distance change
    ltt = t_tdb.light_travel_time(coord)
    
    # BJD_TDB = TDB + LTT correction
    bjd_tdb = t_tdb + ltt
    
    return {
        'bjd_tdb': bjd_tdb.value,
        'tdb': t_tdb.value,
        'ltt_correction_days': ltt.to(u.day).value,
        'ltt_correction_sec': ltt.to(u.second).value,
    }


def get_earth_position(jd):
    """Get Earth position in barycentric coordinate system"""
    t = Time(jd, format='jd', scale='utc')
    earth = get_body('earth', t)
    # Transform to ICRS (barycentric) and get cartesian coordinates
    earth_bary = earth.transform_to('icrs')
    cart = earth_bary.cartesian
    return {
        'x_au': cart.x.to(u.AU).value,
        'y_au': cart.y.to(u.AU).value,
        'z_au': cart.z.to(u.AU).value,
    }

def get_barycentric_rv(jd, ra_hours, dec_deg):
    """Calculate barycentric radial velocity correction for a given target."""
    
    ra_deg = ra_hours * 15
    dt = 1.0 / 86400  # 1 second
    
    t0 = Time(jd, format='jd', scale='utc')
    t_minus = Time(jd - dt, format='jd', scale='utc')
    t_plus = Time(jd + dt, format='jd', scale='utc')
    
    # Get Earth positions in ICRS (barycentric)
    em = get_body('earth', t_minus).transform_to('icrs').cartesian
    ep = get_body('earth', t_plus).transform_to('icrs').cartesian
    
    # Velocity in km/s
    vx = ((ep.x - em.x) / (2 * dt)).to(u.km/u.s).value
    vy = ((ep.y - em.y) / (2 * dt)).to(u.km/u.s).value
    vz = ((ep.z - em.z) / (2 * dt)).to(u.km/u.s).value
    
    # Line of sight unit vector
    cos_dec = np.cos(np.radians(dec_deg))
    sin_dec = np.sin(np.radians(dec_deg))
    cos_ra = np.cos(np.radians(ra_deg))
    sin_ra = np.sin(np.radians(ra_deg))
    
    lx = cos_dec * cos_ra
    ly = cos_dec * sin_ra
    lz = sin_dec
    
    # Barycentric RV
    rv_bary = vx * lx + vy * ly + vz * lz
    
    # Relativistic correction
    v_mag = np.sqrt(vx**2 + vy**2 + vz**2)
    beta = v_mag / const.c.to(u.km/u.s).value
    rel_corr = 0.5 * beta**2 * rv_bary
    
    return {
        'rv_barycentric_km_s': round(float(rv_bary), 4),
        'relativistic_correction_km_s': round(float(rel_corr), 6),
        'total_correction_km_s': round(float(rv_bary + rel_corr), 4),
        'earth_orbital_velocity_km_s': round(float(v_mag), 4),
    }


if __name__ == '__main__':
    print("Starting JD→BJD Converter Service...")
    print("Please visit: http://localhost:5000")
    app.run(host='0.0.0.0', port=5000, debug=True)