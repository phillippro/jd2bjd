"""
Validation script: Compare our JD→BJD converter with Astropy (TEMPO2 algorithm)
"""

from app import jd_to_bjd
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u

def validate():
    """Run validation tests"""
    test_cases = [
        ("Tau Ceti", 1.73447, -15.9375, 768.5, 1726, 981),
        ("Proxima Centauri", 14.48979, -62.6784, 768.5, -3777, 765),
        ("Barnard's Star", 17.962, +04.693, 546.9, -798, 10328),
        ("Sirius", 6.752, -16.716, 379.21, -546, -1223),
        ("Vega", 18.616, +38.784, 130.23, 201, 287),
    ]
    
    jd = 2450000.123
    location = EarthLocation(lat=19.82, lon=-155.48, height=4207*u.m)
    
    print("=== Validation: Our App vs Astropy (TEMPO2) ===\n")
    print(f"JD: {jd}, Observatory: Mauna Kea\n")
    
    max_diff = 0
    for name, ra_h, dec_d, plx, pmra, pmdec in test_cases:
        our = jd_to_bjd(jd, ra_h, dec_d, "Mauna Kea", parallax=plx, pmra=pmra, pmdec=pmdec)
        
        t = Time(jd, format='jd', scale='utc', location=location)
        coord = SkyCoord(ra=ra_h*u.hourangle, dec=dec_d*u.deg,
                        distance=(1000/plx)*u.pc if plx > 0 else None,
                        pm_ra_cosdec=pmra*u.mas/u.yr if pmra != 0 else 0*u.mas/u.yr,
                        pm_dec=pmdec*u.mas/u.yr if pmdec != 0 else 0*u.mas/u.yr)
        ltt = t.tdb.light_travel_time(coord)
        bjd_astropy = t.tdb.value + ltt.to(u.day).value
        
        diff_ms = abs(our['bjd_tdb'] - bjd_astropy) * 86400 * 1000
        max_diff = max(max_diff, diff_ms)
        print(f"{name}: diff = {diff_ms:.3f} ms")
    
    print(f"\n✓ Max difference: {max_diff:.3f} ms")
    return max_diff < 1.0  # Pass if < 1 ms

if __name__ == "__main__":
    validate()
