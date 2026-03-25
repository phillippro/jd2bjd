"""
JD → BJD Converter Service
High-precision astrometry calculations using Astropy

Consider:
1. Earth orbital motion (barycenter correction)
2. Stellar parallax and proper motion (optional)
3. Light-travel time effect
"""

from flask import Flask, render_template_string, request, jsonify
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, get_body
from astropy import units as u
import math

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
    return {
        'x_au': cart.x.to(u.AU).value,
        'y_au': cart.y.to(u.AU).value,
        'z_au': cart.z.to(u.AU).value,
    }


@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)


@app.route('/convert', methods=['POST'])
def convert():
    data = request.json
    
    try:
        # Parse input
        jd = float(data.get('jd', 0))
        ra = data.get('ra', '0')
        dec = data.get('dec', '0')
        obs_name = data.get('observatory', 'Shanghai')
        parallax = float(data.get('parallax', 0))
        pmra = float(data.get('pmra', 0))
        pmdec = float(data.get('pmdec', 0))
        
        if jd == 0:
            return jsonify({'error': 'Please enter a valid JD'})
        
        # Convert
        result = jd_to_bjd(jd, ra, dec, obs_name, parallax, pmra, pmdec)
        
        # Get Earth position (for demonstration)
        earth_pos = get_earth_position(jd)
        
        return jsonify({
            'success': True,
            'input': {
                'jd': jd,
                'ra': ra,
                'dec': dec,
                'observatory': obs_name,
                'parallax': parallax,
                'pmra': pmra,
                'pmdec': pmdec,
            },
            'result': {
                'bjd_tdb': result['bjd_tdb'],
                'tdb': result['tdb'],
                'ltt_correction_days': result['ltt_correction_days'],
                'ltt_correction_sec': result['ltt_correction_sec'],
            },
            'earth_position_au': earth_pos,
        })
        
    except Exception as e:
        return jsonify({'error': str(e)})


@app.route('/batch', methods=['POST'])
def batch_convert():
    """Batch convert multiple JD values"""
    data = request.json
    
    try:
        jd_list = data.get('jd_list', [])
        ra = data.get('ra', '0')
        dec = data.get('dec', '0')
        obs_name = data.get('observatory', 'Shanghai')
        
        results = []
        for jd in jd_list:
            result = jd_to_bjd(float(jd), ra, dec, obs_name)
            results.append({
                'jd': jd,
                'bjd_tdb': result['bjd_tdb'],
                'ltt_correction_sec': result['ltt_correction_sec'],
            })
        
        return jsonify({'success': True, 'results': results})
        
    except Exception as e:
        return jsonify({'error': str(e)})


HTML_TEMPLATE = '''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>JD → BJD Converter</title>
    <style>
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
            min-height: 100vh;
            padding: 20px;
            color: #eee;
        }
        .container {
            max-width: 800px;
            margin: 0 auto;
            background: rgba(255,255,255,0.05);
            border-radius: 16px;
            padding: 30px;
            backdrop-filter: blur(10px);
        }
        h1 {
            text-align: center;
            margin-bottom: 10px;
            background: linear-gradient(90deg, #00d4ff, #7b2ff7);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        .subtitle {
            text-align: center;
            color: #888;
            margin-bottom: 30px;
            font-size: 14px;
        }
        .form-group {
            margin-bottom: 20px;
        }
        label {
            display: block;
            margin-bottom: 8px;
            color: #aaa;
            font-size: 14px;
        }
        input, select {
            width: 100%;
            padding: 12px;
            background: rgba(255,255,255,0.1);
            border: 1px solid rgba(255,255,255,0.1);
            border-radius: 8px;
            color: #fff;
            font-size: 16px;
        }
        input:focus, select:focus {
            outline: none;
            border-color: #00d4ff;
        }
        .row {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
        }
        .btn {
            width: 100%;
            padding: 14px;
            background: linear-gradient(90deg, #00d4ff, #7b2ff7);
            border: none;
            border-radius: 8px;
            color: white;
            font-size: 16px;
            font-weight: bold;
            cursor: pointer;
            transition: transform 0.2s;
        }
        .btn:hover {
            transform: scale(1.02);
        }
        .btn:active {
            transform: scale(0.98);
        }
        .result {
            margin-top: 30px;
            padding: 20px;
            background: rgba(0,212,255,0.1);
            border-radius: 12px;
            border: 1px solid rgba(0,212,255,0.3);
            display: none;
        }
        .result.show { display: block; }
        .result h3 { color: #00d4ff; margin-bottom: 15px; }
        .result-item {
            display: flex;
            justify-content: space-between;
            padding: 10px 0;
            border-bottom: 1px solid rgba(255,255,255,0.1);
        }
        .result-item:last-child { border-bottom: none; }
        .result-label { color: #aaa; }
        .result-value { 
            color: #00d4ff; 
            font-family: 'Monaco', monospace;
            font-size: 18px;
        }
        .info {
            margin-top: 30px;
            padding: 20px;
            background: rgba(255,255,255,0.03);
            border-radius: 12px;
            font-size: 13px;
            color: #888;
            line-height: 1.6;
        }
        .info h4 { color: #aaa; margin-bottom: 10px; }
        .loading {
            text-align: center;
            padding: 20px;
            display: none;
        }
        .loading.show { display: block; }
    </style>
</head>
<body>
    <div class="container">
        <h1>🔭 JD → BJD Converter</h1>
        <p class="subtitle">Julian Date (JD) → Barycentric Julian Date (BJD-TDB)</p>
        
        <div class="row">
            <div class="form-group">
                <label>JD (Julian Date)</label>
                <input type="number" id="jd" placeholder="e.g., 2460000.5" step="0.0001">
            </div>
            <div class="form-group">
                <label>Observatory</label>
                <select id="observatory">
                    <option value="Shanghai">Shanghai Observatory</option>
                    <option value="Mauna Kea">Mauna Kea (Hawaii)</option>
                    <option value="Paranal">Paranal (Chile)</option>
                    <option value="La Silla">La Silla (Chile)</option>
                    <option value="McDonald">McDonald (Texas)</option>
                    <option value="Siding Spring">Siding Spring (Australia)</option>
                    <option value="Cerro Tololo">Cerro Tololo (Chile)</option>
                </select>
            </div>
        </div>
        
        <div class="row">
            <div class="form-group">
                <label>Right Ascension RA (degrees or HH:MM:SS)</label>
                <input type="text" id="ra" placeholder="e.g., 83.633 or 05:34:32" value="83.633">
            </div>
            <div class="form-group">
                <label>Declination Dec (degrees or DD:MM:SS)</label>
                <input type="text" id="dec" placeholder="e.g., -5.391 or -05:23:28" value="-5.391">
            </div>
        </div>
        
        <div class="row">
            <div class="form-group">
                <label>Parallax (mas)</label>
                <input type="number" id="parallax" placeholder="e.g., 76.29" step="0.01" value="0">
            </div>
            <div class="form-group">
                <label>Radial Velocity (km/s)</label>
                <input type="number" id="rv" placeholder="e.g., 5.2" step="0.1" value="0">
            </div>
        </div>
        
        <div class="row">
            <div class="form-group">
                <label>Proper Motion RA (mas/yr)</label>
                <input type="number" id="pmra" placeholder="e.g., 3.5" step="0.01" value="0">
            </div>
            <div class="form-group">
                <label>Proper Motion Dec (mas/yr)</label>
                <input type="number" id="pmdec" placeholder="e.g., -5.8" step="0.01" value="0">
            </div>
        </div>
        
        <button class="btn" onclick="convert()">Convert</button>
        
        <div class="loading" id="loading">Calculating...</div>
        
        <div class="result" id="result">
            <h3>🎯 Conversion Results</h3>
            <div class="result-item">
                <span class="result-label">Input JD (UTC)</span>
                <span class="result-value" id="res-jd"></span>
            </div>
            <div class="result-item">
                <span class="result-label">BJD-TDB (Barycenter)</span>
                <span class="result-value" id="res-bjd"></span>
            </div>
            <div class="result-item">
                <span class="result-label">TDB</span>
                <span class="result-value" id="res-tdb"></span>
            </div>
            <div class="result-item">
                <span class="result-label">Light-travel Time (days)</span>
                <span class="result-value" id="res-ltt-d"></span>
            </div>
            <div class="result-item">
                <span class="result-label">Light-travel Time (sec)</span>
                <span class="result-value" id="res-ltt-s"></span>
            </div>
            <div class="result-item">
                <span class="result-label">Earth Position X (AU)</span>
                <span class="result-value" id="res-ex"></span>
            </div>
            <div class="result-item">
                <span class="result-label">Earth Position Y (AU)</span>
                <span class="result-value" id="res-ey"></span>
            </div>
            <div class="result-item">
                <span class="result-label">Earth Position Z (AU)</span>
                <span class="result-value" id="res-ez"></span>
            </div>
        </div>
        
        <div class="info">
            <h4>📖 Explanation</h4>
            <p>• <b>JD</b>: Julian Date, days since noon January 1, 4713 BC</p>
            <p>• <b>BJD-TDB</b>: Barycentric Julian Date, referenced to solar system barycenter</p>
            <p>• <b>Light-travel Time Correction</b>: Time difference due to Earth-barycenter distance, ~±16 minutes</p>
            <p>• <b>Parallax</b>: Enter parallax for more accurate stellar distance correction</p>
            <p>• <b>Proper Motion</b>: Stellar proper motion, important for observations spanning long periods</p>
        </div>
    </div>
    
    <script>
        async function convert() {
            const jd = document.getElementById('jd').value;
            if (!jd) {
                alert('Please enter JD');
                return;
            }
            
            document.getElementById('loading').classList.add('show');
            document.getElementById('result').classList.remove('show');
            
            try {
                const response = await fetch('/convert', {
                    method: 'POST',
                    headers: {'Content-Type': 'application/json'},
                    body: JSON.stringify({
                        jd: parseFloat(jd),
                        ra: document.getElementById('ra').value,
                        dec: document.getElementById('dec').value,
                        observatory: document.getElementById('observatory').value,
                        parallax: parseFloat(document.getElementById('parallax').value) || 0,
                        pmra: parseFloat(document.getElementById('pmra').value) || 0,
                        pmdec: parseFloat(document.getElementById('pmdec').value) || 0,
                    })
                });
                
                const data = await response.json();
                
                if (data.error) {
                    alert(data.error);
                    return;
                }
                
                document.getElementById('res-jd').textContent = data.input.jd.toFixed(6);
                document.getElementById('res-bjd').textContent = data.result.bjd_tdb.toFixed(6);
                document.getElementById('res-tdb').textContent = data.result.tdb.toFixed(6);
                document.getElementById('res-ltt-d').textContent = data.result.ltt_correction_days.toFixed(6);
                document.getElementById('res-ltt-s').textContent = data.result.ltt_correction_sec.toFixed(3);
                document.getElementById('res-ex').textContent = data.earth_position_au.x_au.toFixed(4);
                document.getElementById('res-ey').textContent = data.earth_position_au.y_au.toFixed(4);
                document.getElementById('res-ez').textContent = data.earth_position_au.z_au.toFixed(4);
                
                document.getElementById('result').classList.add('show');
                
            } catch (e) {
                alert('Error: ' + e);
            } finally {
                document.getElementById('loading').classList.remove('show');
            }
        }
        
        // Default to current time
        document.getElementById('jd').value = Math.floor(Date.now()/86400000 + 2440587.5);
    </script>
</body>
</html>
'''


if __name__ == '__main__':
    print("Starting JD→BJD Converter Service...")
    print("Please visit: http://localhost:5000")
    app.run(host='0.0.0.0', port=5000, debug=True)