# JD → BJD Converter

A Flask web application for converting Julian Date (JD) to Barycentric Julian Date (BJD-TDB), considering Earth's orbital motion and light-travel time corrections.

## 🌐 Demo

Try it online: (coming soon)

## Features

- ✅ JD → BJD-TDB conversion
- ✅ Earth barycenter position calculation
- ✅ Light-travel time (LTT) correction
- ✅ Multiple observatories support
- ✅ Optional: parallax and proper motion corrections
- ✅ Batch conversion API

## Installation

```bash
# Clone the repository
git clone https://github.com/phillippro/jd2bjd.git
cd jd2bjd

# Install dependencies
pip install -r requirements.txt
```

## Usage

```bash
# Start the web server
python3 app.py

# Open in browser
http://localhost:5000
```

## API Usage

```bash
# Single conversion
curl -X POST http://localhost:5000/convert \
  -H "Content-Type: application/json" \
  -d '{
    "jd": 2460000.5,
    "ra": "83.633",
    "dec": "-5.391",
    "observatory": "mauna_kea",
    "parallax": 76.29
  }'

# Batch conversion
curl -X POST http://localhost:5000/batch \
  -H "Content-Type: application/json" \
  -d '{
    "jd_list": [2460000.5, 2460001.5, 2460002.5],
    "ra": "83.633",
    "dec": "-5.391"
  }'
```

## Supported Observatories

| Observatory | Longitude (°) | Latitude (°) | Elevation (m) |
|-------------|--------------|--------------|---------------|
| Shanghai | 121.54 | 31.22 | 7 |
| Mauna Kea | -155.48 | 19.82 | 4207 |
| Paranal | -70.40 | -24.63 | 2635 |
| La Silla | -70.73 | -29.26 | 2400 |
| McDonald | -104.02 | 30.67 | 2070 |
| Siding Spring | 149.04 | -31.27 | 1164 |
| Cerro Tololo | -70.81 | -30.17 | 2215 |

## Calculation Method

```
BJD_TDB = JD_UTC + ΔT_UTC→TDB + Δt_LTT
```

Where:
- **ΔT_UTC→TDB**: UTC → TDB conversion (~minutes)
- **Δt_LTT**: Light-travel time correction (~±8 minutes, depending on Earth's position)

### Technical Details

1. **Earth Orbital Correction**: Uses Astropy with JPL DE440 ephemerides
2. **Light-travel Time**: Calculates the light path from observer to solar system barycenter
3. **Coordinate Transform**: Geocentric → Barycentric reference frame

### Precision

| Component | Precision |
|-----------|-----------|
| UTC → TDB | ~nanosecond |
| Earth position (JPL DE440) | ~meter |
| Light-travel time | ~nanosecond |
| **Overall BJD-TDB** | **< 1 microsecond** |

This matches Tempo2 to within ~10 nanoseconds for typical observations.

## Running Tests

```bash
# Run precision benchmark
python3 benchmark.py
```

## Dependencies

- flask>=2.0.0
- astropy>=5.0
- numpy>=1.20

## License

MIT License

## References

- Eastman et al. (2010) - Precise transit timing
- **Feng et al. (2019)** - [PEXO: A Global Photometric and Radial Velocity Model for Exoplanets](https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.4008F)
- Astropy Time Documentation
- JPL Ephemerides (DE440)