import pygmt
import numpy as np
import os
from pyrocko import util, model, io, trace
from pyrocko import orthodrome as od
import pyrocko.moment_tensor as pmt

# ──────────────────────────────────────────────
#  PATHS
# ──────────────────────────────────────────────
workdir     = '../'
catdir      = os.path.join(workdir, 'CAT')
metadatadir = os.path.join(workdir, 'META_DATA')

# ──────────────────────────────────────────────
#  SWITCHES
# ──────────────────────────────────────────────

# COORDINATES: Gulf map (near) or far map
switch_coord_far = False

# MOMENT TENSOR: deviatoric (True) or double-couple (False)
switch_deviatoric = True

# TIMESTAMPS: print event date next to focal sphere (True/False)
switch_timestamps = False

# ELLIPSES: plot location-uncertainty ellipses (True/False)
switch_ellipses = True

# ──────────────────────────────────────────────
#  MAP REGION
# ──────────────────────────────────────────────
if switch_coord_far:
    minlon, maxlon = 14.055, 14.19
    minlat, maxlat = 40.76,  40.87
    map_name = 'far'
else:
    minlon, maxlon = 14.07,  14.175
    minlat, maxlat = 40.775, 40.855
    map_name = 'gulf'

region     = [minlon, maxlon, minlat, maxlat]
projection = "M6i"   # Mercator, 6-inch width

# ──────────────────────────────────────────────
#  HELPER: build uncertainty ellipse in lon/lat
# ──────────────────────────────────────────────
def build_ellipse(lat0, lon0, east_shift_16, east_shift_84,
                  nord_shift_16, nord_shift_84, n_points=180):
    """
    Return (lon_pts, lat_pts) arrays for the uncertainty ellipse.

    The ellipse centre is the mean of the 16th–84th percentile shifts;
    semi-axes are half the spread along each direction.
    The focal sphere sits at (lat0, lon0), which may fall outside the ellipse.

    Conversion metres → lon/lat uses od.ne_to_latlon().
    """
    center_n = (nord_shift_16 + nord_shift_84) / 2.0   # north offset of centre [m]
    center_e = (east_shift_16 + east_shift_84) / 2.0   # east  offset of centre [m]

    semi_e = (east_shift_84 - east_shift_16) / 2.0     # semi-axis along east   [m]
    semi_n = (nord_shift_84 - nord_shift_16) / 2.0     # semi-axis along north  [m]

    theta    = np.linspace(0.0, 2.0 * np.pi, n_points)
    north_pts = center_n + semi_n * np.sin(theta)
    east_pts  = center_e + semi_e * np.cos(theta)

    lat_pts, lon_pts = od.ne_to_latlon(lat0, lon0, north_pts, east_pts)
    return lon_pts, lat_pts


# ──────────────────────────────────────────────
#  LOAD UNCERTAINTY PARAMETERS  (if needed)
# ──────────────────────────────────────────────
# Expected .txt format (space or comma separated, one event per line):
#   nome_evento  lat_event  lon_event  east_shift_16  east_shift_84  nord_shift_16  nord_shift_84
#
# Units: lat/lon in degrees, shifts in metres.

ellipse_params = {}   # dict: ev.name -> dict of shift values

if switch_ellipses:
    ellipse_filename='uncertainty_ellipses'                 #CHANGE
    ellipse_file = os.path.join(catdir, ellipse_filename + '.txt')
    with open(ellipse_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # accept both comma-separated and whitespace-separated
            toks = line.replace(',', ' ').split()
            if len(toks) < 7:
                print(f"[WARNING] skipping malformed line: {line}")
                continue
            ev_name        = toks[0]
            lat_ev         = float(toks[1])
            lon_ev         = float(toks[2])
            east_shift_16  = float(toks[3])
            east_shift_84  = float(toks[4])
            nord_shift_16  = float(toks[5])
            nord_shift_84  = float(toks[6])
            ellipse_params[ev_name] = dict(
                lat_ev        = lat_ev,
                lon_ev        = lon_ev,
                east_shift_16 = east_shift_16,
                east_shift_84 = east_shift_84,
                nord_shift_16 = nord_shift_16,
                nord_shift_84 = nord_shift_84,
            )
    print(f"[INFO] loaded ellipse parameters for {len(ellipse_params)} events")

# ──────────────────────────────────────────────
#  BUILD FIGURE
# ──────────────────────────────────────────────
fig = pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="ddd.xxF")

fig.basemap(region=region, projection=projection,
            frame='a0.05', map_scale='x3c/-0.7c+w3')

# Topography (1 arc-second shaded relief)
topo_data = pygmt.datasets.load_earth_relief(resolution="01s", region=region)
fig.grdimage(grid=topo_data, region=region, projection=projection,
             shading="+a45+ne0.5", cmap="gray")

# Coastlines
fig.coast(shorelines="1/0.5p,black", resolution="f", water="#EBEBEE")

# ──────────────────────────────────────────────
#  LOAD FOCAL MECHANISM CATALOGUE
# ──────────────────────────────────────────────
filename   = 'catalogue_flegrei_composite_MT_LF_std_reloc_best'     #CHANGE
events_name = os.path.join(catdir, filename + '.pf')
fm_events   = model.load_events(events_name)

# ──────────────────────────────────────────────
#  PLOT UNCERTAINTY ELLIPSES  (below focal spheres → plotted first)
# ──────────────────────────────────────────────
if switch_ellipses:
    for ev in fm_events:
        if ev.name not in ellipse_params:
            continue
        p = ellipse_params[ev.name]
        lon_ell, lat_ell = build_ellipse(
            lat0          = p['lat_ev'],
            lon0          = p['lon_ev'],
            east_shift_16 = p['east_shift_16'],
            east_shift_84 = p['east_shift_84'],
            nord_shift_16 = p['nord_shift_16'],
            nord_shift_84 = p['nord_shift_84'],
        )
        fig.plot(x=lon_ell, y=lat_ell,
                 pen="0.7p,gray30,dashed",
                 region=region, projection=projection)

# ──────────────────────────────────────────────
#  PLOT FOCAL MECHANISMS
# ──────────────────────────────────────────────
for ev in fm_events:

    if switch_deviatoric:
        msix = pmt.to6(ev.moment_tensor.m_up_south_east())
        moment_tensor_par = {
            "mrr": msix[0] * 10**7,
            "mtt": msix[1] * 10**7,
            "mff": msix[2] * 10**7,
            "mrt": msix[3] * 10**7,
            "mrf": msix[4] * 10**7,
            "mtf": msix[5] * 10**7,
            "exponent": 1,
        }
        MT_white = True
        if MT_white:
            fig.meca(spec=moment_tensor_par, convention='mt',
                     longitude=ev.lon, latitude=ev.lat, depth=ev.depth,
                     scale="1.2c",
                     compressionfill="white", extensionfill="white",
                     pen="1p,black", outline="1p,black")
        else:
            fig.meca(spec=moment_tensor_par, convention='mt',
                     longitude=ev.lon, latitude=ev.lat, depth=ev.depth,
                     scale="0.8c",
                     compressionfill="#BD2025", extensionfill="white",
                     pen="0.5p,gray30,solid")
    else:
        moment_tensor_par = {
            "strike":    ev.moment_tensor.strike1,
            "dip":       ev.moment_tensor.dip1,
            "rake":      ev.moment_tensor.rake1,
            "magnitude": ev.magnitude,
        }
        fig.meca(spec=moment_tensor_par,
                 longitude=ev.lon, latitude=ev.lat, depth=ev.depth,
                 scale="0.8c",
                 compressionfill="#BD2025", extensionfill="white",
                 pen="0.5p,gray30,solid")
                 # blue: #0066cc   red: #BD2025

    if switch_timestamps:
        name    = ev.name.split('_')[1:]
        name_ev = (f"{name[0]}-{name[1]}-{name[2]}"
                   f"_{name[3]}:{name[4]}:{name[5]}")
        fig.text(text=name_ev, x=ev.lon, y=ev.lat + 0.0006,
                 font="5p,Helvetica,black", justify="CM")

# ──────────────────────────────────────────────
#  PLOT STATIONS
# ──────────────────────────────────────────────
latsta, lonsta, namsta = [], [], []
with open(metadatadir + '/stations_flegrei_INGV_final.pf', 'r') as f:
    for line in f:
        if line[0] == ' ':
            continue
        toks = line.split()
        latsta.append(float(toks[1]))
        lonsta.append(float(toks[2]))
        namsta.append(toks[0].split('.')[1])

latsta = np.array(latsta)
lonsta = np.array(lonsta)

fig.plot(x=lonsta, y=latsta, style="t0.3",
         fill="#FFCC4E", pen="black")

# ──────────────────────────────────────────────
#  SAVE
# ──────────────────────────────────────────────
fig.show()

suffix = 'deviatoric' if switch_deviatoric else 'dc'
if switch_ellipses:
     suffix += '_ellipses'
out_path = f'../PLOTS/MAPS/{filename}_{suffix}_{map_name}.pdf'
fig.savefig(out_path)
print(f"[INFO] figure saved → {out_path}")
