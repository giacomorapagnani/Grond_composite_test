import pygmt
import numpy as np
import os
from pyrocko import util, model, io, trace, gmtpy
import pyrocko.moment_tensor as pmt
import pyrocko.orthodrome as od

workdir = '../'
catdir = os.path.join(workdir, 'CAT')
metadatadir = os.path.join(workdir, 'META_DATA')

# ─────────────────────────────────────────────────────────────────────────────
# UNCERTAINTY ELLIPSES — style parameters (easy to modify)
# ─────────────────────────────────────────────────────────────────────────────
ELLIPSE_PEN             = "0.8p,gray30,dashed"   # pen: thickness, color, style
ELLIPSE_N_POINTS        = 400                    # number of points along the ellipse
ELLIPSE_FILE_NAME       = 'uncertainty_ellipses.txt'  # path to the ellipse parameter file
# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# LOAD UNCERTAINTY ELLIPSES
# columns: event_name lat_event lon_event east_shift_16 east_shift_84
#          nord_shift_16 nord_shift_84   (shifts in metres)
# ─────────────────────────────────────────────────────────────────────────────
def load_ellipses(filepath):
    """Return a dict  {event_name: dict_of_params}."""
    ellipses = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            toks = line.split()
            if len(toks) < 7:
                continue
            name = toks[0]
            ellipses[name] = {
                'lat_ref'  : float(toks[1]),
                'lon_ref'  : float(toks[2]),
                'e16'      : float(toks[3]),
                'e84'      : float(toks[4]),
                'n16'      : float(toks[5]),
                'n84'      : float(toks[6]),
            }
    return ellipses


def make_ellipse_coords(lat_ref, lon_ref, e16, e84, n16, n84, n_pts=200):
    """
    Build an array of (lon, lat) points tracing an axis-aligned ellipse.

    The ellipse axes are derived from the 16th–84th percentile shifts
    (in metres, east/north):
        centre_e = (e16 + e84) / 2
        centre_n = (n16 + n84) / 2
        semi_a   = (e84 - e16) / 2   ← east–west half-axis
        semi_b   = (n84 - n16) / 2   ← north–south half-axis

    od.ne_to_latlon converts (north_m, east_m) offsets from the reference
    point to geographic coordinates.
    """
    centre_e = (e16 + e84) / 2.0
    centre_n = (n16 + n84) / 2.0
    semi_a   = (e84 - e16) / 2.0   # east–west
    semi_b   = (n84 - n16) / 2.0   # north–south

    theta = np.linspace(0, 2 * np.pi, n_pts)
    east_pts  = centre_e + semi_a * np.cos(theta)   # metres
    north_pts = centre_n + semi_b * np.sin(theta)   # metres

    lats, lons = od.ne_to_latlon(lat_ref, lon_ref, north_pts, east_pts)
    return lons, lats


# ─────────────────────────────────────────────────────────────────────────────
# MAP EXTENT SWITCH
# ─────────────────────────────────────────────────────────────────────────────
##########################################
############## SWITCH ##############
##########################################
switch_coord_far = False

if switch_coord_far:
    # FAR COORD
    minlon, maxlon = 14.055, 14.19
    minlat, maxlat = 40.76,  40.87
    map_name = 'far'
else:
    # GULF COORD (NEAR)
    minlon, maxlon = 14.07,  14.175
    minlat, maxlat = 40.775, 40.855
    map_name = 'gulf'

region     = [minlon, maxlon, minlat, maxlat]
projection = "M6i"

# ─────────────────────────────────────────────────────────────────────────────
# BASE MAP
# ─────────────────────────────────────────────────────────────────────────────
fig = pygmt.Figure()
pygmt.config(FORMAT_GEO_MAP="ddd.xxF")

fig.basemap(region=region, projection=projection, frame='a0.05',
            map_scale='x3c/-0.7c+w3')

topo_data = pygmt.datasets.load_earth_relief(resolution="01s", region=region)
fig.grdimage(grid=topo_data, region=region, projection=projection,
             shading="+a45+ne0.5", cmap="gray")
fig.coast(shorelines="1/0.5p,black", resolution="f", water="#EBEBEE")

# ─────────────────────────────────────────────────────────────────────────────
# LOAD EVENTS AND ELLIPSE PARAMETERS
# ─────────────────────────────────────────────────────────────────────────────
filename   = 'catalogue_flegrei_composite_MT_LF_std_reloc_best'     #CHANGE
events_name = os.path.join(catdir, filename + '.pf')
fm_events  = model.load_events(events_name)

ellipse_filepath = os.path.join(catdir, ELLIPSE_FILE_NAME)
ellipses = load_ellipses(ellipse_filepath)

# ─────────────────────────────────────────────────────────────────────────────
# SWITCHES
# ─────────────────────────────────────────────────────────────────────────────
##########################################
############## SWITCH ##############
##########################################
switch_deviatoric = True    # TRUE → deviatoric MT;  FALSE → DC only
switch_timestamps = False   # TRUE → print date next to each beach ball

# ─────────────────────────────────────────────────────────────────────────────
# PLOT UNCERTAINTY ELLIPSES  (drawn first → below the beach balls)
# ─────────────────────────────────────────────────────────────────────────────
for ev in fm_events:
    if ev.name not in ellipses:
        print(f"[WARNING] no ellipse found for event: {ev.name}")
        continue

    ep = ellipses[ev.name]
    ell_lons, ell_lats = make_ellipse_coords(
        ep['lat_ref'], ep['lon_ref'],
        ep['e16'], ep['e84'],
        ep['n16'], ep['n84'],
        n_pts=ELLIPSE_N_POINTS
    )

    fig.plot(x=ell_lons, y=ell_lats, pen=ELLIPSE_PEN,
             region=region, projection=projection)

# ─────────────────────────────────────────────────────────────────────────────
# PLOT FOCAL MECHANISMS
# ─────────────────────────────────────────────────────────────────────────────
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
            "exponent": 1
        }

        MT_white = True
        if MT_white:
            fig.meca(spec=moment_tensor_par, convention='mt',
                     longitude=ev.lon, latitude=ev.lat, depth=ev.depth,
                     scale="1.2c", compressionfill="white",
                     extensionfill="white", pen="1p,black", outline="1p,black")
        else:
            fig.meca(spec=moment_tensor_par, convention='mt',
                     longitude=ev.lon, latitude=ev.lat, depth=ev.depth,
                     scale="0.8c", compressionfill="#BD2025",
                     extensionfill="white", pen="0.5p,gray30,solid")
    else:
        moment_tensor_par = {
            "strike"   : ev.moment_tensor.strike1,
            "dip"      : ev.moment_tensor.dip1,
            "rake"     : ev.moment_tensor.rake1,
            "magnitude": ev.magnitude
        }
        fig.meca(spec=moment_tensor_par,
                 longitude=ev.lon, latitude=ev.lat, depth=ev.depth,
                 scale="0.8c", compressionfill="#BD2025",
                 extensionfill="white", pen="0.5p,gray30,solid")

    if switch_timestamps:
        name     = ev.name.split('_')[1:]
        name_ev  = f"{name[0]}-{name[1]}-{name[2]}_{name[3]}:{name[4]}:{name[5]}"
        fig.text(text=name_ev, x=ev.lon, y=ev.lat + 0.0006,
                 font="5p,Helvetica,black", justify="CM")

# ─────────────────────────────────────────────────────────────────────────────
# STATIONS
# ─────────────────────────────────────────────────────────────────────────────
latsta, lonsta, namsta = [], [], []
with open(os.path.join(metadatadir, 'stations_flegrei_INGV_final.pf'), 'r') as f:
    for line in f:
        if line[0] == ' ':
            continue
        toks = line.split()
        latsta.append(float(toks[1]))
        lonsta.append(float(toks[2]))
        namsta.append(toks[0].split('.')[1])

latsta = np.array(latsta)
lonsta = np.array(lonsta)

fig.plot(x=lonsta, y=latsta, style="t0.3", fill="#FFCC4E", pen="black")
# fig.text(x=lonsta+0.005, y=latsta+0.002, text=namsta,
#          justify='BR', font='5p', fill="#FFCC4E")

# ─────────────────────────────────────────────────────────────────────────────
# SAVE & SHOW
# ─────────────────────────────────────────────────────────────────────────────
fig.show()
suffix = 'deviatoric_ellipse' if switch_deviatoric else 'dc_ellipse'
fig.savefig(f'../PLOTS/MAPS/{filename}_{suffix}_{map_name}.pdf')
