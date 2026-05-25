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
ELLIPSE_PEN        = "2.0p,gray10,2p_2p"   # pen: thickness, color:gray10 / #BD2025, style
ELLIPSE_N_POINTS   = 400                    # number of points along the ellipse
ELLIPSE_FILE_NAME       = 'uncertainty_ellipses_old.txt'  # path to the ellipse parameter file
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

    The 4 shifts (metres) define the cardinal extremes of the ellipse.
    Each extreme is converted to geographic coordinates independently via
    od.ne_to_latlon(); lat_ref/lon_ref is only the local origin and may
    fall outside the ellipse.

    Returns: (lons, lats, semi_a_m, semi_b_m)
        semi_a_m, semi_b_m  — half-axes in metres (for area computation)
    """
    # Convert the 4 cardinal extreme points to geographic coords
    lat_w, lon_w = od.ne_to_latlon(lat_ref, lon_ref, 0.0,  e16)   # west extreme
    lat_e, lon_e = od.ne_to_latlon(lat_ref, lon_ref, 0.0,  e84)   # east extreme
    lat_s, lon_s = od.ne_to_latlon(lat_ref, lon_ref, n16,  0.0)   # south extreme
    lat_n, lon_n = od.ne_to_latlon(lat_ref, lon_ref, n84,  0.0)   # north extreme

    # Ellipse centre = midpoint of the two axis pairs
    lat_c = (lat_s + lat_n) / 2.0
    lon_c = (lon_w + lon_e) / 2.0

    # Semi-axes in degrees
    semi_a_deg = abs(lon_e - lon_w) / 2.0   # east–west
    semi_b_deg = abs(lat_n - lat_s) / 2.0   # north–south

    # Semi-axes in metres (kept in original units for area)
    semi_a_m = (e84 - e16) / 2.0
    semi_b_m = (n84 - n16) / 2.0

    # Parametric ellipse in geographic coords
    theta     = np.linspace(0, 2 * np.pi, n_pts)
    ell_lons  = lon_c + semi_a_deg * np.cos(theta)
    ell_lats  = lat_c + semi_b_deg * np.sin(theta)

    return ell_lons, ell_lats, semi_a_m, semi_b_m


def make_ellipse_coords_from_axes(lon_c, lat_c, semi_a_deg, semi_b_deg, n_pts=200):
    """Build ellipse coords given a geographic centre and semi-axes in degrees."""
    theta    = np.linspace(0, 2 * np.pi, n_pts)
    ell_lons = lon_c + semi_a_deg * np.cos(theta)
    ell_lats = lat_c + semi_b_deg * np.sin(theta)
    return ell_lons, ell_lats


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
filename   = 'catalogue_flegrei_MT_VLP_reloc'     #CHANGE: catalogue_flegrei_composite_MT_LF_std_reloc_best
                                                            # catalogue_flegrei_MT_VLP_reloc
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
# also collect semi-axes (metres) for mean ellipse
# ─────────────────────────────────────────────────────────────────────────────
all_semi_a_m = []   # east–west half-axes in metres
all_semi_b_m = []   # north–south half-axes in metres

for ev in fm_events:
    if ev.name not in ellipses:
        print(f"[WARNING] no ellipse found for event: {ev.name}")
        continue

    ep = ellipses[ev.name]
    ell_lons, ell_lats, semi_a_m, semi_b_m = make_ellipse_coords(
        ep['lat_ref'], ep['lon_ref'],
        ep['e16'], ep['e84'],
        ep['n16'], ep['n84'],
        n_pts=ELLIPSE_N_POINTS
    )

    all_semi_a_m.append(semi_a_m)
    all_semi_b_m.append(semi_b_m)

    fig.plot(x=ell_lons, y=ell_lats, pen=ELLIPSE_PEN,
             region=region, projection=projection)

# ─────────────────────────────────────────────────────────────────────────────
# MEAN ELLIPSE — plotted in the bottom-left corner of the map
# ─────────────────────────────────────────────────────────────────────────────
if all_semi_a_m:
    mean_semi_a_m = np.mean(all_semi_a_m)   # metres
    mean_semi_b_m = np.mean(all_semi_b_m)   # metres
    mean_area_m2  = np.pi * mean_semi_a_m * mean_semi_b_m

    print(f"[INFO] mean semi-axis east–west : {mean_semi_a_m:.1f} m")
    print(f"[INFO] mean semi-axis north–south: {mean_semi_b_m:.1f} m")
    print(f"[INFO] mean ellipse area         : {mean_area_m2/1e6:.4f} km²")

    # Position the mean ellipse in the bottom-left corner (in degrees)
    MEAN_ELLIPSE_OFFSET_LON = 0.015   # degrees from minlon
    MEAN_ELLIPSE_OFFSET_LAT = 0.015   # degrees from minlat

    # Convert mean semi-axes from metres to degrees at the anchor point
    anchor_lat = minlat + MEAN_ELLIPSE_OFFSET_LAT
    anchor_lon = minlon + MEAN_ELLIPSE_OFFSET_LON

    # 1 degree latitude ≈ 111320 m  (constant)
    # 1 degree longitude ≈ 111320 * cos(lat) m
    deg_per_m_lat = 1.0 / 111320.0
    deg_per_m_lon = 1.0 / (111320.0 * np.cos(np.radians(anchor_lat)))

    mean_semi_a_deg = mean_semi_a_m * deg_per_m_lon   # east–west
    mean_semi_b_deg = mean_semi_b_m * deg_per_m_lat   # north–south

    mean_lons, mean_lats = make_ellipse_coords_from_axes(
        anchor_lon, anchor_lat,
        mean_semi_a_deg, mean_semi_b_deg,
        n_pts=ELLIPSE_N_POINTS
    )

    fig.plot(x=mean_lons, y=mean_lats, pen=ELLIPSE_PEN,
             region=region, projection=projection)

    # Label below the mean ellipse
    fig.text(
        text=f"mean ellipse  ({mean_area_m2/1e6:.3f} km\u00b2)",
        x=anchor_lon,
        y=anchor_lat - mean_semi_b_deg - 0.002,
        font="5p,Helvetica,black",
        justify="CM"
    )

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
suffix = 'deviatoric_ellipse_mean' if switch_deviatoric else 'dc_ellipse_mean'
fig.savefig(f'../PLOTS/MAPS/{filename}_{suffix}_{map_name}.pdf')
