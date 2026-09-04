"""
Map of Campi Flegrei focal mechanisms.

Layers, drawn bottom → top:
    topography → uncertainty ellipses → VT circles → beach balls
    → magnitude legend → stations → depth colorbar

Everything meant to be tuned lives in the USER SETTINGS block below; nothing
below the "END OF USER SETTINGS" line needs to be touched for routine changes.
Run from inside `codes/` (all paths are relative to it).
"""

import os

import numpy as np
import pygmt
import rioxarray
import xarray as xr
from pyrocko import model
import pyrocko.moment_tensor as pmt
import pyrocko.orthodrome as od


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║                            USER SETTINGS                                  ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# ═══════════════════════════════════════════════════════════════
#  1 · PATHS AND INPUT CATALOGUE
# ═══════════════════════════════════════════════════════════════
workdir      = '../'
catdir       = os.path.join(workdir, 'CAT')
metadatadir  = os.path.join(workdir, 'META_DATA')
plotdir      = os.path.join(workdir, 'PLOTS', 'MAPS')

# Catalogue whose moment tensors are drawn as beach balls (no .pf extension).
#   catalogue_flegrei_composite_MT_LF_std_reloc_best_FLIPPED
#   catalogue_flegrei_composite_MT_LF_std_reloc_best
#   catalogue_flegrei_MT_VT
filename = 'catalogue_flegrei_composite_MT_LF_std_reloc_best_FLIPPED'

STATION_FILE = 'stations_flegrei_INGV_final.pf'   # inside metadatadir

# ═══════════════════════════════════════════════════════════════
#  2 · MAP AREA AND PROJECTION
# ═══════════════════════════════════════════════════════════════
#   True  → tight zoom on Pozzuoli
#   False → whole gulf
switch_pozzuoli = False

MAP_WIDTH_IN = 6          # Mercator map width, inches (drives the whole layout)

# ═══════════════════════════════════════════════════════════════
#  3 · TOPOGRAPHY
# ═══════════════════════════════════════════════════════════════
#   'local'  → Tinitaly 10 m DEM (land-only, high-res)
#   'pygmt'  → GMT built-in 01s relief (~30 m, includes bathymetry)
topo_source = 'local'
tif_path    = '/Users/giaco/UNI/PhD_CODE/QGIS/topo_flegrei/w45090_s10.tif'
water_color = "#A9B0B8"   # sea / NaN colour — change freely

# Topography colour palette (tones + hillshading).
#   'gray'   → classic black-and-white
#   'brown'  → soft cream → warm sienna
#   'green'  → light mint → olive
#   'sand'   → pale ivory → dark tan
# Any GMT two-colour gradient also works, e.g. 'lightyellow,chocolate'.
topo_color = 'gray'

# ═══════════════════════════════════════════════════════════════
#  4 · FOCAL MECHANISMS (BEACH BALLS)
# ═══════════════════════════════════════════════════════════════
#   True  → plot the full deviatoric tensor
#   False → plot the double-couple nodal plane 1 (strike/dip/rake)
switch_deviatoric = True

#   'depth_color' → colour gradient by depth (palette set by meca_color)
#   'fixed'       → uniform fill set by fixed_color
color_mode  = 'fixed'
fixed_color = '#BD2025'

# Beach-ball palette, only used when color_mode == 'depth_color'.
#   'gray' | 'blue' | 'red' | 'green' | 'brown' | 'purple' | 'black'
# A custom (light_rgb, dark_rgb) tuple also works,
# e.g. meca_color = ((255, 235, 205), (139, 69, 19))
meca_color = 'red'

# ═══════════════════════════════════════════════════════════════
#  5 · VT EVENT CIRCLES
#  One circle per event of an extra catalogue, radius scaling with magnitude.
#  Meant for the VT seismicity plotted underneath the beach balls.
# ═══════════════════════════════════════════════════════════════
switch_vt_circles = True                          # master on/off switch

VT_FILE_NAME = 'catalogue_flegrei_MT_VT.pf'       # inside catdir

# Circle colour — any GMT colour: 'black', 'gray30', '#BD2025', '255/0/0' …
VT_COLOR = 'gray90'

#   'filled'  → solid disc, no border
#   'outline' → empty circle, coloured border only
#   'both'    → solid disc + border
VT_MODE = 'both'

VT_PEN_WIDTH = '0.8p'      # border thickness ('outline' and 'both')
VT_PEN_COLOR = 'black'     # border colour used by 'both' ('outline' uses VT_COLOR)
VT_ALPHA     = 0           # 0 = opaque, 100 = invisible

# Radius vs magnitude.
#   'meca'   → same diameter a beach ball of that Mw would get, so circles and
#              beach balls are directly comparable in size
#   'linear' → diameter_cm = VT_SIZE_BASE_CM + VT_SIZE_SLOPE_CM * Mw
VT_SIZE_MODE   = 'meca'
VT_SIZE_FACTOR = 0.40      # global multiplier: 1.0 = exactly beach-ball size
VT_SIZE_BASE_CM  = 0.05    # 'linear' mode only
VT_SIZE_SLOPE_CM = 0.30    # 'linear' mode only

#   False → circles below the beach balls (default, keeps mechanisms readable)
#   True  → circles on top of the beach balls
VT_ON_TOP = False

# ═══════════════════════════════════════════════════════════════
#  6 · UNCERTAINTY ELLIPSES
#  Drawn BEFORE the beach balls, so they stay underneath.
# ═══════════════════════════════════════════════════════════════
#   'outline' → dashed border only
#   'filled'  → semi-transparent solid fill
#   'both'    → transparent fill + dashed border
#   'none'    → ellipses switched off
ELLIPSE_MODE = 'none'

# Any GMT colour ('gray40', 'black', '#BD2025', '255/0/0' …).
# Special value 'depth' → each ellipse takes the depth colour of its beach ball.
ELLIPSE_COLOR = '#BD2025'

# Fill transparency: 0 = opaque, 100 = invisible. 70–85 works well when many
# ellipses overlap.
ELLIPSE_ALPHA = 85

ELLIPSE_PEN_WIDTH = '0.8p'      # border thickness
ELLIPSE_DASH      = '4p_3p'     # dash pattern: 'solid' for a continuous line
ELLIPSE_N_POINTS  = 400         # points used to trace each ellipse
ELLIPSE_FILE_NAME = 'uncertainty_ellipses.txt'    # inside catdir

# Reference ("mean") ellipse drawn in the bottom-left corner as a size legend.
switch_mean_ellipse     = True
MEAN_ELLIPSE_OFFSET_LON = 0.015   # degrees from minlon
MEAN_ELLIPSE_OFFSET_LAT = 0.015   # degrees from minlat

# ═══════════════════════════════════════════════════════════════
#  7 · ANNOTATIONS
# ═══════════════════════════════════════════════════════════════
switch_timestamps = True    # date/time label next to every beach ball
switch_sta_names  = False    # station code next to every station triangle

# ═══════════════════════════════════════════════════════════════
#  8 · STATIONS
# ═══════════════════════════════════════════════════════════════
STA_STYLE = 't0.3c'        # GMT symbol: t = triangle, 0.3c = size
STA_FILL  = '#FFCC4E'
STA_PEN   = '0.6p,black'

# ═══════════════════════════════════════════════════════════════
#  9 · LEGEND (top-left box)
#  Two independent magnitude scales, stacked, plus the station symbol:
#    · beach balls (VLP) — circles sized exactly as the beach balls on the map
#    · VT circles        — circles sized by vt_diameter_cm()
#  The two scales are unrelated on purpose: each row shows the size its own
#  symbol really gets on the map.
# ═══════════════════════════════════════════════════════════════
LEG_SHOW_MECA = True                      # magnitude scale of the beach balls
LEG_MECA_TITLE = 'VLP'
ref_mags_meca  = [1.0, 2.0, 3.0]

LEG_SHOW_VT   = True                      # magnitude scale of the VT circles
LEG_VT_TITLE  = 'VT'                      # (skipped when switch_vt_circles is off)
ref_mags_vt   = [1.0, 2.0, 3.0, 4.0]

LEG_SHOW_STATION  = True                  # station triangle + label row
LEG_STATION_LABEL = 'stations'

LEG_PAD        = 0.16    # cm — inner padding of the box
LEG_GAP        = 0.20    # cm — clear space between two circles
LEG_MIN_STEP   = 0.62    # cm — min centre-to-centre distance (keeps labels apart)
LEG_TITLE_H    = 0.34    # cm — row reserved for each scale title
LEG_LABEL_H    = 0.32    # cm — row reserved for the labels under the circles
LEG_STATION_H  = 0.40    # cm — row reserved for the station symbol
LEG_SYM_GAP    = 0.10    # cm — space between station symbol and its label
LEG_MARGIN     = 0.10    # cm — distance from the map frame
LEG_FONT_TITLE = "5p,Helvetica-Bold,black"
LEG_FONT_LABEL = "4.5p,Helvetica,black"

# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║                        END OF USER SETTINGS                               ║
# ╚═══════════════════════════════════════════════════════════════════════════╝


# ═══════════════════════════════════════════════════════════════
#  DERIVED SETTINGS AND COLOUR PRESETS
# ═══════════════════════════════════════════════════════════════
os.makedirs(plotdir, exist_ok=True)

if switch_pozzuoli:
    minlon, maxlon = 14.07, 14.175
    minlat, maxlat = 40.79,  40.845
    map_name = 'pozzuoli'
else:
    minlon, maxlon = 14.07, 14.175
    minlat, maxlat = 40.775, 40.855
    map_name = 'gulf'

region     = [minlon, maxlon, minlat, maxlat]
projection = f"M{MAP_WIDTH_IN}i"
map_w_cm   = MAP_WIDTH_IN * 2.54

_TOPO_PRESETS = {
    'gray':  ('gray',               True),   # (GMT cmap, reverse)
    'brown': ('linen,sienna',       False),  # cream → warm brown
    'green': ('honeydew,olivedrab', False),  # light mint → olive
    'sand':  ('ivory,tan',          False),  # pale ivory → warm tan
}

_MECA_PRESETS = {
    'gray':   ((205, 205, 205), (50,  50,  50)),    # shallow → deep
    'blue':   ((198, 219, 239), (8,   81,  156)),
    'red':    ((252, 187, 161), (165, 15,  21)),
    'green':  ((199, 233, 192), (0,   109, 44)),
    'brown':  ((253, 208, 162), (140, 81,  10)),
    'purple': ((218, 218, 235), (84,  39,  143)),
    'black':  ((255, 255, 255), (0,   0,   0)),
}

MECA_EXPONENT = 1     # exponent field of the 'mt' spec — see meca_diameter_cm()


# ═══════════════════════════════════════════════════════════════
#  LOAD CATALOGUES
# ═══════════════════════════════════════════════════════════════
fm_events = model.load_events(os.path.join(catdir, filename + '.pf'))
fm_events = [ev for ev in fm_events if ev.moment_tensor is not None]

depths_km = np.array([ev.depth / 1000 for ev in fm_events])
depth_min = max(0.0, depths_km.min() - 1.0)
depth_max = depths_km.max() + 1.0

vt_events = []
if switch_vt_circles:
    vt_events = model.load_events(os.path.join(catdir, VT_FILE_NAME))
    vt_events = [ev for ev in vt_events if ev.magnitude is not None]
    # big first, so small circles stay visible on top of large ones
    vt_events.sort(key=lambda ev: ev.magnitude, reverse=True)


# ═══════════════════════════════════════════════════════════════
#  DEPTH → COLOUR  (pure Python, independent of GMT CPT state)
#  shallow → light tone, deep → dark tone (palette set by meca_color)
# ═══════════════════════════════════════════════════════════════
def depth_to_color(depth_km):
    light, dark = meca_color if isinstance(meca_color, tuple) \
        else _MECA_PRESETS.get(meca_color, _MECA_PRESETS['gray'])
    t = np.clip((depth_km - depth_min) / (depth_max - depth_min), 0.0, 1.0)
    rgb = [round(light[i] + t * (dark[i] - light[i])) for i in range(3)]
    return f"{rgb[0]}/{rgb[1]}/{rgb[2]}"


# ═══════════════════════════════════════════════════════════════
#  MAGNITUDE → GMT `scale` PARAMETER
#  NOTE: this is NOT the size drawn on paper. GMT's -S option treats
#  `scale` as the diameter of a magnitude-5 beach ball and shrinks it
#  linearly with magnitude — see meca_diameter_cm() below.
# ═══════════════════════════════════════════════════════════════
def mag_to_size(mag):
    return 1.0 + 0.10 * mag


# ═══════════════════════════════════════════════════════════════
#  MOMENT TENSOR → GMT meca SPEC
# ═══════════════════════════════════════════════════════════════
def meca_spec(ev):
    """Return (spec_dict, convention) for fig.meca()."""
    if switch_deviatoric:
        msix = pmt.to6(ev.moment_tensor.m_up_south_east())
        spec = {
            "mrr": msix[0] * 1e7, "mtt": msix[1] * 1e7, "mff": msix[2] * 1e7,
            "mrt": msix[3] * 1e7, "mrf": msix[4] * 1e7, "mtf": msix[5] * 1e7,
            "exponent": MECA_EXPONENT,
        }
        return spec, 'mt'

    spec = {
        "strike":    ev.moment_tensor.strike1,
        "dip":       ev.moment_tensor.dip1,
        "rake":      ev.moment_tensor.rake1,
        "magnitude": ev.magnitude,
    }
    return spec, 'aki'


# ═══════════════════════════════════════════════════════════════
#  MAGNITUDE → DIAMETER ACTUALLY DRAWN ON PAPER (cm)
#
#  GMT draws a beach ball whose diameter is  scale * Mw / 5 , i.e. `scale`
#  is the diameter at Mw 5 and the size is linear in magnitude (verified by
#  rendering and measuring: aki, scale=1c → 0.203/0.402/0.605/0.804 cm for
#  Mw 1/2/3/4).
#
#  With convention='mt' GMT does not read a magnitude, it derives one from
#  the tensor it is handed:  Mw_gmt = (log10(M0[dyne-cm]) - 16.1) / 1.5.
#  The spec passes full dyne-cm values together with exponent=1, so GMT sees
#  a moment 10x the true one and every beach ball comes out ~0.635 magnitude
#  units too big (a constant +0.13*scale cm — uniform, so relative sizes stay
#  correct). This function reproduces that chain exactly, so the legend shows
#  the size a given Mw really gets on the map.
# ═══════════════════════════════════════════════════════════════
def meca_diameter_cm(mag):
    if switch_deviatoric:
        m0_dyne_cm = pmt.magnitude_to_moment(mag) * 1e7 * 10 ** MECA_EXPONENT
        mag_gmt = (np.log10(m0_dyne_cm) - 16.1) / 1.5
    else:
        mag_gmt = mag
    return mag_to_size(mag) * mag_gmt / 5.0


# ═══════════════════════════════════════════════════════════════
#  VT CIRCLES
# ═══════════════════════════════════════════════════════════════
def vt_diameter_cm(mag):
    """Circle diameter on paper (cm) for a VT event of magnitude `mag`."""
    if VT_SIZE_MODE == 'meca':
        diam = meca_diameter_cm(mag)
    else:                                    # 'linear'
        diam = VT_SIZE_BASE_CM + VT_SIZE_SLOPE_CM * mag
    return max(VT_SIZE_FACTOR * diam, 0.01)  # never let a circle vanish


def draw_vt_circles(events):
    """One circle per event, radius scaling with magnitude (VT_* settings)."""
    fill = None if VT_MODE == 'outline' else VT_COLOR
    pen = None
    if VT_MODE == 'outline':
        pen = f"{VT_PEN_WIDTH},{VT_COLOR}"
    elif VT_MODE == 'both':
        pen = f"{VT_PEN_WIDTH},{VT_PEN_COLOR}"

    for ev in events:
        fig.plot(
            x=ev.lon, y=ev.lat,
            style=f"c{vt_diameter_cm(ev.magnitude):.3f}c",
            fill=fill, pen=pen, transparency=VT_ALPHA,
            region=region, projection=projection,
        )

    if events:
        mags = [ev.magnitude for ev in events]
        print(f"[INFO] VT circles plotted          : {len(events)} "
              f"(Mw {min(mags):.2f} – {max(mags):.2f}, "
              f"{vt_diameter_cm(min(mags)):.2f} – "
              f"{vt_diameter_cm(max(mags)):.2f} cm)")
    else:
        print(f"[WARNING] no VT event with a magnitude in {VT_FILE_NAME}")


# ═══════════════════════════════════════════════════════════════
#  UNCERTAINTY ELLIPSES
#  File columns: event_name lat_ref lon_ref east_16 east_84 north_16 north_84
#                (shifts in metres, relative to lat_ref/lon_ref)
# ═══════════════════════════════════════════════════════════════
def load_ellipses(filepath):
    """Return {event_name: {lat_ref, lon_ref, e16, e84, n16, n84}}."""
    ellipses = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            toks = line.split()
            if len(toks) < 7:
                continue
            ellipses[toks[0]] = {
                'lat_ref': float(toks[1]),
                'lon_ref': float(toks[2]),
                'e16':     float(toks[3]),
                'e84':     float(toks[4]),
                'n16':     float(toks[5]),
                'n84':     float(toks[6]),
            }
    return ellipses


def make_ellipse_coords(lat_ref, lon_ref, e16, e84, n16, n84, n_pts=400):
    """
    Trace an axis-aligned ellipse from the 16th/84th percentile shifts.

    The 4 shifts (metres) are the cardinal extremes of the ellipse. Each one
    is converted to geographic coordinates independently through
    od.ne_to_latlon(); lat_ref/lon_ref is only the local origin and may well
    fall outside the ellipse — the ellipse centre is the midpoint of the two
    axis pairs, NOT the catalogue hypocentre.

    Returns: (lons, lats, semi_a_m, semi_b_m)  — half-axes in metres.
    """
    lat_w, lon_w = od.ne_to_latlon(lat_ref, lon_ref, 0.0, e16)   # west extreme
    lat_e, lon_e = od.ne_to_latlon(lat_ref, lon_ref, 0.0, e84)   # east extreme
    lat_s, lon_s = od.ne_to_latlon(lat_ref, lon_ref, n16, 0.0)   # south extreme
    lat_n, lon_n = od.ne_to_latlon(lat_ref, lon_ref, n84, 0.0)   # north extreme

    lat_c = (lat_s + lat_n) / 2.0
    lon_c = (lon_w + lon_e) / 2.0

    semi_a_deg = abs(lon_e - lon_w) / 2.0    # east–west, degrees
    semi_b_deg = abs(lat_n - lat_s) / 2.0    # north–south, degrees

    semi_a_m = (e84 - e16) / 2.0             # east–west, metres
    semi_b_m = (n84 - n16) / 2.0             # north–south, metres

    theta    = np.linspace(0.0, 2.0 * np.pi, n_pts)
    ell_lons = lon_c + semi_a_deg * np.cos(theta)
    ell_lats = lat_c + semi_b_deg * np.sin(theta)

    return ell_lons, ell_lats, semi_a_m, semi_b_m


def draw_ellipse(lons, lats, color):
    """Draw one ellipse honouring ELLIPSE_MODE / ELLIPSE_ALPHA."""
    if ELLIPSE_MODE in ('filled', 'both'):
        fig.plot(x=lons, y=lats, fill=color, close=True,
                 transparency=ELLIPSE_ALPHA,
                 region=region, projection=projection)

    if ELLIPSE_MODE in ('outline', 'both'):
        fig.plot(x=lons, y=lats, close=True,
                 pen=f"{ELLIPSE_PEN_WIDTH},{color},{ELLIPSE_DASH}",
                 region=region, projection=projection)


def ellipse_color_for(ev):
    """Resolve ELLIPSE_COLOR for a given event ('depth' → beach-ball colour)."""
    if ELLIPSE_COLOR == 'depth':
        return depth_to_color(ev.depth / 1000) if color_mode == 'depth_color' \
            else fixed_color
    return ELLIPSE_COLOR


# ═══════════════════════════════════════════════════════════════
#  FIGURE SETUP
# ═══════════════════════════════════════════════════════════════
fig = pygmt.Figure()
pygmt.config(
    FORMAT_GEO_MAP     = "ddd.xxF",
    FONT_ANNOT_PRIMARY = "9p,Helvetica,black",
    MAP_FRAME_PEN      = "0p,white@100",
)
fig.basemap(region=region, projection=projection,
            frame='a0.05', map_scale='x3c/-0.7c+w3')

# ═══════════════════════════════════════════════════════════════
#  TOPOGRAPHY
# ═══════════════════════════════════════════════════════════════
if topo_source == 'local':
    pygmt.config(COLOR_NAN=water_color)   # NaN (sea) cells → water_color

    buf     = 0.02
    raw     = rioxarray.open_rasterio(tif_path, masked=True).squeeze()
    raw_wgs = raw.rio.reproject('EPSG:4326')
    clipped = raw_wgs.rio.clip_box(
        minx=minlon - buf, miny=minlat - buf,
        maxx=maxlon + buf, maxy=maxlat + buf,
    )
    topo = xr.DataArray(
        data=clipped.values.astype('float32'),
        coords={'lat': clipped.y.values, 'lon': clipped.x.values},
        dims=['lat', 'lon'],
    )
    _cmap_name, _cmap_rev = _TOPO_PRESETS.get(topo_color, (topo_color, False))
    pygmt.makecpt(cmap=_cmap_name, series=[0, 600, 10], reverse=_cmap_rev)
    fig.grdimage(grid=topo, cmap=True, shading="+a315+ne0.6",
                 region=region, projection=projection)
    pygmt.config(COLOR_NAN="white")       # reset so other elements are unaffected

else:   # identical to original script 11
    topo = pygmt.datasets.load_earth_relief(resolution="01s", region=region)
    _cmap_name, _cmap_rev = _TOPO_PRESETS.get(topo_color, (topo_color, False))
    pygmt.makecpt(cmap=_cmap_name, series=[-200, 900, 10], reverse=_cmap_rev)
    fig.grdimage(grid=topo, region=region, projection=projection,
                 shading="+a45+ne0.5", cmap=True)
    fig.coast(shorelines="1/0.5p,black", resolution="f", water="#EBEBEE")

# ═══════════════════════════════════════════════════════════════
#  UNCERTAINTY ELLIPSES  (drawn first → below everything else)
# ═══════════════════════════════════════════════════════════════
all_semi_a_m, all_semi_b_m = [], []
missing = []

if ELLIPSE_MODE != 'none':
    ellipses = load_ellipses(os.path.join(catdir, ELLIPSE_FILE_NAME))

    for ev in fm_events:
        ep = ellipses.get(ev.name)
        if ep is None:
            missing.append(ev.name)
            continue

        ell_lons, ell_lats, semi_a_m, semi_b_m = make_ellipse_coords(
            ep['lat_ref'], ep['lon_ref'],
            ep['e16'], ep['e84'], ep['n16'], ep['n84'],
            n_pts=ELLIPSE_N_POINTS,
        )
        all_semi_a_m.append(semi_a_m)
        all_semi_b_m.append(semi_b_m)

        draw_ellipse(ell_lons, ell_lats, ellipse_color_for(ev))

    if missing:
        print(f"[WARNING] no ellipse for {len(missing)}/{len(fm_events)} events:")
        for name in missing:
            print(f"          - {name}")

    if all_semi_a_m:
        mean_a = float(np.mean(all_semi_a_m))
        mean_b = float(np.mean(all_semi_b_m))
        mean_area_km2 = np.pi * mean_a * mean_b / 1e6
        print(f"[INFO] ellipses plotted            : {len(all_semi_a_m)}")
        print(f"[INFO] mean semi-axis east–west    : {mean_a:.1f} m")
        print(f"[INFO] mean semi-axis north–south  : {mean_b:.1f} m")
        print(f"[INFO] mean ellipse area           : {mean_area_km2:.4f} km²")

        if switch_mean_ellipse:
            anchor_lat = minlat + MEAN_ELLIPSE_OFFSET_LAT
            anchor_lon = minlon + MEAN_ELLIPSE_OFFSET_LON

            deg_per_m_lat = 1.0 / 111320.0
            deg_per_m_lon = 1.0 / (111320.0 * np.cos(np.radians(anchor_lat)))

            theta = np.linspace(0.0, 2.0 * np.pi, ELLIPSE_N_POINTS)
            mean_lons = anchor_lon + mean_a * deg_per_m_lon * np.cos(theta)
            mean_lats = anchor_lat + mean_b * deg_per_m_lat * np.sin(theta)

            ref_color = 'gray40' if ELLIPSE_COLOR == 'depth' else ELLIPSE_COLOR
            draw_ellipse(mean_lons, mean_lats, ref_color)

            fig.text(
                text=f"mean ellipse ({mean_area_km2:.3f} km²)",
                x=anchor_lon,
                y=anchor_lat - mean_b * deg_per_m_lat - 0.002,
                font="5p,Helvetica,black", justify="CM",
            )

# ═══════════════════════════════════════════════════════════════
#  VT CIRCLES  (below the beach balls unless VT_ON_TOP)
# ═══════════════════════════════════════════════════════════════
if switch_vt_circles and not VT_ON_TOP:
    draw_vt_circles(vt_events)

# ═══════════════════════════════════════════════════════════════
#  FOCAL MECHANISMS
#  compressionfill is computed from depth in Python — no GMT CPT needed.
# ═══════════════════════════════════════════════════════════════
for ev in fm_events:
    size = mag_to_size(ev.magnitude)
    spec, convention = meca_spec(ev)
    fill = depth_to_color(ev.depth / 1000) if color_mode == 'depth_color' else fixed_color

    fig.meca(
        spec=spec, convention=convention,
        longitude=ev.lon, latitude=ev.lat, depth=ev.depth / 1000,
        scale=f"{size:.2f}c",
        compressionfill=fill, extensionfill="white",
        pen="0.5p,gray30,solid",
    )

    if switch_timestamps:
        parts = ev.name.split('_')[1:]
        label = f"{parts[0]}-{parts[1]}-{parts[2]} {parts[3]}:{parts[4]}:{parts[5]}"
        fig.text(text=label, x=ev.lon, y=ev.lat + 0.001,
                 font="4p,Helvetica,black", justify="CM")

if switch_vt_circles and VT_ON_TOP:
    draw_vt_circles(vt_events)

# ═══════════════════════════════════════════════════════════════
#  LEGEND — compact box, top-left
#  Stacked magnitude scales (beach balls and/or VT circles), each one a title
#  row + a row of circles + a row of labels, then the station symbol. Every
#  circle has the diameter its own symbol really gets on the map. The layout
#  is worked out in cm and converted to degrees at the end, so the box always
#  fits its content.
# ═══════════════════════════════════════════════════════════════
# cm → degrees. Mercator: 1 cm of latitude spans cos(lat) times the
# longitude span of 1 cm.
deg_lon_per_cm = (maxlon - minlon) / map_w_cm
deg_lat_per_cm = deg_lon_per_cm * np.cos(np.radians((minlat + maxlat) / 2))


def text_width_cm(text, font):
    """Rough width of a label (cm), from the point size in a GMT font string."""
    pt = float(font.split('p,')[0])
    return len(text) * 0.5 * pt * 2.54 / 72.0


def symbol_size_cm(style):
    """Size of a GMT symbol string such as 't0.3c' → 0.3 cm."""
    try:
        return float(style[1:].rstrip('c'))
    except ValueError:
        return 0.3


def mag_scale(title, mags, diam_fn, fill, pen):
    """
    Lay out one magnitude scale.

    Returns a dict with the circle diameters, their centres in cm from the
    left edge of the row, the row width and the total height of the block.
    """
    diams = [diam_fn(m) for m in mags]

    xs_cm, x_cm = [], diams[0] / 2.0
    for i, d in enumerate(diams):
        if i > 0:
            x_cm += max(diams[i - 1] / 2.0 + LEG_GAP + d / 2.0, LEG_MIN_STEP)
        xs_cm.append(x_cm)

    return {
        'title':  title,
        'mags':   mags,
        'diams':  diams,
        'xs_cm':  xs_cm,
        'fill':   fill,
        'pen':    pen,
        'w_cm':   xs_cm[-1] + diams[-1] / 2.0,
        'h_cm':   LEG_TITLE_H + max(diams) + LEG_LABEL_H,
    }


# --- collect the scales to draw -------------------------------------------
scales = []

if LEG_SHOW_MECA:
    meca_fill = depth_to_color((depth_min + depth_max) / 2) \
        if color_mode == 'depth_color' else fixed_color
    scales.append(mag_scale(LEG_MECA_TITLE, ref_mags_meca, meca_diameter_cm,
                            meca_fill, "0.4p,black"))

if LEG_SHOW_VT and switch_vt_circles:
    vt_leg_fill = None if VT_MODE == 'outline' else VT_COLOR
    if VT_MODE == 'outline':
        vt_leg_pen = f"{VT_PEN_WIDTH},{VT_COLOR}"
    elif VT_MODE == 'both':
        vt_leg_pen = f"{VT_PEN_WIDTH},{VT_PEN_COLOR}"
    else:
        vt_leg_pen = None
    scales.append(mag_scale(LEG_VT_TITLE, ref_mags_vt, vt_diameter_cm,
                            vt_leg_fill, vt_leg_pen))

# --- box geometry (cm) ----------------------------------------------------
sta_size_cm  = symbol_size_cm(STA_STYLE)
sta_row_w_cm = 0.0
if LEG_SHOW_STATION:
    sta_row_w_cm = (sta_size_cm + LEG_SYM_GAP
                    + text_width_cm(LEG_STATION_LABEL, LEG_FONT_LABEL))

if scales or LEG_SHOW_STATION:
    box_w_cm = max([sc['w_cm'] for sc in scales] + [sta_row_w_cm]) + 2 * LEG_PAD
    box_h_cm = (sum(sc['h_cm'] for sc in scales)
                + (LEG_STATION_H if LEG_SHOW_STATION else 0.0) + 2 * LEG_PAD)

    box_x0 = minlon + LEG_MARGIN * deg_lon_per_cm
    box_y1 = maxlat - LEG_MARGIN * deg_lat_per_cm      # top edge
    box_y0 = box_y1 - box_h_cm * deg_lat_per_cm        # bottom edge
    box_w  = box_w_cm * deg_lon_per_cm

    fig.plot(
        x=[box_x0, box_x0 + box_w, box_x0 + box_w, box_x0],
        y=[box_y0, box_y0, box_y1, box_y1],
        close=True, pen="0.4p,gray50", fill="white@20",
    )

    def leg_y(depth_cm):
        """Latitude of a point `depth_cm` below the top edge of the box."""
        return box_y1 - depth_cm * deg_lat_per_cm

    # --- one block per magnitude scale, stacked top to bottom -------------
    cur_cm = LEG_PAD
    for sc in scales:
        max_d = max(sc['diams'])
        y_title = leg_y(cur_cm + LEG_TITLE_H / 2.0)
        y_circ  = leg_y(cur_cm + LEG_TITLE_H + max_d / 2.0)
        y_label = leg_y(cur_cm + LEG_TITLE_H + max_d + LEG_LABEL_H / 2.0)

        fig.text(text=sc['title'], x=box_x0 + box_w / 2.0, y=y_title,
                 font=LEG_FONT_TITLE, justify="CM")

        row_x0_cm = (box_w_cm - sc['w_cm']) / 2.0      # centre the row
        for mag, d_cm, cx_cm in zip(sc['mags'], sc['diams'], sc['xs_cm']):
            xc = box_x0 + (row_x0_cm + cx_cm) * deg_lon_per_cm
            fig.plot(x=xc, y=y_circ, style=f"c{d_cm:.3f}c",
                     fill=sc['fill'], pen=sc['pen'])
            fig.text(text=f"{mag:g} Mw", x=xc, y=y_label,
                     font=LEG_FONT_LABEL, justify="CM")

        cur_cm += sc['h_cm']

    # --- station symbol + label -------------------------------------------
    if LEG_SHOW_STATION:
        y_sta     = leg_y(cur_cm + LEG_STATION_H / 2.0)
        sta_x0_cm = (box_w_cm - sta_row_w_cm) / 2.0
        fig.plot(x=box_x0 + (sta_x0_cm + sta_size_cm / 2.0) * deg_lon_per_cm,
                 y=y_sta, style=STA_STYLE, fill=STA_FILL, pen=STA_PEN)
        fig.text(text=LEG_STATION_LABEL,
                 x=box_x0 + (sta_x0_cm + sta_size_cm + LEG_SYM_GAP) * deg_lon_per_cm,
                 y=y_sta, font=LEG_FONT_LABEL, justify="LM")

# ═══════════════════════════════════════════════════════════════
#  STATIONS
# ═══════════════════════════════════════════════════════════════
latsta, lonsta, namsta = [], [], []
with open(os.path.join(metadatadir, STATION_FILE)) as f:
    for line in f:
        if line[0] == ' ':
            continue
        toks = line.split()
        latsta.append(float(toks[1]))
        lonsta.append(float(toks[2]))
        namsta.append(toks[0].split('.')[1])
latsta = np.array(latsta)
lonsta = np.array(lonsta)

fig.plot(x=lonsta, y=latsta, style=STA_STYLE, fill=STA_FILL, pen=STA_PEN)
if switch_sta_names:
    fig.text(x=lonsta + 0.002, y=latsta + 0.001, text=namsta,
             justify='BL', font='5p,Helvetica,black')

# ═══════════════════════════════════════════════════════════════
#  DEPTH COLORBAR — inside map, bottom-centre
#  Strategy:
#    1. shift_origin to bottom-centre of map interior
#    2. semi-transparent white box (background)
#    3. basemap ONCE with full frame — ticks extend outward (SNwe)
#    4. gradient rectangles fill the interior (ticks stay visible outside)
#    5. shift_origin restore
# ═══════════════════════════════════════════════════════════════
if color_mode == 'depth_color':
    cb_w   = 5.0    # colorbar width  (cm)
    cb_h   = 0.45   # colorbar strip height (cm)
    n_grad = 200    # gradient steps

    x_shift = (map_w_cm - cb_w) / 2    # left edge of colorbar from map origin

    # y_shift: distance from map bottom frame to bottom of colorbar strip
    # 0.5c up from the bottom puts it inside the map with room for ticks below
    y_shift = 0.5

    fig.shift_origin(xshift=f"{x_shift:.2f}c", yshift=f"{y_shift:.2f}c")

    # Semi-transparent background box (pad = 0.15c on all sides)
    pad = 0.15

    # ONE basemap call: sets coordinates AND draws frame (ticks go outward).
    # 'S' = annotated bottom axis, 'nwe' = plain lines (uppercase on more than
    # one side would repeat the "Depth (km)" label).
    fig.basemap(
        region=[depth_min, depth_max, 0, 1],
        projection=f"X{cb_w}c/{cb_h}c",
        frame=["Snwe", "a1f0.5", "x+lDepth (km)"],
    )

    # Background box (drawn after basemap so it's in the right coord system,
    # then gradient is drawn on top)
    pad_x = pad * (depth_max - depth_min) / cb_w
    fig.plot(
        x=[depth_min - pad_x, depth_max + pad_x, depth_max + pad_x, depth_min - pad_x],
        y=[-pad / cb_h, -pad / cb_h, 1 + pad / cb_h, 1 + pad / cb_h],
        close=True, fill="white@40", pen="0.4p,gray50",
    )

    # Gradient rectangles (light → dark, shallow → deep).
    # No `pen` at all: "0p" is not "no line", it is a hairline in the current
    # pen colour and turns the whole bar into black stripes. Slices overlap by
    # one step so no hairline gaps show through either.
    dz = (depth_max - depth_min) / n_grad
    for i in range(n_grad):
        z0 = depth_min + i * dz
        z1 = min(z0 + 1.5 * dz, depth_max)
        fig.plot(
            x=[z0, z1, z1, z0], y=[0, 0, 1, 1],
            close=True, fill=depth_to_color(z0),
        )

    # Restore origin
    fig.shift_origin(xshift=f"{-x_shift:.2f}c", yshift=f"{-y_shift:.2f}c")

# ═══════════════════════════════════════════════════════════════
#  SAVE
# ═══════════════════════════════════════════════════════════════
suffix   = 'deviatoric' if switch_deviatoric else 'dc'
ell_tag  = '' if ELLIPSE_MODE == 'none' else f"_ell-{ELLIPSE_MODE}"
vt_tag   = '_vt' if switch_vt_circles else ''
outpath  = os.path.join(plotdir, f"{filename}_{suffix}{ell_tag}{vt_tag}_{map_name}.pdf")
fig.show()
fig.savefig(outpath)
print(f"Saved: {outpath}")
