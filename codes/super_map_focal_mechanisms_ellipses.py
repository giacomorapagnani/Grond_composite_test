import pygmt
import numpy as np
import os
import xarray as xr
import rioxarray
from pyrocko import model
import pyrocko.moment_tensor as pmt
import pyrocko.orthodrome as od

workdir = '../'
catdir = os.path.join(workdir, 'CAT')
metadatadir = os.path.join(workdir, 'META_DATA')
plotdir = os.path.join(workdir, 'PLOTS', 'MAPS')
os.makedirs(plotdir, exist_ok=True)

# ═══════════════════════════════════════════════════════════════
#  SETTINGS
# ═══════════════════════════════════════════════════════════════

filename = 'catalogue_flegrei_composite_MT_LF_std_reloc_best_FLIPPED'
# catalogue_flegrei_composite_MT_LF_std_reloc_best
# catalogue_flegrei_MT_VT

#   'local'  → Tinitaly 10 m DEM (land-only, high-res)
#   'pygmt'  → GMT built-in 01s relief (~30 m, includes bathymetry)
topo_source = 'local'
tif_path    = '/Users/giaco/UNI/PhD_CODE/QGIS/topo_flegrei/w45090_s10.tif'
water_color = "#A9B0B8"   # sea / NaN colour — change freely

# Topography colour palette (tones + hillshading).
# Built-in presets — pick one or write a custom "color1,color2" string.
#   'gray'   → classic black-and-white
#   'brown'  → soft cream → warm sienna
#   'green'  → light mint → olive
#   'sand'   → pale ivory → dark tan
# You can also pass any GMT two-color gradient, e.g. 'lightyellow,chocolate'
topo_color = 'gray'

_TOPO_PRESETS = {
    'gray':  ('gray',               True),   # (GMT cmap, reverse)
    'brown': ('linen,sienna',       False),  # cream → warm brown
    'green': ('honeydew,olivedrab', False),  # light mint → olive
    'sand':  ('ivory,tan',          False),  # pale ivory → warm tan
}

#   'depth_color' → colour gradient by depth (palette set by meca_color)
#   'fixed'       → uniform fill set by fixed_color
color_mode  = 'depth_color'
fixed_color = '#606060'

# Focal mechanism colour palette (only used when color_mode == 'depth_color').
#   'gray' | 'blue' | 'red' | 'green' | 'brown' | 'purple'
# Or pass a custom (light_rgb, dark_rgb) tuple, e.g. meca_color = ((255,235,205),(139,69,19))
meca_color = 'red'

_MECA_PRESETS = {
    'gray':   ((205, 205, 205), (50,  50,  50)),    # shallow → deep
    'blue':   ((198, 219, 239), (8,   81,  156)),
    'red':    ((252, 187, 161), (165, 15,  21)),
    'green':  ((199, 233, 192), (0,   109, 44)),
    'brown':  ((253, 208, 162), (140, 81,  10)),
    'purple': ((218, 218, 235), (84,  39,  143)),
    'black':  ((255, 255, 255), (0,   0,   0)),
}

switch_deviatoric = True
switch_timestamps = False
switch_sta_names  = False

switch_pozzuoli = False
if switch_pozzuoli:
    minlon, maxlon = 14.07, 14.175
    minlat, maxlat = 40.79,  40.845
    map_name = 'pozzuoli'
else:
    minlon, maxlon = 14.07, 14.175
    minlat, maxlat = 40.775, 40.855
    map_name = 'gulf'

region     = [minlon, maxlon, minlat, maxlat]
projection = "M6i"

# ═══════════════════════════════════════════════════════════════
#  UNCERTAINTY ELLIPSES — STYLE SETTINGS
#  Ellipses are drawn BEFORE the beach balls, so they stay underneath.
# ═══════════════════════════════════════════════════════════════

# Display mode:
#   'outline' → dashed border only            (option 1)
#   'filled'  → semi-transparent solid fill   (option 2)
#   'both'    → transparent fill + dashed border
#   'none'    → ellipses switched off
ELLIPSE_MODE = 'filled'

# ───────────────────────────────────────────────────────────────
#  ↓↓↓  ELLIPSE COLOUR — CHANGE IT HERE  ↓↓↓
#  Any GMT colour: 'gray40', 'black', '#BD2025', '255/0/0' …
#  Special value 'depth' → each ellipse takes the same depth colour
#                          as its own beach ball.
ELLIPSE_COLOR = 'gray40'
#  ↑↑↑  ELLIPSE COLOUR — CHANGE IT HERE  ↑↑↑
# ───────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────
#  ↓↓↓  FILL TRANSPARENCY — CHANGE IT HERE  ↓↓↓
#  0   = fully opaque       (overlaps hide each other)
#  100 = fully transparent  (invisible)
#  70–85 works well when many ellipses overlap.
ELLIPSE_ALPHA = 80
#  ↑↑↑  FILL TRANSPARENCY — CHANGE IT HERE  ↑↑↑
# ───────────────────────────────────────────────────────────────

ELLIPSE_PEN_WIDTH = '0.8p'      # border thickness
ELLIPSE_DASH      = '4p_3p'     # dash pattern: 'solid' for a continuous line
ELLIPSE_N_POINTS  = 400         # points used to trace each ellipse
ELLIPSE_FILE_NAME = 'uncertainty_ellipses.txt'

# Reference ("mean") ellipse drawn in the bottom-left corner as a size legend.
switch_mean_ellipse     = True
MEAN_ELLIPSE_OFFSET_LON = 0.015   # degrees from minlon
MEAN_ELLIPSE_OFFSET_LAT = 0.015   # degrees from minlat


# ═══════════════════════════════════════════════════════════════
#  LOAD CATALOGUE
# ═══════════════════════════════════════════════════════════════
fm_events = model.load_events(os.path.join(catdir, filename + '.pf'))
fm_events = [ev for ev in fm_events if ev.moment_tensor is not None]

depths_km = np.array([ev.depth / 1000 for ev in fm_events])
depth_min = max(0.0, depths_km.min() - 1.0)
depth_max = depths_km.max() + 1.0


# ═══════════════════════════════════════════════════════════════
#  LOAD UNCERTAINTY ELLIPSES
#  columns: event_name lat_ref lon_ref east_16 east_84 north_16 north_84
#           (shifts in metres, relative to lat_ref/lon_ref)
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


# ═══════════════════════════════════════════════════════════════
#  DEPTH → COLOUR  (pure Python, independent of GMT CPT state)
#  shallow → light tone,  deep → dark tone  (palette set by meca_color)
# ═══════════════════════════════════════════════════════════════
def depth_to_color(depth_km):
    light, dark = meca_color if isinstance(meca_color, tuple) \
        else _MECA_PRESETS.get(meca_color, _MECA_PRESETS['gray'])
    t = np.clip((depth_km - depth_min) / (depth_max - depth_min), 0.0, 1.0)
    rgb = [round(light[i] + t * (dark[i] - light[i])) for i in range(3)]
    return f"{rgb[0]}/{rgb[1]}/{rgb[2]}"


def ellipse_color_for(ev):
    """Resolve ELLIPSE_COLOR for a given event ('depth' → beach-ball colour)."""
    if ELLIPSE_COLOR == 'depth':
        return depth_to_color(ev.depth / 1000) if color_mode == 'depth_color' \
            else fixed_color
    return ELLIPSE_COLOR


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
MECA_EXPONENT = 1     # exponent field of the 'mt' spec — see meca_diameter_cm()


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
#  UNCERTAINTY ELLIPSES  (drawn first → below the beach balls)
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

# ═══════════════════════════════════════════════════════════════
#  MAGNITUDE LEGEND  — compact horizontal, top-left
#  Circle diameters = meca_diameter_cm(), i.e. exactly the size the beach
#  balls get on the map. The whole layout is worked out in cm and converted
#  to degrees at the end, so the box always fits its content.
# ═══════════════════════════════════════════════════════════════
ref_mags = [1.0, 2.0, 3.0]

LEG_PAD      = 0.16    # cm — inner padding of the box
LEG_GAP      = 0.20    # cm — clear space between two circles
LEG_MIN_STEP = 0.62    # cm — min centre-to-centre distance (keeps labels apart)
LEG_TITLE_H  = 0.34    # cm — row reserved for the "Mag" title
LEG_LABEL_H  = 0.32    # cm — row reserved for the labels under the circles
LEG_FONT_TITLE = "5p,Helvetica-Bold,black"
LEG_FONT_LABEL = "4.5p,Helvetica,black"
LEG_MARGIN     = 0.10  # cm — distance from the map frame

# cm → degrees. Mercator: 1 cm of latitude spans cos(lat) times the
# longitude span of 1 cm.
map_w_cm       = 6 * 2.54                       # projection is M6i
deg_lon_per_cm = (maxlon - minlon) / map_w_cm
deg_lat_per_cm = deg_lon_per_cm * np.cos(np.radians((minlat + maxlat) / 2))

ref_diams = [meca_diameter_cm(m) for m in ref_mags]

# circle centres, in cm from the inner left edge of the box
xs_cm = []
x_cm  = ref_diams[0] / 2.0
for i, d in enumerate(ref_diams):
    if i > 0:
        x_cm += max(ref_diams[i - 1] / 2.0 + LEG_GAP + d / 2.0, LEG_MIN_STEP)
    xs_cm.append(x_cm)

box_w_cm = xs_cm[-1] + ref_diams[-1] / 2.0 + 2 * LEG_PAD
box_h_cm = LEG_TITLE_H + max(ref_diams) + LEG_LABEL_H + 2 * LEG_PAD

box_x0 = minlon + LEG_MARGIN * deg_lon_per_cm
box_y1 = maxlat - LEG_MARGIN * deg_lat_per_cm          # top edge
box_y0 = box_y1 - box_h_cm * deg_lat_per_cm            # bottom edge
box_w  = box_w_cm * deg_lon_per_cm

fig.plot(
    x=[box_x0, box_x0 + box_w, box_x0 + box_w, box_x0],
    y=[box_y0, box_y0, box_y1, box_y1],
    close=True, pen="0.4p,gray50", fill="white@20",
)

y_title = box_y1 - (LEG_PAD + LEG_TITLE_H / 2.0) * deg_lat_per_cm
y_circ  = box_y1 - (LEG_PAD + LEG_TITLE_H + max(ref_diams) / 2.0) * deg_lat_per_cm
y_label = box_y0 + (LEG_PAD + LEG_LABEL_H / 2.0) * deg_lat_per_cm

fig.text(text="Mag", x=box_x0 + box_w / 2.0, y=y_title,
         font=LEG_FONT_TITLE, justify="CM")

legend_fill = depth_to_color((depth_min + depth_max) / 2) \
    if color_mode == 'depth_color' else fixed_color

for rm, d_cm, cx_cm in zip(ref_mags, ref_diams, xs_cm):
    xc = box_x0 + (LEG_PAD + cx_cm) * deg_lon_per_cm
    fig.plot(x=xc, y=y_circ, style=f"c{d_cm:.3f}c",
             fill=legend_fill, pen="0.4p,black")
    fig.text(text=f"{int(rm)} Mw", x=xc, y=y_label,
             font=LEG_FONT_LABEL, justify="CM")

# ═══════════════════════════════════════════════════════════════
#  STATIONS
# ═══════════════════════════════════════════════════════════════
latsta, lonsta, namsta = [], [], []
with open(os.path.join(metadatadir, 'stations_flegrei_INGV_final.pf')) as f:
    for line in f:
        if line[0] == ' ':
            continue
        toks = line.split()
        latsta.append(float(toks[1]))
        lonsta.append(float(toks[2]))
        namsta.append(toks[0].split('.')[1])
latsta = np.array(latsta)
lonsta = np.array(lonsta)

fig.plot(x=lonsta, y=latsta, style="t0.3c", fill="#FFCC4E", pen="0.6p,black")
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

    # map_w_cm (M6i = 15.24 cm) is defined in the legend block above
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
outpath  = os.path.join(plotdir, f"{filename}_{suffix}{ell_tag}_{map_name}.pdf")
fig.show()
fig.savefig(outpath)
print(f"Saved: {outpath}")
