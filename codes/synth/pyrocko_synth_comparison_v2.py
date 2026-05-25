# COMPARE SYNTHETIC TRACES (VT AND VLP SOLUTIONS) WITH OBSERVED TRACES
# Two figures, N component only.
#
# Figure 1 — 2x2 subplot
#   rows: VT (top), VLP (bottom)
#   cols: HF 0.5–2.0 Hz (left), LF 0.075–0.125 Hz (right)
#   sharex per column, sharey per row
#
# Figure 2 — 4x1 subplot, LF only
#   row 1: VT synth
#   row 2: VLP synth
#   row 3: VT+VLP sum
#   row 4: Recorded
#   sharex and sharey global, ylim from max amplitude across all 4 traces

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime

from pyrocko import gf, trace, model, io
from pyrocko.gf import LocalEngine, Target, MTSource

# ─── Directories ──────────────────────────────────────────────────────────────
catdir        = '../../CAT/synth'
metadatadir   = '../../META_DATA'
datadir       = '../../DATA_RESPONSE'

# ─── Load metadata ────────────────────────────────────────────────────────────
stations_name = os.path.join(metadatadir, 'stations_flegrei_INGV_final.pf')
stations      = model.load_stations(stations_name)

catname_VLP   = os.path.join(catdir, 'catalogue_flegrei_test_VLP.pf')
events_VLP    = model.load_events(catname_VLP)

catname_VT    = os.path.join(catdir, 'catalogue_flegrei_test_VT.pf')
events_VT     = model.load_events(catname_VT)

# ─── Select station ───────────────────────────────────────────────────────────
s_name = 'CMIS'
st = next((s for s in stations if s.station == s_name), None)
if st is None:
    raise ValueError(f'Station {s_name} not found')

# ─── Select event ─────────────────────────────────────────────────────────────
e_name = 'flegrei_2023_06_11_06_44_25'

ev_VT = next((e for e in events_VT if e.name == e_name), None)
if ev_VT is None:
    raise ValueError(f'Event {e_name} not found in {catname_VT}')

ev_VLP = next((e for e in events_VLP if e.name == e_name), None)
if ev_VLP is None:
    raise ValueError(f'Event {e_name} not found in {catname_VLP}')

# ─── GF engine & targets ──────────────────────────────────────────────────────
store_id = 'campiflegrei_near_0_dist'
engine   = LocalEngine(store_superdirs=['../../GF_STORES'])

channel_codes = 'ENZ'

targets_VT = [
    Target(lat=st.lat, lon=st.lon, store_id=store_id,
           codes=('', st.station, 'VT', ch))
    for ch in channel_codes]

targets_VLP = [
    Target(lat=st.lat, lon=st.lon, store_id=store_id,
           codes=('', st.station, 'VLP', ch))
    for ch in channel_codes]

# ─── Sources ──────────────────────────────────────────────────────────────────
source_mt_VT  = MTSource.from_pyrocko_event(ev_VT)

source_mt_VLP = MTSource.from_pyrocko_event(ev_VLP)
source_mt_VLP.stf = gf.ResonatorSTF(25., frequency=0.114)

# ─── Compute synthetics ───────────────────────────────────────────────────────
response_VT  = engine.process(source_mt_VT,  targets_VT)
response_VLP = engine.process(source_mt_VLP, targets_VLP)

synthetic_traces_VT  = response_VT.pyrocko_traces()
synthetic_traces_VLP = response_VLP.pyrocko_traces()

# ─── Mirror-padding helper ────────────────────────────────────────────────────
def mirror_pad(tr, pad_seconds):
    """
    Extend a trace with mirror-reflected padding to reduce edge artefacts
    when applying recursive IIR filters.  Pure numpy + pyrocko.
    """
    dt   = tr.deltat
    npad = int(round(pad_seconds / dt))
    y    = tr.get_ydata()
    # reflect first/last npad samples (clip to available length)
    n    = min(npad, len(y) - 1)
    head = y[1:n+1][::-1]          # mirror of the beginning
    tail = y[-(n+1):-1][::-1]      # mirror of the end
    y_padded = np.concatenate((head, y, tail))
    tr_padded = tr.copy()
    tr_padded.set_ydata(y_padded)
    tr_padded.tmin = tr.tmin - n * dt
    return tr_padded, n             # also return the number of samples added


def pad_and_detrend(trs, pad_seconds=240.):
    """Apply mirror padding and remove linear trend (mean + slope)."""
    padded = []
    npad_list = []
    for tr in trs:
        # remove linear trend before padding to avoid slope artefacts
        y = tr.get_ydata().astype(float)
        x = np.arange(len(y))
        p = np.polyfit(x, y, 1)
        y -= np.polyval(p, x)
        tr_detr = tr.copy()
        tr_detr.set_ydata(y)
        tr_p, n = mirror_pad(tr_detr, pad_seconds)
        padded.append(tr_p)
        npad_list.append(n)
    return padded, npad_list


def unpad(tr, npad):
    """Remove the mirror-padding samples after filtering."""
    y   = tr.get_ydata()
    dt  = tr.deltat
    y_c = y[npad: len(y) - npad]
    tr_out = tr.copy()
    tr_out.set_ydata(y_c)
    tr_out.tmin = tr.tmin + npad * dt
    return tr_out


# ─── Pad synthetics ───────────────────────────────────────────────────────────
pad_sec = 300.    # generous padding — will be stripped after filtering

trs_VT_pad,  npad_VT  = pad_and_detrend(synthetic_traces_VT,  pad_sec)
trs_VLP_pad, npad_VLP = pad_and_detrend(synthetic_traces_VLP, pad_sec)

# ─── Build VT+VLP sum (N component, index 1 in ENZ) ──────────────────────────
# We work in the padded domain so that the sum is also filtered correctly later.
# Align on tmin of VLP (usually the longer one).
def sum_traces(tr_VLP, tr_VT):
    """Sum two traces; they may have different tmin — align sample-exactly."""
    dt       = tr_VLP.deltat
    tshift   = int(round((tr_VT.tmin - tr_VLP.tmin) / dt))
    y_vlp    = tr_VLP.get_ydata().copy()
    y_vt     = tr_VT.get_ydata()
    i0       = max(tshift, 0)
    i0_vt    = max(-tshift, 0)
    n_copy   = min(len(y_vt) - i0_vt, len(y_vlp) - i0)
    if n_copy > 0:
        y_vlp[i0: i0 + n_copy] += y_vt[i0_vt: i0_vt + n_copy]
    tr_sum = tr_VLP.copy()
    tr_sum.set_ydata(y_vlp)
    return tr_sum


# N component is index 1 (E=0, N=1, Z=2)
N_IDX = 1

tr_VT_N_pad  = trs_VT_pad[N_IDX]
tr_VLP_N_pad = trs_VLP_pad[N_IDX]
tr_sum_N_pad = sum_traces(tr_VLP_N_pad, tr_VT_N_pad)

# ─── Load observed traces ─────────────────────────────────────────────────────
dir_name  = os.path.join(datadir, e_name)
file_name = os.path.join(dir_name, e_name + '.mseed')
obs_all   = io.load(file_name)
obs_N     = next(
    (tr for tr in obs_all if tr.station == s_name and tr.channel.endswith('N')),
    None)
if obs_N is None:
    raise ValueError(f'No N-component trace found for station {s_name}')
obs_N.location = 'Recorded'

# ─── Filter + chop helper ─────────────────────────────────────────────────────
def filter_and_chop(tr, fmin, fmax, npad, o_t, chop_before, chop_after,
                    remove_pad=True):
    """
    Filter a (possibly padded) trace with 4th-order Butterworth bandpass,
    strip the mirror padding, then chop around the origin time.
    Returns a plain numpy array (ydata) and a corresponding time array (seconds
    relative to o_t).
    """
    tmp = tr.copy()
    tmp.bandpass(4, fmin, fmax)       # pyrocko Butterworth bandpass

    if remove_pad:
        tmp = unpad(tmp, npad)

    # chop
    tmp.chop(o_t - chop_before, o_t + chop_after)

    y   = tmp.get_ydata().astype(float)
    # time axis in seconds relative to origin time
    t0  = tmp.tmin - o_t
    tax = t0 + np.arange(len(y)) * tmp.deltat
    return tax, y


# Frequency bands
HF = (0.5,  2.0)
LF = (0.075, 0.125)

# Chop windows
chop_HF = (5,  15)     # seconds before / after origin
chop_LF = (30, 150)

o_t = ev_VT.time       # origin time (float, UTC seconds)

# ─── Compute all filtered/chopped traces ──────────────────────────────────────
# Synthetics use mirror-padded versions; obs uses its own data (pad it too)
obs_N_pad, npad_obs = pad_and_detrend([obs_N], pad_sec)
obs_N_pad = obs_N_pad[0]
npad_obs  = npad_obs[0]

# VT
t_VT_HF, y_VT_HF = filter_and_chop(tr_VT_N_pad,  *HF, npad_VT[N_IDX],  o_t, *chop_HF)
t_VT_LF, y_VT_LF = filter_and_chop(tr_VT_N_pad,  *LF, npad_VT[N_IDX],  o_t, *chop_LF)

# VLP
t_VLP_HF, y_VLP_HF = filter_and_chop(tr_VLP_N_pad, *HF, npad_VLP[N_IDX], o_t, *chop_HF)
t_VLP_LF, y_VLP_LF = filter_and_chop(tr_VLP_N_pad, *LF, npad_VLP[N_IDX], o_t, *chop_LF)

# SUM
t_sum_LF, y_sum_LF = filter_and_chop(tr_sum_N_pad,  *LF, npad_VLP[N_IDX], o_t, *chop_LF)

# OBS
t_obs_LF, y_obs_LF = filter_and_chop(obs_N_pad, *LF, npad_obs, o_t, *chop_LF)

# ─── Convert time axis to datetime for nicer x labels ─────────────────────────
def to_dt(t_rel):
    """Convert relative-seconds array to list of datetime objects."""
    base = datetime.datetime.utcfromtimestamp(o_t)
    return [base + datetime.timedelta(seconds=float(s)) for s in t_rel]

dt_VT_HF  = to_dt(t_VT_HF)
dt_VT_LF  = to_dt(t_VT_LF)
dt_VLP_HF = to_dt(t_VLP_HF)
dt_VLP_LF = to_dt(t_VLP_LF)
dt_sum_LF = to_dt(t_sum_LF)
dt_obs_LF = to_dt(t_obs_LF)

# ─── Figure 1 — 2×2 ──────────────────────────────────────────────────────────
# sharex per column (same freq band), sharey per row (same source type)
fig1, axs1 = plt.subplots(
    2, 2,
    figsize=(14, 7),
    sharex='col',
    sharey='row')

col_titles = [f'HF  {HF[0]}–{HF[1]} Hz', f'LF  {LF[0]}–{LF[1]} Hz']
row_labels  = ['VT synthetic', 'VLP synthetic']
colors_fig1 = ['#BD2025', '#2563EB']   # red for VT, blue for VLP

datasets = [
    [(dt_VT_HF,  y_VT_HF),  (dt_VT_LF,  y_VT_LF)],
    [(dt_VLP_HF, y_VLP_HF), (dt_VLP_LF, y_VLP_LF)],
]

for row in range(2):
    for col in range(2):
        ax  = axs1[row, col]
        dt_, y_ = datasets[row][col]
        ax.plot(dt_, y_, color=colors_fig1[row], linewidth=1.5)
        ax.axvline(datetime.datetime.utcfromtimestamp(o_t),
                   color='k', linewidth=0.8, linestyle='--', alpha=0.5)
        ax.grid(True, alpha=0.3)
        ax.set_ylabel('Displacement [m]', fontsize=9)
        if row == 0:
            ax.set_title(col_titles[col], fontsize=11, fontweight='bold')
        if row == 1:
            ax.set_xlabel('Time (UTC)', fontsize=9)
        # row label on the left side
        if col == 0:
            ax.annotate(row_labels[row], xy=(0, 0.5),
                        xycoords='axes fraction',
                        xytext=(-60, 0), textcoords='offset points',
                        va='center', ha='center', fontsize=10,
                        rotation=90, fontweight='bold',
                        color=colors_fig1[row])
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha='right')

fig1.suptitle(f'{e_name}  |  Station: {s_name}  |  N component',
              fontsize=12, fontweight='bold', y=1.01)
fig1.tight_layout()

# ─── Figure 2 — 4×1, LF only, shared x and y ─────────────────────────────────
traces_fig2 = [
    (dt_VT_LF,  y_VT_LF,  'VT synthetic',  '#BD2025'),
    (dt_VLP_LF, y_VLP_LF, 'VLP synthetic', '#2563EB'),
    (dt_sum_LF, y_sum_LF, 'VT + VLP',      '#FF7400'),
    (dt_obs_LF, y_obs_LF, 'Recorded',       '#22863A'),
]

# Global ylim: symmetric around 0, based on the max absolute amplitude
global_ymax = max(np.max(np.abs(y)) for _, y, _, _ in traces_fig2) * 1.15

# Global xlim: widest time range among the 4 traces
global_xmin = min(dt[0]  for dt, _, _, _ in traces_fig2)
global_xmax = max(dt[-1] for dt, _, _, _ in traces_fig2)

fig2, axs2 = plt.subplots(4, 1, figsize=(12, 12), sharex=True, sharey=True)

for i, (dt_, y_, label, color) in enumerate(traces_fig2):
    ax = axs2[i]
    ax.plot(dt_, y_, color=color, linewidth=1.5, label=label)
    ax.axvline(datetime.datetime.utcfromtimestamp(o_t),
               color='k', linewidth=0.8, linestyle='--', alpha=0.5,
               label='Origin time')
    ax.set_ylim(-global_ymax, global_ymax)
    ax.set_xlim(global_xmin, global_xmax)
    ax.set_ylabel('Displacement [m]', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=9)
    if i < 3:
        plt.setp(ax.xaxis.get_majorticklabels(), visible=False)
    else:
        ax.set_xlabel('Time (UTC)', fontsize=9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha='right')

fig2.suptitle(
    f'{e_name}  |  Station: {s_name}  |  N component  |  '
    f'LF  {LF[0]}–{LF[1]} Hz',
    fontsize=12, fontweight='bold')
fig2.tight_layout()

plt.show()
