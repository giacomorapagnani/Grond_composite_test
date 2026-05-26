# COMPARE SYNTHETIC TRACES (VT AND VLP) WITH OBSERVED TRACES
# N component only — 2 figures
#
# Figure 1 — 2x2:  rows=source (VT,VLP), cols=freq (HF,LF)
#   sharex per column, sharey per row
#
# Figure 2 — 4x1, LF only:  VT / VLP / VT+VLP / Recorded
#   sharex and sharey global (ylim from max amplitude)

import os
import numpy as np
import matplotlib.pyplot as plt
import datetime

from pyrocko import gf, trace, model, io
from pyrocko.gf import LocalEngine, Target, MTSource

# ── dirs ──────────────────────────────────────────────────────────────────────
catdir      = '../../CAT/synth'
metadatadir = '../../META_DATA'
datadir     = '../../DATA_RESPONSE'

# ── load metadata ─────────────────────────────────────────────────────────────
stations_name = os.path.join(metadatadir, 'stations_flegrei_INGV_final.pf')
stations      = model.load_stations(stations_name)

catname_VLP = os.path.join(catdir, 'catalogue_flegrei_test_VLP.pf')
events_VLP  = model.load_events(catname_VLP)

catname_VT = os.path.join(catdir, 'catalogue_flegrei_test_VT.pf')
events_VT  = model.load_events(catname_VT)

# ── select station ────────────────────────────────────────────────────────────
s_name = 'CMIS'                         # CHANGE
st = False
for s in stations:
    if s.station == s_name:
        st = s
if not st:
    print(f'Error: station {s_name} not found')

# ── select event ──────────────────────────────────────────────────────────────
e_name = 'flegrei_2023_06_11_06_44_25'  # CHANGE

ev_VT = False
for e in events_VT:
    if e.name == e_name:
        ev_VT = e
if not ev_VT:
    print(f'Error: event {e_name} not found in {catname_VT}')

ev_VLP = False
for e in events_VLP:
    if e.name == e_name:
        ev_VLP = e
if not ev_VLP:
    print(f'Error: event {e_name} not found in {catname_VLP}')

# ── GF engine & targets ───────────────────────────────────────────────────────
store_id = 'campiflegrei_near_0_dist'   # CHANGE
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

# ── sources ───────────────────────────────────────────────────────────────────
source_mt_VT  = MTSource.from_pyrocko_event(ev_VT)
source_mt_VLP = MTSource.from_pyrocko_event(ev_VLP)
source_mt_VLP.stf = gf.ResonatorSTF(25., frequency=0.114)

# ── compute synthetics ────────────────────────────────────────────────────────
response_VT  = engine.process(source_mt_VT,  targets_VT)
response_VLP = engine.process(source_mt_VLP, targets_VLP)

synthetic_traces_VT  = response_VT.pyrocko_traces()
synthetic_traces_VLP = response_VLP.pyrocko_traces()

# ── chop edges then add constant padding (same as original) ──────────────────
trs_VT  = [x.chop(x.tmin+1, x.tmax-1) for x in synthetic_traces_VT]
trs_VLP = [x.chop(x.tmin+1, x.tmax-1) for x in synthetic_traces_VLP]

tlen = 180.
newtrs_VT = []
for tr in trs_VT:
    newtr = tr.copy()
    ydata = newtr.get_ydata()
    ydata = ydata - np.mean(ydata)
    first, last = ydata[0], ydata[-1]
    npts  = int(tlen / newtr.deltat)
    ydata = np.concatenate((np.ones(npts)*first, ydata, np.ones(npts)*last))
    newtr.ydata = ydata
    newtr.shift(-tlen)
    newtr.tmax += 2*tlen
    newtrs_VT.append(newtr)

tlen = 240.
newtrs_VLP = []
for tr in trs_VLP:
    newtr = tr.copy()
    ydata = newtr.get_ydata()
    first, last = ydata[0], ydata[-1]
    npts  = int(tlen / newtr.deltat)
    ydata = np.concatenate((np.ones(npts)*first, ydata, np.ones(npts)*last))
    newtr.ydata = ydata
    newtr.shift(-tlen)
    newtr.tmax += 2*tlen
    newtrs_VLP.append(newtr)

# ── sum VT + VLP ──────────────────────────────────────────────────────────────
trs_sum = []
channels = ['E', 'N', 'Z']
for n, ch in enumerate(channels):
    dt     = newtrs_VLP[n].deltat
    tmin   = newtrs_VLP[n].tmin
    tshift = int((newtrs_VT[n].tmin - newtrs_VLP[n].tmin) / dt)
    len_tr_VT = len(newtrs_VT[n].get_ydata())
    tr1    = newtrs_VLP[n].get_ydata()
    tr2    = newtrs_VT[n].get_ydata()
    trsum  = tr1.copy()
    trsum[tshift:tshift+len_tr_VT] += tr2
    trs_sum.append(trace.Trace(
        station=s_name, channel=ch, location='VT+VLP',
        deltat=dt, tmin=tmin, ydata=trsum))

# ── load observed traces ──────────────────────────────────────────────────────
dir_name  = os.path.join(datadir, e_name)
file_name = os.path.join(dir_name, e_name + '.mseed')
obs_trs_all = io.load(file_name)
obs_trs = [tr for tr in obs_trs_all if tr.station == s_name]
for tr in obs_trs:
    tr.location = 'Recorded'

# ── select N component (index 1 in ENZ order) ────────────────────────────────
N = 1   # E=0, N=1, Z=2         #CHANGE

tr_VT_N   = newtrs_VT[N]
tr_VLP_N  = newtrs_VLP[N]
tr_sum_N  = trs_sum[N]
# observed: find the trace whose channel ends with 'N'
tr_obs_N  = next((tr for tr in obs_trs if tr.channel.endswith('N')), None)
if tr_obs_N is None:
    print(f'Error: no N-component observed trace found for station {s_name}')

# ── filter + chop helper (same logic as original) ────────────────────────────
def filter_chop_plot(ax, tr, fq, o_t, color, label):
    """
    Chop first (with extra margin), filter, chop to final window,
    shift tmin to 23:00:00, check array lengths, plot.
    Returns (eq_dates, yax) so the caller can inspect amplitudes if needed.
    """
    tmp = tr.copy()

    if fq == [0.5, 2]:
        chop1, chop2 = 5, 15
        taper_margin = 10.      # extra seconds around final window before filtering
    elif fq == [0.075, 0.125]:
        chop1, chop2 = 30, 150
        taper_margin = 60.      # longer margin for low-frequency filter
    else:
        print('Error: frequency range not recognised')
        chop1, chop2 = 30, 150
        taper_margin = 60.

    # 1) chop with margin so the filter sees real signal, not padding edges
    tmp.chop(o_t - chop1 - taper_margin, o_t + chop2 + taper_margin)
    # 2) filter on the margined trace
    tmp.lowpass(4,  fq[1])
    tmp.highpass(4, fq[0])
    # 3) chop to the final window
    tmp.chop(o_t - chop1, o_t + chop2)

    new_tmin = 23*60*60.
    new_tmax = new_tmin + (tmp.tmax - tmp.tmin)
    tmp.tmin = new_tmin
    tmp.tmax = new_tmax

    tax  = np.arange(new_tmin, new_tmax, tmp.deltat)
    yax  = tmp.get_ydata()
    len_t, len_y = len(tax), len(yax)
    if len_t < len_y:
        eq_dates = [datetime.datetime.fromtimestamp(t) for t in tax]
        yax = yax[:len_t]
        print(f'Warning: {label} tax shorter than ydata ({len_t} vs {len_y})')
    elif len_t > len_y:
        eq_dates = [datetime.datetime.fromtimestamp(t) for t in tax[:len_y]]
        print(f'Warning: {label} tax longer than ydata ({len_t} vs {len_y})')
    else:
        eq_dates = [datetime.datetime.fromtimestamp(t) for t in tax]

    ax.plot(eq_dates, yax, color=color, linewidth=1.5, label=label)
    ax.grid(True)
    ax.legend(loc=1)
    ax.set_ylabel('Displacement [m]')

    return eq_dates, yax


o_t = ev_VT.time   # origin time in seconds (from .pf catalogue)

HF = [0.5,   2.0  ]
LF = [0.075, 0.125]

colors = {
    'VT':       '#BD2025',
    'VLP':      '#FFCC4E', #2563EB
    'VT+VLP':   '#FF7400',
    'Recorded': '#22863A',
}

# ════════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — 2x2
#   row 0: VT   | col 0: HF  col 1: LF
#   row 1: VLP  | col 0: HF  col 1: LF
#   sharex='col', sharey='row'
# ════════════════════════════════════════════════════════════════════════════════
fig1, axs1 = plt.subplots(2, 2, figsize=(14, 7),
                           sharex='col', sharey='col')

# row 0 — VT
filter_chop_plot(axs1[0,0], tr_VT_N,  HF, o_t, colors['VT'],  'VT')
filter_chop_plot(axs1[0,1], tr_VT_N,  LF, o_t, colors['VT'],  'VT')
# row 1 — VLP
filter_chop_plot(axs1[1,0], tr_VLP_N, HF, o_t, colors['VLP'], 'VLP')
filter_chop_plot(axs1[1,1], tr_VLP_N, LF, o_t, colors['VLP'], 'VLP')

axs1[0,0].set_title(f'HF  {HF[0]}–{HF[1]} Hz', fontsize=12, fontweight='bold')
axs1[0,1].set_title(f'LF  {LF[0]}–{LF[1]} Hz', fontsize=12, fontweight='bold')
for row, lbl in enumerate(['VT synthetic', 'VLP synthetic']):
    axs1[row,0].set_ylabel(f'{lbl}\nDisplacement [m]')
for col in range(2):
    axs1[1,col].set_xlabel('Time')

fig1.suptitle(f'{e_name}  |  Station: {s_name}  |  N component',
              fontsize=11, fontweight='bold')
fig1.tight_layout()

# ════════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — 4x1, LF only
#   row 0: VT   row 1: VLP   row 2: VT+VLP   row 3: Recorded
#   sharex and sharey global; ylim = ± max amplitude × 1.15
# ════════════════════════════════════════════════════════════════════════════════
fig2, axs2 = plt.subplots(4, 1, figsize=(12, 12),
                           sharex=True, sharey=True)

traces_fig2 = [
    (tr_VT_N,   colors['VT'],       'VT'),
    (tr_VLP_N,  colors['VLP'],      'VLP'),
    (tr_sum_N,  colors['VT+VLP'],   'VT+VLP'),
    (tr_obs_N,  colors['Recorded'], 'Recorded'),
]

# first pass: collect all ydata to compute global ylim
all_y = []
for tr, col, lbl in traces_fig2:
    tmp = tr.copy()
    tmp.chop(o_t - 30 - 60., o_t + 150 + 60.)
    tmp.lowpass(4,  LF[1])
    tmp.highpass(4, LF[0])
    tmp.chop(o_t - 30, o_t + 150)
    all_y.append(tmp.get_ydata())

global_ymax = max(np.max(np.abs(y)) for y in all_y) * 1.15

# second pass: plot
for i, (tr, col, lbl) in enumerate(traces_fig2):
    filter_chop_plot(axs2[i], tr, LF, o_t, col, lbl)
    axs2[i].set_ylim(-global_ymax, global_ymax)
    if i < 3:
        axs2[i].set_xlabel('')
axs2[3].set_xlabel('Time')

fig2.suptitle(
    f'{e_name}  |  Station: {s_name}  |  N component  |  '
    f'LF  {LF[0]}–{LF[1]} Hz',
    fontsize=11, fontweight='bold')
fig2.tight_layout()

plt.show()