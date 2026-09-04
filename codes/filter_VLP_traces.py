#!/usr/bin/env python3
'''
Filtra e taglia le tracce VLP degli eventi del catalogo e le salva in DATA_FILT.

Per ogni evento NON commentato di CAT/catalogue_flegrei_VLP.pf:
    DATA/<event_name>/*.mseed  ->  DATA_FILT/<event_name>/<event_name>.mseed

Ordine delle operazioni (per traccia):
    1.  degap/deoverlap dei segmenti dello stesso canale
    2.  scarto le tracce che non coprono tutta la finestra richiesta
    3.  demean + detrend lineare sulla traccia intera (10 minuti)
    4.  bandpass causale: lowpass(FMAX) + highpass(FMIN), Butterworth ORDER poli
    5.  downsample a TARGET_DELTAT (anti-alias incluso)
    6.  taglio finale: [OT-T_PRE , OT+T_POST]

Le tracce restano in counts (nessuna deconvoluzione della risposta strumentale).
'''

import os
import glob
import numpy as num

from pyrocko import io, model, trace, util

# ----------------------------------------------------------------------------
# PARAMETRI
# ----------------------------------------------------------------------------

# banda del filtro [Hz]
FMIN = 0.075
FMAX = 0.125
ORDER = 4                 # poli del Butterworth (causale, come in VLP_trace_similarity)

# finestra di taglio rispetto all'origin time [s]
T_PRE = 60.
T_POST = 180.

# ricampionamento: intervallo di campionamento finale [s] (0.05 = 20 Hz)
# metti None per mantenere il dt nativo
TARGET_DELTAT = 0.05

# gap massimo (in campioni) richiuso per interpolazione
MAX_GAP_SAMPLES = 5

# stazioni da tenere; metti None per tenere tutte quelle presenti in DATA
STATIONS = [
    'CAAM', 'CBAC', 'CFMN', 'CMIS', 'CMSN', 'CNIS',
    'COLB', 'CPIS', 'CPOZ', 'CQUE', 'CSFT', 'CSOB', 'CSTH']

# se True, scarta le stazioni elencate in BLACKLIST_VLP/<event_name>.black
USE_BLACKLIST = False

# ----------------------------------------------------------------------------
# PERCORSI
# ----------------------------------------------------------------------------

workdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')

catname = os.path.join(workdir, 'CAT', 'catalogue_flegrei_VLP.pf')
datadir = os.path.join(workdir, 'DATA')
blackdir = os.path.join(workdir, 'BLACKLIST_VLP')
outdir = os.path.join(workdir, 'DATA_FILT')


def load_event_traces(ev_dir):
    '''Carica tutti i .mseed di una cartella evento.

    Gestisce entrambi i layout presenti in DATA: un unico file per evento
    oppure un file per stazione.
    '''
    trs = []
    for fn in sorted(glob.glob(os.path.join(ev_dir, '*.mseed'))):
        trs.extend(io.load(fn))

    return trs


def load_blacklist(event_name):
    '''Legge le stazioni da escludere per un evento (una per riga).'''
    fn = os.path.join(blackdir, event_name + '.black')
    if not os.path.isfile(fn):
        return set()

    with open(fn, 'r') as f:
        return set(line.strip() for line in f if line.strip())


def detrend(tr):
    '''Rimuove media e trend lineare (in place).'''
    y = tr.ydata.astype(num.float64)
    x = num.arange(y.size, dtype=num.float64)
    slope, offset = num.polyfit(x, y, 1)
    tr.set_ydata(y - (slope*x + offset))


def process_trace(tr, ev_time):
    '''Applica la catena di processing a una traccia. Ritorna None se scartata.'''

    tcut_min = ev_time - T_PRE
    tcut_max = ev_time + T_POST
    tol = tr.deltat

    # 2. la traccia deve coprire tutta la finestra finale
    if tr.tmin > tcut_min + tol or tr.tmax < tcut_max - tol:
        return None

    tr = tr.copy()

    # 3. demean + detrend sulla traccia intera
    detrend(tr)

    # 4. bandpass causale, prima lowpass poi highpass
    tr.lowpass(ORDER, FMAX, demean=False)
    tr.highpass(ORDER, FMIN, demean=False)

    # 5. downsample (pyrocko applica da solo il filtro anti-alias)
    if TARGET_DELTAT is not None and tr.deltat < TARGET_DELTAT:
        tr.downsample_to(TARGET_DELTAT, snap=True, demean=False)

    # 6. taglio finale
    tr = tr.chop(tcut_min, tcut_max, inplace=False)

    # float32: dimezza i file, precisione ampiamente sufficiente sui counts
    tr.set_ydata(tr.ydata.astype(num.float32))

    return tr


def main():
    events = model.load_events(catname)   # le righe con '#' sono gia' saltate
    print('eventi nel catalogo (non commentati): %d' % len(events))
    print('banda %.3f - %.3f Hz | finestra OT-%.0f s / OT+%.0f s\n'
          % (FMIN, FMAX, T_PRE, T_POST))

    for ev in events:
        ev_dir = os.path.join(datadir, ev.name)
        if not os.path.isdir(ev_dir):
            print('%-32s cartella dati mancante' % ev.name)
            continue

        trs = load_event_traces(ev_dir)
        if not trs:
            print('%-32s nessun .mseed trovato' % ev.name)
            continue

        # 1. richiude i piccoli gap e risolve le sovrapposizioni
        trs = trace.degapper(
            trs, maxgap=MAX_GAP_SAMPLES, fillmethod='interpolate')

        blacklist = load_blacklist(ev.name) if USE_BLACKLIST else set()
        keep = set(STATIONS) if STATIONS is not None else None

        trs_out = []
        nslc_in = set()
        nslc_out = set()
        for tr in trs:
            if keep is not None and tr.station not in keep:
                continue

            if tr.station in blacklist:
                continue

            nslc_in.add(tr.nslc_id)
            tr_proc = process_trace(tr, ev.time)
            if tr_proc is not None:
                nslc_out.add(tr.nslc_id)
                trs_out.append(tr_proc)

        # canali per cui nessun segmento copre tutta la finestra
        skipped = ['.'.join(k) for k in nslc_in - nslc_out]

        if not trs_out:
            print('%-32s nessuna traccia utilizzabile' % ev.name)
            continue

        ev_outdir = os.path.join(outdir, ev.name)
        util.ensuredirs(os.path.join(ev_outdir, 'dummy'))
        fn_out = os.path.join(ev_outdir, ev.name + '.mseed')
        io.save(trs_out, fn_out)

        sta_out = sorted(set(tr.station for tr in trs_out))
        msg = '%-32s %2d staz, %3d tracce -> %s' % (
            ev.name, len(sta_out), len(trs_out), os.path.relpath(fn_out, workdir))
        if keep is not None:
            missing = sorted(keep - set(sta_out))
            if missing:
                msg += '   [assenti: %s]' % ','.join(missing)
        if skipped:
            msg += '   [scartate: %s]' % ', '.join(sorted(skipped))
        print(msg)


if __name__ == '__main__':
    main()
