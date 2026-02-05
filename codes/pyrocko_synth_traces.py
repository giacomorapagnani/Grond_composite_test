# GENERATE SYNTHETIC TRACES USING MT SOLUTIONS (FOR VT AND VLP EVENTS)

# libs
import os

from pyrocko import gf,trace,model,util,io
from pyrocko.gf import LocalEngine, Target, DCSource, ws,MTSource
from pyrocko.gui.marker import PhaseMarker

import numpy as np
import matplotlib.pyplot as plt
import datetime

# dirs
catdir='../CAT'
metadatadir='../META_DATA'

stations_name=os.path.join(metadatadir,'stations_flegrei_INGV_final.pf')
stations=model.load_stations(stations_name)

catname_VLP=os.path.join(catdir,'catalogue_flegrei_VLP.pf')
events_VLP = model.load_events(catname_VLP)

catname_VT=os.path.join(catdir,'catalogue_flegrei_VT.pf')
events_VT = model.load_events(catname_VT)

# Create trace list with all traces
trs_VT_synth=[]
trs_VLP_synth=[]
trs_VT_VLP_synth=[]

# Select station
for s in stations:
    st=s
    print('synt at station:',st.station)
    
    # Select event
    e_name='flegrei_2023_06_11_06_44_25'

    ev_VT=False
    for e in events_VT:
        if e.name==e_name:
            ev_VT=e
    if not ev_VT:
        print(f'Error: event {e_name} not found in {catname_VT}')

    ev_VLP=False
    for e in events_VLP:
        if e.name==e_name:
            ev_VLP=e
    if not ev_VLP:
        print(f'Error: event {e_name} not found in {catname_VLP}')

    # Synth
    store_id = 'campiflegrei_near'

    engine = LocalEngine(store_superdirs=['../GF_STORES'])

    channel_codes = ['HHE','HHN','HHZ']
    targets_VT = [
        Target(
            lat=st.lat,
            lon=st.lon,
            store_id=store_id,
            codes=('IV', st.station, '', channel_code))
        for channel_code in channel_codes]

    targets_VLP = [
        Target(
            lat=st.lat,
            lon=st.lon,
            store_id=store_id,
            codes=('IV', st.station, '', channel_code))
        for channel_code in channel_codes]

    # MT source representation.

    source_mt_VT = MTSource.from_pyrocko_event(ev_VT)

    source_mt_VLP = MTSource.from_pyrocko_event(ev_VLP)
    source_mt_VLP.stf=gf.ResonatorSTF(30., frequency=0.114)

    # return a pyrocko.gf.Reponse object.
    response_VT = engine.process(source_mt_VT, targets_VT)
    response_VLP = engine.process(source_mt_VLP, targets_VLP)

    # list of the requested traces:
    synthetic_traces_VT = response_VT.pyrocko_traces()
    synthetic_traces_VLP = response_VLP.pyrocko_traces()

    a=synthetic_traces_VT[0].tmax-synthetic_traces_VT[0].tmin
    b=synthetic_traces_VLP[0].tmax-synthetic_traces_VLP[0].tmin
    print(f'BEFORE\nlenght of VT:{a}s\nlengthof VLP: {b}')
    # chop heads and tails
    trs_VT= [ x.chop(x.tmin+1,x.tmax-1) for x in synthetic_traces_VT ]
    trs_VLP= [ x.chop(x.tmin+1,x.tmax-1) for x in synthetic_traces_VLP ]
    # add longer heads and tails
    tlen = 180.
    newtrs_VT = []
    for tr in trs_VT:
        newtr = tr.copy()
        ydata = newtr.get_ydata()
        ydata= ydata - np.mean(ydata)
        first, last = ydata[0], ydata[-1]
        npts = int(tlen/newtr.deltat)
        ydata = np.concatenate( (np.ones(npts) * first, ydata, np.ones(npts) * last) )
        newtr.ydata = ydata
        newtr.shift(-tlen)
        newtr.tmax+= 2* tlen
        newtrs_VT.append(newtr)

    tlen = 240.
    newtrs_VLP = []
    for tr in trs_VLP:
        newtr = tr.copy()
        ydata = newtr.get_ydata()
        ydata= ydata - np.mean(ydata)
        first, last = ydata[0], ydata[-1]
        npts = int(tlen/newtr.deltat)
        ydata = np.concatenate( (np.ones(npts) * first, ydata, np.ones(npts) * last) )
        newtr.ydata = ydata
        newtr.shift(-tlen)
        newtr.tmax+= 2* tlen
        newtrs_VLP.append(newtr)
    
    a=newtrs_VT[0].tmax-newtrs_VT[0].tmin
    b=newtrs_VLP[0].tmax-newtrs_VLP[0].tmin
    print(f'AFTER\nlenght of VT:{a}s\nlengthof VLP: {b}')

    # sum VLP and VT traces
    trs_sum = []
    channels=['HHE','HHN','HHZ']
    for n,ch in enumerate(channels):
        dt=newtrs_VLP[n].deltat
        tmin = newtrs_VLP[n].tmin
        tshift= int( (newtrs_VT[n].tmin - newtrs_VLP[n].tmin) / dt )
        len_tr_VT= len (newtrs_VT[n].get_ydata())
        tr1=newtrs_VLP[n].get_ydata()
        tr2=newtrs_VT[n].get_ydata()
        trsum=tr1.copy()
        trsum[tshift:tshift+len_tr_VT] += tr2
        trs_sum.append(trace.Trace(
                    network='IV', station=st.station, channel=ch,location='', deltat=dt, tmin=tmin, ydata=trsum))

    #trs.extend(trs_VT)     # non chopped
    #trs.extend(trs_VLP)    # non chopped
    trs_VT_synth.extend(newtrs_VT)   # chopped
    trs_VLP_synth.extend(newtrs_VLP)  # chopped
    trs_VT_VLP_synth.extend(trs_sum)     # sum

# snuffler
trace.snuffle(trs_VLP_synth)

# save synth traces (watch out for the 'location' parameter: max 2 letters)
io.save(trs_VT_synth, '../DATA_synth/VT_flegrei_2023_06_11_06_44_25/VT_flegrei_2023_06_11_06_44_25.mseed')
io.save(trs_VLP_synth, '../DATA_synth/VLP_flegrei_2023_06_11_06_44_25/VLP_flegrei_2023_06_11_06_44_25.mseed')
io.save(trs_VT_VLP_synth, '../DATA_synth/VT+VLP_flegrei_2023_06_11_06_44_25/VT+VLP_flegrei_2023_06_11_06_44_25.mseed')