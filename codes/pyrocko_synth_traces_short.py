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
    print('\nSynthetic traces for station:',st.station)
    
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
    store_id = 'campiflegrei_near_0_dist'               ###CHANGE###

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
    synth_tr_VT = response_VT.pyrocko_traces()
    synth_tr_VLP = response_VLP.pyrocko_traces()

    # calculate lenght
    a=synth_tr_VT[0].tmax-synth_tr_VT[0].tmin
    b=synth_tr_VLP[0].tmax-synth_tr_VLP[0].tmin
    #print(f'lenght of VT:{a}s\nlength of VLP: {b}')
    
    # cut first and last element
    #trs_VT= [ x.chop(x.tmin+1,x.tmax-1) for x in synthetic_traces_VT ]
    #trs_VLP= [ x.chop(x.tmin+1,x.tmax-1) for x in synthetic_traces_VLP ]
    # add longer heads and tails
    
    # sum VLP and VT traces
    trs_sum = []
    channels=['HHE','HHN','HHZ']
    for n,ch in enumerate(channels):
        # dt
        dt=synth_tr_VLP[n].deltat
        # tmin
        tmin_vlp = synth_tr_VLP[n].tmin
        tmin_vt  = synth_tr_VT[n].tmin
        if tmin_vlp <= tmin_vt:
            tmin = tmin_vlp
            # shift between tmin VLP & VT
            tshift= int( round( (tmin_vt - tmin_vlp) / dt ) )

            tr_vlp=synth_tr_VLP[n].get_ydata()
            tr_vt=synth_tr_VT[n].get_ydata()
            trsum=tr_vlp.copy()

            len_tr_VT= len (tr_vt)
            trsum[tshift:tshift+len_tr_VT] += tr_vt
        else: 
            tmin=tmin_vt
        
            tshift= int( round( ( tmin_vlp - tmin_vt) / dt ) )

            tr_vlp=synth_tr_VLP[n].get_ydata()
            tr_vt=synth_tr_VT[n].get_ydata()

            # trace's first part
            trsum=tr_vt[:tshift].copy()
            # traces' second part (shared)
            lenght_shared_tr=len(tr_vt)-tshift
            tmp=tr_vt[tshift:].copy() + tr_vlp[0:lenght_shared_tr].copy()
            trsum=np.concatenate( (trsum,tmp) )
            # trace's third part
            tmp=tr_vlp[lenght_shared_tr:].copy()
            trsum=np.concatenate( (trsum,tmp) )
            #print(f'trs_sum len: {len(trsum)}, tr_vlp len: {len(tr_vlp)}, tshift: {tshift}')

        trs_sum.append(trace.Trace(
                    network='IV', station=st.station, channel=ch,location='', deltat=dt, tmin=tmin, ydata=trsum))

        
    trs_VT_synth.extend(synth_tr_VT)   
    trs_VLP_synth.extend(synth_tr_VLP) 
    trs_VT_VLP_synth.extend(trs_sum)     # sum

# snuffler
#trace.snuffle(trs_VLP_synth)

# save synth traces (watch out for the 'location' parameter: max 2 letters)
io.save(trs_VT_synth, '../DATA/VT_flegrei_2023_06_11_06_44_25/VT_short_flegrei_2023_06_11_06_44_25.mseed')
io.save(trs_VLP_synth, '../DATA/VLP_flegrei_2023_06_11_06_44_25/VLP_short_flegrei_2023_06_11_06_44_25.mseed')
io.save(trs_VT_VLP_synth, '../DATA/VT+VLP_flegrei_2023_06_11_06_44_25/VT+VLP_short_flegrei_2023_06_11_06_44_25.mseed')