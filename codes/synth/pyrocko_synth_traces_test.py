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
catdir='../../CAT'
metadatadir='../../META_DATA'

stations_name=os.path.join(metadatadir,'stations_flegrei_INGV_final.pf')               ###CHANGE###
stations=model.load_stations(stations_name)

catname_VLP=os.path.join(catdir,'catalogue_flegrei_VLP_test.pf')               ###CHANGE###
events_VLP = model.load_events(catname_VLP)


# Create trace list with all traces
trs_VLP_synth=[]

# Select station
for s in stations:
    st=s
    print('\nSynthetic traces for station:',st.station)
    
    # Select event
    e_name='test_2023_06_11_06_44_25'

    ev_VLP=False
    for e in events_VLP:
        if e.name==e_name:
            ev_VLP=e
    if not ev_VLP:
        print(f'Error: event {e_name} not found in {catname_VLP}')

    # Synth
    store_id = 'campiflegrei_test_10hz'          ###CHANGE###

    engine = LocalEngine(store_superdirs=['../../GF_STORES'])

    channel_codes = ['HHE','HHN','HHZ']
    targets_VLP = [
        Target(
            lat=st.lat,
            lon=st.lon,
            store_id=store_id,
            codes=('IV', st.station, '', channel_code))
        for channel_code in channel_codes]

    # MT source representation.


    source_mt_VLP = MTSource.from_pyrocko_event(ev_VLP)
    source_mt_VLP.stf=gf.ResonatorSTF(25., frequency=0.114)               ###CHANGE###

    # return a pyrocko.gf.Reponse object.
    response_VLP = engine.process(source_mt_VLP, targets_VLP)

    # list of the requested traces:
    synth_tr_VLP = response_VLP.pyrocko_traces()

    # calculate lenght
    #b=synth_tr_VLP[0].tmax-synth_tr_VLP[0].tmin
    #print(f'lenght of VT:{a}s\nlength of VLP: {b}')
    
    # cut first and last element
    #trs_VT= [ x.chop(x.tmin+1,x.tmax-1) for x in synthetic_traces_VT ]
    #trs_VLP= [ x.chop(x.tmin+1,x.tmax-1) for x in synthetic_traces_VLP ]
    # add longer heads and tails
    
    trs_VLP_synth.extend(synth_tr_VLP) 

# snuffler
#trace.snuffle(trs_VLP_synth)

# save synth traces (watch out for the 'location' parameter: max 2 letters)
io.save(trs_VLP_synth, '../../DATA/SYNTH/VLP_synth_10hz.mseed')