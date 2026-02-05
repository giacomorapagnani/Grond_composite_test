#!/usr/bin/env python3
"""
Script per creare un file StationXML con risposta piatta che matcha
ESATTAMENTE le tracce presenti in un file mseed.

Questo risolve il problema "No matching response information found"
quando si lavora con dati sintetici.

Uso:
    python create_flat_xml_from_mseed.py input.mseed output.xml
"""

import sys
from obspy import read
import xml.etree.ElementTree as ET
from datetime import datetime

def create_flat_xml_from_mseed(mseed_file, output_xml):
    """
    Legge un file mseed e crea un XML con risposta piatta per tutte le tracce.
    """
    
    print(f"Lettura file mseed: {mseed_file}")
    st = read(mseed_file)
    
    print(f"Trovate {len(st)} tracce")
    
    # Crea root element
    root = ET.Element('{http://www.fdsn.org/xml/station/1}FDSNStationXML')
    root.set('schemaVersion', '1.2')
    
    # Header
    source = ET.SubElement(root, '{http://www.fdsn.org/xml/station/1}Source')
    source.text = "Generated for flat response"
    
    sender = ET.SubElement(root, '{http://www.fdsn.org/xml/station/1}Sender')
    sender.text = "Python Script"
    
    created = ET.SubElement(root, '{http://www.fdsn.org/xml/station/1}Created')
    created.text = datetime.utcnow().isoformat() + "Z"
    
    # Raggruppa per network
    networks = {}
    for tr in st:
        net_code = tr.stats.network
        sta_code = tr.stats.station
        loc_code = tr.stats.location
        cha_code = tr.stats.channel
        
        if net_code not in networks:
            networks[net_code] = {}
        if sta_code not in networks[net_code]:
            networks[net_code][sta_code] = {}
        if loc_code not in networks[net_code][sta_code]:
            networks[net_code][sta_code][loc_code] = []
        
        # Aggiungi canale se non già presente
        if cha_code not in networks[net_code][sta_code][loc_code]:
            networks[net_code][sta_code][loc_code].append(cha_code)
    
    # Crea XML structure
    for net_code, stations in networks.items():
        network = ET.SubElement(root, '{http://www.fdsn.org/xml/station/1}Network')
        network.set('code', net_code)
        network.set('startDate', '1900-01-01T00:00:00.000000Z')
        network.set('restrictedStatus', 'open')
        
        desc = ET.SubElement(network, '{http://www.fdsn.org/xml/station/1}Description')
        desc.text = f"Network {net_code}"
        
        for sta_code, locations in stations.items():
            station = ET.SubElement(network, '{http://www.fdsn.org/xml/station/1}Station')
            station.set('code', sta_code)
            station.set('startDate', '1900-01-01T00:00:00.000000Z')
            station.set('restrictedStatus', 'open')
            
            # Coordinate fittizie (puoi modificarle)
            lat = ET.SubElement(station, '{http://www.fdsn.org/xml/station/1}Latitude')
            lat.set('unit', 'DEGREES')
            lat.text = '40.8'
            
            lon = ET.SubElement(station, '{http://www.fdsn.org/xml/station/1}Longitude')
            lon.set('unit', 'DEGREES')
            lon.text = '14.1'
            
            elev = ET.SubElement(station, '{http://www.fdsn.org/xml/station/1}Elevation')
            elev.text = '0.0'
            
            site = ET.SubElement(station, '{http://www.fdsn.org/xml/station/1}Site')
            site_name = ET.SubElement(site, '{http://www.fdsn.org/xml/station/1}Name')
            site_name.text = f"Station {sta_code}"
            
            for loc_code, channels in locations.items():
                for cha_code in channels:
                    channel = ET.SubElement(station, '{http://www.fdsn.org/xml/station/1}Channel')
                    channel.set('code', cha_code)
                    channel.set('locationCode', loc_code)
                    channel.set('startDate', '1900-01-01T00:00:00.000000Z')
                    channel.set('endDate', '2099-12-31T23:59:59.000000Z')
                    channel.set('restrictedStatus', 'open')
                    
                    ch_lat = ET.SubElement(channel, '{http://www.fdsn.org/xml/station/1}Latitude')
                    ch_lat.set('unit', 'DEGREES')
                    ch_lat.text = '40.8'
                    
                    ch_lon = ET.SubElement(channel, '{http://www.fdsn.org/xml/station/1}Longitude')
                    ch_lon.set('unit', 'DEGREES')
                    ch_lon.text = '14.1'
                    
                    ch_elev = ET.SubElement(channel, '{http://www.fdsn.org/xml/station/1}Elevation')
                    ch_elev.text = '0.0'
                    
                    ch_depth = ET.SubElement(channel, '{http://www.fdsn.org/xml/station/1}Depth')
                    ch_depth.text = '0.0'
                    
                    # Azimuth/Dip basati sul channel code
                    azimuth = ET.SubElement(channel, '{http://www.fdsn.org/xml/station/1}Azimuth')
                    azimuth.set('unit', 'DEGREES')
                    if cha_code.endswith('E'):
                        azimuth.text = '90.0'
                    elif cha_code.endswith('N'):
                        azimuth.text = '0.0'
                    else:  # Z
                        azimuth.text = '0.0'
                    
                    dip = ET.SubElement(channel, '{http://www.fdsn.org/xml/station/1}Dip')
                    dip.set('unit', 'DEGREES')
                    if cha_code.endswith('Z'):
                        dip.text = '-90.0'
                    else:
                        dip.text = '0.0'
                    
                    # Sample rate
                    sample_rate = ET.SubElement(channel, '{http://www.fdsn.org/xml/station/1}SampleRate')
                    sample_rate.set('unit', 'SAMPLES/S')
                    sample_rate.text = '100.0'
                    
                    # RISPOSTA PIATTA
                    response = ET.SubElement(channel, '{http://www.fdsn.org/xml/station/1}Response')
                    
                    # InstrumentSensitivity
                    inst_sens = ET.SubElement(response, '{http://www.fdsn.org/xml/station/1}InstrumentSensitivity')
                    
                    sens_value = ET.SubElement(inst_sens, '{http://www.fdsn.org/xml/station/1}Value')
                    sens_value.text = '1.0'
                    
                    sens_freq = ET.SubElement(inst_sens, '{http://www.fdsn.org/xml/station/1}Frequency')
                    sens_freq.text = '1.0'
                    
                    sens_input_units = ET.SubElement(inst_sens, '{http://www.fdsn.org/xml/station/1}InputUnits')
                    sens_input_name = ET.SubElement(sens_input_units, '{http://www.fdsn.org/xml/station/1}Name')
                    sens_input_name.text = 'M/S'
                    
                    sens_output_units = ET.SubElement(inst_sens, '{http://www.fdsn.org/xml/station/1}OutputUnits')
                    sens_output_name = ET.SubElement(sens_output_units, '{http://www.fdsn.org/xml/station/1}Name')
                    sens_output_name.text = 'COUNTS'
                    
                    # Stage 1 con polo e zero che si annullano
                    stage = ET.SubElement(response, '{http://www.fdsn.org/xml/station/1}Stage')
                    stage.set('number', '1')
                    
                    pz = ET.SubElement(stage, '{http://www.fdsn.org/xml/station/1}PolesZeros')
                    
                    pz_input_units = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}InputUnits')
                    pz_input_name = ET.SubElement(pz_input_units, '{http://www.fdsn.org/xml/station/1}Name')
                    pz_input_name.text = 'M/S'
                    
                    pz_output_units = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}OutputUnits')
                    pz_output_name = ET.SubElement(pz_output_units, '{http://www.fdsn.org/xml/station/1}Name')
                    pz_output_name.text = 'COUNTS'
                    
                    pz_type = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}PzTransferFunctionType')
                    pz_type.text = 'LAPLACE (RADIANS/SECOND)'
                    
                    norm_factor = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}NormalizationFactor')
                    norm_factor.text = '1.0'
                    
                    norm_freq = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}NormalizationFrequency')
                    norm_freq.set('unit', 'HERTZ')
                    norm_freq.text = '1.0'
                    
                    # Zero
                    zero = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}Zero')
                    zero.set('number', '0')
                    zero_real = ET.SubElement(zero, '{http://www.fdsn.org/xml/station/1}Real')
                    zero_real.text = '-10.0'
                    zero_imag = ET.SubElement(zero, '{http://www.fdsn.org/xml/station/1}Imaginary')
                    zero_imag.text = '0.0'
                    
                    # Pole (stesso valore dello zero!)
                    pole = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}Pole')
                    pole.set('number', '0')
                    pole_real = ET.SubElement(pole, '{http://www.fdsn.org/xml/station/1}Real')
                    pole_real.text = '-10.0'
                    pole_imag = ET.SubElement(pole, '{http://www.fdsn.org/xml/station/1}Imaginary')
                    pole_imag.text = '0.0'
                    
                    # StageGain
                    stage_gain = ET.SubElement(stage, '{http://www.fdsn.org/xml/station/1}StageGain')
                    gain_value = ET.SubElement(stage_gain, '{http://www.fdsn.org/xml/station/1}Value')
                    gain_value.text = '1.0'
                    gain_freq = ET.SubElement(stage_gain, '{http://www.fdsn.org/xml/station/1}Frequency')
                    gain_freq.text = '1.0'
                    
                    print(f"  Creato: {net_code}.{sta_code}.{loc_code}.{cha_code}")
    
    # Registra namespace
    ET.register_namespace('', 'http://www.fdsn.org/xml/station/1')
    
    # Salva
    tree = ET.ElementTree(root)
    tree.write(output_xml, encoding='UTF-8', xml_declaration=True)
    
    print(f"\n✓ File salvato: {output_xml}")
    print(f"  Networks: {len(networks)}")
    total_channels = sum(len(ch) for sta in networks.values() 
                        for loc in sta.values() for ch in loc)
    print(f"  Canali totali: {total_channels}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python create_flat_xml_from_mseed.py <input.mseed> <output.xml>")
        sys.exit(1)
    
    mseed_file = sys.argv[1]
    output_xml = sys.argv[2]
    
    create_flat_xml_from_mseed(mseed_file, output_xml)
