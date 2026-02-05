#!/usr/bin/env python3
"""
Script per creare una risposta strumentale PIATTA compatibile con ObsPy.
La risposta piatta ha poli e zeri che si annullano, risultando in risposta unitaria.

Autore: Generato per sismologo
"""

import xml.etree.ElementTree as ET
from copy import deepcopy

def create_flat_response_obspy_compatible(input_file, output_file):
    """
    Crea un file StationXML con risposta strumentale piatta ma COMPLETA.
    
    La risposta è configurata con:
    - Sensitivity = 1.0
    - Stage 1: un polo e uno zero che si annullano (risposta = 1)
    - Stage 2+: rimossi o con gain = 1.0
    
    Questo è compatibile con ObsPy e non modifica il segnale.
    """
    
    namespace = {'fdsn': 'http://www.fdsn.org/xml/station/1'}
    ET.register_namespace('', 'http://www.fdsn.org/xml/station/1')
    
    print(f"Lettura file: {input_file}")
    tree = ET.parse(input_file)
    root = tree.getroot()
    
    channels_modified = 0
    
    for network in root.findall('.//fdsn:Network', namespace):
        for station in network.findall('.//fdsn:Station', namespace):
            for channel in station.findall('.//fdsn:Channel', namespace):
                response = channel.find('fdsn:Response', namespace)
                
                if response is not None:
                    # 1. Modifica InstrumentSensitivity
                    instrument_sensitivity = response.find('fdsn:InstrumentSensitivity', namespace)
                    if instrument_sensitivity is not None:
                        value_elem = instrument_sensitivity.find('fdsn:Value', namespace)
                        if value_elem is not None:
                            value_elem.text = "1.0"
                    
                    # 2. Trova tutti gli stage esistenti
                    stages = response.findall('fdsn:Stage', namespace)
                    
                    # 3. Rimuovi tutti gli stage esistenti
                    for stage in stages:
                        response.remove(stage)
                    
                    # 4. Crea UN SOLO stage con risposta piatta
                    # Questo evita problemi con decimation senza filtri
                    stage1 = ET.SubElement(response, '{http://www.fdsn.org/xml/station/1}Stage')
                    stage1.set('number', '1')
                    
                    # Crea PolesZeros con un polo e uno zero che si annullano
                    # Risultato: H(s) = (s + a) / (s + a) = 1
                    pz = ET.SubElement(stage1, '{http://www.fdsn.org/xml/station/1}PolesZeros')
                    
                    input_units = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}InputUnits')
                    n_input = ET.SubElement(input_units, '{http://www.fdsn.org/xml/station/1}Name')
                    n_input.text = "M/S"
                    desc_input = ET.SubElement(input_units, '{http://www.fdsn.org/xml/station/1}Description')
                    desc_input.text = "Velocity in Meters per Second"
                    
                    output_units = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}OutputUnits')
                    n_output = ET.SubElement(output_units, '{http://www.fdsn.org/xml/station/1}Name')
                    n_output.text = "COUNTS"
                    desc_output = ET.SubElement(output_units, '{http://www.fdsn.org/xml/station/1}Description')
                    desc_output.text = "Digital Counts"
                    
                    pz_type = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}PzTransferFunctionType')
                    pz_type.text = "LAPLACE (RADIANS/SECOND)"
                    
                    norm_factor = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}NormalizationFactor')
                    norm_factor.text = "1.0"
                    
                    norm_freq = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}NormalizationFrequency')
                    norm_freq.set('unit', 'HERTZ')
                    norm_freq.text = "1.0"
                    
                    # Zero at s = -10
                    zero = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}Zero')
                    zero.set('number', '0')
                    real_z = ET.SubElement(zero, '{http://www.fdsn.org/xml/station/1}Real')
                    real_z.text = "-10.0"
                    imag_z = ET.SubElement(zero, '{http://www.fdsn.org/xml/station/1}Imaginary')
                    imag_z.text = "0.0"
                    
                    # Pole at s = -10 (si annulla con lo zero!)
                    pole = ET.SubElement(pz, '{http://www.fdsn.org/xml/station/1}Pole')
                    pole.set('number', '0')
                    real_p = ET.SubElement(pole, '{http://www.fdsn.org/xml/station/1}Real')
                    real_p.text = "-10.0"
                    imag_p = ET.SubElement(pole, '{http://www.fdsn.org/xml/station/1}Imaginary')
                    imag_p.text = "0.0"
                    
                    # StageGain = 1.0
                    stage_gain = ET.SubElement(stage1, '{http://www.fdsn.org/xml/station/1}StageGain')
                    gain_value = ET.SubElement(stage_gain, '{http://www.fdsn.org/xml/station/1}Value')
                    gain_value.text = "1.0"
                    gain_freq = ET.SubElement(stage_gain, '{http://www.fdsn.org/xml/station/1}Frequency')
                    gain_freq.text = "1.0"
                    
                    channels_modified += 1
    
    print(f"Scrittura file: {output_file}")
    tree.write(output_file, encoding='UTF-8', xml_declaration=True)
    
    print(f"\n✓ Completato!")
    print(f"  - Canali modificati: {channels_modified}")
    print(f"  - Sensibilità: 1.0")
    print(f"  - Configurazione: 1 polo + 1 zero che si annullano")
    print(f"  - Risultato: Risposta H(f) = 1.0 (piatta)")
    print(f"  - File salvato: {output_file}")


if __name__ == "__main__":
    input_file = "/mnt/user-data/uploads/1770043123243_stations_flegrei_INGV_final.xml"
    output_file = "/mnt/user-data/outputs/stations_flegrei_FLAT_RESPONSE_v2.xml"
    
    create_flat_response_obspy_compatible(input_file, output_file)
    
    print("\n" + "="*70)
    print("File generato e compatibile con ObsPy!")
    print("La risposta è piatta: deconvolvendo non modificherai il segnale.")
    print("="*70)
