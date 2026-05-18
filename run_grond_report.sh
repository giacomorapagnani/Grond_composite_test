#!/bin/bash

events=(
  "flegrei_2018_09_18_21_36_41"
  "flegrei_2023_06_11_06_44_25"
  "flegrei_2023_09_07_17_45_28"
  "flegrei_2023_09_26_07_10_29"
  "flegrei_2023_10_02_20_08_26"
  "flegrei_2024_06_18_01_58_24"
  "flegrei_2024_07_26_11_46_21"
  "flegrei_2024_08_30_19_23_15"
  "flegrei_2025_03_13_00_25_02"
  "flegrei_2025_03_14_18_44_10"
)

echo "Composite STD parameters LF reports" >> errors.log
echo "----------------------------" >> errors.log

for ev in "${events[@]}"; do
  gronf_name=$(echo "$ev" | sed 's/flegrei_/flegrei./; s/_/./g; s/$/.composite.LF.std.gronf/')
  screen -dmS "$ev" bash -c "
    source ../venvs/grond_composite_0/bin/activate
    if [ \$? -ne 0 ]; then echo '[ERROR] $ev: attivazione venv fallita' >> errors.log; exit 1; fi

    grond report runs/cmt_composite_LF_std_${ev}.grun/
    if [ \$? -ne 0 ]; then echo '[ERROR] $ev: grond report fallito' >> errors.log; exit 1; fi

    echo '[OK] $ev: completato' >> errors.log
  "
  echo "Avviato: $ev"
done