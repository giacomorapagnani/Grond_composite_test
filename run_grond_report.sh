#!/bin/bash

events=(
  "flegrei_2024_06_08_01_52_04"
  "flegrei_2024_04_27_03_44_56"
  "flegrei_2024_05_22_06_28_00"
  "flegrei_2025_02_16_14_30_02"
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