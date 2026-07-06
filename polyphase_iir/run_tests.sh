#!/bin/bash

# Configuration format: "order,stages,frequency"
configs=(
    "4,2,0.125"
    "8,4,0.03125"
    "12,4,0.02"
    "16,8,0.015625"
    "16,32,0.0078125"
)

compilers=("g++")

for cc in "${compilers[@]}"; do
  for cfg in "${configs[@]}"; do
    IFS=',' read -r order stages freq <<< "$cfg"
    id="${order}_${stages}"

    python run_tests.py \
      --ref "data/reference_${id}.json" \
      --results "data/results_${id}_${cc}.json" \
      --compiler "$cc" \
      --samples 48000 \
      --benchmark \
      --design "butterworth,${order},${freq},${stages}"
  done
done
