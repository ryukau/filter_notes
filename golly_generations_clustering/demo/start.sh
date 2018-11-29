#!/bin/bash
if [ ! -d "img" ]; then mkdir img; fi
if [ ! -d "data" ]; then mkdir data; fi

python3 mfcc.py $1
# python3 spectrogram.py $1
python3 plot_kmeans.py

python3 tsne.py
python3 plot_tsne.py
