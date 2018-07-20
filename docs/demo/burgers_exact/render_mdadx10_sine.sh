#!/bin/sh
python3 mdadx10_sine.py
ffmpeg -y -framerate 60 -i img/mdadx10_sine_%04d.png -i mdadx10_sine.wav -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -c:a libfdk_aac -b:a 192k mdadx10_sine.mp4
