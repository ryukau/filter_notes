#!/bin/sh
python3 burgers_exact.py
ffmpeg -y -framerate 60 -i img/burgers_%04d.png -i burgers_exact.wav -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -c:a libfdk_aac -b:a 192k burgers_exact.mp4
