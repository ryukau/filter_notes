#!/bin/sh
python3 render_movie.py
ffmpeg -y -framerate 60 -i img/out%04d.png -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p heat_wave1d.mp4
