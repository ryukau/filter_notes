#!/bin/bash

mkdir -p bin
mkdir -p snd

g++ -std=c++17 -lsndfile -O3 -o bin/test test.cpp
./bin/test

python3 plot.py
