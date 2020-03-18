#!/bin/bash

mkdir -p bin
mkdir -p snd

g++ -std=c++17 -lsndfile -O3 -o bin/test test.cpp
./bin/test

g++ -std=c++17 -lsndfile -O3 -o bin/expAD expAD.cpp
./bin/expAD

python3 plot.py
