#!/bin/bash

mkdir -p bin
mkdir -p snd

g++ -std=c++17 -lsndfile -O3 -o bin/expADSR test_exponential_envelope.cpp
./bin/expADSR

g++ -std=c++17 -lsndfile -O3 -o bin/P_expADSR test_P_exponential_envelope.cpp
./bin/P_expADSR

g++ -std=c++17 -lsndfile -O3 -o bin/expAD expAD.cpp
./bin/expAD

python3 plot.py
