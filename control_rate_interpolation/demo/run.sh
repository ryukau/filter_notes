#!/bin/bash
mkdir -p bin
mkdir -p snd

g++ -lsndfile -O3 -o bin/naive naive.cpp
./bin/naive

g++ -lsndfile -O3 -o bin/linterp linterp.cpp
./bin/linterp

g++ -lsndfile -O3 -o bin/smoother smoother.cpp
./bin/smoother

echo bench_linterp
g++ -lsndfile -O3 -o bin/bench_linterp bench_linterp.cpp
./bin/bench_linterp

echo bench_smoother
g++ -lsndfile -O3 -o bin/bench_smoother bench_smoother.cpp
./bin/bench_smoother
