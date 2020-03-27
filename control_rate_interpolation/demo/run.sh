#!/bin/bash
mkdir -p bin
mkdir -p snd

g++ -std=c++17 -lsndfile -O3 -o bin/naive naive.cpp
./bin/naive

g++ -std=c++17 -lsndfile -O3 -o bin/linterp linterp.cpp
./bin/linterp

g++ -std=c++17 -lsndfile -O3 -o bin/smoother smoother.cpp
./bin/smoother

g++ -std=c++17 -lsndfile -O3 -o bin/pcontroller pcontroller.cpp
./bin/pcontroller

g++ -std=c++17 -lsndfile -O3 -o bin/rate_limiter rate_limiter.cpp
./bin/rate_limiter

echo bench_linterp
g++ -std=c++17 -lsndfile -O3 -o bin/bench_linterp bench_linterp.cpp
./bin/bench_linterp

echo bench_smoother
g++ -std=c++17 -lsndfile -O3 -o bin/bench_smoother bench_smoother.cpp
./bin/bench_smoother
