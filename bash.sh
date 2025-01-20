#!/bin/bash

g++ -O3 -march=armv8.4-a+simd -ffast-math BSM.cxx -o BSM
# g++ -O BSM.cxx -o BSM
./BSM 100000 100
./BSM 1000000 100
./BSM 10000000 100
./BSM 100000000 100