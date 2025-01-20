#!/bin/bash

g++ -Ofast -march=native BSM.cxx -o BSM -fopenmp -I/opt/armpl/include -L/opt/armpl/lib -larmpl
# g++ -O BSM.cxx -o BSM
./BSM 100000 100
./BSM 1000000 100
./BSM 10000000 100
./BSM 100000000 100
