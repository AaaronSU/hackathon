#!/bin/bash


g++ -O BSM2.cxx -o BSM2
./BSM2 100000 10
./BSM2 1000000 10
./BSM2 10000000 10
./BSM2 100000000 10

g++ -Ofast -march=native BSM2p.cxx -o BSM2p -fopenmp 
./BSM2p 100000 10
./BSM2p 1000000 10
./BSM2p 10000000 10
./BSM2p 100000000 10

g++ -Ofast -march=native BSM.cxx -o BSM -fopenmp 
# -I/opt/armpl/include -L/opt/armpl/lib -larmpl

./BSM 100000 10
./BSM 1000000 10
./BSM 10000000 10
./BSM 100000000 10