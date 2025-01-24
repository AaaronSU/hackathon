#!/bin/bash

g++ opt1.cxx -o kk1 -Ofast -mcpu=neoverse-v2 -funroll-loops -ftree-vectorize -flto -ffp-contract=fast -fopenmp -DNDEBUG -g0 -lm -fomit-frame-pointer -I/tools/acfl/24.10/armpl-24.10.1_AmazonLinux-2_gcc/include -L/tools/acfl/24.10/armpl-24.10.1_AmazonLinux-2_gcc/lib -ltrng4 -lgsl -lgslcblas -lm -larmpl


export PATH=/tools/acfl/24.10/arm-linux-compiler-24.10.1_AmazonLinux-2/bin:$PATH
armclang++ -mcpu=neoverse-512tvb -O3 -fopenmp -funroll-loops -fvectorize -ffinite-math-only -funsafe-math-optimizations -fno-math-errno -finline-functions -armpl -lamath -lm -g -fno-omit-frame-pointer opt1.cxx -o kk2  -I/usr/include -L/usr/lib64 -lgsl

# ../maqao.aarch64.2.21.1/bin/
maqao OV -R1 -c="g++.json" -xp=report_g++ --replace --executable=./kk1
# ../maqao.aarch64.2.21.1/bin/
maqao OV -R1 -c="arm.json" -xp=report_arm --replace --executable=./kk2