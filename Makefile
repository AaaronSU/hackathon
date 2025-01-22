CXX = g++
CXXFLAGS = -O3 -g
LDFLAGS = -lm
OPENMP_FLAGS = -fopenmp

# Lib en dur 
ARMPL_INCLUDE = /tools/acfl/24.10/armpl-24.10.1_AmazonLinux-2_gcc/include
ARMPL_LIB = /tools/acfl/24.10/armpl-24.10.1_AmazonLinux-2_gcc/lib
GSL_INCLUDE = /path/to/gsl/include
GSL_LIB = /path/to/gsl/lib

TARGETS = BSM_armpl BSM_original BSM_sve BSM_gsl

SOURCES_ARMPL = BSM_armpl.cxx
SOURCES_ORIGINAL = BSM_original.cxx
SOURCES_SVE = BSM_sve.cxx
SOURCES_GSL = BSM_gsl.cxx

all: $(TARGETS)

BSM_armpl: $(SOURCES_ARMPL)
	$(CXX) -mcpu=native $(CXXFLAGS) -fopenmp -funroll-loops -ftree-vectorize -finline-functions \
	-I$(ARMPL_INCLUDE) -L$(ARMPL_LIB) -larmpl -larmpl_mp -lamath $(LDFLAGS) $< -o $@

BSM_original: $(SOURCES_ORIGINAL)
	$(CXX) -O BSM_original.cxx -o BSM_original

BSM_sve: $(SOURCES_SVE)
	g++ -O3 -fopenmp -march=native  BSM_sve.cxx -o BSM -fopenmp -lsleef -lm

BSM_gsl: $(SOURCES_GSL)
	$(CXX) $(CXXFLAGS) -I$(GSL_INCLUDE) -L$(GSL_LIB) -lgsl -lgslcblas $(SOURCES_GSL) -o BSM_gsl

clean:
	rm -f $(TARGETS)

.PHONY: all clean
