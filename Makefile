CXX = g++
CXXFLAGS = -O3 -g
LDFLAGS = -lm
OPENMP_FLAGS = -fopenmp

# Lib en dur 
ARMPL_INCLUDE = /tools/acfl/24.10/armpl-24.10.1_AmazonLinux-2_gcc/include
ARMPL_LIB = /tools/acfl/24.10/armpl-24.10.1_AmazonLinux-2_gcc/lib
GSL_INCLUDE = /usr/include
GSL_LIB = /usr/lib64
SLEEF_INCLUDE = /usr/local/include
SLEEF_LIB = /usr/local/lib64

TARGETS = BSM_armpl BSM_original BSM_sve BSM_gsl

SOURCES_ARMPL = BSM_armpl.cxx
SOURCES_ORIGINAL = BSM2_original.cxx
SOURCES_SVE = BSM_sve.cxx
SOURCES_GSL = BSM_gsl.cxx

all: $(TARGETS)

BSM_armpl: $(SOURCES_ARMPL)
	$(CXX) -mcpu=native $(CXXFLAGS) $(OPENMP_FLAGS) -funroll-loops -ftree-vectorize -finline-functions \
	-I$(ARMPL_INCLUDE) -L$(ARMPL_LIB) -larmpl -larmpl_mp -lamath $(LDFLAGS) $< -o $@

BSM_original: $(SOURCES_ORIGINAL)
	$(CXX) -O $< -o $@

BSM_sve: $(SOURCES_SVE)
	$(CXX) -O3 -march=native $(OPENMP_FLAGS) -I$(SLEEF_INCLUDE) -L$(SLEEF_LIB) $(LDFLAGS) $< -lsleef -o $@

BSM_gsl: $(SOURCES_GSL)
	$(CXX) $(CXXFLAGS) -I$(GSL_INCLUDE) -L$(GSL_LIB) -lgsl -lgslcblas $< -o $@

clean:
	rm -f $(TARGETS)

.PHONY: all clean
