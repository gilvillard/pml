
INC_DIR = ../../include
LIB_DIR = ../../lib

OBJS = $(wildcard $(LIB_DIR)/*)

SRC = $(wildcard time*cpp)
EXEC = $(patsubst %.cpp, %, $(SRC))
DATA = $(patsubst %.cpp, %.dat, $(SRC))
COMPILE = $(patsubst %.cpp, %.exe, $(patsubst time_%, %, $(SRC)))
RUN = $(patsubst %.cpp, %.dat, $(patsubst time_%, %, $(SRC)))

CXX = g++
CXXFLAGS = -Wall -std=c++11 -O3 -march=native -I$(INC_DIR)
LIBS = -L/dsk/l1/vneiger/lib -lntl -lpthread -lgmp -lm

LBCXXFLAGS = -DHAVE_CONFIG_H -O2 -Wall -g -DNDEBUG -U_LB_DEBUG -DDISABLE_COMMENTATOR -fopenmp -fabi-version=6 -I/usr/local/include -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -mavx2 -mfma
LBLIBS = -lntl -L/usr/local/lib -lopenblas -lpthread -lgfortran -lgivaro -lgmp -lgmpxx -fopenmp

.PHONY: clean

clean:
	rm -f $(EXEC) $(DATA)

all:	$(COMPILE)

run:	$(RUN)

ntlx_linbox: ntlx_linbox.cpp
	$(CXX) -o ntlx_linbox ntlx_linbox.cpp $(OBJS) $(CXXFLAGS) $(LBCXXFLAGS) $(LIBS) $(LBLIBS)

%.exe:
	if test -e time_$*.cpp; then $(CXX) -o time_$* time_$*.cpp $(OBJS) $(CXXFLAGS) $(LIBS); fi

%.dat: 
	if test -e time_$*; then ./time_$* > time_$*.dat; fi
