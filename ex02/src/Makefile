SHELL=/bin/sh

include ./make.inc

.PHONY: clean cleanall

all: waveguide

HDA_INC = ./
NHEADERS =
HEADERS = $(addprefix $(HDA_INC), $(NHEADERS))
OBJECTS = kernels.o utils.o main.o

#asmcode:
#	$(RZ_CXX) $(RZ_COMP) $(RZ_CXXFLAGS) -c relax_line.cc -S -fsource-asm -o relax_line.s

waveguide: $(HEADERS) $(OBJECTS)
	$(RZ_CXX) $(RZ_COMP) $(RZ_LDFLAGS) $(OBJECTS) -o $@

%.o : %.c
	$(RZ_CC) $(RZ_CFLAGS) -c $<

%.o : %.cc
	$(RZ_CXX) $(RZ_COMP) $(RZ_CXXFLAGS) -I$(HDA_INC) -c $<

clean:
	rm -f *.o
cleanall: clean
	rm -f waveguide
