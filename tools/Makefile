# ===================================================================
# Makefile for COMMON - Programs
# 
# created 01/02/1999    EF
#
# last changed : 01/02/1999
#
# =================================================================== 

PLATFORM=$(HOSTTYPE)
# C++ compilation
#CC=g++
#EXT=++

# C compilation
CC=gcc
EXT=

CFLAGS=-O -g -Wall -D_GNU_SOURCE

apps :
	$(MAKE) -C . SV="`date` `hostname -s` `whoami`" $(PGNS)

OBJDIR=../obj_$(PLATFORM)
LDFLAGS= -L $(OBJDIR) -lcommon
DSV=-DSRC_VERSION='"$(SV)"'

default : 
	@echo "*** Please use top makefile"

.PHONY : clean pu

include Pgnlist

fft1d : fft1d.c ../obj_$(PLATFORM)/libcommon.a
	$(CC) $(DSV) $(CFLAGS) $(XFLAG) -I../include $< -o $@ $(LDFLAGS) -lfftw -lm

% : %.c ../obj_$(PLATFORM)/libcommon.a
	$(CC) $(DSV) $(CFLAGS) $(XFLAG) -I../include $< -o $@ $(LDFLAGS) -lm

clean :
	rm -f $(PGNS)

pu:
	rm -f *~ *% *.bak
