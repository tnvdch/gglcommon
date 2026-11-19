# ===================================================================
# Makefile for COMMON - Programs
# 
# submakefile which is invoqued from top makefile
#
# created 09/06/98      MC - EF
#
# last changed : 18/12/98 
#
# =================================================================== 

VERBOSE=0
ifeq ($(VERBOSE),1)
V=
else
V=@
endif

OBJDIR=../obj
LDFLAGS= -L $(OBJDIR) -lcommon${EXT} -lm
LIBVER='"'[ GGL Common Library $(LIBVERSION) ]'"'
DSV=-DSRC_VERSION='"$(SV)"'

default : 
	@echo "*** Please use top makefile"

.PHONY : clean pu

include Pgnlist

RPC_FILES=$(foreach file,$(PGNS),$(REPOSITORY)/apps/$(file).c,v )

% : %.c ${OBJDIR}/libcommon.a
	$(V)$(CC) $(DSV) $(CFLAGS) -I../include $< -o $@ $(LDFLAGS)

clean :
	$(V)rm -f $(PGNS)

pu:
	$(V)rm -f *~ *% *.bak
