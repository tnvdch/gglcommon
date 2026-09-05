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
	$(CC) $(DSV) $(CFLAGS) -I../include $< -o $@ $(LDFLAGS)

coro :
	-co -r $(RPC_FILES) $(REPOSITORY)/apps/Makefile,v $(REPOSITORY)/apps/Pgnlist,v

co :
	-co -l $(RPC_FILES) $(REPOSITORY)/apps/Makefile,v $(REPOSITORY)/apps/Pgnlist,v

ci :
	-ci -u $(RPC_FILES) $(REPOSITORY)/apps/Makefile,v $(REPOSITORY)/apps/Pgnlist,v

clean :
	rm -f $(PGNS)

pu:
	rm -f *~ *% *.bak
