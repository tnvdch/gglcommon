# ===================================================================
# TOP Makefile for COMMON LIB - Programs
# 
# submakefile which is invoqued from top makefile
#
# created 09/06/98      MC - EF
#
# see Log at end of file
#
# =================================================================== 

VERBOSE=0
ifeq ($(VERBOSE),1)
V=
else
V=@
endif


include ./VERSION

COMPILINFO=compiled `date` by `whoami` at `uname -n`

# What compiler to use and which flags to set.
# extension used in filenames 
# 
# C++ compilation
#CC=g++
#EXT=++

# C compilation
CC=gcc
EXT=

# To enable Debugging Checks in matcalc, define DEVEL_MAT_DEBUG !
CFLAGS=-O -g -Wall -DLIBCOMP="\"[ $(COMPILINFO) ]"\" -DLIBVER='"[ GGL Common Library '$(LIBVERSION)' ]"' -D_GNU_SOURCE -DDEVEL_MAT_DEBUG -fPIC

ifeq ($(M32),1)
CFLAGS += -m32
endif

#
# where the RCS Repository is
#

REPOSITORY=/mnt/molnya00/linux/gglutils/src/RCS/common
INST_DIR=/mnt/molnya00/linux/gglutils

#
# Where to install things 
#
INST_LIB_DIR=${INST_DIR}/lib_${PLATFORM}
INST_INC_DIR=${INST_DIR}/include
INST_PGN_DIR=${INST_DIR}/bin_${PLATFORM}

export

default :
	@echo "use with :"
	@echo "  dep       make dependencies"
	@echo "  lib       make library"
	@echo "  instlib   install library to: $(INST_LIB_DIR)"
	@echo "  apps      make applications"
	@echo "  instapps  install applications to: $(INST_PGN_DIR)"
	@echo "  install   install library and applications"
	@echo "  cleandep  remove dependency file"
	@echo "  cleanlib  clean library"
	@echo "  cleanapps clean applications"
	@echo "  cleandoc  clean documentation"
	@echo "  clean     clean library, applications and documentation"
	@echo "  doc       create documentation from sgml"
	@echo "  coro      check out (read only)"
	@echo "  co        check out (with lock)"
	@echo "  ci        check in (of a previous co)"

.PHONY:  cleanlib cleanapps cleandep clean pu instlib instapps install

#------------------------------  check-in check-out

co_lib : 
	$(V)touch libsrc/.depend
	$(V)$(MAKE) -C libsrc co

coro_lib :
	$(V)touch libsrc/.depend
	$(V)$(MAKE) -C libsrc coro

ci_lib :
	$(V)touch libsrc/.depend
	$(V)$(MAKE) -C libsrc ci

co_doc : 
	$(V)$(MAKE) -C Documentation/sgml co

coro_doc :
	$(V)$(MAKE) -C Documentation/sgml coro

ci_doc :
	$(V)$(MAKE) -C Documentation/sgml ci

co_apps :
	$(V)$(MAKE) -C appssrc co

coro_apps :
	$(V)$(MAKE) -C appssrc coro

ci_apps:
	$(V)$(MAKE) -C appssrc ci

#------------------------------  make liberary

include libsrc/Objlist

dep :
	$(V)$(MAKE) -C libsrc dep

lib : dep
	$(V)$(MAKE) -C libsrc lib

instlib : dep lib
	$(V)$(MAKE) -C libsrc install


#------------------------------  make apps

include appssrc/Pgnlist

% :
	$(V)$(MAKE) -C appssrc SV="`date` `hostname -s` `whoami`" $@

apps : lib
	$(V)$(MAKE) -C appssrc SV="`date` `hostname -s` `whoami`" $(PGNS)

instapps :  apps 
	$(V)if [ ! -d $(INST_PGN_DIR) ] ; then mkdir $(INST_PGN_DIR) ; fi
	$(V)$(foreach pgn, $(PGNS),install -m 0771 appssrc/$(pgn) $(INST_PGN_DIR)/$(pgn) && ) echo ""

install : instlib instapps

#------------------------------  create documentation

doc :
	$(V)if [ ! -d Documentation/bin ] ; then mkdir Documentation/bin ; fi
	$(V)cd Documentation/bin; tar xf $(REPOSITORY)/Documentation/bin/pgn.tar; makedoc

#------------------------------  cleaning

cleandep :
	$(V)$(MAKE) -C libsrc cleandep VERBOSE=$(VERBOSE)

cleanlib :
	$(V)$(MAKE) -C libsrc clean VERBOSE=$(VERBOSE)

cleanapps :
	$(V)$(MAKE) -C appssrc clean VERBOSE=$(VERBOSE)

cleandoc :
	$(V)rm -f Documentation/{latex,html}/*
	$(V)rm -f Documentation/sgml/DOC*.sgml

clean : cleanlib cleanapps cleandoc

pu : 
	$(V)rm -f *~ *.bak *%
	$(V)$(MAKE) -C appssrc pu
	$(V)$(MAKE) -C libsrc pu

#
# $Log: Makefile,v $
# Revision 1.98  2001/08/31 12:09:06  favey
# cygwin adaptions
#
# Revision 1.97  2001/02/08 09:36:43  cocard
# *** empty log message ***
#
# Revision 1.96  2001/02/08 09:31:51  cocard
# *** empty log message ***
#
# Revision 1.95  2001/01/31 17:02:08  cocard
# *** empty log message ***
#
# Revision 1.94  2001/01/31 16:48:51  cocard
# *** empty log message ***
#
# Revision 1.93  2000/12/21 16:16:17  favey
# change interpolation method of geoid progs to bilinear
#
# Revision 1.92  2000/12/11 16:09:51  favey
# new options, exponent and z
#
# Revision 1.91  2000/08/17 09:53:51  favey
# gmt_writegrid: allow comment
#
# Revision 1.90  2000/08/15 16:39:11  favey
# mapras now working
#
# Revision 1.89  2000/08/12 10:54:11  favey
# bugfix in free chain
#
# Revision 1.88  2000/07/26 08:45:12  favey
# bugfix int vs. char
#
# Revision 1.87  2000/07/25 13:19:13  favey
# introducing OS_AIX
#
# Revision 1.86  2000/07/06 16:54:58  cocard
# *** empty log message ***
#
# Revision 1.85  2000/05/12 09:23:11  favey
# DEVEL_MAT_DEBUG
#
# Revision 1.84  2000/04/17 16:29:34  favey
# by Marcel Sutter
#
# Revision 1.83  2000/04/17 15:54:25  cocard
# ..
#
# Revision 1.82  2000/03/28 17:09:04  cocard
# *** empty log message ***
#
# Revision 1.81  2000/03/17 10:10:14  favey
# bugfix memory bursting
#
# Revision 1.80  2000/03/09 16:17:07  ypeter
# - vitacc : Fully revised. More options, more error checks
#
# Revision 1.79  2000/03/06 17:02:27  favey
# GMT colortable and grid output
#
# Revision 1.78  2000/02/25 13:02:28  ypeter
# Added applic. vitacc
# Modified function call to m_polynom_approx (added weighted obs.)
# Moved m_polynom_interpol from winterpol to matcalc
#
# Revision 1.77  2000/02/22 13:24:07  ypeter
# tconv: added option -N
#
# Revision 1.76  1999/12/02 14:10:23  favey
# *** empty log message ***
#
# Revision 1.75  1999/11/08 15:10:16  favey
# mkgrd: new option -f, -D
#
# Revision 1.74  1999/11/03 07:53:05  ypeter
# global.h: remove definitions of CMAX, cmax, CMIN and cmin
#
# Revision 1.73  1999/11/02 16:02:34  favey
# fixed bug memory leak in tconv
#
# Revision 1.72  1999/10/14 11:20:05  favey
# bugfixes
#
# Revision 1.71  1999/10/12 12:52:08  favey
# new option -M to mkgrd
#
# Revision 1.70  1999/10/11 15:47:22  favey
# beta release of mkgrd
#
# Revision 1.69  1999/10/07 16:36:43  favey
# new unit and pgn xyz2grd, mkgrd
#
# Revision 1.68  1999/10/07 14:52:31  cocard
# *** empty log message ***
#
# Revision 1.67  1999/10/07 09:07:21  favey
# bugfix in gps_to_dat
#
# Revision 1.66  1999/09/28 13:10:30  favey
# new function gps_to_dat
#
# Revision 1.65  1999/09/03 14:53:22  cocard
# *** empty log message ***
#
# Revision 1.64  1999/08/18 08:58:43  cocard
# *** empty log message ***
#
# Revision 1.63  1999/08/17 13:12:49  cocard
# *** empty log message ***
#
# Revision 1.62  1999/08/17 12:59:23  cocard
# *** empty log message ***
#
# Revision 1.61  1999/08/06 16:08:34  cocard
# *** empty log message ***
#
# Revision 1.60  1999/07/29 13:22:58  cocard
# *** empty log message ***
#
# Revision 1.59  1999/07/06 09:38:54  ypeter
# *** empty log message ***
#
# Revision 1.58  1999/07/06 07:54:10  favey
# mjdnow bugfix microsecond instead of millisecond
#
# Revision 1.57  1999/06/30 15:35:50  favey
# new app mapras
#
# Revision 1.56  1999/06/30 12:17:47  favey
# new function iszipped()
#
# Revision 1.55  1999/06/21 12:23:35  favey
# *** empty log message ***
#
# Revision 1.54  1999/06/21 09:48:21  favey
# added file name to struct for easier logging
#
# Revision 1.53  1999/06/10 15:45:43  cocard
# *** empty log message ***
#
# Revision 1.52  1999/06/03 14:59:47  cocard
# *** empty log message ***
#
# Revision 1.51  1999/06/03 10:10:00  cocard
# *** empty log message ***
#
# Revision 1.50  1999/06/03 09:47:29  cocard
# *** empty log message ***
#
# Revision 1.49  1999/04/19 09:51:26  cocard
# *** empty log message ***
#
# Revision 1.48  1999/04/08 13:33:11  favey
# bugfix in time_index_data_window
#
# Revision 1.47  1999/04/08 12:32:44  cocard
# *** empty log message ***
#
# Revision 1.46  1999/04/07 11:46:01  cocard
# *** empty log message ***
#
# Revision 1.45  1999/04/07 08:00:38  cocard
# changed to molnya00
#
# Revision 1.44  1999/03/26 14:53:34  favey
# molnya
#
# Revision 1.43  1999/03/26 14:52:45  favey
# molnya00
#
# Revision 1.42  1999/03/08 17:25:28  favey
# time string argument conversion added in cmdline.c
#
# Revision 1.41  1999/02/13 09:36:23  favey
# added winterpol
#
# Revision 1.40  1999/02/10 18:28:54  favey
# changed to make commonlib compile on solaris
#
# Revision 1.39  1999/02/10 17:28:10  favey
# changed PI to M_PI for portability
#
# Revision 1.38  1999/02/01 17:03:29  favey
# mat3 changed to sed version
#
# Revision 1.37  1999/01/26 17:19:54  favey
# angles.h moved out of swiss.h
#
# Revision 1.36  1999/01/25 16:04:11  favey
# added put value routines in gridop
#
#
