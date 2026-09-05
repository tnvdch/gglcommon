#
# Makefile for the GGL Common Library
#
# We inherit the following from the top Makefile:
#
# COMMONDIR
APPDIR=../Documentation
# LIBVERSIONMAJOR
# LIBVERSIONMINOR
LIBVER="\"GGL Common Library Version $(LIBVERSIONMAJOR).$(LIBVERSIONMINOR)\""
# PLATFORM
# CC
# CFLAGS
# LDFLAGS

THISPART=sgml
#
# The object binaries we will create
#
OBJECTS=errproc parse inout window \
        cmdline mat3 quatern matcalc gridop ggl_time angles \
        trafo crdtrans swiss chain intersect closestpair

DREPOSITORY=$(REPOSITORY)/Documentation

##############################################################
#
# End of development configuration
#
##############################################################

#
# Where we will put our objects and headers
#

RCSFILES=$(foreach thename,$(OBJECTS),$(DREPOSITORY)/$(THISPART)/$(thename).sgml,v )

#
# make binaries
#
lib :
	@echo "Please run \`make' in the distribution root directory"

co :
	co -l $(DREPOSITORY)/$(THISPART)/Makefile,v
	co -l $(DREPOSITORY)/$(THISPART)/comabstract.sgml,v
	co -l $(RCSFILES)

coro :
	co -r $(DREPOSITORY)/$(THISPART)/Makefile,v
	co -r $(DREPOSITORY)/$(THISPART)/comabstract.sgml,v
	co -r $(RCSFILES)

ci :
	ci -u $(DREPOSITORY)/$(THISPART)/Makefile,v
	ci -u $(DREPOSITORY)/$(THISPART)/comabstract.sgml,v
	ci -u $(RCSFILES)

