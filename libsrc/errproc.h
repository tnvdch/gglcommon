/**
 **        FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/errproc.h,v $
 **
 ** DESCRIPTION: This hack I made, because I did not find the
 **              two useful variables `program_invocation_name'
 **              and `program_invocation_short_name' to do
 **              my own error procrssing on my local machine.
 **
 **   COPYRIGHT: Copyright (C) 1996 by Etienne Favey
 ** $Id: errproc.h,v 1.6 2000/07/25 13:18:46 favey Exp $
 **
 **  DISCLAIMER:   This program is free software; you can redistribute
 **              it and/or modify it under the terms of the GNU General
 **              Public License as published by the Free Software
 **              Foundation; either version 2 of the License, or
 **              (at your option) any later version.
 **                This program is distributed in the hope that it will
 **              be useful, but WITHOUT ANY WARRANTY; without even the
 **              implied warranty of MERCHANTABILITY or FITNESS FOR A
 **              PARTICULAR PURPOSE. See the GNU General Public License
 **              for more details.
 **                You should have received a copy of the GNU General
 **              Public License along with this program; if not, write
 **              to the Free Software Foundation, Inc., 675 Mass Ave,
 **              Cambridge, MA 02139, USA.
 **
 **      AUTHOR: Etienne Favey <eti@favey.lugs.ch>
 **              Zurich, Switzerland
 **
 **     CHANGES: dd-mm-yyyy -- description
 **              20-10-1996 -- written and tested
 **              19-08-1997 -- added error codes
 ** $Log: errproc.h,v $
 ** Revision 1.6  2000/07/25 13:18:46  favey
 ** introducing OS_AIX
 **
 ** Revision 1.5  1998/11/20 08:18:01  favey
 ** new function gputerr
 **
 **
 **/

/* I like to use each identifier just once ;-) */
#ifndef fifi_ERR_PROC
#define fifi_ERR_PROC

/* for better readability for former Pascal users */
typedef int boolean;
#ifndef OS_AIX
/* #define true (1) */
/* #define false (0) */
#define TRUE (1)
#define FALSE (0)
#endif

/*
 *  error codes to make source code reading easier
 *  we start at 500 so we do not overlay the errno error codes
 */

#define GGLERR_INOUT 501		/* error in file input/output */
#define GGLERR_INVAL 502		/* invalid format used */
#define GGLERR_NOMEM 503		/* no more memory available */

#define GGLWRN_OVWRI 601		/* overwriting old values */

/*
 *  Errornumbers:
 *  700 to 799 are reserved for local use in libfunctions
 */

/*
 *  extend global variables for verbosity
 * strange problems in c++ --> removed in c++ version
 */
#ifndef CPP
char *program_invocation_name;        /* dup of the same as argv[0] */
char *program_invocation_short_name;  /* above w/o preceeding path  */
#endif
/*
 *  this funcion is defined in errproc.c
 */
void                              /* Initialize the above variables */
errinit (const char *program_name);     /*   with the appropriate values. */


/*
 *  this funcion is defined in errproc.c
 */
void puterr (
	const char *gglpart,	/* name of library object */
	const char *routine)	/* name of function */
;
void gputerr (
	const char *gglpart,	/* name of library object */
	const char *routine)	/* name of function */
;
#endif







