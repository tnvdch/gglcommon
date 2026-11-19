/**
 **        FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/errproc.c,v $
 **
 ** DESCRIPTION: This hack I made, because I did not find the
 **              two useful variables `program_invocation_name'
 **              and `program_invocation_short_name' to do
 **              my own error processing on my local machine.
 **
 **   COPYRIGHT: Copyright (C) 1996 by Etienne Favey
 ** $Id: errproc.c,v 1.8 2001/06/08 12:32:23 cocard Exp favey $
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
 **              20-10-1996 -- written and tested :-)
 ** $Log: /GPS/SOFTWARE/common_devel/libsrc/errproc.c $
 * 
 * 1     31.08.01 14:51 Etienne
 ** Revision 1.8  2001/06/08 12:32:23  cocard
 ** *** empty log message ***
 **
 ** Revision 1.7  1999/09/08 09:12:18  favey
 ** copyright info added
 **
 ** Revision 1.6  1998/11/20 08:18:01  favey
 ** new function gputerr
 **
 ** Revision 1.5  1998/11/16 17:58:06  favey
 ** changed syntax of puterr
 **
 ** Revision 1.4  1998/10/30 14:22:07  favey
 ** changed config.h to global.h
 **
 ** Revision 1.3  1998/05/29 16:56:53  favey
 ** *** empty log message ***
 **
 ** Revision 1.2  1998/05/04 07:02:14  favey
 ** *** empty log message ***
 **
 ** Revision 1.1  1998/04/22 12:31:10  favey
 ** Initial revision
 **
 ** Revision 0.3  1997/12/24 12:10:20  favey
 ** followup of V. 0.29
 **
 **
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
/*
 *  unix specific headers
 */
#include <syslog.h>
/*
 *  our own headers
 */
#include "global.h"
#include "errproc.h"

/*
 *  this funcion is declared in errproc.h
 *  it should receive argv[0] as parameter
 */
void                                    /* Initialize the variables below */
errinit (const char *program_name)      /*   with the appropriate values. */
{
#if 0
  extern char *program_invocation_name;        /* dup of the same as argv[0] */
  extern char *program_invocation_short_name;  /* above w/o preceeding path  */
  char *argument;
  char *token, *final;

  argument = (char *) malloc(strlen(program_name)+2);
  /*
   *  let's be independent of what a user does later with argv[0]
   */
  program_invocation_name = strdup (program_name);

  /*
   *  We have to duplicate things, because strtok(3) modifies
   *  the string passed to it! We want to leave argv[0] untouched!
   */
  strcpy (argument, program_name);
  /*
   *  add a sentinel ...
   */
  strncat (argument, "/", 1);
  /*
   *  ... initialize final with a first call to strtok(3) ...
   */
  final = strtok(argument, "/");
  /*
   *  ... repeat until the last word is found, and strtok(3) returns NULL
   */
  while ( (token = strtok(NULL, "/")) )
    final = token;
  /*
   *  `argument', thus `final' as well, is locally allocated
   *  and will be discarded, so we duplicate again
   */ 
  program_invocation_short_name = strdup (final);
  free (argument);
#endif
}

/*
 *  this funcion is declared in errproc.h
 *  purpose: write an error message to stderr or to a logfile
 *  or to the syslog facility
 *
 *  PORTING: Be sure to port this function if you are using a 
 *  non unixish platform
 */
void puterr (
	const char *gglpart,	/* name of library object */
	const char *routine)	/* name of function */
{
	fprintf (stderr, "GGL Common Library -- %s\n",
		 (gglpart == NULL)?"general":gglpart);
	fprintf (stderr, "%s: ", (routine == NULL)?"nofunc":routine); 
}


void gputerr (
	const char *gglpart,	/* name of library object */
	const char *routine)	/* name of function */
{
	fprintf (stderr, "GGL GPSkin Library -- %s\n",
		 (gglpart == NULL)?"general":gglpart);
	fprintf (stderr, "%s: ", (routine == NULL)?"nofunc":routine); 
}


