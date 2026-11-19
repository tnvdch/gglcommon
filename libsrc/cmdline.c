/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/cmdline.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: command line parsing made easy
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: cmdline.c,v 1.9 2001/06/08 12:32:23 cocard Exp favey $
 *
 *  DISCLAIMER:   This program is free software; you can redistribute
 *              it and/or modify it under the terms of the GNU General
 *              Public License as published by the Free Software
 *              Foundation; either version 2 of the License, or
 *              (at your option) any later version.
 *                This program is distributed in the hope that it will
 *              be useful, but WITHOUT ANY WARRANTY; without even the
 *              implied warranty of MERCHANTABILITY or FITNESS FOR A
 *              PARTICULAR PURPOSE. See the GNU General Public License
 *              for more details.
 *                You should have received a copy of the GNU General
 *              Public License along with this program; if not, write
 *              to the Free Software Foundation, Inc., 675 Mass Ave,
 *              Cambridge, MA 02139, USA.
 *
 *   AUTHORS: Favey Etienne   <favey@geod.ethz.ch>
 *
 *   HISTORY: 
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/cmdline.c $
 * 
 * 3     30.01.02 16:44 Etienne
 * added incremental option
 * 
 * 2     12.12.01 15:37 Etienne
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.9  2001/06/08 12:32:23  cocard
 * *** empty log message ***
 *
 * Revision 1.8  2000/07/26 08:44:52  favey
 * bugfix int vs. char
 *
 * Revision 1.7  1999/09/08 09:32:31  favey
 * copyright info added
 *
 * Revision 1.6  1999/03/08 17:24:14  favey
 * added new type CMDL_TIME
 * automatic conversion of time string argument to double
 *
 * Revision 1.5  1998/10/30 14:22:07  favey
 * changed config.h to global.h
 *
 * Revision 1.4  1998/05/05 08:18:21  favey
 * *** empty log message ***
 *
 * Revision 1.3  1998/05/05 08:05:57  favey
 * added description line
 *
 * Revision 1.2  1998/05/04 10:05:36  favey
 * *** empty log message ***
 *
 * Revision 1.1  1998/05/04 07:02:14  favey
 * Initial revision
 *
 *
 *
 */

#include <errno.h>	/*  the allround headers */
#include <stdlib.h>	/*  the allround headers */
#include <stdio.h>	/*  crucial for file input and output   */
#include <unistd.h>	/*  to parse the command line */
#include <ctype.h>	/*  character handling */
#include <string.h>	/*  string handling */
#include <math.h>	/*  mathematics */


#include "global.h"
#include "errproc.h"
#include "ggl_time.h"
#include "cmdline.h"



#define GGL_MAXOPTS 512





/*******************************************************************/
/*
 *  register commandline option
 */
void print_usage_info (
	const int verb,
	struct mainoptions o[],
	const int numelem,
	const char *description,
	const char *nonopt)
{
	extern char *program_invocation_short_name;
	int numopt;
	int wrlen;

	fprintf (stderr, "%s --- %s\n\n", 
		 program_invocation_short_name, description);
	fprintf (stderr, "usage:\n   %s [-h]" , program_invocation_short_name);
	wrlen = 8 + strlen(program_invocation_short_name);
	for (numopt=0;numopt<numelem;numopt++) {
		if (wrlen > 55) {
			wrlen = 6;
			fprintf (stderr, "\n      ");
		}
		if (o[numopt].valid) {
			wrlen += 5 + strlen (o[numopt].shortdesc);
			fprintf (stderr," [-%c%s]", o[numopt].optchar,
				 o[numopt].shortdesc);
		}
	}
	fprintf (stderr, "\n      %s", nonopt);
	fprintf (stderr, "\n\n");
	
	if (verb != FALSE) for (numopt=0;numopt<numelem;numopt++) {
	  if (o[numopt].valid) {
		 fprintf (stderr," -%c     %s", o[numopt].optchar,
					 o[numopt].longdesc);
		 switch (o[numopt].argtype) {
		 case 1: // int
			fprintf (stderr, ", default: %d", *((int *)o[numopt].argument));
			break;
		 case 2: // double
			fprintf (stderr, ", default: %15.9f", *((double *)o[numopt].argument));
			break;
		 default:
		   ;
		 }
		 fprintf (stderr, "\n");
	  }
	} else
		fprintf (stderr,
			 "  -h    Prints more detailed help on each option\n");
	fprintf (stderr, "\n");
}




/*******************************************************************/
/*
 *  fill argument and flag
 */
static int eval_cmdline_opt(
	struct mainoptions *os,
	const char *oa)
{
	extern char *program_invocation_short_name;
	double tmpt;
	char *ptmp = NULL;

	if (os->argument != NULL) switch (os->argtype) {
	case 0: /* char * */
		/* we overwrite the default string leaving rubbish */
		*((char **) os->argument) = strdup (oa);
		break;
	case 1: /* int */
		*((int *) os->argument) = strtol(oa, &ptmp, 10);
		if (ptmp == oa)
		{
			fprintf (stderr, 
				 "%s: Option `-%c' requires integer argument.\n",
				 program_invocation_short_name, os->optchar);
			return -1;
		}
		break;
	case 2: /* double */
		*((double *) os->argument) = strtod(oa, &ptmp);
		if (ptmp == oa)
		{
			fprintf (stderr, 
				 "%s: Option `-%c' requires floating point argument.\n",
				 program_invocation_short_name, os->optchar);
			return -1;
		}
		break;
	case 3: /* time converted to double */
		tmpt = sreadtime2 (oa);
		if (tmpt < 0.0)
		{
			fprintf (stderr, 
				 "%s: Option `-%c' requires time string argument.\n",
				 program_invocation_short_name, os->optchar);
			return -1;
		}
		*((double *) os->argument) = tmpt;
		break;
	default:
		fprintf (stderr, "%s: DEBUG: Unknown option type %d.\n", 
			 program_invocation_short_name, os->argtype);
		return -2;
	}

	if (os->flag != NULL) {
	  if (os->argtype == CMDL_INCR) {
		 *(os->flag) += os->valseen;
	  } else {
		 *(os->flag) = os->valseen;
	  }
	}

	return 0;
}


/*******************************************************************/
/*
 *  register commandline option
 */
int process_cmdline_opt (
	const int argcount,
	char **argvect,
	struct mainoptions o[],
	const int numelem,
	const char *description,
	const char *nonopt)
{
	extern char *program_invocation_short_name;
	
        int numopt;
	char shortopts[GGL_MAXOPTS];
	int so_len = 1;
	int arg_char;
	int err;

	shortopts[0] = 'h';
	opterr = 0;

	for (numopt=0;numopt<numelem;numopt++) {
		if (o[numopt].optchar != '\0') {
			o[numopt].valid = TRUE;
			shortopts[so_len++] = o[numopt].optchar;
			if (o[numopt].argument != NULL)
				shortopts[so_len++] = ':';
			if (o[numopt].argument == NULL &&
			    o[numopt].flag     == NULL   )
				o[numopt].valid = FALSE;
		} else {
			o[numopt].valid = FALSE;
		}
	}
	shortopts[so_len] = '\0';

	while ((arg_char = getopt (argcount, argvect, shortopts)) != -1) {
		if (arg_char == '?') {   /* oops, an unknown option */
			if (isprint (optopt))
				fprintf (stderr, 
					 "%s: Unknown option `-%c'.\n",
					 program_invocation_short_name,
					 optopt);
			else
				fprintf (stderr,
				"%s: Unknown option character `\\x%x'.\n",
					 program_invocation_short_name, 
					 optopt);
			return -1;
		}

		if (arg_char == 'h') {   /* -h   the help option */
			print_usage_info (TRUE, o, numelem, description,nonopt);
			return -3;
		}

		for (numopt=0;numopt<numelem;numopt++)
			if ((o[numopt].valid) && o[numopt].optchar==arg_char) {
				err = eval_cmdline_opt(&(o[numopt]), optarg);
				if (err != 0) return err;
			}
	}

	return optind;
}
