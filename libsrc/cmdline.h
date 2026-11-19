/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/cmdline.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: command line parsing made easy
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: cmdline.h,v 1.5 1999/09/08 09:32:31 favey Exp favey $
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
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/cmdline.h $
 * 
 * 2     30.01.02 16:44 Etienne
 * added incremental option
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.5  1999/09/08 09:32:31  favey
 * copyright info added
 *
 * Revision 1.4  1999/03/08 17:24:14  favey
 * added new type CMDL_TIME
 *
 * Revision 1.3  1998/05/13 13:52:21  favey
 * *** empty log message ***
 *
 * Revision 1.2  1998/05/05 08:06:11  favey
 * added description line
 *
 * Revision 1.1  1998/05/04 07:02:34  favey
 * Initial revision
 *
 *
 *
 */

#ifndef GGL_CMDLINE
#define GGL_CMDLINE

#ifdef OS_CYGWIN
#include <getopt.h>
#endif

#define CMDL_CHAR   (0)
#define CMDL_INT    (1)
#define CMDL_DOUBLE (2)
#define CMDL_TIME   (3)
#define CMDL_INCR   (4)


/*******************************************************************/
/*
 *  The structure to hold all the valid program options
 *
 *  valid	a flag for internal checks, please do not set.
 *  optchar	short option, '\0' if no short option
 *  argument	a pointer where to store the agument, NULL no arg
 *  argtype	0 = char * ; 1 = int ; 2 = double ; 3 = time ;
 *  flag	variable to set, when option is seen, NULL no flag
 *  valseen	the value to set flag to when option was seen
 *  shortdesc	usage info : command line element: e.g. "[-i <input_file>]"
 *  longdesc	usage info : help line: e.g." -i   data input file"
 *  
 *  you must provide either `argument' or `flag' or both.
 */
struct mainoptions {
	int valid;
	int optchar;
	void *argument;
	int argtype;
	int *flag;
	int valseen;
	char *shortdesc;
	char *longdesc;
} ;


/*******************************************************************/
/*
 *  print commandline option information
 */
void print_usage_info (
	const int verb,
	struct mainoptions o[],
	const int numelem,
        const char *description,
	const char *nonopt);




/*******************************************************************/
/*
 *  register commandline options
 */
int process_cmdline_opt (
	const int argcount,
	char **argvect,
	struct mainoptions o[],
	const int numelem,
        const char *description,
	const char *nonopt);


#endif
