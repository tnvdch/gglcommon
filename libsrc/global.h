/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/global.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: main definitions valid for entire library compilation
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: global.h,v 1.10 2001/07/03 08:47:13 cocard Exp favey $
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
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch>
 *            Favey Etienne   <favey@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *                  1987  --  created
 *            19-03-1997  --  pascal version to port
 *            07-04-1997  --  c port
 *            21-04-1998  --  restructured makefiles
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/global.h $
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.10  2001/07/03 08:47:13  cocard
 * *** empty log message ***
 *
 * Revision 1.9  2001/06/08 12:32:23  cocard
 * *** empty log message ***
 *
 * Revision 1.8  2000/07/25 13:18:46  favey
 * introducing OS_AIX
 *
 * Revision 1.7  1999/11/03 07:52:04  ypeter
 * Remove definition of CMAX, cmax, CMIN and cmin (unused until now)
 * Definition conflict of CMIN with /usr/include/sys/ttydefaults.h
 *
 * Revision 1.6  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.5  1999/02/10 18:03:59  favey
 * changes to make commonlib compile on Solaris
 *
 * Revision 1.4  1999/01/07 12:57:32  cocard
 * SQR added
 *
 * Revision 1.3  1998/11/03 18:00:51  favey
 * added NEARLYZERO, removed from matcalc.h
 *
 * Revision 1.2  1998/10/30 16:40:34  favey
 * added RCS verbosity
 *
 *            
 *
 */

/*
 *  Standard Libraries -- math functions
 *  please compile your programs with -lm flag!
 */
#include <math.h>

#ifndef GGL_COMMON
#define GGL_COMMON

/*
 **********************************************************************
 * just for those people who prefer to use PI instead of M_PI
 */
#ifndef PI
#define PI M_PI
#endif

/*
 **********************************************************************
 *
 *  Start of configurable section
 *
 **********************************************************************
 *
 *  Do you want the lib functions to produce error messages?
 */
#define WITH_ERRORS
/*
 *  Can you rely on GNU and POSIX specific routines?
 */
#define OS_LINUX
/* possible options: OS_LINUX OS_AIX OS_SOLARIS OS_CYGWIN */

/*
 *  Do you need additional debugging information?
 */
/* #define WITH_DEBUG */


/*
 **********************************************************************
 *
 *  End of configurable section
 *
 **********************************************************************
 *
 *  some miscellanous definitions (I should remove them sooner or later)
 */

#define sqr(x) ((x)*(x))
#define SQR(x) ((x)*(x))

#define sign(x) ((x<0)?-1:+1)
#define SIGN(x) ((x<0)?-1:+1)

//#define cmax(a,b) ((a<b)?b:a)
//#define CMAX(a,b) ((a<b)?b:a)

//#define cmin(a,b) ((a<b)?a:b)
//#define CMIN(a,b) ((a<b)?a:b)

#define trunc2(x) floor(x)
#define TRUNC2(x) floor(x)

#define NEARLYZERO (1.0e-15)

#define MAXLINELENGTH (1024)

#if defined(OS_SOLARIS) || defined(OS_AIX) || defined(OS_CYGWIN)
int isblank(int c);
#endif

#endif
