/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/global.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: main definitions valid for entire library compilation
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: global.c,v 1.7 2000/07/25 13:18:46 favey Exp favey $
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
 *            30-10-1998  --  created global
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/global.c $
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.7  2000/07/25 13:18:46  favey
 * introducing OS_AIX
 *
 * Revision 1.6  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.5  1999/09/03 14:53:12  cocard
 * *** empty log message ***
 *
 * Revision 1.4  1999/02/10 18:03:59  favey
 * changes to make commonlib compile on Solaris
 *
 * Revision 1.3  1998/10/30 16:40:34  favey
 * added RCS verbosity
 *
 *            
 *
 */

#include "global.h"

/*
 *  Programs should use these string in version information to show
 *  the library version (common_version) they are linked with.
 *  common_compileinfo contains info on the compilation time 
 *  of the common library.  Makes debugging easier!
 *
 */
char *common_version = LIBVER;
char *common_compileinfo = LIBCOMP;

#if defined(OS_SOLARIS) || defined(OS_AIX) || defined(OS_CYGWIN)
int isblank(int c)
{
	if (c == ' ')  return 1;
	if (c == '\t') return 1;
	return 0;
}
#endif
