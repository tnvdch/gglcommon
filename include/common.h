/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/global.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: main definitions valid for entire library compilation
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: global.h,v 1.11 2001/08/31 12:08:24 favey Exp $
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
 * $Log: global.h,v $
 * Revision 1.11  2001/08/31 12:08:24  favey
 * cygwin adaptions
 *
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
/**
 **        FILE: angles.h
 **
 ** DESCRIPTION: Functions to transform coordinates
 **              from the bessel ellipsoid to the
 **              Swiss military coordinates and the
 **              vice versa.
 **              Constants and ellipsoid data for calculations
 **              with the swiss projection system.
 **
 **              This header file contains the declarations for
 **              swiss.c and angles.c
 **
 **   COPYRIGHT: Copyright (C) 1996 by Etienne Favey
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
 **     VERSION: 4
 **
 **     CHANGES: dd-mm-yyyy -- description
 **              03-12-1996 -- first release (Version 1.05)
 **              18-12-1996 -- bugfix xyz2blh iteration criterion
 **              10-04-1997 -- adjusted to use with GGL coordtr.c
 **                            taken from public version 1.08,
 **                            renamed blh to llh, introduced
 **                            struct ellipsoid and vectors
 **              07-05-1997 -- added admsstr(), ellsupport, det_ell()
 **              14-01-1998 -- added struct frameconv, select_ellipsoid(),
 **                            select_frameconv()
 **              20-01-1998 -- split swiss.c into swiss.c and trafo.c
 **              13-05-1998 -- parameter for ch_magnetic_deviation changed
 **
 **/
#ifndef ANGLES
#define ANGLES


/* ----------------------------------------------------------------
 * for angles.c
 */

struct hms {
    int deg;
    int min;
    double sec;
  } ;

struct Angle {
  double rad;
  double grad;
  double deg;
  struct hms dms;
} ;

struct Angle afromrad (double angrad);
struct Angle afromdeg (double angdeg);
struct Angle afromgrad(double anggrad);

struct Angle afromdmsverb(int degrees, int minutes, double seconds);
struct Angle afromdms(struct hms mydms);
struct Angle afromstr(char *s);

char * admsstr (const struct Angle a, const int mode);
/*
 * mode: 0 -> 46 30 25.643346789
 *       1 -> 46d30m25.643346789
 *       2 -> 46\t30\t25.643346789
 */
#endif

/**
 **        FILE: trafo.h
 **
 ** DESCRIPTION: Functions to transform coordinates
 **              from the bessel ellipsoid to the
 **              Swiss military coordinates and the
 **              vice versa.
 **              Constants and ellipsoid data for calculations
 **              with the swiss projection system.
 **
 **              This file contains the declarations for trafo.c
 **
 **   COPYRIGHT: Copyright (C) 1996 by Etienne Favey
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
 **     VERSION: 4
 **
 **     CHANGES: dd-mm-yyyy -- description
 **              03-12-1996 -- first release (Version 1.05)
 **              18-12-1996 -- bugfix xyz2blh iteration criterion
 **              10-04-1997 -- adjusted to use with GGL coordtr.c
 **                            taken from public version 1.08,
 **                            renamed blh to llh, introduced
 **                            struct ellipsoid and vectors
 **              07-05-1997 -- added admsstr(), ellsupport, det_ell()
 **              14-01-1998 -- added struct frameconv, select_ellipsoid(),
 **                            select_frameconv()
 **              20-01-1998 -- split swiss.c into swiss.c and trafo.c
 **
 **/
#ifndef TRAFO
#define TRAFO


/*
 *  Frame Conversion Structure to simplify datum shifts
 */
struct frameread {
	int fwd;	/* sign for forward or reverse conversion */
	double ref_epoch;
	double dx[3];	/* translation in [m] */
	double scale;	/* scale - 1 in [ppm] */
	double rot[3];	/* small rotation elements [rad] */
	double ddx[3];
	double dscale;	/* rates of the above per year */
	double drot[3];
} ;


struct frameconv {
	int fwd;
	double scale;
	double rm[3][3];
	double dx[3];
} ;


/*
 * Ellipsoid patameters for several frequently used frames
 */
struct ellipsoid {
  double a;
  double e2;
} ;

/*
 * conversion of cartesian coordinates to spherical
 */
void cartesian_to_spherical(
			    double xyz[3], 
			    double *lat, 
			    double *lon, 
			    double *radius);

/*
 * conversion of spherical to cartesian coordinates
 */
void spherical_to_cartesian(
			    double lat, 
			    double lon, 
			    double radius, 
			    double xyz[3]);

/*
 *  geocentric to ellipsoidal
 *  parameters have to be passed and are returned in rad
 */ 
int xyz_to_llh(
	       const struct ellipsoid ell,
	       const double xyz[3],
	       double *lat,
	       double *lon,
	       double *h);


/*
 *  ellipsoidal to geocentric
 *  parameters have to be passed and are returned in rad
 */
int llh_to_xyz(
	       const struct ellipsoid ell,
	       const double lat,
	       const double lon,
	       const double h,
	       double xyz[3]);


/*
 *  select ellipsoid from configuration file
 *  returns pointer to struct or NULL upon failure
 *  free struct with free(3) after finishing using it
 */
struct ellipsoid * select_ellipsoid (
	const char *ellstr)
;


/*
 *  transform arbitrary frame given source epoch
 */
int frame_conv (
	struct frameconv *fcp,
	double source[3],
	double target[3])
;


/*
 *  select frame conversion parameters from file
 *  returns pointer to struct or NULL upon failure
 *  free struct with free(3) after finishing using it
 */
struct frameconv * select_frameconv (
	const char *source,
	const double inepoch,
	const char *target,
	const double outepoch)
;

/*
 *  select form types for input
 *  see trafo config file for more information
 */
int select_formtype (
	const char *formtype,
	int ft_type[3],
	int ft_cols[3])
;

/*
 *  select geoid from trafo config file
 *  returns pointer to struct or NULL upon failure
 *  free struct with free(3) after finishing using it
 */
struct Grid * select_geoid (
	const char *geoidstr,
	int g_ftyp[3])
;

#endif

/*
 *      FILE: chain
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: routines for keeping data in a chain structure
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            20-05-1998  --  created
 */

/*
 * structure of a chain :
 *
 *  nil <- head <--> node1 <--> node2 ....  <--> last node  -> nil
 *          |          |          |                  |  
 *         nil       data1      data2              data
 *
 * for cyclic chains : node1 and last node are connected 
 * (but head->next still points to node1)
 */

/*
 * main definition of node
 */
#ifndef CHAIN
#define CHAIN

typedef int (*comparefunction)(void *, void *);

struct Node
{
  struct Node *prev;  /* pointer to previous node */
  struct Node *next;  /* pointer to next node */ 
  void *data;         /* void pointer to data */
};

struct Node * chain_init(int *count);           /* initialize a new chain */
struct Node * chain_newelement(void * data);    /* allocates memory for a new node */
/* free memory for an element */
void chain_freeelement (struct Node * p, int removedata);
void chain_free(struct Node *chain, int removedata, int count); /* frees the  entire chain */

/* inserts  the elememt in the chain (after chain) */
struct Node * chain_insert(struct Node *chain, struct Node *element, int *count);

/* inserts the element as last member in the chain */
struct Node * chain_append(struct Node *chain, struct Node *element, int *count);

/* inserts the  element at the correct position using compare,
 *  return values of compare :
 *      positive  . if element1 > element2
 *      zero        if elememt1 = element2
 *      negative    if element1 < element2
 */

struct Node * chain_insertsorted (struct Node *chain, struct Node *element, 
				 comparefunction compare, int *count);

/* remove actual element from chain */
struct Node * chain_remove (struct Node * chain, int *count);

struct Node * chain_head(struct Node *chain);  /* returns the pointer to the head of the chain */
struct Node * chain_last(struct Node *chain);  /* returns the pointer to the last node */
struct Node * chain_next(struct Node *chain);  /* next node */
struct Node * chain_prev(struct Node *chain);  /* previous node */

/*
  makes the chain cyclic by closing last with first element
  head is returned must be saved for freeing the chain
 */
struct Node * chain_cyclic(struct Node *chain);

struct Node * chain_acyclic(struct Node *head);   /*   makes the cyclic chain acylcic  */

/* converts a chain to an array (with allocation of memory for array of pointers) */
void * chain_to_array(struct Node * chain, int count);

/*  converts an array to a chain (with allocation of memory for chain)  */
struct Node * array_to_chain(void **array, int count);

void array_free(void **array, int removedata, int count);  /*  frees the array  */

/* sorts the elements from node1 to node2 using quicksort algorithm */
void chain_quicksort1(struct Node *node1,  struct Node *node2, comparefunction compare);

/* sorts the entire chain and returns head of sorted chain (uses chain_quicksort1)  */
struct Node * chain_quicksort2(struct Node *chain,  comparefunction compare);

#endif
/****************************************************************************

Routines:	routines necessary to run the programs coneccted with
		the robot calibration project

Copyright:      (c) 1995 Geodesy and Geodynamics Lab,
                     ETH Zurich

Author:		Etienne Favey & Andreas Schlatter
Tutor:		Hans-Gerd Maas

Address:	Institute of Geodesy and Photogrammetry
		ETH - Hoenggerberg
		CH - 8093 Zurich

Creation Date:	11/1995

Description:		
	
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
****************************************************************************/


/* used in closest pair, source: clop.c */
#define DISTANCE(a,b) sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y))

typedef struct		/* used in closest pair, source: closestpair.c */
{ 
	long int pnr;
	double x;
	double y;
} coord_2d;

double	clop();		/* closest pair algorhithm, source: closestpair.c */

int countinbox();	/* find in box a set of points, source: findinbox.c */

/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/cmdline.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: command line parsing made easy
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: cmdline.h,v 1.7 2002/01/30 16:07:10 favey Exp $
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
 * $Log: cmdline.h,v $
 * Revision 1.7  2002/01/30 16:07:10  favey
 * cmdline: added output of default values and incremental option flag
 *
 * Revision 1.6  2001/08/31 12:08:24  favey
 * cygwin adaptions
 *
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
/*
 *      FILE: crdtrans
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations on 3x3 matrix and 3x1 vectors
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *
 *
 */

#ifndef GGL_CRDTRANS
#define GGL_CRDTRANS


double mod2pi (double x)
;

double center_angle( const double lat1,
		     const double lon1,
		     const double lat2,
		     const double lon2)
;

double ellipse_circumference (const double a,
			      const double ex)
;

double ellipse_approxdist (const double lat1,
			   const double lon1,
			   const double lat2,
			   const double lon2,
			   const struct ellipsoid ell)
;

int xyz_to_azel (double xyz1[3],
		 double xyz2[3],
		 struct ellipsoid ell,
		 double *dist,
		 double *elev,
		 double *azi)
;

int top_to_xyz (double lat,
		double lon,
		double dtop[3],
		double dxyz[3])
;

int top_to_xyz2 (double xyz[3],
		 double dtop[3],
		 double dxyz[3],
		 struct ellipsoid ell)
;

int azel_to_xyz (double xyz[3],
		 double azi,
		 double ele,
		 double dxyz[3],
		 struct ellipsoid ell)
;

int xyz_to_top (double lat,
		double lon,
		double dxyz[3],
		double dtop[3])
;

int xyz_to_top2 (double xyz[3],
		 double dxyz[3],
		 double dtop[3],
		 struct ellipsoid ell)
;

#endif
/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/d1linear.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: Linear one dimensional interpolation
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: d1linear.h,v 1.3 1999/09/08 09:32:31 favey Exp $
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
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            28-04-1998  --  created
 * $Log: d1linear.h,v $
 * Revision 1.3  1999/09/08 09:32:31  favey
 * copyright info added
 *
 * Revision 1.2  1998/04/28 16:02:11  favey
 * *** empty log message ***
 *
 *
 *
 */

#ifndef COMMON_D1LINEAR
#define COMMON_D1LINEAR

/*
 *  Simple linear interpolation. 
 */
int d1ipol_linear (
	const double x1,
	const double x2,
	const double y1,
	const double y2,
	const double xip,
	double *yres)
;


#endif
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







/*
 *      FILE: ggl_time
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: time operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: ggl_time.h,v 1.5 1999/09/28 13:10:03 favey Exp $
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
 *   AUTHORS: Cocard Marc     (pascal version)
 *            Favey Etienne   (c port)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *                  1987  --  created
 *            19-03-1997  --  pascal version to port
 *            07-04-1997  --  c port
 *            
 *
 */
#ifndef GGL_TIME
#define GGL_TIME

struct ttime {
  int year;
  int month;
  int day;
  int hour;
  int min;
  double sec;
  double mjd;
} ;


/*
 *   System time 0 is 00.00.00 January 1, 1970 UTC on POSIX systems
 */
#define SYSMJD0 (40587.0)

/*
 *  JulDat0 + mjd = JulDat
 */
#define JULDAT0  (2400000.5)

/*
 *   mjd for gpsweek = 0 and sec = 0
 */
#define GPSMJD0  (44244.0)

#define leapyear(year) (!((year%4!=0) || (year%100==0 && year%400!=0)) )

double dat_to_mjd (const int year, const int month, const int day,
		   const int hour, const int min, const double sec);

int datmjd (struct ttime *t);

int mjd_to_dat (double mjd,
		int *year, int *month, int *day,
		int *hour, int *min, double *sec);

int mjddat (struct ttime *t);

double gps_to_mjd (double sec,int week);

int mjd_to_gps (double mjd, double *sec, int *week);

/* returns secofday, struct ttime *t is filled with correct data */
/* routine addresses problem of insufficient double resolution of mjd */
double gps_to_dat (struct ttime *t, const int week, const double sec);

double yearday_to_mjd (int year,int daynr);

int mjd_to_yearday (double mjd, int *year, int *daynr);

int dayofyear1(struct ttime *t);
int dayofyear2(double mjd);

double mjdnow();

int systimenow (struct ttime *t);
int localtimenow (struct ttime *t);

char * strtime1 (struct ttime *t, int i);
char * strtime2 (double mjd, int i);
char * strtime3 (struct ttime *t);

int    readtime1  (FILE *fpin, struct ttime *t);
int    sreadtime1 (const char *instring, struct ttime *t);
double readtime2  (FILE *fpin);
double sreadtime2 (const char *instring);

double GMST (double mjd);

double secofday (const char *mjd);

#endif
/*
 *      FILE: gridop
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: grid operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            31-08-1997  --  created
 *
 */

#ifndef GGL_GRIDOP
#define GGL_GRIDOP

/* this is for making it compile on Solaris */
#ifndef NAN
#define NAN HUGE_VAL
#endif

/* environment variables for GMT rasters */
#define GMT_RASTERINFO_FILE "RASTERINFO"
#define GMT_GRID_DIR "GMT_GRIDDIR"

/* surfer specific */
#define SURFERMASK (1.70141e38)
#define SURFERMASK_STR " 1.70141e38"

struct Grid {
  char *grdfilename;     /* name of grid file */
  FILE *fp;              /* file pointer to gridfile */
  char *id;              /* identification string */
  long int nx, ny;       /* number of columns and rows */
  double xmin, xmax;     /* x-range */
  double ymin, ymax;     /* y-range */
  double zmin, zmax;     /* y-range */
  double dx, dy;         /* spacing */
  
  double offset, scale;  /* zvalue = offset + data*scale */
  double mask;           /* mask value */

  char datatype;         /* u = unsigned char */
			 /* c = char */
			 /* d = unsigned short int */
			 /* i = short int */
			 /* f = float */
			 /* e = double */      

  unsigned char datasize; /* size of corresponding datatype */
  long int d_offset;      /* offset to first element */
  long int d_row;         /* length of row */
  char *data;             /* pointer to data buffer */
};

/*--------------------------------------------------------------------------------- 
 *  GRID : initialize basic griddata structure without memory allocation for data 
 */
struct Grid * g_init0();

/*--------------------------------------------------------------------------------- 
 *  GRID : initialize memory for griddata depending on datatype  
 */
void * g_getmem4data(struct Grid * g);

/* ---------------------------------------------------------------------------------
 * GRID : set d_offset and d_row according to modus
 *   modus 0 = upside up   --> first element in databuffer = data(xmin, ymin)
 *   modus 1 = upside down --> first element in databuffer = data(xmin, ymax)
 */
void g_setdataparameter(struct Grid *g, int modus);

/*--------------------------------------------------------------------------------- 
 *  GRID : free grid 
 */
void g_free(struct Grid * g);

/*--------------------------------------------------------------------------------- 
 *  GRID : return offset in bytes for grid point (ix, iy)  
 */
long int g_index(struct Grid *g, const long int ix, const long int iy);

/*--------------------------------------------------------------------------------- 
 *  GRID : write grid info to file (for debugging purposes)
 */
void g_printinfo(FILE *fp, struct Grid *g);

/* ---------------------------------------------------------------------------------
 * GRID : scale limits  by factor scalex, scaley resp.
 */
void g_scalelimits(struct Grid *g, double xscale, double yscale);

/*--------------------------------------------------------------------------------- 
 *  GRID : return  value as a double from grid at position (ix,iy) 
 */
double g_getval(struct Grid *g, long int offset);

/*--------------------------------------------------------------------------------- 
 *  GRID : return  value as a double from grid at position (ix,iy) 
 */
double g_getvalue(struct Grid *g, const long int ix, const long int iy);

/* ---------------------------------------------------------------------------------
 * GRID : put  value as a double to grid at position (ix,iy) 
 */
int g_putval(struct Grid *g, long int offset, const double val);

/* ---------------------------------------------------------------------------------
 * GRID : put  value as a double to grid at position (ix,iy) 
 */
int g_putvalue(struct Grid *g, const long int ix, const long int iy, const double v);

/* ---------------------------------------------------------------------------------
 * GRID : looks for minimal and maximal values which are stored in g->zmin, g->zmax
 */
void g_minmax(struct Grid *g);

/*---------------------------------------------------------------------------------
 *  GRID : return the 4 corner values and x- and y-values for point 0
 *         ordering of the      y
 *         points 0,1,2,3 :     |
 *                              3 - 2
 *                              | x |
 *                              0 - 1 -> x
 */
int g_getcorner(struct Grid *g, const double x, const double y, 
		double z[4], double * x0, double *y0);

/*--------------------------------------------------------------------------------- 
 *  GRID : return  a vector xyz[3] with x-, y-, and z-value from grid at position (ix,iy) 
*/
int g_getpoint (struct Grid *g, const long int ix, 
                                const long int iy, double xyz[3]);

/* ---------------------------------------------------------------------------------
 * GRID : return the maximal value of the four neighbor
 * only 1 grid point needed
 */
int g_maxof4neighbor(struct Grid *g, const double x, const double y, double *zvalue);

/* ---------------------------------------------------------------------------------
 * GRID : interpolation for point with coordinates (x,y) as weighted average over the
 *        4 corner points (weight = 1/sqr(distance))
 *
 *        g_inter_simple : new version with additional parameter
 *             minpoint = minimal number of corners needed (ranging from 1 - 4)
 */
int g_inter_simple(struct Grid *g, const double x, const double y, int minpoint, double *zvalue);

/*---------------------------------------------------------------------------------
 *  GRID : bilinear interpolation for point with coordinates (x,y)
 */ 
int g_inter_bilinear (struct Grid *g, const double x, const double y, double *zvalue);

/*---------------------------------------------------------------------------------
 *  GMT : returns pointer to start of next substring and length in *n
 *       (may be moved to inout of libcommon)
 */
char * next_substring(char *line, int *n);

/*---------------------------------------------------------------------------------
 *  GMT : decodes a R-option string  (example : -R0/360/-90/90)
 */
int R_option( char * s, double *xmin, double *xmax, double *ymin, double *ymax);


/* ---------------------------------------------------------------------------------
 *  GMT : decodes an I-option string  (example : -I15m)
 */
int I_option( char * s, double *dx, double *dy);

/* ---------------------------------------------------------------------------------
 *  GMT : get info from grdraster.info file
 */
struct Grid * gmt_rasterinfo(int grdnr);

/* ---------------------------------------------------------------------------------
 *  GMT : open rasterfile for reading
 */
int gmt_openrasterfile(struct Grid *g);

/* ---------------------------------------------------------------------------------
 *  GMT : close rasterfile
 */
void gmt_closerasterfile(struct Grid *g);

/* ---------------------------------------------------------------------------------
 *  GMT : get memory for data and read whole raster data from file
 */
int gmt_readgrid(struct Grid *g);

/* ---------------------------------------------------------------------------------
 *  GMT : initialize subgrid and read data from file
 */
struct Grid * gmt_readsubgrid(struct Grid *oldgrid, 
			      double xmin, double xmax, double ymin, double ymax );

/* ---------------------------------------------------------------------------------
 *  SURFER : read surfer grid file ASCII (grid is initialized in routine)
 */
struct Grid *surf_readgrid(const char *grdfname);

/* ---------------------------------------------------------------------------------
 *  SURFER : write surfer grid to file ASCII (format string of 1 number i.e."10.3lf")
 */
int surf_writegrid(const char *grdfname, struct Grid *g, const char *format);

/***************************************************************************/
/*  GMT: write GMT grid to file (uses xyz2grd of GMT)
 *       unit: will be set for all three elements in -D option of xyz2grd
 */
int gmt_writegrid (const char *grdfn, struct Grid *g,
		   const char *unit, const char *comment, const int verbose);

/***************************************************************************/
/*   write grid to ascii file
 *   format entire line with format string. e.g. "%10.3f %10.3f %10.3f\n"
 */
int g_writeasciigrid (const char *grdfn, struct Grid *g, const char *format);

#endif
/*
 *      FILE: inout
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *            16-09-1997  --  added parsing functions
 *            
 *
 */
#ifndef GGL_INOUT
#define GGL_INOUT


/* 
 * function for scanning 1 line from a data file 
 * (used in read_colfile)
 */
typedef void * (*readlinefunction)(char *, void *);

/*
 *  Struct to hold char columns
 */
struct coordcol {
        int count;
        char *col[7];
} ;

struct scoordcol {
	int count;
	int slen;
	char **col;
} ;

int read_comment (FILE *fpin, const int idchar)
;

/*
 *  Overread comment lines, assuming next char read is first of a line.
 *  RETURN VALUE: zero on success, EOF upon end of file.
 */
int overread_comment (FILE *fpin,	/* input stream */
		      FILE *fpout,	/* output: if NULL, lines are discarded */
		      const int commentchar);

char * read_option (FILE *fpin, const int idchar);

char * read_soption (FILE *fpin, const int idchar);

/*
 *  The following functions will have to be reviewed when porting
 *  the library to another plattform or using another compiler 
 *  than the GNU C Compiler (gcc)!
 */
#if defined(OS_LINUX) || defined(OS_AIX) || defined(OS_CYGWIN)
int read_delim (FILE *fpin, int delimiter, char **lineptr);

int read_line (FILE *fpin, char **lineptr);
int read_line_static (FILE *fpin, char *s, const int size);

/* 
 * read column file 
 */
void ** read_colfile (FILE *fp,       int *count, void *param, readlinefunction readline);
void ** read_colfile2(char *filename, int *count, void *param, readlinefunction readline);


/* for use with returned int of %n on sscanf */
int strpos_withblanks(const char *s, const int count);

int erase_file (const char *filename);

/*
 *  uses the access function from unistd.h
 */
int check_file_name (const char *filename,  int mode);
FILE * open_text_file (const char *fname, const char *mode);
FILE * open_binary_file (const char *fname, const char *mode);

/*
 * check, whether the file fname is a zipped file by looking at its magic
 * - if the file is zipped, the function closes it, 
 *   and a pipe on the unzipped file is passed back.
 * - if the file is unzipped, the filepointer is left untouched
 * ispipe: TRUE is zipped, FALSE if unzipped, EOF on error
 */
FILE * openzipped (const char *fname, int *ispipe);

#endif /* OS_LINUX */ 

/*
 * returns TRUE or FALSE
 * lines starting with # and blank lines are considered as comments
 */
int iscomment(char *line);

/*
 * skipping col values in string, returning offset 
 */
int col_offset(char * s, int col);


/*
 *  this function modifies the string passed to it!
 *  the extension including its leading point will
 *  be removed (in fact: the point is replaced by a
 *  null character)
 */
void remove_extension (char *filename);

/*
 *  this function returns a new allocated string
 *  containing the extension including the leading point
 */
char * get_extension (const char *filename);

/*
 *  this function returns a new allocated string
 *  containing the filename with a new extension
 */
char * add_extension (
	const char *oldname,
	const char *ext );

/*
 *  this function returns a new allocated string
 *  containing the filename excluding the leading path
 */
char * get_fname_wopath (const char *filename);

/*
 *  Beware, this function will produce a segmentation fault
 *  if the string is not null terminated
 */
void uppercase (char *s);

/*
 *  Beware, this function will produce a segmentation fault
 *  if the string is not null terminated
 */
void lowercase (char *s);

/*
 *  Beware, this function will produce a segmentation fault
 *  if the string is not null terminated
 */
void replacechar (char *s, int oldc, int newc);

/*
 *  extract columns from an input line
 */
int extract_values (
	const char *line,
	char **prefix,
	const unsigned int prefix_size,
	struct coordcol *infix,
	char **suffix,
	const unsigned long linecount);

/*
 *  extract columns from an input line
 */
int extract_values_static (
	const char *line,
	char prefix[],
	const int prefix_size,
	struct scoordcol *infix,
	char suffix[],
	const int max_suffix_size,
	const unsigned long linecount);

/*
 *  strip directory and suffix from filenames
 *  returns a string 
 */ 
char * base_name(char *filename, char *ext);



/*
************************************************************
* parse.c
************************************************************
*/
#define ASCII_LINE 4
#define ASCII_BLOCK 5
#define BINARY_BLOCK 1

/* byte of type meaning */

#define BLOCK 1
/* unused 2 */
#define ASCII 4
/* unused 8 */

#define SKIP_NONE 0
#define SKIP_CHAR 1
#define SKIP_FIELD 2

int parse_line (
  FILE *fpin,
  const char *comment_delimiters,
  const char *column_delimiters,
  const int skip,
  const int skip_type,
  char **irrelevant,
  char **datastring)
;

char * get_fileformat (
  FILE *fpin,
  const int type,
  const char *comment_delimiters,
  const char *column_delimiters)
;

#endif


/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sunras.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: I/O of SUN raster graphics files
 *
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: sunras.h,v 1.4 2001/02/21 14:00:19 cocard Exp $
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
 *   AUTHORS: Marc Troller (VTB Sommer 1999)
 *            Favey Etienne <favey@geod.ethz.ch>
 *
 *   HISTORY: 
 * $Log: sunras.h,v $
 * Revision 1.4  2001/02/21 14:00:19  cocard
 * *** empty log message ***
 *
 * Revision 1.3  2000/03/06 17:01:42  favey
 * GMT colortable and grid output
 *
 * Revision 1.2  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.1  1999/06/30 15:34:56  favey
 * Initial revision
 *
 */
struct rasterfile {
  unsigned int ras_magic;              /* magic number */
  int     ras_width;              /* width (pixels) of image */
  int     ras_height;             /* height (pixels) of image */
  int     ras_depth;              /* depth (1, 8, or 24 bits) of pixel */
  int     ras_length;             /* length (bytes) of image */
  int     ras_type;               /* type of file; see RT_* below */
  int     ras_maptype;            /* type of colormap; see RMT_* below */
  int     ras_maplength;          /* length (bytes) of following map */
  /* color map follows for ras_maplength bytes, followed by image */
};

#define RT_BYTE_ENCODED 2       /* Run-length compression of bytes */
#define RMT_NONE        0       /* ras_maplength is expected to be 0 */

int read_sunras_header (FILE *fpras, struct rasterfile *header);
int write_sunras_header (FILE *fpras, struct rasterfile *header, const int inv);

int read_colormap (FILE *fpras, unsigned int colmapsize,  unsigned char *colmap);
int write_colormap (FILE *fpras,unsigned int colmapsize,  unsigned char *colmap);
int write_cptfile (FILE *fpcpt, unsigned int colmapsize, unsigned char *colmap);

void write_sunras_info ( struct rasterfile rh);
/*
 *      FILE: intersect
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: routines for geometrical intersections
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            15-05-1998  --  created
 *
 */

#ifndef INTERSECT
#define INTERSECT

/*
 * calculates the coordinates of a point on a line
 */
void line_coord
(const double x[3],           /* input:  coordinates of point */
 const double e[3],           /* input:  unary direction vector */
 const double dist,           /* input:  distance to new point */
 double xnew[3])              /* output: coordinates of new point */
;

/*
 * calculates the intersection of a line and an ellipsoid 
 * defined by its three axis
 * return value = #count of intersections (0,1, or 2)
 */
int intersection_line_ell3
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double axis[3],        /* input:  all 3 major axis  */ 
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
;

/*
 * calculates the intersection of a line and a ellipoid of rotation
 * return value = #count of intersections (0,1, or 2)
 * derivated from intersection_line_ell3
*/
int intersection_line_ell
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double major_axis,     /* input:  major axis */ 
 const double ex2,            /* input:  square value of excentricity */
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
;

/*
 * calculates the intersection of a line and a sphere
 * return value = #count of intersections (0,1, or 2)
 * derivated from intersection_line_ell3
 */
int intersection_line_sphere
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double radius,         /* input:  radius of sphere */ 
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
;

/*
 * calculates the intersection of a line with a surface with
 * constant longitude (plane surface)
 * return value = #count of intersections (-1, 0 or 1)
 */
int intersection_line_longitude
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double longitude,      /* input:  longitude */
 double *dist)                /* output: distance to intersection */
;

/*
 * calculates the intersection of a line with a surface with
 * constant elliposidal latitude (conical surface).
 * return value = #count of intersections (-1, 0, 1, or 2)
 */
int intersection_line_latitude
(const double x[3],           /* input:  point  coordinates */
 const double e[3],           /* input:  unary direction vector */
 const double latitude,       /* input:  latitude */
 const double major_axis,     /* input:  major axis */ 
 const double ex2,            /* input:  excentricity of ellipsoid squared */
 double *dist1,               /* output: distance from point to first intersection */
 double *dist2)               /* output: distance from point to second intersection */
;

/*
 * 2D routine :
 * This function calculates the intersection of two circles given
 * by the coordinates  of the center and the radius. The resulting
 * intersection points are returned in res1 and res2. The return
 * value of the function specifies the number of valid inter-
 * sections (-1, 0, 1, 2) where -1 means circles are identical /
 * 0 = no intersection / 1 = one tangential point / 2 = 2 inter-
 * section points.
 */
int intersection_circles (
  const double center1[2],  /* input:  coordinates of center of the first circle */
  const double radius1,     /* input:  radius of the first circle */
  const double center2[2],  /* input:  coordinates of center of the second circle */
  const double radius2,     /* input:  radius of the second circle */
  double res1[2],           /* output: coordinates of intersection 1 */
  double res2[2])           /* output: coordinates of intersection 2 */
;

int intersection_line_segment(
  double origin[2], 
  double azi, 
  double p1[2], 
  double p2[2], 
  double ipoint[2])
;

#endif
/*  DO NOT EDIT THIS FILE, EDIT sed_matcalc INSTEAD!!! */
/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sed_matcalc.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: sed_matcalc.h,v 1.8 2000/05/12 09:22:53 favey Exp $
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
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch> (pascal version)
 *            Favey Etienne   <favey@geod.ethz.ch>  (c port)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *                  1987  --  created
 *            19-03-1997  --  pascal version to port
 *            07-04-1997  --  c port
 *            11-12-1997  --  procedure m_write debugged  
 *                            procedure m_writefmt added
 *                            procedure m_writefilled modified and completed (MC)
 *            01-07-1998  --  m_separate_sym added (MC)
 *            02-11-1998  --  introduced SED functionality
 *
 * $Log: sed_matcalc.h,v $
 * Revision 1.8  2000/05/12 09:22:53  favey
 * changed devel version to DEVEL_MAT_DEBUG
 *
 * Revision 1.7  2000/04/17 15:54:16  cocard
 * *** empty log message ***
 *
 * Revision 1.6  2000/02/25 12:57:58  ypeter
 * m_polynom_approx: modified function call, added possible weightedobs.
 * m_polynom_interpol: moved from winterpol to matcalc
 *
 * Revision 1.5  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.4  1999/07/29 13:22:28  cocard
 * $new function m_resize added (MC)
 *
 * Revision 1.3  1999/04/19 09:51:06  cocard
 * *** empty log message ***
 *
 * Revision 1.2  1999/04/07 11:45:32  cocard
 * read/write binary debugged
 *
 * Revision 1.1  1998/11/03 18:00:51  favey
 * Initial revision
 *
 *            
 *   CAVEATS: contains error messages exiting main()
 *            when compiled with DEVEL_MAT_DEBUG defined
 *            returning EXIT_FAILURE
 *
 */

#ifndef GGL_iMATCALC
#define GGL_iMATCALC

struct int_matrix {
	char *name;
	long int **mat;	/* normal access via mat[][] */
	long int *fast;	/* allows fast access to matrix */
	int row;
	int col;
	int dim_row;
	int dim_col;
};


/* writes some error on stderr */
void im_error (
	const char *procname,
	struct int_matrix *m1, 
	struct int_matrix *m2,
	struct int_matrix *m3,
	int *row, 
	int *col);

/* initializes a matrix structure */
struct int_matrix * im_init (int row, int col);

/* initializes the name of a matrix */
int im_name (struct int_matrix *m, const char *newname);

/* initializes a matrix structure with name */
struct int_matrix * im_init_name (int row, int col, char *newname);

/* frees the structure matrix */
void im_free (struct int_matrix *m);


/* resizes the matrix without new memory allocation 
 * sets all values to zero if set_zero == TRUE
 */
int im_resize (struct int_matrix *m, 
		   int row, int col, int set_zero);


/* fill one element at position (row, col) with value */
int im_fill (struct int_matrix *m,
		 int row, int col, long int value);

/* fill all elements of row  with value */
int im_fill_row (struct int_matrix *m,
		     int row, long int value);

/* fill all elements of column col  with value */
int im_fill_col (struct int_matrix *m,
		     int col, long int value);

/* copies upper part to lower part */
int im_fillsym (struct int_matrix *m);

/* fills matrix completely with zero */
int im_zero (struct int_matrix *m);

/* fills one row of matrix with zero */
int im_zero_row (struct int_matrix *m,	int row);

/* makes the matrix to a unit matrix */
int im_unit (struct int_matrix *m);

/* fills all elements of matrix with value */
int im_set (struct int_matrix *m, long int value);

/* fills all elements with random values in the range 0 - fac */
int im_random (struct int_matrix *m, long int fac);

/* deep copy from source to target */
int im_copy (struct int_matrix *source,
		 struct int_matrix *target);

/* ????????? */
int im_pushup (struct int_matrix *m);
int im_pushdown (struct int_matrix *m);


int im_max_col (
	struct int_matrix *m,
	const int col,
	long int *max,
	int *row)
;
int im_min_col (
	struct int_matrix *m,
	const int col,
	long int *min,
	int *row)
;
int im_max (
	struct int_matrix *m,
	long int *max,
	int *row,
	int *col)
;
int im_min (
	struct int_matrix *m,
	long int *min,
	int *row,
	int *col)
;
int im_absmax (
	struct int_matrix *m,
	long int *max,
	int *row,
	int *col)
;
int im_write (
	FILE *fpout,
	struct int_matrix *m,
	const int fmt1,
	const int fmt2,
	int n,
	const int tr,
	const int title)
;
int im_writefmt (
	FILE *fpout,
	struct int_matrix *m,
	const char * fmt,
	int n,
	const int tr,
	const int title)
;
int im_write_filled (
	FILE *fpout,
	struct int_matrix *m,
	const int tr,
	long int nearlyzero)
;
int im_writebin (
	FILE *fpout,
	struct int_matrix *m)
;
int im_writebin2 (
	char * fname, 
	struct int_matrix *m)
;
int im_readbin (
	FILE *fpin,
	struct int_matrix **m)
;
int im_readbin2 (
	char * fname, 
	struct int_matrix **m)
;
int im_readascii(
	char *fname, 
	struct int_matrix **m, 
	int row2read,
	int col2read)
;

int im_add (
	struct int_matrix *m1,
	struct int_matrix *m2,
	struct int_matrix *m_res)
;
int im_sub (
	struct int_matrix *m1,
	struct int_matrix *m2,
	struct int_matrix *m_res)
;
int im_mult (
	struct int_matrix *m1,
	const int m1tr,
	struct int_matrix *m2,
	const int m2tr,
	struct int_matrix *m_res)
;
int im_multx (
	struct int_matrix *m,
	long int scalar,
	struct int_matrix *m_res)
;
long int im_trace (struct int_matrix *m)
;
int im_trans (
	struct int_matrix *m,
	struct int_matrix *m_res)
;
int im_atpa1 (
	struct int_matrix *a,
	struct int_matrix *p,
	struct int_matrix *atpa)
;
int im_det (
        struct int_matrix *m, 
        long int *det)
;

int im_pivot (
	struct int_matrix *m,
	int row,
	int col)
;
int im_inv (struct int_matrix *m)
;
int im_cholesky_decomposition (struct int_matrix *m)
;
int im_choleskyinv (struct int_matrix *m)
;
long int im_sym_check (struct int_matrix *m)
;
int im_add_block (
	struct int_matrix *m,
	struct int_matrix *block,
	int row,
	int col)
;
int im_pick_block (
	struct int_matrix *m,
	struct int_matrix *block,
	int row,
	int col,
	int drow,
	int dcol)
;
int im_del_row (
	struct int_matrix *m,
	int row)
;
int im_del_col (
	struct int_matrix *m,
	int col)
;
int im_polynom_approx (
	struct int_matrix *data,
	int col_x,
	int col_y,
        int col_w,
	const int col_app,
	const int col_res,
	const int row_beg,
	const int row_end,
	const int degree,
	long int x0,
	long int *dummy,
	struct int_matrix *coeff,
	struct int_matrix *qxx,
	long int *s0)
;

long int im_polynom_interpol (
	struct int_matrix *coeff,
        int num_coeff,
	const long int dx)
;

int im_separate_sym(
	struct int_matrix *m, 
	struct int_matrix *sym, 
	struct int_matrix *asym)
;
#endif
/*  DO NOT EDIT THIS FILE, EDIT sed_matcalc INSTEAD!!! */
/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sed_matcalc.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: sed_matcalc.h,v 1.8 2000/05/12 09:22:53 favey Exp $
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
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch> (pascal version)
 *            Favey Etienne   <favey@geod.ethz.ch>  (c port)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *                  1987  --  created
 *            19-03-1997  --  pascal version to port
 *            07-04-1997  --  c port
 *            11-12-1997  --  procedure m_write debugged  
 *                            procedure m_writefmt added
 *                            procedure m_writefilled modified and completed (MC)
 *            01-07-1998  --  m_separate_sym added (MC)
 *            02-11-1998  --  introduced SED functionality
 *
 * $Log: sed_matcalc.h,v $
 * Revision 1.8  2000/05/12 09:22:53  favey
 * changed devel version to DEVEL_MAT_DEBUG
 *
 * Revision 1.7  2000/04/17 15:54:16  cocard
 * *** empty log message ***
 *
 * Revision 1.6  2000/02/25 12:57:58  ypeter
 * m_polynom_approx: modified function call, added possible weightedobs.
 * m_polynom_interpol: moved from winterpol to matcalc
 *
 * Revision 1.5  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.4  1999/07/29 13:22:28  cocard
 * $new function m_resize added (MC)
 *
 * Revision 1.3  1999/04/19 09:51:06  cocard
 * *** empty log message ***
 *
 * Revision 1.2  1999/04/07 11:45:32  cocard
 * read/write binary debugged
 *
 * Revision 1.1  1998/11/03 18:00:51  favey
 * Initial revision
 *
 *            
 *   CAVEATS: contains error messages exiting main()
 *            when compiled with DEVEL_MAT_DEBUG defined
 *            returning EXIT_FAILURE
 *
 */

#ifndef GGL_sMATCALC
#define GGL_sMATCALC

struct sint_matrix {
	char *name;
	short int **mat;	/* normal access via mat[][] */
	short int *fast;	/* allows fast access to matrix */
	int row;
	int col;
	int dim_row;
	int dim_col;
};


/* writes some error on stderr */
void sm_error (
	const char *procname,
	struct sint_matrix *m1, 
	struct sint_matrix *m2,
	struct sint_matrix *m3,
	int *row, 
	int *col);

/* initializes a matrix structure */
struct sint_matrix * sm_init (int row, int col);

/* initializes the name of a matrix */
int sm_name (struct sint_matrix *m, const char *newname);

/* initializes a matrix structure with name */
struct sint_matrix * sm_init_name (int row, int col, char *newname);

/* frees the structure matrix */
void sm_free (struct sint_matrix *m);


/* resizes the matrix without new memory allocation 
 * sets all values to zero if set_zero == TRUE
 */
int sm_resize (struct sint_matrix *m, 
		   int row, int col, int set_zero);


/* fill one element at position (row, col) with value */
int sm_fill (struct sint_matrix *m,
		 int row, int col, short int value);

/* fill all elements of row  with value */
int sm_fill_row (struct sint_matrix *m,
		     int row, short int value);

/* fill all elements of column col  with value */
int sm_fill_col (struct sint_matrix *m,
		     int col, short int value);

/* copies upper part to lower part */
int sm_fillsym (struct sint_matrix *m);

/* fills matrix completely with zero */
int sm_zero (struct sint_matrix *m);

/* fills one row of matrix with zero */
int sm_zero_row (struct sint_matrix *m,	int row);

/* makes the matrix to a unit matrix */
int sm_unit (struct sint_matrix *m);

/* fills all elements of matrix with value */
int sm_set (struct sint_matrix *m, short int value);

/* fills all elements with random values in the range 0 - fac */
int sm_random (struct sint_matrix *m, short int fac);

/* deep copy from source to target */
int sm_copy (struct sint_matrix *source,
		 struct sint_matrix *target);

/* ????????? */
int sm_pushup (struct sint_matrix *m);
int sm_pushdown (struct sint_matrix *m);


int sm_max_col (
	struct sint_matrix *m,
	const int col,
	short int *max,
	int *row)
;
int sm_min_col (
	struct sint_matrix *m,
	const int col,
	short int *min,
	int *row)
;
int sm_max (
	struct sint_matrix *m,
	short int *max,
	int *row,
	int *col)
;
int sm_min (
	struct sint_matrix *m,
	short int *min,
	int *row,
	int *col)
;
int sm_absmax (
	struct sint_matrix *m,
	short int *max,
	int *row,
	int *col)
;
int sm_write (
	FILE *fpout,
	struct sint_matrix *m,
	const int fmt1,
	const int fmt2,
	int n,
	const int tr,
	const int title)
;
int sm_writefmt (
	FILE *fpout,
	struct sint_matrix *m,
	const char * fmt,
	int n,
	const int tr,
	const int title)
;
int sm_write_filled (
	FILE *fpout,
	struct sint_matrix *m,
	const int tr,
	short int nearlyzero)
;
int sm_writebin (
	FILE *fpout,
	struct sint_matrix *m)
;
int sm_writebin2 (
	char * fname, 
	struct sint_matrix *m)
;
int sm_readbin (
	FILE *fpin,
	struct sint_matrix **m)
;
int sm_readbin2 (
	char * fname, 
	struct sint_matrix **m)
;
int sm_readascii(
	char *fname, 
	struct sint_matrix **m, 
	int row2read,
	int col2read)
;

int sm_add (
	struct sint_matrix *m1,
	struct sint_matrix *m2,
	struct sint_matrix *m_res)
;
int sm_sub (
	struct sint_matrix *m1,
	struct sint_matrix *m2,
	struct sint_matrix *m_res)
;
int sm_mult (
	struct sint_matrix *m1,
	const int m1tr,
	struct sint_matrix *m2,
	const int m2tr,
	struct sint_matrix *m_res)
;
int sm_multx (
	struct sint_matrix *m,
	short int scalar,
	struct sint_matrix *m_res)
;
short int sm_trace (struct sint_matrix *m)
;
int sm_trans (
	struct sint_matrix *m,
	struct sint_matrix *m_res)
;
int sm_atpa1 (
	struct sint_matrix *a,
	struct sint_matrix *p,
	struct sint_matrix *atpa)
;
int sm_det (
        struct sint_matrix *m, 
        short int *det)
;

int sm_pivot (
	struct sint_matrix *m,
	int row,
	int col)
;
int sm_inv (struct sint_matrix *m)
;
int sm_cholesky_decomposition (struct sint_matrix *m)
;
int sm_choleskyinv (struct sint_matrix *m)
;
short int sm_sym_check (struct sint_matrix *m)
;
int sm_add_block (
	struct sint_matrix *m,
	struct sint_matrix *block,
	int row,
	int col)
;
int sm_pick_block (
	struct sint_matrix *m,
	struct sint_matrix *block,
	int row,
	int col,
	int drow,
	int dcol)
;
int sm_del_row (
	struct sint_matrix *m,
	int row)
;
int sm_del_col (
	struct sint_matrix *m,
	int col)
;
int sm_polynom_approx (
	struct sint_matrix *data,
	int col_x,
	int col_y,
        int col_w,
	const int col_app,
	const int col_res,
	const int row_beg,
	const int row_end,
	const int degree,
	short int x0,
	short int *dummy,
	struct sint_matrix *coeff,
	struct sint_matrix *qxx,
	short int *s0)
;

short int sm_polynom_interpol (
	struct sint_matrix *coeff,
        int num_coeff,
	const short int dx)
;

int sm_separate_sym(
	struct sint_matrix *m, 
	struct sint_matrix *sym, 
	struct sint_matrix *asym)
;
#endif
/*  DO NOT EDIT THIS FILE, EDIT sed_matcalc INSTEAD!!! */
/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sed_matcalc.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: sed_matcalc.h,v 1.8 2000/05/12 09:22:53 favey Exp $
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
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch> (pascal version)
 *            Favey Etienne   <favey@geod.ethz.ch>  (c port)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *                  1987  --  created
 *            19-03-1997  --  pascal version to port
 *            07-04-1997  --  c port
 *            11-12-1997  --  procedure m_write debugged  
 *                            procedure m_writefmt added
 *                            procedure m_writefilled modified and completed (MC)
 *            01-07-1998  --  m_separate_sym added (MC)
 *            02-11-1998  --  introduced SED functionality
 *
 * $Log: sed_matcalc.h,v $
 * Revision 1.8  2000/05/12 09:22:53  favey
 * changed devel version to DEVEL_MAT_DEBUG
 *
 * Revision 1.7  2000/04/17 15:54:16  cocard
 * *** empty log message ***
 *
 * Revision 1.6  2000/02/25 12:57:58  ypeter
 * m_polynom_approx: modified function call, added possible weightedobs.
 * m_polynom_interpol: moved from winterpol to matcalc
 *
 * Revision 1.5  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.4  1999/07/29 13:22:28  cocard
 * $new function m_resize added (MC)
 *
 * Revision 1.3  1999/04/19 09:51:06  cocard
 * *** empty log message ***
 *
 * Revision 1.2  1999/04/07 11:45:32  cocard
 * read/write binary debugged
 *
 * Revision 1.1  1998/11/03 18:00:51  favey
 * Initial revision
 *
 *            
 *   CAVEATS: contains error messages exiting main()
 *            when compiled with DEVEL_MAT_DEBUG defined
 *            returning EXIT_FAILURE
 *
 */

#ifndef GGL_lMATCALC
#define GGL_lMATCALC

struct long_matrix {
	char *name;
	long double **mat;	/* normal access via mat[][] */
	long double *fast;	/* allows fast access to matrix */
	int row;
	int col;
	int dim_row;
	int dim_col;
};


/* writes some error on stderr */
void lm_error (
	const char *procname,
	struct long_matrix *m1, 
	struct long_matrix *m2,
	struct long_matrix *m3,
	int *row, 
	int *col);

/* initializes a matrix structure */
struct long_matrix * lm_init (int row, int col);

/* initializes the name of a matrix */
int lm_name (struct long_matrix *m, const char *newname);

/* initializes a matrix structure with name */
struct long_matrix * lm_init_name (int row, int col, char *newname);

/* frees the structure matrix */
void lm_free (struct long_matrix *m);


/* resizes the matrix without new memory allocation 
 * sets all values to zero if set_zero == TRUE
 */
int lm_resize (struct long_matrix *m, 
		   int row, int col, int set_zero);


/* fill one element at position (row, col) with value */
int lm_fill (struct long_matrix *m,
		 int row, int col, long double value);

/* fill all elements of row  with value */
int lm_fill_row (struct long_matrix *m,
		     int row, long double value);

/* fill all elements of column col  with value */
int lm_fill_col (struct long_matrix *m,
		     int col, long double value);

/* copies upper part to lower part */
int lm_fillsym (struct long_matrix *m);

/* fills matrix completely with zero */
int lm_zero (struct long_matrix *m);

/* fills one row of matrix with zero */
int lm_zero_row (struct long_matrix *m,	int row);

/* makes the matrix to a unit matrix */
int lm_unit (struct long_matrix *m);

/* fills all elements of matrix with value */
int lm_set (struct long_matrix *m, long double value);

/* fills all elements with random values in the range 0 - fac */
int lm_random (struct long_matrix *m, long double fac);

/* deep copy from source to target */
int lm_copy (struct long_matrix *source,
		 struct long_matrix *target);

/* ????????? */
int lm_pushup (struct long_matrix *m);
int lm_pushdown (struct long_matrix *m);


int lm_max_col (
	struct long_matrix *m,
	const int col,
	long double *max,
	int *row)
;
int lm_min_col (
	struct long_matrix *m,
	const int col,
	long double *min,
	int *row)
;
int lm_max (
	struct long_matrix *m,
	long double *max,
	int *row,
	int *col)
;
int lm_min (
	struct long_matrix *m,
	long double *min,
	int *row,
	int *col)
;
int lm_absmax (
	struct long_matrix *m,
	long double *max,
	int *row,
	int *col)
;
int lm_write (
	FILE *fpout,
	struct long_matrix *m,
	const int fmt1,
	const int fmt2,
	int n,
	const int tr,
	const int title)
;
int lm_writefmt (
	FILE *fpout,
	struct long_matrix *m,
	const char * fmt,
	int n,
	const int tr,
	const int title)
;
int lm_write_filled (
	FILE *fpout,
	struct long_matrix *m,
	const int tr,
	long double nearlyzero)
;
int lm_writebin (
	FILE *fpout,
	struct long_matrix *m)
;
int lm_writebin2 (
	char * fname, 
	struct long_matrix *m)
;
int lm_readbin (
	FILE *fpin,
	struct long_matrix **m)
;
int lm_readbin2 (
	char * fname, 
	struct long_matrix **m)
;
int lm_readascii(
	char *fname, 
	struct long_matrix **m, 
	int row2read,
	int col2read)
;

int lm_add (
	struct long_matrix *m1,
	struct long_matrix *m2,
	struct long_matrix *m_res)
;
int lm_sub (
	struct long_matrix *m1,
	struct long_matrix *m2,
	struct long_matrix *m_res)
;
int lm_mult (
	struct long_matrix *m1,
	const int m1tr,
	struct long_matrix *m2,
	const int m2tr,
	struct long_matrix *m_res)
;
int lm_multx (
	struct long_matrix *m,
	long double scalar,
	struct long_matrix *m_res)
;
long double lm_trace (struct long_matrix *m)
;
int lm_trans (
	struct long_matrix *m,
	struct long_matrix *m_res)
;
int lm_atpa1 (
	struct long_matrix *a,
	struct long_matrix *p,
	struct long_matrix *atpa)
;
int lm_det (
        struct long_matrix *m, 
        long double *det)
;

int lm_pivot (
	struct long_matrix *m,
	int row,
	int col)
;
int lm_inv (struct long_matrix *m)
;
int lm_cholesky_decomposition (struct long_matrix *m)
;
int lm_choleskyinv (struct long_matrix *m)
;
long double lm_sym_check (struct long_matrix *m)
;
int lm_add_block (
	struct long_matrix *m,
	struct long_matrix *block,
	int row,
	int col)
;
int lm_pick_block (
	struct long_matrix *m,
	struct long_matrix *block,
	int row,
	int col,
	int drow,
	int dcol)
;
int lm_del_row (
	struct long_matrix *m,
	int row)
;
int lm_del_col (
	struct long_matrix *m,
	int col)
;
int lm_polynom_approx (
	struct long_matrix *data,
	int col_x,
	int col_y,
        int col_w,
	const int col_app,
	const int col_res,
	const int row_beg,
	const int row_end,
	const int degree,
	long double x0,
	long double *dummy,
	struct long_matrix *coeff,
	struct long_matrix *qxx,
	long double *s0)
;

long double lm_polynom_interpol (
	struct long_matrix *coeff,
        int num_coeff,
	const long double dx)
;

int lm_separate_sym(
	struct long_matrix *m, 
	struct long_matrix *sym, 
	struct long_matrix *asym)
;
#endif
/*  DO NOT EDIT THIS FILE, EDIT sed_matcalc INSTEAD!!! */
/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sed_matcalc.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: sed_matcalc.h,v 1.8 2000/05/12 09:22:53 favey Exp $
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
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch> (pascal version)
 *            Favey Etienne   <favey@geod.ethz.ch>  (c port)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *                  1987  --  created
 *            19-03-1997  --  pascal version to port
 *            07-04-1997  --  c port
 *            11-12-1997  --  procedure m_write debugged  
 *                            procedure m_writefmt added
 *                            procedure m_writefilled modified and completed (MC)
 *            01-07-1998  --  m_separate_sym added (MC)
 *            02-11-1998  --  introduced SED functionality
 *
 * $Log: sed_matcalc.h,v $
 * Revision 1.8  2000/05/12 09:22:53  favey
 * changed devel version to DEVEL_MAT_DEBUG
 *
 * Revision 1.7  2000/04/17 15:54:16  cocard
 * *** empty log message ***
 *
 * Revision 1.6  2000/02/25 12:57:58  ypeter
 * m_polynom_approx: modified function call, added possible weightedobs.
 * m_polynom_interpol: moved from winterpol to matcalc
 *
 * Revision 1.5  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.4  1999/07/29 13:22:28  cocard
 * $new function m_resize added (MC)
 *
 * Revision 1.3  1999/04/19 09:51:06  cocard
 * *** empty log message ***
 *
 * Revision 1.2  1999/04/07 11:45:32  cocard
 * read/write binary debugged
 *
 * Revision 1.1  1998/11/03 18:00:51  favey
 * Initial revision
 *
 *            
 *   CAVEATS: contains error messages exiting main()
 *            when compiled with DEVEL_MAT_DEBUG defined
 *            returning EXIT_FAILURE
 *
 */

#ifndef GGL_MATCALC
#define GGL_MATCALC

struct matrix {
	char *name;
	double **mat;	/* normal access via mat[][] */
	double *fast;	/* allows fast access to matrix */
	int row;
	int col;
	int dim_row;
	int dim_col;
};


/* writes some error on stderr */
void m_error (
	const char *procname,
	struct matrix *m1, 
	struct matrix *m2,
	struct matrix *m3,
	int *row, 
	int *col);

/* initializes a matrix structure */
struct matrix * m_init (int row, int col);

/* initializes the name of a matrix */
int m_name (struct matrix *m, const char *newname);

/* initializes a matrix structure with name */
struct matrix * m_init_name (int row, int col, char *newname);

/* frees the structure matrix */
void m_free (struct matrix *m);


/* resizes the matrix without new memory allocation 
 * sets all values to zero if set_zero == TRUE
 */
int m_resize (struct matrix *m, 
		   int row, int col, int set_zero);


/* fill one element at position (row, col) with value */
int m_fill (struct matrix *m,
		 int row, int col, double value);

/* fill all elements of row  with value */
int m_fill_row (struct matrix *m,
		     int row, double value);

/* fill all elements of column col  with value */
int m_fill_col (struct matrix *m,
		     int col, double value);

/* copies upper part to lower part */
int m_fillsym (struct matrix *m);

/* fills matrix completely with zero */
int m_zero (struct matrix *m);

/* fills one row of matrix with zero */
int m_zero_row (struct matrix *m,	int row);

/* makes the matrix to a unit matrix */
int m_unit (struct matrix *m);

/* fills all elements of matrix with value */
int m_set (struct matrix *m, double value);

/* fills all elements with random values in the range 0 - fac */
int m_random (struct matrix *m, double fac);

/* deep copy from source to target */
int m_copy (struct matrix *source,
		 struct matrix *target);

/* ????????? */
int m_pushup (struct matrix *m);
int m_pushdown (struct matrix *m);


int m_max_col (
	struct matrix *m,
	const int col,
	double *max,
	int *row)
;
int m_min_col (
	struct matrix *m,
	const int col,
	double *min,
	int *row)
;
int m_max (
	struct matrix *m,
	double *max,
	int *row,
	int *col)
;
int m_min (
	struct matrix *m,
	double *min,
	int *row,
	int *col)
;
int m_absmax (
	struct matrix *m,
	double *max,
	int *row,
	int *col)
;
int m_write (
	FILE *fpout,
	struct matrix *m,
	const int fmt1,
	const int fmt2,
	int n,
	const int tr,
	const int title)
;
int m_writefmt (
	FILE *fpout,
	struct matrix *m,
	const char * fmt,
	int n,
	const int tr,
	const int title)
;
int m_write_filled (
	FILE *fpout,
	struct matrix *m,
	const int tr,
	double nearlyzero)
;
int m_writebin (
	FILE *fpout,
	struct matrix *m)
;
int m_writebin2 (
	char * fname, 
	struct matrix *m)
;
int m_readbin (
	FILE *fpin,
	struct matrix **m)
;
int m_readbin2 (
	char * fname, 
	struct matrix **m)
;
int m_readascii(
	char *fname, 
	struct matrix **m, 
	int row2read,
	int col2read)
;

int m_add (
	struct matrix *m1,
	struct matrix *m2,
	struct matrix *m_res)
;
int m_sub (
	struct matrix *m1,
	struct matrix *m2,
	struct matrix *m_res)
;
int m_mult (
	struct matrix *m1,
	const int m1tr,
	struct matrix *m2,
	const int m2tr,
	struct matrix *m_res)
;
int m_multx (
	struct matrix *m,
	double scalar,
	struct matrix *m_res)
;
double m_trace (struct matrix *m)
;
int m_trans (
	struct matrix *m,
	struct matrix *m_res)
;
int m_atpa1 (
	struct matrix *a,
	struct matrix *p,
	struct matrix *atpa)
;
int m_det (
        struct matrix *m, 
        double *det)
;

int m_pivot (
	struct matrix *m,
	int row,
	int col)
;
int m_inv (struct matrix *m)
;
int m_cholesky_decomposition (struct matrix *m)
;
int m_choleskyinv (struct matrix *m)
;
double m_sym_check (struct matrix *m)
;
int m_add_block (
	struct matrix *m,
	struct matrix *block,
	int row,
	int col)
;
int m_pick_block (
	struct matrix *m,
	struct matrix *block,
	int row,
	int col,
	int drow,
	int dcol)
;
int m_del_row (
	struct matrix *m,
	int row)
;
int m_del_col (
	struct matrix *m,
	int col)
;
int m_polynom_approx (
	struct matrix *data,
	int col_x,
	int col_y,
        int col_w,
	const int col_app,
	const int col_res,
	const int row_beg,
	const int row_end,
	const int degree,
	double x0,
	double *dummy,
	struct matrix *coeff,
	struct matrix *qxx,
	double *s0)
;

double m_polynom_interpol (
	struct matrix *coeff,
        int num_coeff,
	const double dx)
;

int m_separate_sym(
	struct matrix *m, 
	struct matrix *sym, 
	struct matrix *asym)
;
#endif
/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/nfmatcalc.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations (non free stuff taken from numerical recipies)
 *
 * $Id: nfmatcalc.h,v 1.3 1998/11/02 10:42:44 favey Exp $
 *
 *   AUTHORS: Cocard Marc 
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            23-06-1998  --  created
 * $Log: nfmatcalc.h,v $
 * Revision 1.3  1998/11/02 10:42:44  favey
 * added RCS verbosity
 *
 *
 *   CAVEATS: 
 *
 */

#ifndef NF_MATCALC
#define NF_MATCALC
#ifndef GGL_MATCALC
*** ERROR: you must include matcalc.h first, nfmatcalc relies on it!
#endif
/*
 *  JACOBI : calculates Eigenvalues and Eigenvectors from a symmetric matrix a
 *
 * m_jacobi1 : input  : matrix a is conserved
 *             output : eigenvalues (eval) stored in a vector (matrix 1 x DIM)
 *             output : evectors (evec) in form of a rotation matrix
 *             output : total number of individual rotations during  iteration
 *             return value = error
 *
 * m_jacobi2 : input/output : matrix a is overwritten. New matrix a is diagonal
 *                            and contains the eigenvalues (A replaced by D so 
 *                            that a = ev * d * evT)
 *             output : evectors (evec) in form of a rotation matrix
 *             output : total number of individual rotations during  iteration
 *             return value = error
 */

int m_jacobi1(struct matrix *a, struct matrix *eval, struct matrix *evec, int *nrot);

int m_jacobi2(struct matrix *a, struct matrix *evec, int *nrot);

#endif  /* NF_MATCALC */   
/*  DO NOT EDIT THIS FILE, EDIT sed_mat3 INSTEAD!!! */
/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sed_mat3.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch> (pascal version)
 *            Favey Etienne   <favey@geod.ethz.ch>  (c port)
 * $Id: sed_mat3.h,v 1.4 2001/02/21 13:28:19 cocard Exp $
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *                  1987  --  created
 *            19-03-1997  --  pascal version to port
 *            07-04-1997  --  c port
 *            14-05-1997  --  fixed boundary bug, thanks Mario!
 * $Log: sed_mat3.h,v $
 * Revision 1.4  2001/02/21 13:28:19  cocard
 * *** empty log message ***
 *
 * Revision 1.3  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.2  1999/02/01 17:04:08  favey
 * initial update
 *
 * Revision 1.1  1999/02/01 17:02:54  favey
 * Initial revision
 *           
 *
 */

/*
 *  now our own definitions
 */
#ifndef GGL_MAT3
#define GGL_MAT3

int m3_add            ( double a[3][3], double b[3][3], double c[3][3]);
int m3_sub            ( double a[3][3], double b[3][3], double c[3][3]);
int m3_unit           ( double a[3][3]);
int m3_copy           ( double a[3][3], double c[3][3]);
int m3_trans          ( double a[3][3], double c[3][3]);
int m3_inv            ( double a[3][3]);
void m3_write         ( FILE *fpout, double a[3][3]);

int v3_add            ( double a[3], double b[3], double c[3]);
int v3_sub            ( double a[3], double b[3], double c[3]);
int v3_copy           ( double a[3], double c[3]);
int v3_multx          ( double a[3], double x, double c[3]);
int v3_norm           ( double a[3], double c[3]);
double v3_length ( double a[3]);
double v3_scalar ( double a[3], double b[3]);
int v3_cross          ( double a[3], double b[3], double c[3]);

int m3v3_mult         ( double a[3][3], double b[3], double c[3]);
int m3m3_mult         ( double a[3][3], double b[3][3], double c[3][3]);

int rot_getmat        ( int axis, double angle, double rotmat[3][3]);
int rotate            ( double x1[3], int axis, double angle, double x2[3]);
int rotvec_to_mat     ( double rotvec[3], double angle,double rotmat[3][3]);
int rotmat_to_vec     ( double rotmat[3][3], double rotvec[3], double *angle);

#endif

/*  DO NOT EDIT THIS FILE, EDIT sed_mat3 INSTEAD!!! */
/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/sed_mat3.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: matrix operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *   AUTHORS: Cocard Marc     <cocard@geod.ethz.ch> (pascal version)
 *            Favey Etienne   <favey@geod.ethz.ch>  (c port)
 * $Id: sed_mat3.h,v 1.4 2001/02/21 13:28:19 cocard Exp $
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *                  1987  --  created
 *            19-03-1997  --  pascal version to port
 *            07-04-1997  --  c port
 *            14-05-1997  --  fixed boundary bug, thanks Mario!
 * $Log: sed_mat3.h,v $
 * Revision 1.4  2001/02/21 13:28:19  cocard
 * *** empty log message ***
 *
 * Revision 1.3  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.2  1999/02/01 17:04:08  favey
 * initial update
 *
 * Revision 1.1  1999/02/01 17:02:54  favey
 * Initial revision
 *           
 *
 */

/*
 *  now our own definitions
 */
#ifndef GGL_lMAT3
#define GGL_lMAT3

int lm3_add            ( long double a[3][3], long double b[3][3], long double c[3][3]);
int lm3_sub            ( long double a[3][3], long double b[3][3], long double c[3][3]);
int lm3_unit           ( long double a[3][3]);
int lm3_copy           ( long double a[3][3], long double c[3][3]);
int lm3_trans          ( long double a[3][3], long double c[3][3]);
int lm3_inv            ( long double a[3][3]);
void lm3_write         ( FILE *fpout, long double a[3][3]);

int lv3_add            ( long double a[3], long double b[3], long double c[3]);
int lv3_sub            ( long double a[3], long double b[3], long double c[3]);
int lv3_copy           ( long double a[3], long double c[3]);
int lv3_multx          ( long double a[3], long double x, long double c[3]);
int lv3_norm           ( long double a[3], long double c[3]);
long double lv3_length ( long double a[3]);
long double lv3_scalar ( long double a[3], long double b[3]);
int lv3_cross          ( long double a[3], long double b[3], long double c[3]);

int lm3v3_mult         ( long double a[3][3], long double b[3], long double c[3]);
int lm3m3_mult         ( long double a[3][3], long double b[3][3], long double c[3][3]);

int lrot_getmat        ( int axis, long double angle, long double rotmat[3][3]);
int lrotate            ( long double x1[3], int axis, long double angle, long double x2[3]);
int lrotvec_to_mat     ( long double rotvec[3], long double angle,long double rotmat[3][3]);
int lrotmat_to_vec     ( long double rotmat[3][3], long double rotvec[3], long double *angle);

#endif

/* Declarations for parse.c are in inout.h */
/*
 *      FILE: quatern
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: quaternions
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *   AUTHORS: Cocard Marc <cocard@geod.ethz.ch>    
 *
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            18-09-1997  --  created
 */

/*
 *  now our own definitions
 */

#ifndef GGL_QUATER
#define GGL_QUATER

/* addition of two quaternions : result = q1 + q2 */
int q_add ( double q1[4], double q2[4], double result[4]);

/* multiplication of two quaternions: result = q1 * q2 */
int q_mult( double q1[4], double q2[4], double result[4]);

/* conjugate of quaternion 'q' is put in 'result' */
int q_conjugate( double q[4], double result[4]);

/* multipilcation of quaternion 'q' with scalar 'factor' put in 'result' */
int q_scalarmult ( double q[4], double factor, double result[4]);

/* returns the norm of the quaternion q */
double q_norm( double q[4]);

/* inverse of quaternion 'q' is put in 'result' */
int q_inverse( double q[4], double result[4]);

/* converts the rotation matrix 'rotmat' to the quaternion 'result' */
int rotmat_to_quat( double rotmat[3][3], double result[4]);

/* converts a rotation vector (must be a unit vector) and the associated
 * rotation angle (radian) to the corresponding quaternion 'result' */
int rotvec_to_quat( double rotvec[3], double angle, double result[4]);

/* converts the quaternion 'q' to the corresponding rotation matrix 'result' */
int quat_to_rotmat( double q[4], double result[3][3]);

/* rotates the vector 'x' using a quaternion 'q' to the rotated vector 'result' */
int q_rotate( double x[3], double q[4], double result[3]);

#endif

/**
 **        FILE: swiss.h
 **
 ** DESCRIPTION: Functions to transform coordinates
 **              from the bessel ellipsoid to the
 **              Swiss military coordinates and the
 **              vice versa.
 **              Constants and ellipsoid data for calculations
 **              with the swiss projection system.
 **
 **              This header file contains the declarations for
 **              swiss.c and angles.c
 **
 **   COPYRIGHT: Copyright (C) 1996 by Etienne Favey
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
 **     VERSION: 4
 **
 **     CHANGES: dd-mm-yyyy -- description
 **              03-12-1996 -- first release (Version 1.05)
 **              18-12-1996 -- bugfix xyz2blh iteration criterion
 **              10-04-1997 -- adjusted to use with GGL coordtr.c
 **                            taken from public version 1.08,
 **                            renamed blh to llh, introduced
 **                            struct ellipsoid and vectors
 **              07-05-1997 -- added admsstr(), ellsupport, det_ell()
 **              14-01-1998 -- added struct frameconv, select_ellipsoid(),
 **                            select_frameconv()
 **              20-01-1998 -- split swiss.c into swiss.c and trafo.c
 **              13-05-1998 -- parameter for ch_magnetic_deviation changed
 **
 **/
#ifndef SWISS
#define SWISS

/* ----------------------------------------------------------------
 * for swiss.c
 */

int ch_elltosph(double elllat, double elllon, double *sphlat, double *sphlon);
/*
 *  conform projection from ellipsoid to sphere
 *  parameters have to be passed and are returned in rad
 */

int ch_sphtoell(double sphlat,double sphlon, double *elllat, double *elllon);
/*
 *  conform reverse projection from sphere to ellipsoid
 *  parameters have to be passed and are returned in rad
 */

int ch_scale_etom(double elllat, double sphlat, double pselat,double *scale);
/*
 *  scale factor of projection from ellipsiod to plane
 *  is only dependent on latitude
 */

int ch_sphtopse(double sphlat, double sphlon, double *pselat, double *pselon);
/* 
 *  projection from sphere to pseudo equator system
 *  parameters have to be passed and are returned in rad
 */

int ch_psetosph(double pselat, double pselon, double *sphlat, double *sphlon);
/* 
 *  reverse projection from pseudo to sphere
 *  parameters have to be passed and are returned in rad
 */

int ch_merconv(double sphlat, double sphlon, double *conv);
/*
 *  meridian convergence for swiss projection
 */

int ch_psetomil(double pselat, double pselon, double *milx, double *mily);
/* 
 *  mercator projection pseudo to plane
 *  parameters have to be passed and are returned in rad
 */

int ch_miltopse(double milx, double mily, double *pselat, double *pselon);
/*
 *  reverse mercator projection plane to pseudo
 *  parameters have to be passed and are returned in rad
 */

int ch_elltomil(double elllat, double elllon, double *milx, double *mily);
/*
 *  conform projection from ellipsoid to plane (mil)
 *  parameters have to be passed and are returned in rad
 */


int ch_miltoell(double milx, double mily, double *elllat, double *elllon);
/*
 *  reverse mercator projection plane to ellipsoid
 *  parameters have to be passed and are returned in rad
 */

int ch_e2m(double elllat, double elllon, 
	   double *milx, double *mily, 
	   double *scale, double *merconv);
/*
 *  conform projection from ellipsoid to plane (mil)
 *  parameters have to be passed and are returned in rad
 */

int ch_m2e(double milx, double mily,
	   double *elllat, double *elllon,
	   double *scale, double *merconv);
/* 
 *  reverse mercator projection plane to ellipsoid
 *  parameters have to be passed and are returned in rad
 */

double ch_magnetic_deviation (
  const double milx,
  const double mily,
  const double year);

#endif

/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/mapnr.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: determine map number from swiss coordinates
 *
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: mapnr.h,v 1.2 1999/09/08 09:12:18 favey Exp $
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
 *   AUTHORS: Marc Troller (VTB Sommer 1999)
 *            Favey Etienne <favey@geod.ethz.ch>
 *
 *   HISTORY: 
 * $Log: mapnr.h,v $
 * Revision 1.2  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.1  1999/06/30 15:34:56  favey
 * Initial revision
 *
 */
#define BREITE25 (17500.0)
#define HOEHE25 (12000.0)
#define BREITE100 (4*17500.0)
#define HOEHE100 (4*12000.0)

int nrsheets25 (const double west,const double north,
		const double east,const double south,
		int *sheets, const int size);
void edge_25 (int sheet_nr, double *links, double *oben);
int map_25 (const double y, const double x) ;
int nrsheets100 (const double west,const double north,
		 const double east,const double south,
		 int *sheets, const int size);
void edge_100 (int sheet_nr, double *links, double *oben);
int map_100 (const double y, const double x) ;
/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/mapcut.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: fill target map from L+T pixel maps
 *
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: mapcut.h,v 1.4 2000/08/15 16:38:13 favey Exp $
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
 *   AUTHORS: Marc Troller (VTB Sommer 1999)
 *            Favey Etienne <favey@geod.ethz.ch>
 *
 *   HISTORY: 
 * $Log: mapcut.h,v $
 * Revision 1.4  2000/08/15 16:38:13  favey
 * mapras: bugfix and complete rewrite
 *
 * Revision 1.3  2000/03/06 17:01:42  favey
 * GMT colortable and grid output
 *
 * Revision 1.2  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.1  1999/06/30 15:34:56  favey
 * Initial revision
 *
 */

/* test if all sheets are available in the specified directory */
int maptest (const char *mappath,const int anz_karten,int *sheets);

/* calculate the size of the variable outmap */
int map_size(
	const double west,const double north,
	const double east,const double south,
	const double spacing,
	int *outmap_length,int *outmap_width);

void rrborder (
	int sheet_nr, const int map_type,
	const double west,const double north,
	const double east,const double south,
	double *rrwest, double *rrnorth,
	double *rrsouth, double *rreast,
	double *sheet_west, double *sheet_north);

/*zusammenfuegen der einzelnen Ausschnitte der jeweiligen Karten */
int mappaste (
	const double west,const double north,
	const double east,const double south,
	const char *mappath, const char *outfname,
	const int map_type, const int do_gmt, 
	const int lpmm, const int verbose,
	const int anz_karten, int *sheets);

/*
 *      FILE: tools
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: 
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            31-08-1997  --  created
 *
 */

/*
 * Calculates an intersection of two arcs. The routine
 * gives the same result as intersection_circles from
 * unit intersect. Latter should be used. (note that
 * the return value of arc_intersect (= error) differ from
 * that of intersection_circles (=# of intersections).
 */

int arc_intersect (
  const double p1[3], /* x, y, radius of circle 1 */
  const double p2[3], /* x, y, radius of circle 2 */
  double na[2],       /* x, y of intersection 1 */
  double nb[2])       /* x, y of intersection 2 */
;

/*
 * initializes the first combination. All following
 * can be obtained by using combMofN_next
 */
int combMofN_first (int *comb, const int M)
;

/*
 * replaces comb by the next combination. returns 0 on success
 * and -1 after the last combination was reached. 
 * example :
 * combMofN_init(4,3)  comb = 0 1 2  (return value =  0)
 * combMofN_next(4,3)  comb = 0 1 3  (return value =  0)
 * combMofN_next(4,3)  comb = 0 2 3  (return value =  0)
 * combMofN_next(4,3)  comb = 1 2 3  (return value =  0)
 * combMofN_next(4,3)  comb = 1 2 3  (return value = -1 -> no new combination returned)
 */
int combMofN_next (int *comb, const int M, const int N)
;

/*
 * returns the number of combinations (M of N) as a double 
 */
double combMofN_count(const int M, const int N)
;

/*----------------------------------------------------------------------
 * routines for determining all permutation within an array
 * example for count = 3
 * 1:   0  1  2
 * 2:   1  0  2
 * 3:   0  2  1
 * 4:   2  0  1
 * 5.   1  2  0
 * 6:   2  1  0
 *  routines
 *      permutation_init  
 *      permutation_next 
 *      permutation_free
 *      permutation_ctrlprint
 *
 * ( for good understanding of how to use this routines see
 *   the code of permutation_ctrlprint )
 *---------------------------------------------------------------------- */


/*----------------------------------------------------------------------
 * initialize memory for the array <sequence>,  and for an additional
 * control array <ctrlarray> which is needed in permutation_next
 * n is the dimension of the array
 * fill the array with the first sequence
 */
 int permutation_init(int n, int **sequence, int **ctrlarray);

/*--------------------------------------------------------------------
 * permutation to the next sequence
 */
int permutation_next(int n, int * sequence, int * ctrlarray);

/*--------------------------------------------------------------------
 * returns the number of possible permutations
 * works up to n=12 then overflow
 */
long int permutation_count(int n);

/*--------------------------------------------------------------------
 * free the memory allocated by permutation_init
 */
int permutation_free(int *sequence, int *ctrlarray);

/*--------------------------------------------------------------------
 * routine for writing a complete set of all permutations 
 * of n integer in the range from 0..(n-1) to a file 
 * if addcounter == true the lines are numbered
 */
int permutation_ctrlprint(FILE *fp, int n, int addcounter);


/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/window.h,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: file input and output functions including
 *            error handling
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: window.h,v 1.11 1999/09/08 09:12:18 favey Exp $
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
 * $Log: window.h,v $
 * Revision 1.11  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.10  1999/06/30 12:17:30  favey
 * new function iszipped()
 *
 * Revision 1.9  1999/06/21 12:22:44  favey
 * bugfix in get_fname_wopath
 *
 * Revision 1.8  1999/06/21 09:47:52  favey
 * added file name to struct for easier logging
 *
 * Revision 1.7  1999/06/03 10:09:13  cocard
 * print_data_window: 2 format parameters added
 *
 * Revision 1.6  1999/04/08 13:32:52  favey
 * bugfix in time_index_data_window
 *
 * Revision 1.5  1998/10/06 16:19:02  cocard
 *  function time_index_data_window added
 *
 * Revision 1.4  1998/10/02 18:47:40  cocard
 * major changes to readfunction
 *
 * Revision 1.3  1998/07/21 12:40:31  favey
 * added error function
 *
 * Revision 1.2  1998/05/04 07:02:34  favey
 * *** empty log message ***
 *
 * Revision 1.1  1998/04/29 10:29:13  favey
 * Initial revision
 *
 *
 */

#ifndef GGL_MATCALC
*** ERROR: include matcalc.h prior to window.h !
#endif

#ifndef GGL_WINDOW
#define GGL_WINDOW


/*------------------------------------------------------------------
 *  Structure to hold file pointer, matrix pointer, time values, etc.
 */ 

struct data_window {
	double mid_time;	/* mid time of actual window */
	double tail_time;	/* start time of actual window */
	double head_time;	/* end time of actual window */
	double halfwinsize;	/* deltatime = (endtime - starttime)/2 */
	int bufsize;		/* maximal number of rows for a window */
	int valid_count;	/* actual number of rows */
	FILE *fpdata;		/* file pointer to data file */
	char *fname;		/* file name of data file */
	int ispipe;		/* flag whether fpdata is a pipe */
	struct matrix *m;	/* struct containing the data */
};

/*------------------------------------------------------------------
 *  Declaration of a reading function to be called
 *
 *  This function has to be written by the user. It must fill the 
 * specified row of the matrix. The value in the first column is 
 * considered to be the time.
 *
 *  The return value should be one of the following:
 *       1         when the reading of one line was successful
 *       0         when a file input/output error occurred
 *     EOF         when we reached the end of file
 *
 *  here an example:
 *    int read_data(FILE *fpin, struct matrix *m, int row)
 */


/*------------------------------------------------------------------
 *  initialize windowing structure
 *
 *  The return values may be:
 *       = 0       when the update was successful
 *       < 0       error: see function window_err
 */
int init_data_window (
	const double half_window_size,
	const int buffer_size,
	const int numcols,
	const char *datafname,
	struct data_window *w,
	int (*read_new_element) (FILE *fpin, struct matrix *m, int row));

/*------------------------------------------------------------------
 *  update windowing structure
 * 
 *  Moves window forward over file by at most buffersize values
 *  or until tail and head times 
 *
 *  The return values may be:
 *      >= 0       the amount of lines read
 *       < 0       error: see function window_err
 */
int update_data_window (
	struct data_window *w,
	const double new_time,
	int (*read_new_element)(FILE *fpin, struct matrix *m, int row));

/*------------------------------------------------------------------
 *  eliminate windowing structure
 *  Frees the matrix and close the file, sets all values to zero.
 */
void eliminate_data_window (struct data_window *w);

/*------------------------------------------------------------------
 *  print structure window (for debugging purposes)
 */
void print_data_window (FILE *fpout, struct data_window *w, int fmt1, int fmt2);

/*------------------------------------------------------------------
 * returns index i of row for  t(i) <= t < t(i+1)
 */
int time_index_data_window(struct data_window *w, double t);

/*------------------------------------------------------------------
 *  print appropriate error message
 *
 *    -1: end of file
 *    -2: I/O error or incorrect matrix update by your function
 *    -3: necessary memory could not be allocated
 *    -4: time not in window
 *    -5: buffer too small
 */
void window_err (const int errval, const char *fname);

#endif
/*
 *
 *  COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
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
 *  AUTHORS:  Etienne Favey <favey@geod.ethz.ch>
 *
 */
#ifndef GGL_WINTERPOL
#define GGL_WINTERPOL

#ifndef GGL_WINDOW
*** ERROR: include window.h prior to winterpol.h !
#endif

struct wfileinfo {
	FILE *fp;
	int blocksize;
	int num_cols;
	double half_window_size;
	struct data_window dw;
	int change;
	int init;
	struct matrix *coeff[10];
	struct matrix *qxx[10];
	double sigma[10];
	int (*read_func)(FILE *fpin, struct matrix *m, int row);
} ;



/*******************************************************************/
/*
 *  RETURN VALUES: 0 on success, -1 on error
 */
int init_join (
	struct wfileinfo *fi,
	const char *fname);



/*******************************************************************/
/*
 *  RETURN VALUES:	negative:	file update error
 *			positive:	no decent interpolation scheme
 *			zero:		updated and interpolated
 */
int interpol_join (
	struct wfileinfo *fi,	/* the file we operate on */
	const double thetime,	/* time of interpolation */
	const int num_coeff,	/* how many coefficients */
	const int donorm,	/* norm time between -1 .. +1 */
	const int cols,		/* column 1..cols are interpolated */
	double value[]);	/* value[1] .. value[cols] hold results */

#endif /* GGL_WINTERPOL */
/*
 *      FILE: xyz2grd
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: create grid from arbitrary ground points
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: xyz2grd.h,v 1.8 2000/12/11 16:09:00 favey Exp $
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
 *   AUTHORS: Favey Etienne <favey@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *               10-1999  --  created
 *            
 *
 */
#ifndef XYZ2GRD
#define XYZ2GRD

#include <stdio.h>
/***************************************************************************/

struct blockchains {
  double bxmin;
  double bxmax;
  double bxinc;
  long int bxn;
  double bymin;
  double bymax;
  double byinc;
  long int byn;
  int *count;
  struct Node ***ch; /* doppeltes Array von Node Pointern */
};

/***************************************************************************/
struct blockchains * blockchain_init (const double xmin, const double xmax,
				      const double ymin, const double ymax,
				      const int xn, const int yn);


/***************************************************************************/
void blockchain_free (struct blockchains *ablock);


/***************************************************************************/
int blockpos (struct blockchains *b, const double xp, const double yp,
	      int *col, int *row);

/***************************************************************************/
int blockchain_insert (struct blockchains *b, double newp[3],
			const int col, const int row);


/***************************************************************************
 *  interpolation for point with coordinates (x,y) as weighted 
 *  average over the circle points (weight = 1/sqr(distance))
 */
double circle_ipol (const double gx, const double gy,
		    struct Node *circle, const int pcount, const int skipg,
		    const int exponent, const double eps);


/***************************************************************************/
int xyz2grd (FILE *fpxyz, FILE *fppts, const char *grdfn,
	     const double xmin, const double xmax,
	     const double ymin, const double ymax,
	     const double xinc, const double yinc,
	     const double radius, const int exponent, const double eps,
	     const int xn, const int yn,
	     const int surfer, const int quadrant,
	     const int minpt, const int dodiff,  const int skipg,
	     const int schlauch, const int verbose);




#endif

