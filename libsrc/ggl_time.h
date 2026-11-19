/*
 *      FILE: ggl_time
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: time operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: ggl_time.h,v 1.5 1999/09/28 13:10:03 favey Exp favey $
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

/*
 * MJD for 0.5 days into the year 2000
 */
#define MJD_J2000 (51544.5)

/*
 * number of days in a Julian century
 */
#define JULCENTURY (36525.0)


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
