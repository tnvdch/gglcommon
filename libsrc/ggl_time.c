/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/ggl_time.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: time operations
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: ggl_time.c,v 1.15 2001/02/21 14:00:19 cocard Exp favey $
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
 * $Log: /GPS/SOFTWARE/common_devel/libsrc/ggl_time.c $
 * 
 * 2     10.12.04 16:14 Nobs
 * make compile under darwin
 * 
 * 1     31.08.01 14:51 Etienne
 * Revision 1.15  2001/02/21 14:00:19  cocard
 * *** empty log message ***
 *
 * Revision 1.14  2001/02/21 13:28:19  cocard
 * *** empty log message ***
 *
 * Revision 1.13  2000/07/25 13:18:46  favey
 * introducing OS_AIX
 *
 * Revision 1.12  2000/03/06 17:01:42  favey
 * YP
 *
 * Revision 1.11  1999/10/07 09:07:05  favey
 * bugfix in gps_to_dat
 *
 * Revision 1.10  1999/09/28 13:10:03  favey
 * new function gps_to_dat, addresses double resolution problem
 *
 * Revision 1.9  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.8  1999/07/06 09:37:21  ypeter
 * Changed timenow to systimenow. Add localtimenow.
 *
 * Revision 1.7  1999/07/06 07:53:49  favey
 * bugfix: microseconds instead of milli in mjdnow
 *
 * Revision 1.6  1998/10/30 14:22:07  favey
 * changed config.h to global.h
 *
 * Revision 1.5  1998/09/28 15:35:32  cocard
 * warnings fixed
 *
 * Revision 1.4  1998/08/05 17:24:24  cocard
 * small changes for c++ compilation
 *
 * Revision 1.3  1998/08/04 15:09:28  cocard
 * change in strtime1
 *
 * Revision 1.2  1998/05/15 19:17:07  cocard
 * *** empty log message ***
 *
 * Revision 1.1  1998/04/22 12:34:18  favey
 * Initial revision
 *
 * Revision 0.4  1998/01/21 11:20:02  favey
 * fixed bug in daymonthacc index handling
 *
 * Revision 0.3  1997/12/24 12:10:32  favey
 * followup of V. 0.29
 *
 *
 */

/*
 *  standard libraries
 */
#include <errno.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *  our own headers
 */
#include "global.h"
#if defined(OS_LINUX) || defined(OS_AIX) || defined(OS_CYGWIN) || defined(__APPLE__)
#include <sys/time.h>
#endif
#include "ggl_time.h"

/*
 *  leapyear defined in header
 */
const char monthname[12][4] = {{"jan"},{"feb"},{"mar"},
			       {"apr"},{"may"},{"jun"},
			       {"jul"},{"aug"},{"sep"},
			       {"oct"},{"nov"},{"dec"}};

const int dayinmonth[12] = {31, 28, 31, 30, 31, 30,
			    31, 31, 30, 31, 30, 31};

const int dayinmonthacc[12] = {  0,  31,  59,  90, 120, 151,
			       181, 212, 243, 273, 304, 334};

double dat_to_mjd (const int year,
		   const int month,
		   const int day,
		   const int hour,
		   const int min,
		   const double sec)
{
  const double k = -15.0 + 1720996.5 - 2400000.5;
  int y, m;
  double t;

  if (month <= 2)
    {
      y = year - 1;
      m = month + 12;
    }
  else
    {
      y = year;
      m = month;
    }
  t = (floor(365.25*y) + floor(30.6001*(m+1)) + k + day);
  return (t + ((sec/60.0+min)/60.0+hour)/24.0);
}

int datmjd (struct ttime *t)
{
  t->mjd = dat_to_mjd (t->year, t->month, t->day, t->hour, t->min, t->sec);
  return 0;
}

int mjd_to_dat (double mjd,
		int *year,
		int *month,
		int *day,
		int *hour,
		int *min,
		double *sec)
{
  double t1, t2, t3, t4;
  double dummy;

  t1 = floor(mjd) + 2400001.0;
  t4 = modf (mjd, &dummy);
  *day = (int) floor((t1-1867216.25)/36524.25);
  t2 = t1 + 1 + *day - floor (*day/4);
  t3 = t2 - 1720995.0;
  *year = (int) floor((t3-122.1)/365.25);
  t1 = 365.25 * *year - modf ((365.25* *year), &dummy);
  *month = (int) floor ((t3-t1)/30.6001);
  *sec = t3 - t1 - floor(30.6001* *month) + t4;
  if (*month>13) *month -= 13; else *month -= 1;
  if (*month<=2) *year += 1;
  *day = (int) floor(*sec);
  *sec = modf (*sec, &dummy) * 24; *hour = (int) floor(*sec);
  *sec = modf (*sec, &dummy) * 60; *min  = (int) floor(*sec);
  *sec = modf (*sec, &dummy) * 60;
  return 0;
}

int mjddat (struct ttime *t)
{
  return mjd_to_dat (t->mjd, &(t->year), &(t->month), &(t->day),
		     &(t->hour), &(t->min), &(t->sec));
}

double gps_to_mjd (double sec,
		   int week)
{
  return (GPSMJD0 + week * 7.0 + sec / 86400.0);
}

int mjd_to_gps (double mjd,
		double *sec,
		int *week)
{
  double w;
  double t = mjd - GPSMJD0;
  modf (t/7.0, &w);
  *sec = (t - (*week = (int) w) * 7.0) * 86400.0;
  return 0;
}

/* returns secofday, struct ttime *t is filled with correct data */
/* routine addresses problem of insufficient double resolution of mjd */
double gps_to_dat (struct ttime *t,
		const int week,
		const double sec)
{
  double dummy;
  double secofday = fmod (sec, 86400.0);
  t->mjd = gps_to_mjd (sec, week);
  mjddat (t);

  t->hour= (int) floor (secofday / (60.0*60.0) );
  t->min = (int) fmod (secofday, 60.0*60.0);
  t->sec = fmod (t->min, 60.0) + modf (sec, &dummy);
  t->min = (int) floor (t->min / 60.0);
  return secofday;
}

double yearday_to_mjd (int year,
		       int daynr)
{
  return rint (dat_to_mjd(year,1,1,0,0,0.01)+daynr-1.0);
}

int mjd_to_yearday (double mjd,
		    int *year,
		    int *daynr)
{
  int month, day, hour, min;
  double sec;
  mjd_to_dat(mjd, year, &month, &day, &hour, &min, &sec);
  *daynr = dayinmonthacc[month-1] + day;
  if ( (month > 2) && leapyear(*year) ) *daynr += 1;
  return 0;
}

int dayofyear1(struct ttime *t)
{
  return (int) floor(t->mjd - dat_to_mjd (t->year,1,1,0,0,0.0)) + 1;
}

int dayofyear2(double mjd)
{
  int year, day;
  mjd_to_yearday (mjd, &year, &day);
  return day;
}

#if defined(OS_LINUX) || defined(OS_AIX) || defined(OS_CYGWIN) || defined(__APPLE__)
double mjdnow ()
{
  struct timeval t;
  int errchk;

  errchk=gettimeofday(&t,NULL);
  return ((t.tv_sec+(t.tv_usec/(double)1000000.0)) /(double) 86400.0 + SYSMJD0);
}

#else

double mjdnow ()
{
  return (time(NULL) /(double) 86400.0 + SYSMJD0);
}
#endif /* OS_LINUX */

int systimenow (struct ttime *t)
{
  t->mjd = mjdnow();
  return mjddat (t);
}

#ifndef OS_CYGWIN
int localtimenow (struct ttime *t)
{
  /* mjdnow return the hardware clock value
   * How do we get the right local time incl. daylight saving ??
   * My guess is that we have to substract the 'timezone' from
   * the system time (YP/06-MAR-2000)
   */
  t->mjd = mjdnow();
  tzset();
#ifndef __APPLE__ // FIXME: fix this !
  t->mjd -= (double)(timezone) / 86400.0;
#endif
  return mjddat (t);
}
#endif

char * strtime1 (struct ttime *t,
		 int i)
{
  char format[256], outstr[256];
  int digits;

  digits=(i==0)?2:3;
  sprintf (format, "%%04d %%02d %%02d   %%02d %%02d %%0%d.%df%c", i+digits, i, '\0');
  sprintf (outstr, format, t->year, t->month, t->day, t->hour,
	   t->min, t->sec);
  return strdup (outstr);
}

char * strtime2 (double mjd,
		 int i)
{
  struct ttime t;
  t.mjd = mjd + exp((-i-2)*log(10))/86400.0;
  mjddat(&t);
  return strtime1(&t,i);
}

char * strtime3 (struct ttime *t)
{
  char outstr[256];
  struct tm mytime;

  memset (&mytime, 0, sizeof(struct tm));
  mytime.tm_sec = (int) rint(t->sec);
  mytime.tm_min = t->min;
  mytime.tm_hour = t->hour;
  mytime.tm_mday = t->day;
  mytime.tm_mon = t->month -1;
  mytime.tm_year = t->year - 1900;
  mktime (&mytime); /* normalize contents of tm struct to prevent SIGSEGV */
  strftime (outstr, 128, "%Y-%b-%d   %H:%M", &mytime);
  return strdup (outstr);
}

int readtime1 (FILE *fpin,
	       struct ttime *t)
{
  if (fscanf (fpin, " %d %d %d %d %d %lf", &(t->year), &(t->month),
	      &(t->day), &(t->hour), &(t->min), &(t->sec)) == EOF )
    {
      fprintf (stderr, "GGL Common Library -- readtime1 %s: ", strerror (errno));
      return -1;
    }
  if (t->year < 100) 
    t->year+=(t->year>70)?1900:2000;

  datmjd (t);
  return 0;
}

int sreadtime1 (const char *instring,
	       struct ttime *t)
{
  if (sscanf (instring, " %d %d %d %d %d %lf", &(t->year), &(t->month),
	      &(t->day), &(t->hour), &(t->min), &(t->sec)) < 6 )
    {
      return -1;
    }

  if (t->year < 100) 
    t->year+=(t->year>70)?1900:2000;

  datmjd (t);
  return 0;
}

double readtime2 (FILE *fpin)
{
  struct ttime t;
  if (readtime1 (fpin, &t) < 0) return -1;
  return (t.mjd);
      
}

double sreadtime2 (const char *instring)
{
  struct ttime t;
  if (sreadtime1 (instring, &t) < 0) return -1.0;
  return (t.mjd);
      
}

double GMST (double mjd)
{
  double t1, t2;
  double dummy;
  t1 = (floor(mjd) + (JULDAT0 - 2451545.0))/36525.0;
  t2 = (6.0+41.0/60.0+50.54841/3600.0);
  t2 += (8640184.812866/3600.0)*t1;
  t2 += (0.093104/3600.0)*t1*t1;
  t2 -= (6.2E-6/3600.0)*t1*t1*t1;
  t2 += modf (mjd, &dummy)*24.0*366.2422/365.2422;
  while (t2<0.0) t2 += 24.0;
  while (t2>24.0) t2 -= 24.0;
  return t2;
}

double secofday (const char *mjd)
{
  char *timestring, *token, *endtok;
  double seconds = -1.0;
  double days;

  timestring = strdup (mjd);
  if (timestring != NULL)
  {
	  token = strchr (timestring, '.');
	  if (token != NULL)
	  {
		  days = strtod (token, &endtok);
		  if (token != endtok)
		  {
			  seconds = days * 86400.0;
		  }
	  }
	  free(timestring);
  }
  return seconds;
}
