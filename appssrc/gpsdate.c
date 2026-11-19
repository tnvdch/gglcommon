/*--------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <common.h>

/*  ETHZ-IGP-GGL  *   MC   */
#define COPYRIGHT "(c) 1999 Geodesy and Geodynamics Lab, ETH Zurich"
/*--------------------------------------------------------------------------*/
/* Purpose : Commandline transformation of following different time formats */
/*           (1)   Year  Month  Day                                         */
/*           (2)   Year  DayofYear                                          */
/*           (3)   GPSWeek  DayofWeek                                       */
/*           (4)   modified julian day                                      */
/*--------------------------------------------------------------------------*/
/*  $Id: gpsdate.c,v 1.4 2001/02/21 13:28:26 cocard Exp favey $  */
/*
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
 */
/*                                                                          */
/*  Remarks   :  none                                                       */
/*--------------------------------------------------------------------------*/

void display_version()
{
        extern char *program_invocation_short_name;
        fprintf (stderr,
                 "%s version %s \n  %s\n  %s\n",
                 program_invocation_short_name, SRC_VERSION, LIBVER, COPYRIGHT);
        fprintf (stderr,
                 "  maintainers: E. Favey <favey@geod.ethz.ch>\n\n");
  fprintf (stderr," *  DISCLAIMER:   This program is free software; you can redistribute\n");
  fprintf (stderr," *              it and/or modify it under the terms of the GNU General\n");
  fprintf (stderr," *              Public License as published by the Free Software\n");
  fprintf (stderr," *              Foundation; either version 2 of the License, or\n");
  fprintf (stderr," *              (at your option) any later version.\n");
  fprintf (stderr," *                This program is distributed in the hope that it will\n");
  fprintf (stderr," *              be useful, but WITHOUT ANY WARRANTY; without even the\n");
  fprintf (stderr," *              implied warranty of MERCHANTABILITY or FITNESS FOR A\n");
  fprintf (stderr," *              PARTICULAR PURPOSE. See the GNU General Public License\n");
  fprintf (stderr," *              for more details.\n");
  fprintf (stderr," *                You should have received a copy of the GNU General\n");
  fprintf (stderr," *              Public License along with this program; if not, write\n");
  fprintf (stderr," *              to the Free Software Foundation, Inc., 675 Mass Ave,\n");
  fprintf (stderr," *              Cambridge, MA 02139, USA.\n\n");
}

/*--------------------------------------------------------------------------*/

static void complete_date (
	const int id,
	struct ttime *t,
	int *DayinYear,
	int *DayinWeek,
	int *GPSWeek)
{
	double GPSSec = 0.0;
	switch (id) {

	case 1:
		/* given: t->year, t->month, t->day */
		datmjd(t);
		*DayinYear = dayofyear1(t);
		mjd_to_gps(t->mjd, &GPSSec, GPSWeek);
		*DayinWeek = (long)(GPSSec / 86400L);
		break;

	case 2:
		/* given: t->year, DayinYear */
		t->month = t->day = 1;
		datmjd(t);
		t->mjd += *DayinYear - 1;
		mjddat(t);
		*DayinYear = dayofyear1(t);
		mjd_to_gps(t->mjd, &GPSSec, GPSWeek);
		*DayinWeek = (long)(GPSSec / 86400L);
		break;

	case 3:
		/* given: GPSWeek, DayofWeek */
		GPSSec = *DayinWeek * 86400.0;
		t->mjd = gps_to_mjd(GPSSec, *GPSWeek);
		mjddat(t);
		*DayinYear = dayofyear1(t);
		break;

	case 4:
		/* given: modified julian day t->mjd */
		mjddat(t);
		mjd_to_gps(t->mjd, &GPSSec, GPSWeek);
		*DayinYear = dayofyear1(t);
		*DayinWeek = (long)(GPSSec / 86400L);
		break;

	case 5:
		/* given: nothing */
		systimenow(t);
		mjd_to_gps(t->mjd, &GPSSec, GPSWeek);
		*DayinYear = dayofyear1(t);
		*DayinWeek = (long)(GPSSec / 86400L);
		break;
    }/*case*/

}

static void print_date (
	struct ttime t,
	int DayinYear,
	int DayinWeek,
	int GPSWeek)
{
      printf("\n  Year  Month   Day : %5d%3d%3d\n", t.year, t.month, t.day);
      printf("  Year    DayofYear : %5d%6d\n", t.year, DayinYear);
      printf("  GPSWeek DayofWeek : %5d%6d\n", GPSWeek, DayinWeek);
      printf("          MJD       : %11.0f\n\n", t.mjd);
}

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{

        int ind_read;           /* index of first non option argument */
        int pr_version = FALSE;
        int dogps = FALSE;
        int donow = FALSE;
	struct ttime t;
	int DayinYear, DayinWeek, GPSWeek;

        struct mainoptions myopts[] =                   /* from cmdline.h */
        {
                {0,'V', NULL, 0, &pr_version, TRUE, "",
                 "Version string."},
                {0,'g', NULL, 0, &dogps, TRUE, "", 
                 "Interpret <date> as GPS_week and dayofweek."},
                {0,'n', NULL, 0, &donow, TRUE, "", 
                 "Now! Print dates of today."},
        };
        int optcount = (sizeof(myopts)/sizeof(struct mainoptions));
        char description[] = "calculate GPS dates given one version of them";
        char nonoptarg[] = " [<date>]\n    <date> may be (see also option -g):\n    year month day\n    year dayofyear\n    MJD";

        errinit (argv[0]);                              /* from errproc.h */
        ind_read = process_cmdline_opt (argc, argv, myopts, optcount, 
                                        description, nonoptarg);
        if (pr_version) { display_version(); return EXIT_SUCCESS; }
        if (ind_read == -3) return EXIT_SUCCESS;
	memset(&t, 0, sizeof(struct ttime));
	t.hour = 1;
	switch (argc-ind_read) {
	case 3:
		printf("Interpreting 3 arguments as year, month, day.\n");
		t.year = strtod(argv[ind_read], NULL);
		t.month = strtod(argv[ind_read+1], NULL);
		t.day = strtod(argv[ind_read+2], NULL);
		complete_date(1, &t, &DayinYear, &DayinWeek, &GPSWeek);
		print_date(t, DayinYear, DayinWeek, GPSWeek);
		return EXIT_SUCCESS;
	case 2:
		if (dogps) {
			printf("Interpreting 2 arguments as GPSweek, dayofweek.\n");
			GPSWeek = strtod(argv[ind_read], NULL);
			DayinWeek = strtod(argv[ind_read+1], NULL);
			complete_date(3, &t, &DayinYear, &DayinWeek, &GPSWeek);
		} else {
			printf("Interpreting 2 arguments as year, dayofyear.\n");
			t.year = strtod(argv[ind_read], NULL);
			DayinYear = strtod(argv[ind_read+1], NULL);
			complete_date(2, &t, &DayinYear, &DayinWeek, &GPSWeek);
		}
		print_date(t, DayinYear, DayinWeek, GPSWeek);
		return EXIT_SUCCESS;
	case 1:
		printf("Interpreting 1 argument as MJD.\n");
		t.mjd = strtod(argv[ind_read], NULL);
		complete_date(4, &t, &DayinYear, &DayinWeek, &GPSWeek);
		print_date(t, DayinYear, DayinWeek, GPSWeek);
		return EXIT_SUCCESS;		
	case 0:
		if (donow) {
			printf("Today is:\n");
			complete_date(5, &t, &DayinYear, &DayinWeek, &GPSWeek);
			print_date(t, DayinYear, DayinWeek, GPSWeek);
			return EXIT_SUCCESS;
		}
		/* fall through */
	default:
		print_usage_info(FALSE, myopts, optcount, 
				 description, nonoptarg);
	}

	return EXIT_FAILURE;
}
/*--------------------------------------------------------------------------*/
/* End. */
