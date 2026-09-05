/*
 *      FILE: gpstime.c
 *
 *   PROJECT: GGL Common Applications
 *
 *   PURPOSE:  Interactive transformation of following 
 *             different time formats 
 *                 (1)   Year  Month  Day   
 *                 (2)   Year  DayofYear    
 *                 (3)   GPSWeek  DayofWeek 
 *                 (4)   modified julian day
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
 *   AUTHORS: Cocard Marc 
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            15-08-1995  --  created
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "common.h"

#define COPYRIGHT "(c) 1995 Geodesy and Geodynamics Lab, ETH Zurich"

int main(int argc, char * argv[])
{

  struct ttime t;
  int id, DayinYear, DayinWeek, GPSWeek;
  double GPSSec;

  printf (" *  DISCLAIMER:   This program is free software; you can redistribute\n");
  printf (" *              it and/or modify it under the terms of the GNU General\n");
  printf (" *              Public License as published by the Free Software\n");
  printf (" *              Foundation; either version 2 of the License, or\n");
  printf (" *              (at your option) any later version.\n");
  printf (" *                This program is distributed in the hope that it will\n");
  printf (" *              be useful, but WITHOUT ANY WARRANTY; without even the\n");
  printf (" *              implied warranty of MERCHANTABILITY or FITNESS FOR A\n");
  printf (" *              PARTICULAR PURPOSE. See the GNU General Public License\n");
  printf (" *              for more details.\n");
  printf (" *                You should have received a copy of the GNU General\n");
  printf (" *              Public License along with this program; if not, write\n");
  printf (" *              to the Free Software Foundation, Inc., 675 Mass Ave,\n");
  printf (" *              Cambridge, MA 02139, USA.\n");
  do 
    {
      printf("\n");
      printf("      (1)   Year  Month  Day     (1) \n");
      printf("      (2)   Year  DayofYear      (2) \n");
      printf("      (3)   GPSWeek  DayofWeek   (3) \n");
      printf("      (4)   modified julian day  (4) \n");
      printf("      (5)   Actual time          (5) \n");
      printf("      (0)   -- prog stop --      (0) \n");
      printf("                                    ----> ");
      scanf("%d%*[^\n]", &id);
      getchar();

      putchar('\n');
      memset(&t, 0, sizeof(struct ttime));
      t.hour = 1;

      switch (id) 
	{

	case 1:
	  {
	    printf("      (1)   Year  Month  Day     (1)   ");
	    scanf("%d%d%d%*[^\n]", &t.year, &t.month, &t.day);
	    getchar();
	    datmjd(&t);
	    DayinYear = dayofyear1(&t);
	    mjd_to_gps(t.mjd, &GPSSec, &GPSWeek);
	    DayinWeek = (long)(GPSSec / 86400L);
	    break;
	  }

	case 2:
	  {
	    printf("      (2)   Year  DayofYear      (2)   ");
	    scanf("%d%d%*[^\n]", &t.year, &DayinYear);
	    getchar();
	    t.month = 1;
	    t.day = 1;
	    datmjd(&t);
	    t.mjd += DayinYear - 1;
	    mjddat(&t);
	    DayinYear = dayofyear1(&t);
	    mjd_to_gps(t.mjd, &GPSSec, &GPSWeek);
	    DayinWeek = (long)(GPSSec / 86400L);
	    break;
	  }

	case 3:
	  {
	    printf("      (3)   GPSWeek  DayofWeek   (3)   ");
	    scanf("%d%d%*[^\n]", &GPSWeek, &DayinWeek);
	    getchar();
	    GPSSec = DayinWeek * 86400.0;
	    t.mjd = gps_to_mjd(GPSSec, GPSWeek);
	    mjddat(&t);
	    DayinYear = dayofyear1(&t);
	    break;
	  }

	case 4:
	  {
	    printf("      (4)   modified julian day  (4)   ");
	    scanf("%lg%*[^\n]", &t.mjd);
	    getchar();
	    mjddat(&t);
	    mjd_to_gps(t.mjd, &GPSSec, &GPSWeek);
	    DayinYear = dayofyear1(&t);
	    DayinWeek = (long)(GPSSec / 86400L);
	    break;
	  }

	case 5:
	  {
	    systimenow(&t);
	    mjd_to_gps(t.mjd, &GPSSec, &GPSWeek);
	    DayinYear = dayofyear1(&t);
	    DayinWeek = (long)(GPSSec / 86400L);
	    break;
	  }
    }/*case*/

      if (id != 0) 
	{
	  printf("\n\n      (1)   Year  Month  Day   (1)      : ");
	  printf("%5d%3d%3d\n", t.year, t.month, t.day);
	  printf("      (2)   Year  DayofYear    (2)      : ");
	  printf("%5d%6d\n", t.year, DayinYear);
	  printf("      (3)   GPSWeek  DayofWeek (3)      : ");
	  printf("%5d%6d\n", GPSWeek, DayinWeek);
	  printf("      (4)        MJD           (4)      : ");
	  printf("%11.0f\n\n", t.mjd);
	  printf("        Press Enter to continue \n");
	  scanf("%*[^\n]");
	  getchar();
	}
    } while (id != 0);

  return 0;
}
