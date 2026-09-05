/*
 *      FILE: mapras
 *
 *   PROJECT: GGL Common Applications
 *
 *   PURPOSE: cut a region from the L+T pixelmaps
 *
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
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
 *   AUTHORS: Marc Troller (VTB Sommer 1999)
 *            Favey Etienne
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *      
 *            30-06-1999  --  created
 *
 */

#define COPYRIGHT "(c) 1999 Geodesy and Geodynamics Lab, ETH Zurich"
#define SIZE 1000

#include <errno.h>	/*  the allround headers */
#include <stdlib.h>	/*  the allround headers */
#include <stdio.h>	/*  crucial for file input and output   */
#include <ctype.h>	/*  character handling */
#include <string.h>	/*  string handling */
#include <math.h>	/*  mathematics */

#include <common.h>	/*  our own institute's headers */

/*
 *  Print version information to stderr
 *  You should not need to modify this
 */
void display_version()
{
	extern char *program_invocation_short_name;
	extern char *common_version;
	fprintf (stderr,
		 "%s version %s\n  %s\n  %s\n",
		 program_invocation_short_name, SRC_VERSION, common_version, COPYRIGHT);
	fprintf (stderr,
		 "  maintainers: M. Troller <troller@geod.baug.ethz.ch>\n\n");
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


/*  That's where we do the work 
 */
int main ( int argc, char **argv)
{
  int ind_read, lv, anz_karten;   /* index of first non option argument */
  int sheets[SIZE];
  int err;
  /*  variables to hold arguments from command line parsing
   *  and default settings for command line options 
   */
  int verbose = FALSE;
  int pr_version = FALSE;
  int linespermm = 20;
  int do_gmt = FALSE;
  char *region = NULL;
  char *mappath = NULL;
  double west, east, south, north;
  int map_type = 25;
  
  /*  let's define our command line options
   *  -h is already included in library 
   */
  struct mainoptions myopts[] =			/* from cmdline.h */
  {
    {0,'V', NULL, 0, &pr_version, TRUE, "",
     "Version string."},
    {0,'v', NULL, 0, &verbose, TRUE, "", 
     "Print information while calculating."},
    {0,'l', &linespermm, CMDL_INT, NULL, 0, " <10|20>", 
     "Scanning resolution in lines per millimeter, default=20"},
    {0,'R', &region, CMDL_CHAR, NULL, 0, " <west/east/south/north>",
     "Region to cut in swiss coordinates"},
    {0,'t', &map_type, CMDL_INT, NULL, 0, " <25|100>",
     "Input map is 1:25000 or 1:100000"},
    {0,'G', NULL, 0, &do_gmt, TRUE, "",
     "Output map is GMT colortable and coordinates with gridvalue"},
    {0,'M', &mappath, CMDL_CHAR, NULL, 0, " <path_to_maps>",
     "Directory where the maps are stored"}
  };
  
  int optcount = (sizeof(myopts)/sizeof(struct mainoptions));
  char description[] = "cut pixel map from repository";
  char nonoptarg[] = " <output_filename>";
  extern char *program_invocation_name;		/* from errproc.h */
  extern char *program_invocation_short_name;	/* from errproc.h */
  
  errinit (argv[0]);				/* from errproc.h */
    ind_read = process_cmdline_opt (argc, argv, myopts, optcount, 
				  description, nonoptarg);
  if (pr_version) { display_version(); return EXIT_SUCCESS; }
  if (ind_read == -3) return EXIT_SUCCESS;
  if (ind_read < 0 || argc - ind_read != 1) {
    print_usage_info(FALSE, myopts, optcount, 
		     description, nonoptarg);
    return EXIT_FAILURE;
  }

  /******************************************************************
   *  String default values are filled only down here
   */
  if (mappath == NULL) mappath = strdup("/cdrom");
  if (region == NULL) { 
    print_usage_info(FALSE, myopts, optcount, description, nonoptarg); 
    return EXIT_FAILURE; 
  } else {
    int i, n;
    for (i=0;i<strlen(region);i++)
      if (region[i] == '/') region[i] = ' ';
    n = sscanf (region, "%lf %lf %lf %lf", &west, &east, &south, &north);
    if (n<4 || west > east || north < south) {
      print_usage_info(FALSE, myopts, optcount, description, nonoptarg); 
      return EXIT_FAILURE; 
    }
  }

  /* 	Auflistung der einzelnen Kartenblaetter mit Gesamtanzahl;        */
  switch (map_type) {
  case 25: 
    anz_karten=nrsheets25(west, north, east, south, sheets, SIZE);
    break;
  case 100: 
    anz_karten=nrsheets100(west, north, east, south, sheets, SIZE);
    break;
  default:
    print_usage_info(FALSE, myopts, optcount, description, nonoptarg); 
    return EXIT_FAILURE; 
  }
  if (verbose) {
    fprintf (stderr, "Number of sheets: %d\n", anz_karten);
    for (lv=0;lv<anz_karten;lv++) fprintf (stderr, "%d\n", sheets[lv]);
  }
  
  err = maptest (mappath, anz_karten, sheets);
  if (err != 0) {
    fprintf (stderr, "%s: Some required map missing!\n",
	     program_invocation_short_name);
    return EXIT_FAILURE;
  }

  err = mappaste (west, north, east, south, mappath, argv[ind_read],
		  map_type, do_gmt, linespermm, verbose, anz_karten, sheets);
  if (err != 0) {
    fprintf(stderr, "%s: Error while pasting maps together!\n",
	    program_invocation_short_name);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


