/*
 *      FILE: xyz2xyz
 *
 *   PROJECT: GGL Common Applications
 *
 *   PURPOSE: geodetic coordinate transformation
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
 *   AUTHORS: Favey Etienne   (c version)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            16-01-1998  --  created
 *
 *   CAVEATS: Needs the GGL Common Library!
 *
 *  Please compile with compiler flags "-lm -lcommon" !
 */

#define COPYRIGHT "(c) 1998 Geodesy and Geodynamics Lab, ETH Zurich"

#include <errno.h>	/*  the allround headers */
#include <stdlib.h>	/*  the allround headers */
#include <stdio.h>	/*  crucial for file input and output   */
#include <unistd.h>	/*  to parse the command line */
#include <ctype.h>	/*  character handling */
#include <string.h>	/*  string handling */
#include <math.h>	/*  mathematics */

#include "common.h"	/*  our own institute's headers */

/*
 *  Print version information to stderr
 *  You should not need to modify this
 */
void display_version()
{
	extern char *program_invocation_short_name;
	fprintf (stderr,
		 "%s version %s\n  %s\n  %s\n",
		 program_invocation_short_name, SRC_VERSION, LIBVER, COPYRIGHT);
	fprintf (stderr,
		 "  maintainers: E. Favey <favey@geod.ethz.ch>, ");
	fprintf (stderr,
		 "M. Cocard <cocard@geod.ethz.ch>\n\n");
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


static int validate_col (int xyz_col[3], int xyz_index[3])
{
  int i;
  for (i=0;i<3;i++) xyz_index[i] = -1;
 
  for (i=0;i<3;i++) {
    if (xyz_col[i] ==  1 ) xyz_index[0] = i;
    if (xyz_col[i] ==  2 ) xyz_index[1] = i;
    if (xyz_col[i] ==  3 ) xyz_index[2] = i;
  }
  for (i=0;i<3;i++) {
#ifdef LIB_DEBUG
    fprintf (stderr, "xyz[%d] %d ",i, xyz_index[i]);
#endif
    if (xyz_index[i] < 0 )  /* needed index not found */
      return GGLERR_INVAL;
  }	
  return 0;
}


/* ---------------------------------------------------------------
 *  That's where we do the work 
 */
int main ( int argc, char **argv)
{
  int ind_read;
  /* for use in command line parsing */
  int verbose = FALSE;
  int pr_version = FALSE;
  char *inframestr = NULL;
  double inepoch = 0.0;
  char *outframestr = NULL;
  double outepoch = 0.0;
  char *inft = NULL;
  char *outft = NULL;
  int prefix_char= 0;
  
  /* for use in main */
  char entireline[MAXLINELENGTH];
  struct scoordcol coordstr; 
  char *prefixstr;  
  char *suffixstr;
  int charcount;
  double source[3], target[3];
  struct frameconv *fc;
  int err;
  unsigned long linecount;
  int in_ftyp[3], in_fcol[3];
  int out_ftyp[3], out_fcol[3];
  int in_index[3], out_index[3];
  

  /* 
   *  let's define our command line options
   *  -h is already included in library 
   */
  struct mainoptions myopts[] =			/* from cmdline.h */
  {
    /* a flag type entry */
    {0,'V', NULL, 0, &pr_version, TRUE, "",
     "Version string."},
    /* a flag type entry */
    {0,'v', NULL, 0, &verbose, TRUE, "", 
     "Print information while calculating."},
    {0,'s', &inframestr, CMDL_CHAR, NULL, 0, " <source frame>",
     "Assume this input frame (*)"},
    {0,'S', &inepoch, CMDL_DOUBLE, NULL, 0, " <source epoch>",
     "Assume this input frame epoch"},
    {0,'t', &outframestr, CMDL_CHAR, NULL, 0, " <target frame>",
     "Assume this output frame (*)"},
    {0,'T', &outepoch, CMDL_DOUBLE, NULL, 0, " <target epoch>",
     "Assume this output frame epoch"},
    {0,'p', &prefix_char, CMDL_INT, NULL, 0, " <prefix_char>",
     "Skip <prefix_char> before parsing values"},
/*    {0,'N', &prefix_col, CMDL_INT, NULL, 0, " <prefix_col>",
     "Skip <prefix_col> before parsing values"}, */
    {0,'I', &inft, CMDL_CHAR, NULL, 0, " <input formtype>",
     "Alternate formtype on input. Default: xyz"},
    {0,'O', &outft, CMDL_CHAR, NULL, 0, " <output formtype>",
     "Alternate formtype on output. Default: xyz"}
  };
  int optcount = (sizeof(myopts)/sizeof(struct mainoptions));  
  char description[] = "geodetic datum conversion of geocentric coordinates";
  char nonoptarg[] = "\n reads from standard in and writes to standard out\n (*) This command line switch may not be ommitted";
  
  /* 
   *  proper error cursing with our own command name 
   */
  extern char *program_invocation_name;        /* from errproc.h */
  extern char *program_invocation_short_name;  /* from errproc.h */
  
  errinit (argv[0]);      /* from errproc.h */
  
  /*
   * Now we process all command line options
   */
  ind_read = process_cmdline_opt (argc, argv, myopts, optcount, 
				  description, nonoptarg);
  
  if (pr_version) { display_version(); return EXIT_SUCCESS; }
  if (ind_read == -3) return EXIT_SUCCESS;
  
  if (argc - ind_read!=0 || inframestr == NULL || outframestr == NULL) {
    print_usage_info(FALSE, myopts, optcount, 
		     description, nonoptarg);
    return EXIT_FAILURE;
  }

  /*
   *  String default values are filled only down here
   */
  if (inft == NULL) inft = strdup ("xyz");
  if (outft == NULL) outft = strdup ("xyz");

  if (verbose != FALSE) display_version();

  /*
   *  this is where we do all the action!
   *  read all kinds of stuff from file
   */
  
  err = select_formtype(inft, in_ftyp, in_fcol);
#ifdef LIB_DEBUG
  fprintf (stderr, "%d in typ %d %d %d col %d %d %d\n", err,
	   in_ftyp[0], in_ftyp[1], in_ftyp[2],
	   in_fcol[0], in_fcol[1], in_fcol[2]);
#endif

  if (err != 0) return EXIT_FAILURE;
  coordstr.count = in_fcol[0] + in_fcol[1] + in_fcol[2];
  err = select_formtype(outft, out_ftyp, out_fcol);
#ifdef LIB_DEBUG
  fprintf (stderr, "%d out typ %d %d %d col %d %d %d\n", err,
	   out_ftyp[0], out_ftyp[1], out_ftyp[2],
	   out_fcol[0], out_fcol[1], out_fcol[2]);
#endif

  if (err != 0) return EXIT_FAILURE;
  
  fc = select_frameconv (inframestr, inepoch, outframestr, outepoch);
  if (fc == NULL) return EXIT_FAILURE;
  
  linecount = 0L;
  
  /* validate column types and resolve indices */
  err = validate_col (in_ftyp, in_index);
  if (err == 0)
    err = validate_col (out_ftyp, out_index);
  if (err != 0) {
    fprintf (stderr,"%s: Tell me how to resolve these columns...\n",
	     program_invocation_short_name);
    fprintf (stderr, "%s: I would guess a misconfigured config file or\n",
	     program_invocation_short_name);
    fprintf (stderr, "%s: are you trying something like latitude in meters?\n",
	     program_invocation_short_name);
    return EXIT_FAILURE;
  }
  
  /* from here on NEW simplied version MC (06.12.2001) 
     strange behaviour of old version */
    
  prefixstr = (char *) malloc((prefix_char+1)*sizeof(char));
  memset(prefixstr, '\0', prefix_char+1);
  
  while (read_line_static (stdin, entireline, MAXLINELENGTH) != EOF) {
    linecount++;
    if (iscomment(entireline)) 	{
      fprintf (stdout, "%s", entireline);
      continue;
    }
    
    memcpy(prefixstr, entireline, prefix_char);
    err=sscanf(entireline+prefix_char, "%lf %lf %lf %n\n",
	       &(source[0]), &(source[1]), &(source[2]), &charcount);
    
    if (err != 3) {
      fprintf(stderr,"%s: invalid input on line nr %ld\n",
	      program_invocation_short_name, linecount); 
      return EXIT_FAILURE;  
    }
    
    suffixstr=entireline+prefix_char+charcount;
    while(suffixstr[-1] == ' ') suffixstr--;
    
    //   fprintf(stderr,"%s|%f %f %f|%s",prefixstr, source[0], source[1], source[2], suffixstr);
    
    /* conversion routines */
    frame_conv(fc, source, target);
    fprintf(stdout,"%s", prefixstr);
    fprintf(stdout," %14.4f %14.4f %14.4f",target[0], target[1], target[2]);
    fprintf(stdout,"%s\n",suffixstr);
  }   

  return EXIT_SUCCESS;
}
