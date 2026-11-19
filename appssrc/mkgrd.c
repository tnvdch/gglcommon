/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/apps/mkgrd.c,v $
 *
 *   PROJECT: GGL GPSkin Applications
 *
 *   PURPOSE: create grid from arbitrary points
 *
 * COPYRIGHT: (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 *
 *   AUTHORS: Etienne Favey
 * $Id: mkgrd.c,v 1.12 2000/12/11 16:09:31 favey Exp favey $
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *             7-10-1999  --  created, linear approximation
 * $Log: /GPS/SOFTWARE/common_devel/appssrc/mkgrd.c $
 * 
 * 1     31.08.01 14:49 Etienne
 * Revision 1.12  2000/12/11 16:09:31  favey
 * new options, exponent and z
 *
 * Revision 1.11  2000/03/17 16:07:06  favey
 * new option -S
 *
 * Revision 1.10  1999/11/08 16:23:06  favey
 * minor bugfix
 *
 * Revision 1.9  1999/11/08 15:10:01  favey
 * new option -f, -D
 *
 * Revision 1.8  1999/10/14 11:48:25  favey
 * default changed
 *
 * Revision 1.7  1999/10/14 11:19:51  favey
 * bugfix in array position calculation
 *
 * Revision 1.6  1999/10/12 14:56:20  favey
 * *** empty log message ***
 *
 * Revision 1.5  1999/10/12 14:40:03  favey
 * changed default behaviour of option -b
 *
 * Revision 1.4  1999/10/12 12:51:51  favey
 * new option -M
 * Revision 1.3  1999/10/11 15:45:08  favey
 * beta release
 *
 *
 *   CAVEATS: Needs the GGL Common Library!
 *
 *  Please compile with compiler flags "-lm -lcommon" !
 */

#define COPYRIGHT "(c) 1999 Geodesy and Geodynamics Lab, ETH Zurich"

#include <errno.h>      /*  the allround headers */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include <common.h>

void display_version() 
{
  extern char *program_invocation_short_name;
  extern char *common_version;
  fprintf (stderr, "%s version %s\n  %s\n  %s\n",
	   program_invocation_short_name, SRC_VERSION, 
	   common_version, COPYRIGHT);
  fprintf (stderr, "  maintainers: E. Favey <favey@geod.baug.ethz.ch>\n\n");
}
/***************************************************************************/
static int b_option (char *s, int *xn, int *yn)
{
  char *p1, *p2;
  if (s==NULL) return -1;
  p1=(s[1]=='b')?s+2:s;
  p2=p1; *xn = strtol(p2, &p1, 10); if (p1==p2) return -1;
  if (p1[0] == '\0') *yn=*xn;
  else { p2=p1+1; *yn = strtol(p2, &p1, 10); if (p1==p2) return -1; }
  if (*xn < 0 || *yn < 0) return -1;
  return 0;
}

/***************************************************************************/
static void b_calcoption (const double r,
			 const double xmin, const double xmax,
			 const double ymin, const double ymax,
			 int *xn, int *yn)
{
  double px, py;
  px = floor( (xmax-xmin)/(2.0*r) );
  py = floor( (ymax-ymin)/(2.0*r) );
  *xn = (int) (px<4.0)?4.0:px;
  *yn = (int) (py<4.0)?4.0:py;
}

/***************************************************************************/
int main (int argc, char **argv)
{
  double radius=20.0;
  int xn=5, yn=5, err=0;
  double xmin, xmax, ymin, ymax, xinc, yinc;
  int surfer=FALSE, schlauch=FALSE;
  int quadrant=TRUE, minpt=1;
  int dodiff=FALSE;
  int skipg=FALSE;
  int exponent=2;
  double eps=1.0e-5;

  extern char *program_invocation_name;
  extern char *program_invocation_short_name;
  char description[] = "reads x y z and creates grid";
  char nonoptarg[] = " [<xyz-file>]";
  FILE *fpxyz, *fpgrd, *fppts = NULL;
  int 
    ind_read,	/* index of first non option argument */
    verbose = TRUE,
    version_flag = FALSE;
  char *grdfn = NULL;
  char *ptfn = NULL;
  char *grdrange = NULL;
  char *grdincr = NULL;
  char *bxnyn = NULL;
  struct mainoptions myopts[] =
  {
    {0, 'V', NULL, 0, &version_flag, TRUE, "",
     "version string."}
    , {0, 'q', NULL, 0, &verbose, FALSE, "",
       "run in quiet mode."}
    , {0, 'G', &grdfn, CMDL_CHAR, NULL, 0, " <grd_fname>",
       "grid file to write output to. compulsory option"}
    , {0, 'R', &grdrange, CMDL_CHAR, NULL, 0, " <west/east/north/south>",
       "boundaries of final grid. compulsory option"}
    , {0, 'I', &grdincr, CMDL_CHAR, NULL, 0, " <x-spacing[/y-spacing]>",
       "spacing of final grid. compulsory option"}
    , {0, 'r', &radius, CMDL_DOUBLE, NULL, 0, " <radius>",
       "radius of points around new grid point to include [meters] default: 20."}
    , {0, 'n', &exponent, CMDL_INT, NULL, 0, " <exponent>",
       "exponent n in weight function (see Option -e). default: 2."}
    , {0, 'e', &eps, CMDL_DOUBLE, NULL, 0, " <epsilon>",
       "smoothing: eps in weight function: w=eps^n/(d^n+eps^n), default: 1.0e-10"}
    , {0, 'b', &bxnyn, CMDL_CHAR, NULL, 0, " <xn[/yn]>",
       "number of blocks within range. default: calculated from <radius>"}
    , {0, 'Q', NULL, 0, &quadrant, FALSE, "",
       "turn off quadrant method, use min. points method."}
    , {0, 'M', &minpt, CMDL_INT, &quadrant, FALSE, " <minimum>",
       "use min. points method. default: 1, implies -Q"}
    , {0, 's', NULL, 0, &surfer, TRUE, "",
       "create surfer ASCII grid instead of GMT grid."}
    , {0, 'S', NULL, 0, &schlauch, TRUE, "",
       "create ASCII stream instead of GMT grid."}
    , {0, 'f', &ptfn, CMDL_CHAR, NULL, 0, " <pt_fname>",
       "point file to use instead of grid points as target points."}
    , {0, 'D', NULL, 0, &dodiff, TRUE, "",
       "additionally print differences, used only with -f."}
    , {0, 'z', NULL, 0, &skipg, TRUE, "",
       "use weight zero for r<0.1, used only with -f."}
  };
  int optcount = (sizeof(myopts)/sizeof(struct mainoptions));
      
  errinit (argv[0]);
  ind_read = process_cmdline_opt (argc, argv, myopts, optcount, 
				  description, nonoptarg);
  if (version_flag) { 
    display_version();
    return EXIT_SUCCESS;
  }
  if (ind_read == -3)  return EXIT_SUCCESS;
  if (ind_read <0 || argc - ind_read < 0){
    print_usage_info(FALSE, myopts, optcount, 
		     description, nonoptarg);
    return EXIT_FAILURE;
  }
#define XOR(a,b) ( !(( (a) && !(b) ) || ( !(a) && (b) )) )
  if (grdrange == NULL || XOR(grdincr==NULL, ptfn==NULL) || grdfn == NULL) {
    fprintf (stderr, "Options -R, [ -I | -f ] and -G must be supplied!\n\n");
    print_usage_info(FALSE, myopts, optcount, 
		     description, nonoptarg);
    return EXIT_FAILURE;
  }
#undef XOR
  if (argc - ind_read == 1) {
    fpxyz = open_text_file (argv[ind_read], "r");
    if (!fpxyz) return EXIT_FAILURE;
  } else {
    if (isatty(STDIN_FILENO)) {
      fprintf (stderr, "You don't wanna type in the data, do you?\n\n");
      print_usage_info(FALSE, myopts, optcount, 
		       description, nonoptarg);
      return EXIT_FAILURE;
    }
    fpxyz = stdin;
  }
  if (R_option (grdrange, &xmin, &xmax, &ymin, &ymax) != 0 ){ 
    fprintf (stderr, "Invalid -R option arguments!\n\n");
    print_usage_info(FALSE, myopts, optcount, 
		     description, nonoptarg);
    return EXIT_FAILURE;
  }
  if (ptfn == NULL) {
    if (I_option (grdincr, &xinc, &yinc) != 0) {
      fprintf (stderr, "Invalid -I option arguments!\n\n");
      print_usage_info(FALSE, myopts, optcount, 
		       description, nonoptarg);
      return EXIT_FAILURE;
    }
  } else {
    fppts = open_text_file (ptfn, "r");
    if (!fppts) return EXIT_FAILURE;
  }
  if (bxnyn != NULL && b_option (bxnyn, &xn, &yn) != 0) {
    fprintf (stderr, "Invalid -b option arguments!\n\n");
    print_usage_info(FALSE, myopts, optcount, 
		     description, nonoptarg);
    return EXIT_FAILURE;
  }
  if (bxnyn == NULL) b_calcoption (radius, xmin, xmax, ymin, ymax, &xn, &yn);
  

  /* where the work is done */
  if (!(fpgrd=open_text_file(grdfn,"w"))) return EXIT_FAILURE;
  fclose (fpgrd);
  
  err = xyz2grd (fpxyz, fppts, grdfn, xmin, xmax, ymin, ymax, xinc, yinc,
		 radius, exponent, eps, xn, yn, surfer, quadrant, minpt, 
		 dodiff, skipg, schlauch, verbose);

  if (fpxyz) fclose (fpxyz);
  if (fppts) fclose (fppts);

  if (err != 0) switch (err) {
  case GGLERR_INOUT:
    fprintf (stderr, "%s: Input/output error.\n",
	     program_invocation_short_name);
    break;
  case GGLERR_NOMEM:
    fprintf (stderr, "%s: Memory allocation error.\n",
	     program_invocation_short_name);
    break;
  case GGLERR_INVAL:
    fprintf (stderr, "%s: Invalid input error.\n",
	     program_invocation_short_name);
    break;
  default:
    fprintf (stderr, "%s: Unknown error %d.\n",
	     program_invocation_short_name, err);
    break;
  }

  return err;
}
