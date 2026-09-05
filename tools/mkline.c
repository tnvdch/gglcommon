#define COPYRIGHT "(c) 2000 Geodesy and Geodynamics Lab, ETH Zurich"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <common.h>

void display_version()
{
        extern char *program_invocation_short_name;
        extern char *common_version;
        fprintf (stderr,
                 "%s version %s\n  %s\n  %s\n",
                 program_invocation_short_name, SRC_VERSION, 
                 common_version, COPYRIGHT);
        fprintf (stderr,
                 "  maintainers: E. Favey <favey@geod.ethz.ch>\n\n");
}

int main ( int argc, char **argv)
{
  int ind_read;           /* index of first non option argument */
  
  int verbose = FALSE;
  int pr_version = FALSE;
  double stepping = 25.0;
  double p1[3], p2[3], pt[3];
  double alpha=0.0, de, dn, dh, dist=0.0, vs, hstep=0.0, vstep=0.0;
  int err;
  
  struct mainoptions myopts[] =                   /* from cmdline.h */
  {
    {0,'V', NULL, 0, &pr_version, TRUE, "",
     "Version string."},
    {0,'s', &stepping, CMDL_DOUBLE, NULL, 0, " <stepping>",
     "Every <stepping> along value we want a point. default=25.0"},
    {0,'v', NULL, 0, &verbose, TRUE, "", 
     "Print information while calculating."}
  };
  int optcount = (sizeof(myopts)/sizeof(struct mainoptions));
  char description[] = "create multiple points within two waypoints";
  char nonoptarg[] = " ";
  extern char *program_invocation_name;           /* from errproc.h */
  extern char *program_invocation_short_name;     /* from errproc.h */
  
  errinit (argv[0]);                              /* from errproc.h */
  ind_read = process_cmdline_opt (argc, argv, myopts, optcount, 
				  description, nonoptarg);
  if (pr_version) { display_version(); return EXIT_SUCCESS; }
  if (ind_read == -3) return EXIT_SUCCESS;
  if (ind_read < 0 || isatty(STDIN_FILENO)) {
    print_usage_info(FALSE, myopts, optcount, 
		     description, nonoptarg);
    return EXIT_FAILURE;
  }

  err = scanf (" %lf %lf %lf %lf %lf %lf", 
		&p1[0], &p1[1], &p1[2], &p2[0], &p2[1], &p2[2]);
  if (err < 6) {
    fprintf (stderr, "%s: invalid or no input in stdin\n",
	     program_invocation_short_name);
    return EXIT_FAILURE;
  }
  de = p2[0]-p1[0];
  dn = p2[1]-p1[1];
  dh = p2[2]-p1[2];
  alpha = atan2 (dn, de);
  dist = sqrt(dn*dn + de*de);
  if (dist < NEARLYZERO) {
    fprintf (stderr, "%s: both waypoints too close to each other\n",
	     program_invocation_short_name);
    return EXIT_FAILURE;
  }
  vs = stepping * dh / dist;
  printf ("%12.3f %12.3f %12.3f\n", p1[0], p1[1], p1[2]);
  do {
    vstep += vs; hstep += stepping;
    pt[0]=p1[0] + hstep * cos(alpha);
    pt[1]=p1[1] + hstep * sin(alpha);
    pt[2]=p1[2] + vstep;
/*     printf ("%12.3f %12.3f ", fabs(pt[1]-p2[1]), fabs(pt[0]-p2[0])); */
    printf ("%12.3f %12.3f %12.3f\n", pt[0], pt[1], pt[2]);
  } while (!((fabs(pt[1]-p2[1]) < 1.0*stepping) &&
	   (fabs(pt[0]-p2[0]) < 1.0*stepping)));
  printf ("%12.3f %12.3f %12.3f\n", p2[0], p2[1], p2[2]);

  return EXIT_SUCCESS;
}
