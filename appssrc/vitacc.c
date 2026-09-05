/*
 *       FILE:  $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/apps/vitacc.c,v $
 *
 *    PROJECT:  GGL Common Applications
 *
 *    PURPOSE:  compute the velocitiy and acceleration using a position time serie
 *
 *  COPYRIGHT:  (c) 1999 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * DISCLAIMER:    This program is free software; you can redistribute
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
 * $Id: vitacc.c,v 1.3 2001/08/31 12:08:55 favey Exp $
 *
 *    AUTHORS:  Geiger Alain  (Pascal version), 
                Peter Yannick (c version)
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            21-02-2000  --  C version created
 *
 * $Log: vitacc.c,v $
 * Revision 1.3  2001/08/31 12:08:55  favey
 * cygwin adaptions
 *
 * Revision 1.2  2000/03/09 16:16:26  ypeter
 * Fully revised. More options, more carefully tested, more error checks.
 *
 * Revision 1.1  2000/02/25 13:01:42  ypeter
 * Initial revision
 *
 *
 *   CAVEATS: Needs the GGL Common Library!
 *
 *  Please compile with compiler flags "-lm -lcommon" !
 */

#define COPYRIGHT "(c) 2000 Geodesy and Geodynamics Lab, ETH Zurich"
#define PROG_VERSION 2.0

/* Input types (default = 1) */
#define GENERIC    1
#define GPSKIN     2

/* Output types (default = 1) */
#define FULL       1
#define POS        2
#define VELO       3
#define ACCEL      4

/* Default values for dim, deg and Time factor !!!Not used in help strings!!! */
#define DIM       10
#define DEG        2
#define TFACTOR    1

/* Dimension for the interpolation */
#define INTERDIM   2
#define INTERDEG   1

/* standard headers */
#include <errno.h>     /* the allround headers */
#include <stdlib.h>    /* the allround headers */
#include <stdio.h>     /* crucial for file input and output   */
#include <ctype.h>     /* character handling */
#include <string.h>    /* string handling */
#include <math.h>      /* mathematics */
#include <unistd.h>    /* to delete a file from disk */
#include <sys/stat.h>  /* to examine file attributes */
#include <time.h>      /* standard time and date */
#include <sys/time.h>  /* extended time and date */
#include <pwd.h>       /* User and group info */

/* our own institute's headers */
#include <common.h>

/* Option parameters are saved in a Param Structure */
struct Param { 
  int version_flag;
  int threeD_flag;    /* Compute the accumulated dist in 3D */
  int dim;            /* Window dimenstion */
  int deg;            /* Degree of polynomial approx */
  int weighted;       /* Use weight in computation */
  char *inpfname;     /* Name of input file */
  char *outfname;     /* Name of output file */
  int inpformat;      /* Input format */
  int outformat;      /* Output format */
  int t_factor;       /* Time factor */
  double dt;          /* Interpolation offset */
  double inc;         /* Interpolation increment */
  int interdim;       /* size of the window for the interpolation */
  int interdeg;       /* Degree of Polynom used for the interpolation */
};

/* --------------------------------------------------------------------------
 *  Print version information to stderr 
 */
void display_version()
{
  extern char *program_invocation_short_name;
  extern char *common_version;
  fprintf (stderr,
	   "%s V%3.2f %s\n  %s\n  %s\n",
	   program_invocation_short_name, PROG_VERSION, SRC_VERSION,
	   common_version, COPYRIGHT);
  fprintf (stderr,
	   "  maintainers: Y. Peter <ypeter@geod.baug.ethz.ch>, ");
  fprintf (stderr,
	   "M. Cocard <cocard@geod.baut.ethz.ch\n\n");
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

/* --------------------------------------------------------------------------
 * Read 1 observation from input file 
 */
int read1obs ( FILE *fp,
	       struct Param *param,
	       struct matrix *data,
	       int mpos,
	       long int *count)
{
  extern char *program_invocation_short_name;
  struct ttime t;
  char line[MAXLINELENGTH];
  int i, readresult, nread, nexpected;
  int push = FALSE;

  if (mpos < 0) {
    mpos = data->row-1;
    push = TRUE;
  }

  while((readresult = read_line_static(fp,line,MAXLINELENGTH)) != EOF) {
    (*count)++;
    i=0; while(isblank(line[i])) i++; /* overread comments and empty lines*/
    if (readresult <= i+1 || strchr("#\n\f", (int)line[i])) continue;

    if (push)
      m_pushup(data);
    switch (param->inpformat){
    case GENERIC:      /* FORMAT: time Xpos Ypos Zpos */
      if (param->weighted)  nexpected = 7;
      else                  nexpected = 4;
      nread = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf",
		     &data->mat[mpos][0], &data->mat[mpos][1],
		     &data->mat[mpos][2], &data->mat[mpos][3],
		     &data->mat[mpos][4], &data->mat[mpos][5],
		     &data->mat[mpos][6]);
      break;
    case GPSKIN:      /* FORMAT: GPSKIN output */
      if (param->weighted)  nexpected = 12;
      else                  nexpected = 9;
      nread = sscanf(line, "%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf",
		     &t.year, &t.month, &t.day, &t.hour, &t.min, &t.sec,
		     &data->mat[mpos][1], &data->mat[mpos][2],
		     &data->mat[mpos][3], &data->mat[mpos][4],
		     &data->mat[mpos][5], &data->mat[mpos][6]);
      if (nread < nexpected)
	break;
      datmjd(&t);
      data->mat[mpos][0] = t.mjd;
      break;
    default: 
      fprintf(stderr, "\n%s: Input format %d is invalid. Aborting!\n\n",
	      program_invocation_short_name, param->inpformat);
      return EOF;
    }
    if (nread < nexpected && readresult != EOF) {
      fprintf(stderr, "%s: Error reading file %s at line %ld\n", 
	      program_invocation_short_name, param->inpfname, *count);
      fprintf(stderr, "\t Element expected: %d  read: %d\n",
	      nexpected, nread);
    } else {
      if (param->weighted)
	for (i=4;i<7;i++)
	  data->mat[mpos][i] = 1.0 / sqr(data->mat[mpos][i]); /* Weigth = 1/sqr(sigma) */
      else
	for (i=4;i<7;i++)
	  data->mat[mpos][i] = 1.0;
      return readresult;
    }
  } /* end while */
  return readresult;
}
/* ----------------------------------------
 * write header as comment to result file 
 */
int write_header ( FILE *fp, 
		   struct Param *param)
{
  extern char *program_invocation_short_name;
  int colnr;
  char *host;
  struct ttime *now;
  struct passwd *user;

  /* Allocate memory for user and time info */
  now=(struct ttime*)malloc(sizeof(struct ttime));
  user=(struct passwd*)malloc(sizeof(struct passwd));
  host=(char*)malloc(255*sizeof(char));

  /* get user, host and time info */
#ifndef OS_CYGWIN
  localtimenow(now);
#endif
  user = getpwuid(getuid());
  gethostname(host, 255);

  /* print file header */
  fprintf(fp,"#\n");
  fprintf(fp,"# program %s (V%3.2f) \n# %s\n",
	  program_invocation_short_name, PROG_VERSION, SRC_VERSION);
  fprintf(fp,"#\n");
#ifndef OS_CYGWIN
  fprintf(fp,"# Run date         : %s\n",strtime3(now));
#endif
  fprintf(fp,"# Run by           : %s@%s\n#\n",user->pw_name, host);
  fprintf(fp,"# Input file       : %s\n",  param->inpfname);
  fprintf(fp,"# Output file      : %s\n",  param->outfname);
  fprintf(fp,"# Data window size : %3d\n", param->dim);
  fprintf(fp,"# Polynom degree   : %3d\n", param->deg);
  fprintf(fp,"# Time factor      : %3d\n", param->t_factor);
  fprintf(fp,"# Weighted input   : %3s\n", (param->weighted)?"yes":"no");
  if (param->inc <= 0)
    fprintf(fp,"# Interpolation    : no\n");
  else
    fprintf(fp,"# Interpolation    : yes (inc= %g offset= %g dim= %d deg= %d)\n",
	    param->inc, param->dt, param->interdim, param->interdeg);

  fprintf(fp,"#\n");

  if (param->deg < 2)
    fprintf(fp,"# CAUTION - Polynom degree < 2  (Can not estimate the vel. and accel.!!)\n");
  else if (param->deg < 3)
    fprintf(fp,"# CAUTION - Polynom degree < 3  (Can not estimate the accel.!!)\n");
  if (param->dim - param->deg < 2 && param->weighted )
    fprintf(fp,"# CAUTION - Degree of Freedom is ZERO. Sigma0 a posteriori set to 1!!\n");

  fprintf(fp,"#\n");

  colnr=1;
  fprintf(fp,"# col %2d      : time\n", colnr++);
  switch (param->outformat){
  case FULL:
    fprintf(fp,"# col %2d - %2d : position (X Y Z)\n", colnr, colnr+2); colnr+=3;
    fprintf(fp,"# col %2d - %2d : velocity (X Y Z)\n", colnr, colnr+2); colnr+=3;
    fprintf(fp,"# col %2d - %2d : acceleration (X Y Z)\n", colnr, colnr+2); colnr+=3;
    fprintf(fp,"# col %2d      : integrated %s distance\n", 
	    colnr++, (param->threeD_flag == TRUE)?"3D":"horizontal");
    fprintf(fp,"# col %2d      : azimuth [deg]\n", colnr++);
    fprintf(fp,"# col %2d      : elevation [deg]\n", colnr++);
    break;
  case POS:
    fprintf(fp,"# col %2d - %2d : position (X Y Z)\n", colnr, colnr+2); colnr+=3;
    fprintf(fp,"# col %2d - %2d : sigma pos (X Y Z)\n", colnr, colnr+2); colnr+=3;
    fprintf(fp,"# col %2d      : integrated %s distance\n", 
	    colnr++, (param->threeD_flag == TRUE)?"3D":"horizontal");
    fprintf(fp,"# col %2d      : azimuth [deg]\n", colnr++);
    fprintf(fp,"# col %2d      : elevation [deg]\n", colnr++);
    break;
  case VELO:
    fprintf(fp,"# col %2d - %2d : velocity (X Y Z)\n", colnr, colnr+2); colnr+=3;
    fprintf(fp,"# col %2d - %2d : sigma vel. (X Y Z)\n", colnr, colnr+2); colnr+=3;
    fprintf(fp,"# col %2d      : velocity amount\n", colnr++);
    break;
  case ACCEL:
    fprintf(fp,"# col %2d - %2d : acceleration (X Y Z)\n", colnr, colnr+2); colnr+=3;
    fprintf(fp,"# col %2d - %2d : sigma accel (X Y Z)\n", colnr, colnr+2); colnr+=3;
    fprintf(fp,"# col %2d      : acceleration amount\n", colnr++);
    break;
  default:
    fprintf(stderr, "\n%s: Output format %d is invalid. Aborting!\n\n",
	      program_invocation_short_name, param->outformat);
      return EXIT_FAILURE;
  } /* end switch */
  fprintf(fp,"#-----------------------------------------\n");
  
  return EXIT_SUCCESS;
}


/* --------------------------------------------------------------------------
 * Write to output file 
 */
int write_result ( FILE *fp,
		   struct Param *param,
		   double time,
		   struct matrix *coeff[3],
		   struct matrix *qxx[3],
		   double s0[3],
		   double dist,
		   double azi,
		   double elev)
{
  int i;
  double abs, fact;
  double pos[3], vel[3], acc[3];

  fact = param->t_factor;

  for (i=0; i<3; i++){
    pos[i] = coeff[i]->mat[0][0];             /* polynom at x=0 (x=T0-T0) */
    vel[i] = coeff[i]->mat[1][0]*fact;        /* 1st derivatif of the polynom at x=0 */
    acc[i] = 2*coeff[i]->mat[2][0]*sqr(fact); /* 2nd derivatif of the polynom at x=0 */
  }
  
  switch (param->outformat) {
  case FULL: /* all values, without sigmas */
    fprintf(fp,
	    "%12.5f  %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f %11.3f",
	    time, pos[0],pos[1],pos[2], vel[0],vel[1],vel[2], acc[0],acc[1],acc[2]);
    fprintf(fp, " %11.3f %11.5f %11.5f", dist, azi*180.0/PI, elev*180.0/PI);
    break;
  case POS: /* position, sigma pos, dist, azi elev */
    fprintf(fp,"%12.5f  %11.3f %11.3f %11.3f %8.3f %8.3f %8.3f %11.3f %11.5f %11.5f",
	    time, pos[0],pos[1],pos[2],
	    s0[0]*sqrt(qxx[0]->mat[0][0]),
	    s0[1]*sqrt(qxx[1]->mat[0][0]),
	    s0[2]*sqrt(qxx[2]->mat[0][0]),
	    dist, azi*180.0/PI, elev*180.0/PI);
    break;
  case VELO: /* velocity only */
    abs=sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);
    fprintf(fp, "%12.5f  %11.3f %11.3f %11.3f %8.3f %8.3f %8.3f %11.3f",
	    time, vel[0],vel[1],vel[2],
	    s0[0]*fact*sqrt(qxx[0]->mat[1][1]), 
	    s0[1]*fact*sqrt(qxx[1]->mat[1][1]), 
	    s0[2]*fact*sqrt(qxx[2]->mat[1][1]),
	    abs);
    break;
  case ACCEL: /* acceleration only */
    abs=sqrt(acc[0]*acc[0]+acc[1]*acc[1]+acc[2]*acc[2]);
    fprintf(fp, "%12.5f  %11.3f %11.3f %11.3f %8.3f %8.3f %8.3f %11.3f",
	    time, acc[0],acc[1],acc[2], 
	    s0[0]*fact*2*sqrt(qxx[0]->mat[2][2]), 
	    s0[1]*fact*2*sqrt(qxx[1]->mat[2][2]), 
	    s0[2]*fact*2*sqrt(qxx[2]->mat[2][2]),
	    abs);
    break;
  default:
    fprintf(stderr, "\n%s: Output format %d is invalid. Aborting!\n\n",
	      program_invocation_short_name, param->outformat);
      return EXIT_FAILURE;
  }
  fprintf(fp, "\n");

  return EXIT_SUCCESS;
}

/* --------------------------------------------------------------------------
 * Interpolate 
 */
int interpolate(struct matrix *data,
		struct matrix *source, double time, int deg)
{
  int i;
  struct matrix *coeff;
  struct matrix *Qxx;
  double s0;

  coeff = m_init(source->row, 1);
  Qxx   = m_init(source->row, source->row);

  m_pushup(data);
  data->mat[data->row-1][0]=time;
  for (i=0;i<3;i++) {
    m_polynom_approx(source, 0, i+1, i+4, -1, -1,
		     /* obs 1,2,3, weight 4,5,6, model and res are not saved */
		     0, source->row-1, deg+1, time, NULL, coeff, Qxx, &s0);    
    data->mat[data->row-1][i+1] = coeff->mat[0][0];
    data->mat[data->row-1][i+4] = 1.0/Qxx->mat[0][0];
  }
  m_free(coeff);
  m_free(Qxx);

  return EXIT_SUCCESS;
}

/* --------------------------------------------------------------------------*/
int main ( int argc, char **argv )
{
  FILE *fpin;
  FILE *fpout;

  struct Param param;            /* contains all option parameters */
  struct matrix *last;           /* latest observations read from file */
  struct matrix *data;           /* Matrix dim x 10: Data buffer (T X Y Z Wx Wy Wz Mx My Mz */
  struct matrix *coeffmat[3];    /* polynom coefficient for X, Y, and Z */
  struct matrix *Qxx[3];         /* Covariance of the coeff. for X, Y and Z */
  double s0[3];                  /* Sigma a posteriori for X, Y, Z coeff */
  double oldpos[3];              /* Position at time tn-1 */
  double T;                      /* Time in the center of the data windows */
  double dist = 0;               /* Integrated distance (2D or 3D) */
  double azi = 0, elev = 0;      /* actual azimuth and elevation */
  double Tint = 0;               /* Time for the interpolation */

  /* Help variables */
  double dx[3];               
  long int linecount = 0;
  int firstrun = TRUE;
  int i, err = 1;
  int weight_col = 4;      /* First column with weights in data matrix */
  int ind_read;		   /* Index of first non option argument */

/* let's define our command line options */
  struct mainoptions myopts[] =	   /* from cmdline.h */
  {
    {0,'V', NULL, 0, &param.version_flag, TRUE, "",
     "Version string."},
    {0,'d', &param.dim, CMDL_INT, NULL, 0, " <dimension>",
     "window dimension (default = 10) (#obs. used for fitting the polynom)"},
    {0,'g', &param.deg, CMDL_INT, NULL, 0, " <deg>", 
     "polynomial degree (default = 2)"},
    {0,'f', &param.t_factor, CMDL_INT, NULL, 0, " <factor>",
     "multiplication factor for the time used\n        by the calculaton of the vel. and accel. (default = 1)"},
    {0,'3', NULL, 0, &param.threeD_flag, TRUE, "",
     "integrate the distance in 3D (default = 2D) "},
    {0,'t', &param.inc, CMDL_DOUBLE, NULL, 0, " <xinc>",
     "resample using <xinc> as sampling interval"},
    {0,'T', &param.dt, CMDL_DOUBLE, NULL, 0, " <offset>",
     "offset to apply to <xinc> for the first observation"},
    {0,'W', NULL, 0, &param.weighted, TRUE, "",
     "compute weighted polynomial approx (require sigmas as input)"},
    {0,'D', &param.interdim, CMDL_INT, NULL, 0, " <dimension>",
     "window dimension for the interpolation (default = 2)"},
    {0,'G', &param.interdeg, CMDL_INT, NULL, 0, " <deg>", 
     "polynomial degree for the interpolation (default = 1 (linear))"},
    {0,'i', &param.inpfname, CMDL_CHAR, NULL, 0, " <inpname>",
    "input file name (default = stdin)"},
    {0,'o', &param.outfname, CMDL_CHAR, NULL, 0, " <outname>",
    "output file name (default = stdout)"},
    {0,'I', &param.inpformat, CMDL_INT, NULL, 0, " <format>",
     "Input format (default = 1):\n\
\t   1 = GENERIC: time X Y Z (sX sY sZ)\n\
\t   2 = GPSKIN output: YYYY MM DD hh mm ss.s X Y Z (sX sY sZ)"},
    {0,'O', &param.outformat, CMDL_INT, NULL, 0, " <format>",
     "Output format (default = 1) :\n\
\t   1 = time, posX/Y/Z, velX/Y/Z, accX/Y/Z, dist, azi, elev\n\
\t   2 = time, posX/Y/Z, sigmaX/Y/Z dist, azi, elev\n\
\t   3 = time, velX/Y/Z, sigmaX/Y/Z velabs\n\
\t   4 = time, accX/Y/Z, sigmaX/Y/Z accabs"}
  };
   int optcount = (sizeof(myopts)/sizeof(struct mainoptions));

/* Description of what the program does */
   char description[] = "Compute vel. and accel. from pos. time series";
/* Information about our non option arguments */
   char nonoptarg[] = " [inpname]";
   
/* proper error cursing with our own command name */
   extern char *program_invocation_name;        /* from errproc.h */
   extern char *program_invocation_short_name;  /* from errproc.h */
   errinit (argv[0]);                           /* from errproc.h */
   
/* Initialize param and set default 
 * settings for command line options  
 */
   memset(&param, 0, sizeof(struct Param));
   param.version_flag = FALSE;
   param.threeD_flag = FALSE;
   param.weighted = FALSE;
   param.dim = DIM;
   param.deg = DEG;
   param.inpformat = GENERIC;
   param.outformat = FULL;
   param.t_factor = TFACTOR;
   param.interdim = INTERDIM;
   param.interdeg = INTERDEG;

/* Now we process all command line options */
   ind_read = process_cmdline_opt (argc, argv, myopts, optcount,
                                   description, nonoptarg);

/* let's evaluate if the Version flag was set */
   if (param.version_flag) { 
     display_version();
     return EXIT_SUCCESS;
   }

/* when we got -3, the user called the -h help command */
   if (ind_read == -3)
     return EXIT_SUCCESS;

/* Set weight_col to a negativ value 
 * if we do not use a weighted input */
   if (! param.weighted)
     weight_col = -100;

/* Check the degree of freedoms */
   if (param.dim - param.deg < 1){
     fprintf(stderr,"\n%s: Degree of freedom is smaller than zero!\n",
	     program_invocation_short_name);
     fprintf(stderr,"\t Reduce the polynom degree (-g) or enlarge the window size (-d)!\n\n");
     return EXIT_FAILURE;
   }
   if (param.interdim - param.interdeg < 1) {
     fprintf(stderr,"\n%s: Degree of freedom for interpolation is smaller than zero!\n",
	     program_invocation_short_name);
     fprintf(stderr,"\t Reduce the polynom degree (-G) or enlarge the window size (-D)!\n\n");
     return EXIT_FAILURE;
   }

/* Check the input and output formats */
   if (param.inpformat < 1 || param.inpformat > 2){
     fprintf(stderr, "\n%s: input format %d is invalid. Aborting!\n\n",
	     program_invocation_short_name, param.inpformat);
     return EXIT_FAILURE;
   }
   if (param.outformat < 1 || param.outformat > 4){
     fprintf(stderr, "\n%s: output format %d is invalid. Aborting!\n\n",
	     program_invocation_short_name, param.outformat);
     return EXIT_FAILURE;
   }
     
/* check if we have a non optional arguments => inpname */
   if (argc - ind_read > 0 ) {
     if (param.inpfname != NULL) {
       fprintf (stderr, 
		"\n%s: non optional argument is the inpname. It exclude the option -i!\n\n",
		program_invocation_short_name);
       print_usage_info(FALSE, myopts, optcount, 
                        description, nonoptarg);
       return EXIT_FAILURE;
     } else {
       param.inpfname = strdup(argv[ind_read]);
     }
   }
/* Open input file (inpname or stdin) */
   if (param.inpfname == NULL) {
     fpin = stdin;
     param.inpfname = strdup("(stdin)");
   } else 
     if (!(fpin=open_text_file(param.inpfname, "r"))) 
       return EXIT_FAILURE;

/* Open output file (outname or stdin) */
   if (param.outfname == NULL) {
     fpout = stdout;
     param.outfname = strdup("(stdout)");
   } else
     if (!(fpout=open_text_file(param.outfname, "w"))) 
       return EXIT_FAILURE;

/* Write header to output file */
   if ((err = write_header(fpout, &param)) != EXIT_SUCCESS)
     return err;

/* initialize matrices */
   data = m_init(param.dim,7);
   last = m_init(param.interdim,7);
   for (i=0; i<3;i++){
     /* We need minimum a 3x1 mat and param.dim can be 1 */
     coeffmat[i] = m_init(param.dim+2, 1); 
     Qxx[i]      = m_init(param.dim+2, param.dim+2);
     coeffmat[i]->row = Qxx[i]->row = Qxx[i]->col = param.dim;
   }

/* Fill in the data matrix with initial values */
   if (param.inc != 0) { /* Interpolation */
     for (i=0;i<last->row;i++)
       if ((err = read1obs(fpin, &param, last, i, &linecount)) == EOF)
       return EXIT_SUCCESS;
     Tint = rint((last->mat[0][0]-param.dt)/param.inc)*param.inc + param.dt;
     interpolate(data, last, Tint, param.interdeg);
     for (i=1;i<data->row;i++){
       Tint += param.inc;
       while (last->mat[last->row/2][0] < Tint) {
	 if ((err = read1obs(fpin, &param, last, -1, &linecount)) == EOF)
	   return EXIT_SUCCESS;
       }
       interpolate(data, last, Tint, param.interdeg);
     }
   } else { /* No interpolation */
     for (i=0; i<data->row; i++) 
       if ((err = read1obs(fpin, &param, data, i, &linecount)) == EOF)
	 return EXIT_SUCCESS;
   }
     
/* where the work is done */
   while (err != EOF) {
     T = data->mat[data->row/2][0]; /* Time in the middle of the data window */
     for (i=0; i<3; i++) {
       m_polynom_approx(data, 0, i+1, i+weight_col, -1, -1,
			/* weight 4,5,6, model and res are not saved */
			0, data->row-1, param.deg+1, T, NULL,
			coeffmat[i], Qxx[i], &s0[i]);
       dx[i] = coeffmat[i]->mat[0][0] - oldpos[i];
       oldpos[i] = coeffmat[i]->mat[0][0];
     }

     /* Compute the dist, azi and elev excepted for the first run */
     if (firstrun)
       firstrun = FALSE;
     else { 
       if (param.threeD_flag) 
         dist = dist+sqrt(sqr(dx[0]) + sqr(dx[1] + sqr(dx[2])));
       else 
         dist = dist+sqrt(sqr(dx[0]) + sqr(dx[1]));
       elev = atan2(sqrt(sqr(dx[0]) + sqr(dx[1])), dx[2]);
       azi = PI/2.0 - atan2(dx[0], dx[1]);
       if (azi < 0) azi += 2.0*PI;
     }

     if (param.dim - param.deg < 2 && param.weighted ) /* If n-u = 0, s0 is undefined */
       for (i=0; i<3; i++) s0[i]=1;
     /* Write out the results */
     err = write_result(fpout, &param, T, coeffmat, Qxx, s0, dist, azi, elev);
     if (err != EXIT_SUCCESS)
       return err;

     /* push a new value in the data stack */
     if (param.inc != 0) { /* Interpolation */
       Tint += param.inc;
       while(last->mat[last->row/2][0] < Tint && err != EOF) {
	 err = read1obs(fpin, &param, last, -1, &linecount);
       }
       /* Also interpolate Trailing part of the data */
       if (err == EOF && Tint < last->mat[last->row-1][0])
	 err = FALSE;
       interpolate(data, last, Tint, param.interdeg);
     } else { /* No interpolation */
       err = read1obs(fpin, &param, data, -1, &linecount);
     }
   } /* End while loop */
   
   /* Free memory and close files */
   m_free(data);
   m_free(last);
   for (i=0;i<3;i++) {
     m_free(coeffmat[i]);
     m_free(Qxx[i]);
   }
   fclose(fpin);
   fclose(fpout);
   
   return EXIT_SUCCESS;
}
/* End of file */
