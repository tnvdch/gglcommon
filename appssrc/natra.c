/*
 *      FILE: natra
 *
 *   PROJECT: GGL Common Applications
 *
 *   PURPOSE: 3D helmert coordinate transformation
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
 *            17-05-1997  --  created
 *            30-10-1998  --  all in one file h3d.c + natra.c > natra.c
 *
 *   CAVEATS: Needs the GGL Common Library!
 *
 *  Please compile with compiler flags "-lm -lcommon" !
 */

#define COPYRIGHT "(c) 1997 Geodesy and Geodynamics Lab, ETH Zurich"

/*
 *  We need some standard headers:
 *
 *  the allround headers
 */
#include <errno.h>
#include <stdlib.h>
/*
 *  this one is crucial for file input and output  
 */
#include <stdio.h>
/*
 *  this one to parse the command line
 */
#include <unistd.h>
/*
 *  these to do character and string handling
 */
#include <ctype.h>
#include <string.h>
/*
 *  provide some error handling that is GNU extension
 *  (should be, I do have an older system, so I had to hack it)
 */
#include <math.h>
/*
 *  our own institute's headers
 */
#include "common.h"

/*
 *  This struct holds one point to transform
 */
struct hpoint {
	char name[64];
	double x[3];
} ;

/*
 *  reads coordinates in the format 
 *  point_name      X   Y   Z
 *  %s [max=64]   %lf %lf %lf
 *
 *  point_names appearing twice are taken as control points
 *  the latter is always assumed to be of the target system
 */
int helmert (
	FILE *fpin,
	FILE *fpout,
	FILE *fplog)
{
	extern char *program_invocation_short_name;

	char *oneline;
	struct hpoint cp_source, cp_target, ap;
	double rotmat[3][3];
	double dtrans[3];
	double alpha0=0.0, beta0=0.0, gamma0=0.0;
	double scale[3];
	double so[3] = { 0.0, 0.0, 0.0 };
	double to[3] = { 0.0, 0.0, 0.0 };
	double tmpd;
	int cp = 0;
	int cp2 = 0;
	int i;
	struct matrix *xs=NULL, *t=NULL, *ttmp=NULL, *ttmp2=NULL, *xt=NULL;

	cp_source.name[0] = cp_target.name[0] = ap.name[0] = '%';
	cp_source.name[1] = cp_target.name[1] = ap.name[1] = '\0';

	if (fplog != NULL)
		fprintf (fplog, "Searching controlpoints ...\n");
	for (i=0; i<2; i++) {
		while (read_line (fpin, &oneline) != EOF) {
			if (sscanf (oneline, " %s %lf %lf %lf", ap.name,
				    &ap.x[0], &ap.x[1], &ap.x[2]) < 1) {
				free(oneline);
				continue;
			}
#ifdef WITH_DEBUG
			if (fplog != NULL)
				fprintf (fplog, "i %d / ap %s source %s\n", i, ap.name,
					 cp_source.name);
#endif
			if (ap.name[0] == '$' || ap.name[0] == '#') { /* a comment, lets log it */
				if (fplog != NULL && i == 1)
					fprintf (fplog, "%s", oneline);
				free(oneline); 
				continue; 
			}
			if (ap.name[0] == '%') {
				if (fplog != NULL)
					fprintf (fplog, "Found end of controlpoints %d\n", i);
				free(oneline); 
				break; 
			}
			if (i == 0 && strcmp (cp_source.name, ap.name) == 0) {
				v3_add (so, cp_source.x, so);
				v3_add (to, ap.x, to);
				cp++;
				cp_source.name[0] = '%';
				continue;
			}
			if (i == 1 && strcmp (cp_source.name, ap.name) == 0) {
				int j;
				memcpy (&cp_target, &ap, sizeof(struct hpoint));
				for (j=0;j<3;j++) {
					xs->mat[cp2][j] = cp_source.x[j];
					xt->mat[cp2][j] = cp_target.x[j];
				}
				
				/* we switch left hand to right, be sure to check conversion
				 *  xs->mat[cp2][0] = cp_source.x[1];
				 *  xt->mat[cp2][0] = cp_target.x[0];
				 *  xs->mat[cp2][1] = cp_source.x[0];
				 *  xt->mat[cp2][1] = cp_target.x[1];
				 *  xs->mat[cp2][2] = cp_source.x[2];
				 *  xt->mat[cp2][2] = cp_target.x[2];
				 *  xs->mat[cp2][3] = xt->mat[cp2][3] = 1.0;
				 */
				cp2++;
				if (fplog != NULL) {
					fprintf (fplog, "Found controlpoint: %s\n",
						 cp_source.name);
				}
				cp_source.name[0] = '%';
				continue;
			}
			memcpy (&cp_source, &ap, sizeof(struct hpoint));
		} /* while not EOF */
		rewind (fpin);
		switch (i) {
			int j,k;
		case 0:
			v3_multx (so, 1.0/cp, so);
			v3_multx (to, 1.0/cp, to);
			if (fplog != NULL) {
				fprintf (fplog, "Got %d controlpoints\n", cp);
				fprintf (fplog, "Source Origin: %14.4f %14.4f %14.4f\n",
					 so[0], so[1], so[2]);
				fprintf (fplog, "Target Origin: %14.4f %14.4f %14.4f\n",
					 to[0], to[1], to[2]);
			}
			xs = m_init (cp, 4);
			xt = m_init (cp, 4);
			t = m_init (4, 4);
			ttmp = m_init (4, 4);
			ttmp2 = m_init (4, 4);
			break;
		case 1:
			m_mult (xs, TRUE, xs, FALSE, ttmp);
			m_inv (ttmp);
			m_mult (xs, TRUE, xt, FALSE, ttmp2);
			m_mult (ttmp, FALSE, ttmp2, FALSE, t);
			/* extract translation */
			for (j=0;j<3;j++)
				dtrans[j] = t->mat[3][j];
			/* extract scales */
			for (j=0;j<3;j++)
				scale[j] = sqrt (rotmat[0][j]*rotmat[0][j] +
						 rotmat[1][j]*rotmat[1][j] +
						 rotmat[2][j]*rotmat[2][j]);
			/* orthonormalize */
			rotmat[2][0] = t->mat[0][2] / scale[0];
			rotmat[1][0] = t->mat[0][1] / scale[0];
			rotmat[2][1] = t->mat[1][2] / scale[1];
			rotmat[0][0] = sqrt (1 - rotmat[2][0]*rotmat[2][0] -
					     rotmat[1][0]*rotmat[1][0]);
			rotmat[2][2] = sqrt (1 - rotmat[2][0]*rotmat[2][0] -
					     rotmat[2][1]*rotmat[2][1]);
			rotmat[0][0] = copysign (rotmat[0][0], t->mat[0][0]);
			rotmat[2][2] = copysign (rotmat[2][2], t->mat[2][2]);
			tmpd = 1 - rotmat[2][0]*rotmat[2][0];
			rotmat[0][1] = -(rotmat[2][2]*rotmat[1][0] + 
					 rotmat[2][1]*rotmat[2][0]*rotmat[0][0]) / tmpd;
			rotmat[0][2] =  (rotmat[2][1]*rotmat[1][0] + 
					 rotmat[2][2]*rotmat[2][0]*rotmat[0][0]) / tmpd;
			rotmat[1][1] =  (rotmat[2][2]*rotmat[0][0] + 
					 rotmat[2][1]*rotmat[2][0]*rotmat[1][0]) / tmpd;
			rotmat[1][2] = -(rotmat[2][1]*rotmat[0][0] + 
					 rotmat[2][2]*rotmat[2][0]*rotmat[1][0]) / tmpd;
			if (fplog != NULL)
				fprintf (fplog, "init matrix :\n");
			for (j=0;j<3;j++) {
				for (k=0;k<3;k++) {
					if (fplog != NULL)
						fprintf (fplog, " %14.9f",
							 rotmat[j][k]);
				}
				if (fplog != NULL)
					fprintf (fplog, "\n");	  
			}
			/* extract angles from first approximation */
			/* beta */
			beta0 = asin(rotmat[2][0]) + (rotmat[2][0] < 0.0) ? 2*M_PI : 0.0;
			/* gamma */
			gamma0 = atan(-rotmat[1][0]/rotmat[0][0]);
			if (rotmat[0][0]<0.0) gamma0 += M_PI;
			else 
				if (rotmat[1][0]>0.0) gamma0 += 2*M_PI;
			/* alpha */
			alpha0 = atan(-rotmat[2][1]/rotmat[2][2]);
			if (rotmat[2][2]<0.0) alpha0 += M_PI;
			else 
				if (rotmat[2][1]>0.0) alpha0 += 2*M_PI;
			if (fplog != NULL)
				fprintf (fplog, "init a %14.10g b %14.10g c %14.10g\n",
					 alpha0 *180 / M_PI, beta0 *180/M_PI, gamma0*180/M_PI);
			
			m_free (ttmp);
			m_free (ttmp2);
			break;
		default:
			fprintf (fplog, "fatal: aborting...\n");
			fflush (fplog);
			abort();
		}
	} /* for twice */
	
	
	/*
	 *  now we have got a primary
	 *  we do iterations
	 */
	{
		struct matrix *ma, *aa, *mn, *mdn, *mdi, *mi, *mb, *mw, 
			*baam, *dd;
		double krit = 999.999;
		int itcount = 0;
		int ki, kj;
		double nenner;
		double sigma;
		
		ma = m_init (3, 6);
		aa = m_init (3, 3);
		mb = m_init (3, 6);
		mw = m_init (3, 1);
		mn = m_init (6, 6);
		mdn = m_init (6, 6);
		mi = m_init (6, 1);
		mdi = m_init (6, 1);
		baam = m_init (6, 3);
		dd = m_init (6, 1);
		
		while (krit > 1.0e-10 && itcount++ < 30) {
			if (fplog != NULL) {
				fprintf (fplog, "performing %d iteration...",
					 itcount);
				fflush (fplog);
			}
			/* fix scales */
			scale[0] = scale[1] = scale[2] = 1;
	
			/* fill matrices */
			m_zero (ma);
			m_zero (mb);
			for (ki=0;ki<3;ki++) {
				for (kj=0;kj<3;kj++) {
					ma->mat[ki][kj] = rotmat[ki][kj]*scale[kj];
				}
				ma->mat[ki][ki+3] = -1.0;
				mb->mat[ki][ki+3] = +1.0;
			}
			m_mult (ma, FALSE, ma, TRUE, aa);
			m_inv (aa);
			m_zero (mn);
			m_zero (mi);
			nenner = sqrt (1 - rotmat[2][0]*rotmat[2][0]);
			
			/* operations including controlpoint data */
			for (ki=0;ki<cp;ki++) {
				/* w */
				double vw[3];
				m3v3_mult (rotmat, xs->mat[ki], vw);
				v3_add (vw, dtrans, vw);
				v3_sub (vw, xt->mat[ki], vw);
				for (kj=0;kj<3;kj++) mw->mat[kj][0] = -vw[kj];
				/* B */
				for (kj=0;kj<3;kj++) {
					mb->mat[kj][0] = rotmat[kj][1] * scale[2] * xs->mat[ki][2]
						- rotmat[kj][2] * scale[1] * xs->mat[ki][1];
					mb->mat[kj][1] = -(xt->mat[ki][2] - dtrans[2])*rotmat[kj][0]
						/ nenner;
				}
				mb->mat[2][1] += scale[0]*xs->mat[ki][0] / nenner;
				mb->mat[0][2] =   xt->mat[ki][1] - dtrans[1];
				mb->mat[1][2] = -(xt->mat[ki][0] - dtrans[0]);
				/* N */
				m_mult (mb, TRUE, aa, FALSE, baam);
				m_mult (baam, FALSE, mb, FALSE, mdn);
				m_add (mn, mdn, mn);
				/* I */
				m_mult (baam, FALSE, mw, FALSE, mdi);
				m_add (mi, mdi, mi);
			}
			/* finally solve the whole thing */
			m_inv (mn);
			m_mult (mn, FALSE, mi, FALSE, dd);
			/* iteration kriterion */
			krit = sqrt (dd->mat[3][0]*dd->mat[3][0] +
				     dd->mat[4][0]*dd->mat[4][0] +
				     dd->mat[5][0]*dd->mat[5][0]);
			if (fplog != NULL) fprintf (fplog, " krit: %14.10g", krit);
			/* update approxies */
			for (ki=0;ki<3;ki++)
				dtrans[ki] += dd->mat[ki+3][0];
			alpha0 += dd->mat[0][0];
			alpha0 -= (alpha0>2.0*M_PI) ? 2.0*M_PI : 0.0;
			beta0 += dd->mat[1][0];
			beta0 -= (beta0>2.0*M_PI) ? 2.0*M_PI : 0.0;
			gamma0 += dd->mat[2][0];
			gamma0 -= (gamma0>2.0*M_PI) ? 2.0*M_PI : 0.0;
			if (fplog != NULL)
				fprintf (fplog, " a %14.10g b %14.10g c %14.10g\n",
					 alpha0 *180 / M_PI, beta0 *180/M_PI, gamma0*180/M_PI);
			/* update rotmat */
			rotmat[2][0] = sin (beta0);
			rotmat[1][0] = -cos(beta0)*sin(gamma0);
			rotmat[2][1] = -cos(beta0)*sin(alpha0);
			rotmat[0][0] = cos(beta0)*cos(gamma0);
			rotmat[2][2] = cos(beta0)*cos(alpha0);
			rotmat[0][1] = cos(alpha0)*sin(gamma0) +
				sin(alpha0)*sin(beta0)*cos(gamma0);
			rotmat[0][2] = sin(alpha0)*sin(gamma0) -
				cos(alpha0)*sin(beta0)*cos(gamma0);
			rotmat[1][1] = cos(alpha0)*cos(gamma0) -
				sin(alpha0)*sin(beta0)*sin(gamma0);
			rotmat[1][2] = sin(alpha0)*cos(gamma0) +
				cos(alpha0)*sin(beta0)*sin(gamma0);
			
		} /* while krit */
		
		/* calculate vPv */
		m_zero (mi);
		for (ki=0;ki<cp;ki++) {
			/* w */
			double vw[3];
			m3v3_mult (rotmat, xs->mat[ki], vw);
			v3_add (vw, dtrans, vw);
			v3_sub (vw, xt->mat[ki], vw);
			for (kj=0;kj<3;kj++) mw->mat[kj][0] = -vw[kj];
			
			m_mult (ma, TRUE, aa, FALSE, baam);
			m_mult (baam, FALSE, mw, FALSE, mdi);
			m_add (mi, mdi, mi);
		}
		m_free (mn);
		mn = m_init (1, 1);
		m_mult (mi, TRUE, mi, FALSE, mn);
		sigma = sqrt (mn->mat[0][0]/(cp*3 - 6));
		if (fplog != NULL) fprintf (fplog, "sigma apost: %14.10g\n", sigma);
		
		m_free(ma);
		m_free(aa);
		m_free(mb);
		m_free(mw);
		m_free(mn);
		m_free(mdn);
		m_free(mi);
		m_free(mdi);
		m_free(baam);
		m_free(dd);
	}
	
	{
		int j,k;
		if (fplog != NULL)
			fprintf (fplog, "using matrix :\n");
		for (j=0;j<3;j++) {
			for (k=0;k<3;k++) {
				if (fplog != NULL)
					fprintf (fplog, "    %14.9f",
						 rotmat[j][k]);
			}
			if (fplog != NULL)
				fprintf (fplog, "\n");	  
		}
		fprintf (fplog, "  a %14.10g  b %14.10g  c %14.10g\n",
			 alpha0 *180 / M_PI, beta0 *180/M_PI, gamma0*180/M_PI);
		fprintf (fplog, " dx %14.10g dy %14.10g dz %14.10g\n",
			 dtrans[0], dtrans[1], dtrans[2]);
	}
	
	if (fplog != NULL)
		fprintf (fplog, "Converting all points ...\n");
	while (read_line (fpin, &oneline) != EOF) {
		if (sscanf (oneline, " %s %lf %lf %lf", ap.name,
			    &ap.x[0], &ap.x[1], &ap.x[2]) < 1) {
			fprintf (fpout, "\n");
			continue;
		}
		if (ap.name[0] == '$' || ap.name[0] == '#') { /* a comment, lets print it */
			fprintf (fpout, "%s", oneline);
			free(oneline); 
			continue; 
		}
		if (ap.name[0] == '%') {
			free(oneline);
			fprintf (fpout, "%%\n");
			continue;
		}
		if (strcmp (cp_source.name, ap.name) == 0) {
			fprintf (fpout, " %s %14.4f %14.4f %14.4f checkpoint!\n\n",
				 ap.name, ap.x[0], ap.x[1], ap.x[2]);
			cp_source.name[0] = '%';
		} else {
			memcpy (&cp_source, &ap, sizeof(struct hpoint));
			m3v3_mult (rotmat, cp_source.x, cp_target.x);
			/* you would have to add the scale if.. */
			v3_add (cp_target.x, dtrans, cp_target.x);
			fprintf (fpout, " %s %14.4f %14.4f %14.4f\n", ap.name,
				 cp_target.x[0], cp_target.x[1], cp_target.x[2]);
		}
	} /* while not EOF */
	return cp;
}



/*
 *  Print version information to stderr
 *  You should not need to modify this
 */
void
display_version()
{
   extern char *program_invocation_short_name;
   extern char *common_version;
  fprintf (stderr,
	   "%s version %s\n  %s\n  %s\n",
	   program_invocation_short_name, SRC_VERSION, common_version, COPYRIGHT);
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

/*
 *  Print usage information to stderr
 *  Please modify string to reflect correct usage
 */
void
display_usage (
	       const int verbose
	       )
{
   extern char *program_invocation_short_name;
   fprintf (stderr,
   "usage:\n   %s [-Vhv] { -d || input_file output_file }\n\n" ,
	    program_invocation_short_name);
   if (verbose == 0)
     {
       fprintf (stderr,
		"  -V    Version string\n");
       fprintf (stderr,
		"  -h    This help message\n");
       fprintf (stderr,
		"  -v    verbose mode, tell about conversion progress\n");
       fprintf (stderr,
       "  -d    Run as a pipe, read from stdin, write to stdout\n");
       fprintf (stderr,
       "  <file>   Input file processed\n");
       fprintf (stderr,
 "  point_names appearing twice are taken as control points\n");
       fprintf (stderr,
 "  the latter is always assumed to be of the target system\n\n");
     }
   else
     {
        fprintf (stderr,
		"  -h    Prints more detailed help on each option\n\n");      
     }
}

/*
 *  That's where we do the work 
 */
int
main (
      int argc,        /* Number of arguments         */
      char **argv      /* array of argument strings   */
      )
{
  int
/*    arg_index,         / * index variable for loop       */
    arg_char;          /* argument character to process */
  int use_pipe = 1;
  int verbose = 1;

  extern char *program_invocation_name;        /* from errproc.h */
  extern char *program_invocation_short_name;  /* from errproc.h */

  extern int opterr;      /* from unistd.h: let's do our own cursing */

  errinit (argv[0]);      /* from errproc.h */
  opterr = 0;

  /*
   * Now we process all command line options
   */
  while ((arg_char = getopt (argc, argv, "Vhvd")) != -1)
    switch (arg_char)
      {
      case 'h':             /* -h   the help option */
	display_usage(0);
      case 'V':             /* -v   the version string */
	display_version();
	return EXIT_SUCCESS;
      case 'd':
	use_pipe = 0;
	break;
      case 'v':
	verbose = 0;
	break;
      case '?':             /* oops, an unknown option */
	if (isprint (optopt))
	  fprintf (stderr, "%s: Unknown option `-%c'.\n",
		   program_invocation_short_name, optopt);
	else
	  fprintf (stderr, "%s: Unknown option character `\\x%x'.\n",
		   program_invocation_short_name, optopt);
	display_usage(1);
	return EXIT_FAILURE;
      default:              /* actually, this should never happen */
	abort ();
      }

  /*
   *  Check if we do have all the necessary items
   */
  if (  ((argc - optind) < 2) && use_pipe == 1 )
    {
      display_usage(1);
      return EXIT_SUCCESS;     /* If there were no arguments ... */
    }


  /*
   * Now we parse the non-option arguments
   */
    {
      /* Those variables are only needed locally */
      FILE *fpin, *fpout, *fplog;

      fplog = (verbose == 0) ? stderr : NULL;

      if (use_pipe == 0)
	{
	  fpin = stdin;
	  fpout = stdout;
	}
      else
	{
	  if ( !(fpin = fopen (argv[optind], "r")) )
	    {
	      fprintf (stderr,
		       "%s: Could not open file `%s' for reading: %s\n",
		       program_invocation_short_name,
		       argv[optind],
		       strerror (errno) );
	      return EXIT_FAILURE;
	    }
	  if ( !(fpout = fopen (argv[optind+1], "w")) )
	    {
	      fprintf (stderr,
		       "%s: Could not open file `%s' for writing: %s\n",
		       program_invocation_short_name,
		       argv[optind+1],
		       strerror (errno) );
	      return EXIT_FAILURE;
	    }
	}

      /*
       *  this is where we do all the action!
       *  read all kinds of stuff from file
       */
      helmert (fpin, fpout, fplog);

      /* clean up */
      fclose (fpin);
      fclose (fpout);
    }
  /*
   * We are all done :-)
   * Exit properly!
   */
  return EXIT_SUCCESS;
}





