/**
 **        FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/trafo.c,v $
 **
 ** DESCRIPTION: Functions to transform coordinates
 **              from the bessel ellipsoid to the
 **              Swiss military coordinates and the
 **              vice versa.
 **              Constants and ellipsoid data for calculations
 **              with the swiss projection system.
 **
 **   COPYRIGHT: Copyright (C) 1996 by Etienne Favey
 ** $Id: trafo.c,v 1.10 1999/06/02 16:13:19 cocard Exp favey $
 **
 **  DISCLAIMER:   This program is free software; you can redistribute
 **              it and/or modify it under the terms of the GNU General
 **              Public License as published by the Free Software
 **              Foundation; either version 2 of the License, or
 **              (at your option) any later version.
 **                This program is distributed in the hope that it will
 **              be useful, but WITHOUT ANY WARRANTY; without even the
 **              implied warranty of MERCHANTABILITY or FITNESS FOR A
 **              PARTICULAR PURPOSE. See the GNU General Public License
 **              for more details.
 **                You should have received a copy of the GNU General
 **              Public License along with this program; if not, write
 **              to the Free Software Foundation, Inc., 675 Mass Ave,
 **              Cambridge, MA 02139, USA.
 **
 **      AUTHOR: Etienne Favey <eti@favey.lugs.ch>
 **              Zurich, Switzerland
 **
 **     CHANGES: dd-mm-yyyy -- description
 **              27-08-1996 -- first release (Version 1.01)
 **              03-12-1996 -- included xyz_to_blh, blh_to_xyz
 **              18-12-1996 -- bugfix xyz2blh iteration criterion
 **              10-04-1997 -- adjusted to use with GGL coordtr.c
 **                            taken from public version 1.08,
 **                            renamed blh to llh, introduced
 **                            struct ellipsoid and vectors
 ** $Log: /GPS/SOFTWARE/common_devel/libsrc/trafo.c $
 * 
 * 1     31.08.01 14:51 Etienne
 ** Revision 1.10  1999/06/02 16:13:19  cocard
 ** *** empty log message ***
 **
 ** Revision 1.9  1999/01/26 17:19:19  favey
 ** moved out of swiss.h
 **
 ** Revision 1.8  1999/01/08 10:54:43  cocard
 ** *** empty log message ***
 **
 ** Revision 1.7  1999/01/07 12:57:32  cocard
 ** conversion cartesian - spherical added
 **
 ** Revision 1.6  1998/12/18 10:49:50  favey
 ** no change
 **
 ** Revision 1.5  1998/11/16 17:58:06  favey
 ** improved pattern matching, away from strstr to regex
 **
 ** Revision 1.4  1998/10/30 14:22:07  favey
 ** changed config.h to global.h
 **
 ** Revision 1.3  1998/09/28 15:34:11  cocard
 ** adapted to new gridop.h
 **
 ** Revision 1.2  1998/05/15 09:54:11  favey
 ** *** empty log message ***
 **
 ** Revision 1.1  1998/04/22 12:35:28  favey
 ** Initial revision
 **
 ** Revision 1.4  1998/01/26 13:07:43  favey
 ** fixed bug calculating cc
 **
 ** Revision 1.3  1998/01/22 10:55:12  favey
 ** changed iteration criterion in xyz_to_llh
 **
 ** Revision 1.2  1998/01/20 10:16:34  favey
 ** split swiss.c into trafo.c and swiss.c
 **
 ** Revision 1.1  1998/01/20 08:50:40  favey
 ** Initial revision
 **
 ** Revision 0.4  1998/01/15 17:22:24  favey
 ** new functions select_[<ellipsoid><frameconv><formtype>]
 **
 **
 **/
/*
 *  the standard C headers
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <regex.h>

/*
 *  the GGL Common Library headers
 */
#include "global.h"
#include "errproc.h"
#include "inout.h"
#include "mat3.h"
#include "gridop.h"
#include "angles.h"
#include "trafo.h"

static void compex_error (
			  const int errcode,
			  regex_t *compiled,
			  const char *routine)
{
#ifdef WITH_ERRORS
  size_t length = regerror (errcode, compiled, NULL, 0);
  char *buffer = (char *) malloc (length*sizeof(char));
  regerror (errcode, compiled, buffer, length);
  puterr ("trafo", routine);
  fprintf (stderr, "regex: %s\n", buffer);
  free (buffer);
#endif
}

/*
 * conversion of cartesian coordinates to spherical
 */
void cartesian_to_spherical(double xyz[3], double *lat, double *lon, double *radius)
{
  *radius=v3_length(xyz);
  *lon=atan2(xyz[1],xyz[0]);
  *lat=asin(xyz[2]/(*radius));
}

/*
 * conversion of spherical to cartesian coordinates
 */
void spherical_to_cartesian(double lat, double lon, double radius, double xyz[3])
{
  xyz[0]=radius*cos(lat)*cos(lon);
  xyz[1]=radius*cos(lat)*sin(lon);
  xyz[2]=radius*sin(lat);
}

/*
 *  geocentric to ellipsoidal
 *  parameters have to be passed and are returned in rad
 */
int xyz_to_llh (
		const struct ellipsoid ell,
		const double xyz[3],
		double *lat,
		double *lon,
		double *h)
{
    double a, ex2, b;
  double n,p,fac;
  double h0,b0;
  a = ell.a;
  ex2 = ell.e2;
  int niter = MAX_ITERATIONS;

  /* The calculation fails near the poles. A different approach is used within
   * 0.5mm from the poles. */
  if (xyz[0] < 5e-4 && xyz[1] < 5e-4)
  {
      *lat = sign(xyz[2]) * M_PI / 2;
      *lon = 0; /* FIXME: or rather: undefined (nan)? */
      b = sqrt(a*a*(1-ex2));
      *h = ( xyz[2] - b * sign(xyz[2]) ) * sign(xyz[2]);

#ifdef WITH_ERRORS
      puterr ("trafo", "xyz_to_llh");
      fprintf (stderr, "Warning: longitude is not defined at the poles.\n");
#endif
    return EXIT_SUCCESS;
  }

  *lon = atan2(xyz[1],xyz[0]);
  p = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
  *lat = atan(xyz[2]/(1-ex2)/p);
  *h = 0.0;

  do
    {
      h0 = *h;
      b0 = *lat;
      n = a/sqrt(1-ex2*(sin(*lat)*sin(*lat)));
      fac = 1-ex2*n/(n+*h);
      *lat = atan(xyz[2]/p/fac);
      *h   = p/cos(*lat)-n;
    }
  while ( ((fabs(h0-*h)>1.0e-6) || (fabs(b0-*lat)>1.0e-11)) && --niter );

#ifdef WITH_ERRORS
  if (niter == 0)
	{
	  puterr ("trafo", "xyz_to_llh");
	  fprintf (stderr, "%f/%f/%f did not converge: d_height=%f, d_lat=%f iter=%i\n",
			   xyz[0], xyz[1], xyz[2], fabs(h0-*h), fabs(b0-*lat), niter);
	}
#endif

  return EXIT_SUCCESS;
}

/*
 *  ellipsoidal to geocentric
 *  parameters have to be passed and are returned in rad
 */
int llh_to_xyz (
		const struct ellipsoid ell,
		const double lat,
		const double lon,
		const double h,
		double xyz[3])
{
  double a, ex2;
  double n;
  a = ell.a;
  ex2 = ell.e2;
  
  n = a/sqrt(1-ex2*(sin(lat)*sin(lat)));
  xyz[0] =(n+h)*cos(lat)*cos(lon);
  xyz[1] =(n+h)*cos(lat)*sin(lon);
  xyz[2] =(n*(1-ex2)+h)*sin(lat);
  
  return EXIT_SUCCESS;
}


/*
 *  transform arbitrary frame given source epoch
 */
static struct frameconv * frame_init (
				      const double epoch,
				      struct frameread *fcparam)
{
  struct frameconv *newf;
  double dt, scale, rot[3];
  double tmpmat[3][3];
  int i,err;
  
  err = 0;
  if (fcparam == NULL) return NULL;
  newf = (struct frameconv *) malloc (sizeof(struct frameconv));
  
  if (epoch > 1.0) 
    {
      if (fcparam->ref_epoch > 1.0) 
	{
	  dt = epoch - fcparam->ref_epoch;
	}
      else 
	{
#ifdef WITH_ERRORS
	  char *inoutep;
	  
	  if (newf->fwd > 0) inoutep = strdup ("input");
	  else inoutep = strdup ("output");
	  fprintf (stderr, "GGL Common Library -- swiss\n");
	  fprintf (stderr, "frame_init: %s epoch specified but not needed, ignored\n", inoutep);
	  free (inoutep);
#endif
	  dt = 0.0;
	}
    } 
  else 
    {
      if (fcparam->ref_epoch < 1.0) 
	{
	  dt = 0.0;
	} 
      else 
	{
#ifdef WITH_ERRORS
	  char *inoutep;
	  
	  if (newf->fwd > 0) inoutep = strdup ("input");
	  else inoutep = strdup ("output");
	  fprintf (stderr, "GGL Common Library -- swiss\n");
	  fprintf (stderr,
		   "frame_init: no %s epoch specified but needed\n",
		   inoutep);
	  free (inoutep);
#endif
	  free (newf);
	  return NULL;
	}
    }
  scale = fcparam->scale + dt*fcparam->dscale;
  for (i=0;i<3;i++) 
    {
      newf->dx[i] = fcparam->dx[i] + dt*fcparam->ddx[i];
      rot[i] = fcparam->rot[i] + dt*fcparam->drot[i];
    }
  newf->fwd = fcparam->fwd;
  
  m3_unit (newf->rm);
  for (i=0;i<3;i++) 
    {
      rot_getmat (i+1, rot[i], tmpmat);
      m3m3_mult (tmpmat, newf->rm, newf->rm);
    }
  if (newf->fwd < 0) m3_trans (newf->rm, newf->rm);

  return newf;
}


int frame_conv (
		struct frameconv *fcp,
		double source[3],
		double target[3])
{
  double tmp[3];
  
  if (fcp == NULL) return GGLERR_INVAL;
  if (fcp->fwd > 0) 
    {
      /* fcp->rm is not inverted! */
      m3v3_mult (fcp->rm, source, tmp);
      v3_add (fcp->dx, tmp, target);
    } 
  else 
    {
      /* fcp->rm is inverted by frame_init */
      v3_sub (source, fcp->dx, tmp);
      m3v3_mult (fcp->rm, tmp, target);
    }
  return 0;
}


/*
 *  select ellipsoid from configuration file
 */
struct ellipsoid * select_ellipsoid (
				     const char *ellstr)
{
  struct ellipsoid *theell;
  char *pattern, *line, *conffile;
  regex_t *compex;
  FILE *fpconf;
  int err, lcount;
  double a, e2;
  
  /* determine and read config file */
  conffile = getenv ("TRAFOCONF");
  if (conffile == NULL) 
    {
#ifdef WITH_ERRORS
      puterr ("trafo", "select_ellipsoid");
      fprintf (stderr, "environment variable TRAFOCONF not set.\n");
#endif
      return NULL;
    }
  fpconf = open_text_file (conffile, "r");
  if (!fpconf) return NULL;
  
  /* initialize regex pattern */
  pattern = (char *) malloc (sizeof(char)*(20+strlen(ellstr)));
  if (pattern == NULL)
  {
    fclose(fpconf);
	return NULL;
  }

  compex = (regex_t *) malloc (sizeof(regex_t));
  if (compex == NULL)
  {
    fclose(fpconf);
	free(pattern);
	return NULL;
  }
  sprintf(pattern, "^ellipsoid  *%s ", ellstr);
  if ((err = regcomp (compex, pattern, REG_ICASE))!=0) 
    {
      compex_error (err, compex, "select_ellipsoid");
      fclose(fpconf);
	  regfree(compex);
      free(pattern);
      return NULL;
    }
  theell = NULL;
  lcount = 0;
  
  while ( (err = read_line (fpconf, &line)) != EOF ) 
    {
      if (line == NULL) break;
      lcount++;
      
      err = regexec (compex, line, 0, NULL, 0);
      if (err == REG_NOMATCH) 
	{
	  free (line);
	  continue;
	}
    if (err != 0) 
	{
	  compex_error (err, compex, "select_ellipsoid");
	  free(line);
	  free (pattern);
      regfree (compex);
  	  fclose (fpconf);
	  return NULL;
	}
      err = sscanf (line, "%*s %*s %lf %lf", &a, &e2);
      if (err < 2) 
	{
#ifdef WITH_ERRORS
	  puterr ("trafo", "select_ellipsoid");
	  fprintf (stderr, "invalid config file %s in line %d\n",
		   conffile, lcount);
#endif
	  free (line);
	  continue;
	}
      free (line);
      theell = (struct ellipsoid *) malloc (sizeof (struct ellipsoid));
      theell->a = a;
      theell->e2 = e2;
      break;
    }
  if (theell == NULL) 
    {
#ifdef WITH_ERRORS
      puterr ("trafo", "select_ellipsoid");
      fprintf (stderr, "ellipsoid %s not found in %s\n",
	       ellstr, conffile);
#endif
    }
  
  free (pattern);
  regfree (compex);
  fclose (fpconf);
  return theell;
}


/*
 *  select frame conversion parameters from file
 */
struct frameconv * select_frameconv (
				     const char *source,
				     const double inepoch,
				     const char *target,
				     const double outepoch)
{
  struct frameread *fc;
  struct frameconv *ftr;
  regex_t *compex[4];
  char *pattern[4], *line, *conffile;
  FILE *fpconf;
  int i, err, lcount, set;
  
  /* determine and read config file */
  conffile = getenv ("TRAFOCONF");
  if (conffile == NULL) 
    {
#ifdef WITH_ERRORS
      puterr ("trafo", "select_frameconv");
      fprintf (stderr, "environment variable TRAFOCONF not set.\n");
#endif
      return NULL;
    }
  fpconf = open_text_file (conffile, "r");
  if (!fpconf) return NULL;

  /* initialize regex pattern */
  for (i=0;i<4;i++) 
    {
      pattern[i]=(char *)malloc(sizeof(char)*
				(20+strlen(source)+strlen(target)));
      compex[i] = (regex_t *) malloc (sizeof(regex_t));
    }
  sprintf(pattern[0], "^frame  *%s  *%s ", source, target);
  sprintf(pattern[1], "^rate  *%s  *%s ", source, target);
  sprintf(pattern[2], "^frame  *%s  *%s ", target, source);
  sprintf(pattern[3], "^rate  *%s  *%s ", target, source);
  for (i=0;i<4;i++) 
    if ((err = regcomp (compex[i], pattern[i], REG_ICASE))!=0) {
      compex_error (err, compex[i], "select_frameconv");
      return NULL;
    }
  fc = (struct frameread *) malloc (sizeof(struct frameread));
  memset (fc, 0, sizeof(struct frameread));
  ftr = NULL;
  set = 0;
  lcount = 0;
  
  while ( (err = read_line (fpconf, &line)) != EOF ) 
    {
      double dx[3], scale, rot[3];
      struct Angle ac[3];
      double fepoch;
      int angtyp;
      char element[256];
      int israte = FALSE;
      int fwd = 0;
      if (line == NULL) break;
      lcount++;
      i = 0;
      do 
	{
	  err = regexec (compex[i], line, 0, NULL, 0);
	  if (err != 0 && err != REG_NOMATCH) 
	    {
	      compex_error (err, compex[i], "select_frameconv");
	      return NULL;
	    }
	} while (err != 0 && ++i<4); /* i will NOT increment on match */
      
      if (err != 0) { free (line); continue; }
      switch (i) 
	{
	case 0: israte=FALSE; fwd= 1; break;
	case 1: israte=TRUE;  fwd= 1; break;
	case 2: israte=FALSE; fwd=-1; break;
	case 3: israte=TRUE;  fwd=-1; break;
	default: puterr ("trafo","select_frameconv");
	  fprintf (stderr, "internal error: i=%d",i);
	  return NULL;
	}
#ifdef LIB_DEBUG
      fprintf (stderr, "%s fwd=%d rate=%d i=%d\n", line, fwd, israte, i);
#endif
      if (israte) 
	{
	  err = sscanf (line, 
			"%*s %*s %*s %lf %lf %lf %lf %s %lf %lf %lf",
			&dx[0], &dx[1], &dx[2], &scale,
			element, &rot[0], &rot[1], &rot[2]);
	} 
      else 
	{
	  err = sscanf (line, 
			"%*s %*s %*s %lf %lf %lf %lf %s %lf %lf %lf %lf",
			&dx[0], &dx[1], &dx[2], &scale,
			element, &rot[0], &rot[1], &rot[2], &fepoch);
	}
      free (line);
      if (err < ((israte != FALSE)?8:9)) {
#ifdef WITH_ERRORS
	puterr ("trafo","select_frameconv");
	fprintf (stderr, "invalid config file %s in line %d\n",
		 conffile, lcount);
#endif
	continue;
      }
      
      angtyp = 0;
      if (strstr (element, "deg") != NULL) angtyp = 1;
      if (strstr (element, "''") != NULL) angtyp = 2;
      if (strstr (element, "\"") != NULL) angtyp = 2;
      if (strstr (element, "rad") != NULL) angtyp = 3;
      if (strstr (element, "cc") != NULL) angtyp = 4;
      if (angtyp == 0) {
#ifdef WITH_ERRORS
	puterr ("trafo","select_frameconv");
	fprintf (stderr, "invalid config file %s in line %d\n",
		 conffile, lcount);
#endif
	continue;
      }
      for (i=0;i<3;i++) 
	{
	  switch (angtyp) 
	    {
	    case 1: ac[i] = afromdeg (rot[i]);            break;
	    case 2: ac[i] = afromdeg (rot[i]/60.0/60.0);  break;
	    case 3: ac[i] = afromrad (rot[i]);            break;
	    case 4: ac[i] = afromgrad (rot[i] / 10000.0); break;
	    }
	}
      
      if (israte == FALSE) 
	{
	  fc->scale = scale;
	  fc->ref_epoch = fepoch;
	  fc->fwd = fwd;
	  for (i=0;i<3;i++) 
	    {
	      fc->dx[i] = dx[i];
	      fc->rot[i] = ac[i].rad;
	    }
	  set = 1;
	} 
      else 
	{
	  fc->dscale = scale;
	  if ((fc->fwd - fwd) != 0 || fc->ref_epoch < 1.0) 
	    {
	     /* rate direction does not match frame direction */
	      set = 2;
	    }
	  for (i=0;i<3;i++) 
	    {
	      fc->ddx[i] = dx[i];
	      fc->drot[i] = ac[i].rad;
	    }
	}
    } /* end of while loop */

  switch (set) 
    {
    case 0:
#ifdef WITH_ERRORS
      puterr ("trafo", "select_ellipsoid");
      fprintf (stderr, "frame conversion %s --> %s not found in %s\n",
	       source, target, conffile);
#endif
      break;
    
    case 2:
#ifdef WITH_ERRORS
      puterr ("trafo", "select_ellipsoid");
      fprintf (stderr,
	       "invalid frame conversion specification %s --> %s in %s\n",
	       source, target, conffile);
#endif
      break;

    default:
      if (fc->fwd > 0)
	ftr = frame_init (inepoch, fc);
      else
	ftr = frame_init (outepoch, fc);
    }

  free (fc);
  for (i=0;i<4;i++) { free(pattern[i]); regfree (compex[i]); }
  fclose (fpconf);
  return ftr;
}


static int set_formtype (const char *type)
{
  char *t;
  unsigned int i;
  t = strdup (type);
  for (i=0;i<strlen(t);i++)
    t[i] = tolower(t[i]);
  
  if (strncmp (t, "x", 1) == 0 )       {free (t); return  1;}
  if (strncmp (t, "y", 1) == 0 )       {free (t); return  2;}
  if (strncmp (t, "z", 1) == 0 )       {free (t); return  3;}
  if (strncmp (t, "n", 1) == 0 )       {free (t); return  4;}
  if (strncmp (t, "e", 1) == 0 )       {free (t); return  5;}
  if (strncmp (t, "h", 1) == 0 )       {free (t); return  6;}
  if (strncmp (t, "lat", 3) == 0 )     {free (t); return  7;}
  if (strncmp (t, "lon", 3) == 0 )     {free (t); return  8;}
  if (strncmp (t, "rlat", 4) == 0 )    {free (t); return  9;}
  if (strncmp (t, "rlon", 4) == 0 )    {free (t); return 10;}
  if (strncmp (t, "glat", 4) == 0 )    {free (t); return 11;}
  if (strncmp (t, "glon", 4) == 0 )    {free (t); return 12;}
  if (strncmp (t, "dmslat", 6) == 0 )  {free (t); return 13;}
  if (strncmp (t, "dmslon", 6) == 0 )  {free (t); return 14;}
  
  free (t);
  return 0;
}

static int set_formcols (const int typenr)
{
  if (typenr < 1) return 0;
  if (typenr > 12) return 3;
  else return 1;
  
}

/*
 *  select form types for input
 */
int select_formtype (
		     const char *formtype,
		     int ft_type[3],
		     int ft_cols[3])
{
  char *pattern, *line, *conffile;
  regex_t *compex;
  FILE *fpconf;
  int i, err, lcount;
  
  /* determine and open config file */
  conffile = getenv ("TRAFOCONF");
  if (conffile == NULL) 
    {
#ifdef WITH_ERRORS
      puterr ("trafo", "select_formtype");
      fprintf (stderr, "environment variable TRAFOCONF not set.\n");
#endif
      return GGLERR_INOUT;
    }
  fpconf = open_text_file (conffile, "r");
  if (!fpconf) return GGLERR_INOUT;

  /* initialize regex pattern */
  pattern = (char *) malloc (sizeof(char)*(20+strlen(formtype)));
  compex = (regex_t *) malloc (sizeof(regex_t));
  sprintf(pattern, "^formtype  *%s ", formtype);
  if ((err = regcomp (compex, pattern, REG_ICASE))!=0) 
    {
      compex_error (err, compex, "select_formtype");
      free(pattern);
      return GGLERR_INOUT;
    }
  for (i=0;i<3;i++) 
    {
      ft_type[i] = 0;
      ft_cols[i] = 0;
    }
  lcount = 0;
  
  while ( (err = read_line (fpconf, &line)) != EOF ) 
    {
      char col0[256], col1[256], col2[256];
      if (line == NULL) break;
      lcount++;
      
      err = regexec (compex, line, 0, NULL, 0);
      if (err == REG_NOMATCH) { free (line); continue; }
      if (err != 0) 
	{
	  compex_error (err, compex, "select_formtype");
	  return GGLERR_INOUT;
	}
      err = sscanf (line, "%*s %*s %s %s %s", col0, col1, col2);
      if (err < 3) 
	{
#ifdef WITH_ERRORS
	  puterr ("trafo", "select_formtype");
	  fprintf (stderr, "invalid config file %s in line %d\n",
		   conffile, lcount);
#endif
	  free (line);
	  continue;
	}
      ft_type[0] = set_formtype (col0);
      ft_type[1] = set_formtype (col1);
      ft_type[2] = set_formtype (col2);
      for (i=0;i<3;i++)
	ft_cols[i] = set_formcols (ft_type[i]);
      break;
    }
  err = 0;

  if (ft_type[0] == 0 || ft_type[1] == 0 || ft_type[2] == 0) {
#ifdef WITH_ERRORS
    puterr ("trafo", "select_formtype");
    fprintf (stderr, "no valid formtype `%s' found in %s\n",
	     formtype, conffile);
#endif
    err = GGLERR_INVAL;
  }
  
  
  free (pattern);
  regfree (compex);
  fclose (fpconf);
  return err;
}


/*
 *  select geoid (read full grid into memory)
 */
struct Grid * select_geoid (
			    const char *geoidstr,
			    int g_ftyp[3])
{
  struct Grid *gg;
  char gfname[256];
  regex_t *compex;
  char *pattern, *line, *conffile;
  FILE *fpconf;
  int err, lcount;
  /* determine and read config file */
  conffile = getenv ("TRAFOCONF");
  if (conffile == NULL) 
    {
#ifdef WITH_ERRORS
    puterr ("trafo", "select_geoid");
    fprintf (stderr, "environment variable TRAFOCONF not set.\n");
#endif
    return NULL;
  }
  fpconf = open_text_file (conffile, "r");
  if (!fpconf) return NULL;
  
  /* initialize regex pattern */
  pattern = (char *) malloc (sizeof(char)*(20+strlen(geoidstr)));
  compex = (regex_t *) malloc (sizeof(regex_t));
  sprintf(pattern, "^geoid  *%s ", geoidstr);
  if ((err = regcomp (compex, pattern, REG_ICASE))!=0) 
    {
      compex_error (err, compex, "select_geoid");
      free(pattern);
      return NULL;
    }
  gg = NULL;
  lcount = 0;
  
  while ( (err = read_line (fpconf, &line)) != EOF ) 
    {
      char col0[256], col1[256], col2[256];
      if (line == NULL) break;
      lcount++;
      err = regexec (compex, line, 0, NULL, 0);
      if (err == REG_NOMATCH) {
	free (line);
	continue;
      }
      if (err != 0) 
	{
	  compex_error (err, compex, "select_geoid");
	  free(pattern);
	  return NULL;
	}
      err = sscanf (line, "%*s %*s %s %s %s %s", col0, col1, col2, gfname);
      if (err < 3) 
	{
#ifdef WITH_ERRORS
	  puterr ("trafo", "select_geoid");
	  fprintf (stderr, "invalid config file %s in line %d\n",
		   conffile, lcount);
#endif
	  free (line);
	  continue;
	}
      g_ftyp[0] = set_formtype (col0);
      g_ftyp[1] = set_formtype (col1);
      g_ftyp[2] = set_formtype ("h");
      if (strstr (col2, "dsaa") == NULL) {
#ifdef WITH_ERRORS
	puterr ("trafo", "select_geoid");
	fprintf (stderr, 
		 "invalid config file %s in line %d\nselect_geoid: ",
		 conffile, lcount);
	fprintf (stderr, 
		 "currently only style `dsaa' supported\n");
#endif
      }
      gg = surf_readgrid (gfname);
    }
#ifdef WITH_ERRORS
  if (gg == NULL) 
    {
      puterr ("trafo", "select_geoid");
      fprintf (stderr,"no valid geoid `%s' found in %s\n",
	       geoidstr, conffile);
    }
#endif
  
  free (pattern);
  regfree (compex);
  fclose (fpconf);
  return gg;
}

