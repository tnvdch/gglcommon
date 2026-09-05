/**
 **        FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/angles.c,v $
 **
 ** DESCRIPTION: Functions to convert angles from and to
 **              + radians (2*M_PI), 
 **              + degrees (360o00'00. or 360.0000), 
 **              + grads (400.0000).
 **              Functions to parse a string and convert its
 **              contents to an angle struct.
 **
 **   COPYRIGHT: Copyright (C) 1996 by Etienne Favey
 ** $Id: angles.c,v 1.5 1999/01/26 17:19:19 favey Exp $
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
 **     CHANGES: 
 ** $Log: angles.c,v $
 ** Revision 1.5  1999/01/26 17:19:19  favey
 ** moved out of swiss.h
 **
 ** Revision 1.4  1998/10/30 14:22:07  favey
 ** changed config.h to global.h
 **
 ** Revision 1.3  1998/05/06 08:58:08  cocard
 ** *** empty log message ***
 **
 ** Revision 1.2  1998/05/04 07:05:55  favey
 ** fixed bug in afromdeg (deg was not set)
 **
 ** Revision 1.1  1998/04/22 12:35:00  favey
 ** Initial revision
 **
 ** Revision 0.4  1998/01/20 10:15:25  favey
 ** minor adjustments to splitting swiss.c into trafo.c and swiss.c
 ** our header is now swiss.h instead of trafo.h
 **
 ** Revision 0.3  1997/12/24 12:04:51  favey
 ** release following 0.29 of the common lib
 **
 **
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "global.h"
#include "angles.h"

struct Angle afromrad(double angrad)
{
  double tmp;
  struct Angle a;

  a.rad = angrad;
  a.grad = angrad/M_PI*200.0;
  a.deg = tmp = angrad/M_PI*180.0;
  tmp = 60*(tmp - (a.dms.deg = (int) floor(tmp) ));
  a.dms.sec = 60*(tmp - (a.dms.min = (int) floor(tmp) ));
  return(a);
}

struct Angle afromdeg(double angdeg)
{
  double tmp;
  struct Angle a;

  a.rad = angdeg/180.0*M_PI;
  a.grad = angdeg/180.0*200.0;
  a.deg = tmp = angdeg;
  tmp = 60*(tmp - (a.dms.deg = (int) floor(tmp) ));
  a.dms.sec = 60*(tmp - (a.dms.min = (int) floor(tmp) ));
  return(a);
}

struct Angle afromgrad(double anggrad)
{
  double tmp;
  struct Angle a;

  a.rad = anggrad/200.0*M_PI;
  a.deg = tmp = anggrad/200.0*180.0;
  a.grad = anggrad;
  tmp = 60*(tmp - (a.dms.deg = (int) floor(tmp) ));
  a.dms.sec = 60*(tmp - (a.dms.min = (int) floor(tmp) ));
  return(a);
}

struct Angle afromdmsverb(int degrees, int minutes, double seconds)
{
  double tmp;
  struct Angle a;

  tmp =   (double) (a.dms.sec = seconds)        / 60.0;
  tmp = ( (double) (a.dms.min = minutes) + tmp) / 60.0;
  a.deg = (double) (a.dms.deg = degrees) + copysign(tmp,degrees);
  a.rad = a.deg/180.0*M_PI;
  a.grad = a.deg/180.0*200.0;
  return(a);
}

struct Angle afromdms(struct hms mydms)
{
  return(afromdmsverb(mydms.deg, mydms.min, mydms.sec));
}

struct Angle afromstr(char *s)
{
  struct Angle a;
  char *tok;
    
  if (((strchr(s,'o')!=NULL) || (strchr(s,'d')!=NULL)) &&
      ((strchr(s,'m')!=NULL) || (strchr(s,'\'')!=NULL)))
    {
      /* string is dms */
      tok = strtok(s,"od");
      a.dms.deg = atoi(tok);
      tok = strtok(NULL,"m'");
      a.dms.min = atoi(tok);
      tok = strtok(NULL,"\0");
      a.dms.sec = atof(tok);
      return(afromdms(a.dms));
    }
  else
    {
      /* string is decimal */
      return(afromdeg(atof(s)));
    }
	  
}

char * admsstr (const struct Angle a, const int mode)
  /*
   * mode: 0 -> 46 30 25.643346789
   *       1 -> 46d30m25.643346789
   *       2 -> 46\t30\t25.643346789
   */
{
  char degstr[256];

  switch (mode)
    {
    case 1:
      sprintf (degstr, "%02do%02dm%08.5f", a.dms.deg, a.dms.min, a.dms.sec);
      break;
    case 2:
      sprintf (degstr, "%02d\t%02d\t%012.9f", a.dms.deg, a.dms.min, a.dms.sec);
      break;
    default:
      sprintf (degstr, "%02d %02d %012.9f", a.dms.deg, a.dms.min, a.dms.sec);
    }
  return strdup(degstr);
}
