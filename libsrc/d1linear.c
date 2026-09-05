/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/d1linear.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: Linear one dimensional interpolation
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: d1linear.c,v 1.5 2000/11/14 18:00:04 favey Exp $
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
 *   AUTHORS: Favey Etienne   <favey@geod.ethz.ch>
 *
 *   HISTORY: dd-mm-yyyy  --  description
 *            
 *            28-04-1998  --  created
 * $Log: d1linear.c,v $
 * Revision 1.5  2000/11/14 18:00:04  favey
 * draft of acceleration integration, incomplete
 *
 * Revision 1.4  1999/09/08 09:32:31  favey
 * copyright info added
 *
 * Revision 1.3  1998/10/30 14:22:07  favey
 * changed config.h to global.h
 *
 * Revision 1.2  1998/04/28 16:02:06  favey
 * *** empty log message ***
 *
 *
 *
 */

/*
 *  Standard Libraries
 */
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 *  our special headers
 */
#include "global.h"
#include "errproc.h"
#include "d1linear.h"



/*
 *  Simple linear interpolation. 
 */
int d1ipol_linear (
	const double x1,
	const double x2,
	const double y1,
	const double y2,
	const double xip,
	double *yres)
{
	double dx;

	dx = x2 - x1;
	if (fabs(dx) < 1e-25) return EINVAL;

	*yres = ( (xip-x1)/dx )*(y2-y1) + y1;
	return 0;
}


/* Simple Integration of Stapdown Accelerations */
/*int acc2pos (double attitude[3], double accel_bodyframe[3], 
 *	     double posxyz[3], const struct ellipsoid ell, double deltat)
 */   /* attitude in radian [roll,pitch,yaw] */
     /* acceleration in bodyframe [m/s^2] */
     /* position geocentric xyz [m] in: initial position, out: result */
     /* which ellipsoid do we use here? */
     /* time interval for simple integration [s] */
/*{
 *  double RE=6378815.9;           / * Erdradius Bessel CH * /
 *  double ws=(1.0027379*2*3.1415926535/24/3600);/ * Erdrotation in rad / sec * /
 *  int i;
 *  double c1,c2,c3,s1,s2,s3, DGA[3][3];
 *  double aG[3], lat, lon, hgt;
 *  double tm1, R, VE, VN, VV, dV[3], Vsav[3];
 *  xyz_to_llh (ell, posxyz, &lat, &lon, &hgt);
 *  / * compute rotation matrix * /
 *  c1=cos(attitude[0]);    s1=sin(attitude[0]);
 *  c2=cos(attitude[1]);    s2=sin(attitude[1]);
 *  c3=cos(attitude[2]);    s3=sin(attitude[2]);
 *
 *  DGA[0][0]=c3*c2;  DGA[0][1]=-s3*c2+c3*s2*s1;  DGA[0][2]= s3*s1+c3*s2*c1;
 *  DGA[1][0]=s3*c2;  DGA[1][1]= c3*c1+s3*s2*s1;  DGA[1][2]=-c3*s1+s3*s2*c1;
 *  DGA[2][0]=-s2;    DGA[2][1]= c2*s1;           DGA[2][2]= c2*c1;
 *  / * rotation of acceleration * /
 *  for (i=0;i<3;i++) aG[i]=0.0;
 *  for (i=0;i<3;i++) {
 *    aG[0]+=(DGA[0][i]*accel_bodyframe[i]);
 *    aG[1]+=(DGA[1][i]*accel_bodyframe[i]);
 *    aG[2]+=(DGA[2][i]*accel_bodyframe[i]);
 *  }
 *  / * standard navigation equations * /
 *  R = RE+hgt;
 *  tm1 = 2*ws*sin(lat)+(VE/R)*tan(lat);
 *  c1 = 2*ws*cos(lat);
 *
 *  dV[0] = aG[0]    - (tm1)*VE + VN*VV/R;
 *  dV[1] = aG[1]    + (tm1)*VN + VV*c1 + VE*VN/R;
 *  dV[2] = aG[2] / *+g* / - VE*c1 - (VN*VN + VE*VE)/R;
 *
 *  Vsav[0]=VN;
 *  Vsav[1]=VE;
 *  Vsav[2]=VV;
 *                
 *  VN=Vsav[0]+deltat*dV[0];
 *  VE=Vsav[1]+deltat*dV[1];
 *  VV=Vsav[2]+deltat*dV[2];
 *
 *  VN=Vsav[0]+deltat*aG[0];
 *  VE=Vsav[1]+deltat*aG[1];
 *  VV=Vsav[2]+deltat*(aG[2]/ * +g * /);
 *  
 *  posEsav=0; / * DEBUG incomplete !* /
 *  posNsav=0;
 *  posVsav=0;
 *
 *  posE=posEsav+deltat*VE;
 *  posN=posNsav+deltat*VN;
 *  posV=posVsav+deltat*VV;
 *
 *  return 0;
 *}
 */



