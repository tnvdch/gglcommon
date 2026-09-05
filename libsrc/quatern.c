/*
 *      FILE: $Source: /mnt/molnya00/linux/gglutils/src/RCS/common/lib/quatern.c,v $
 *
 *   PROJECT: GGL Common Library
 *
 *   PURPOSE: quaternions 
 *
 * COPYRIGHT: (c) 1997 Geodesy and Geodynamics Lab,
 *                     ETH Zurich
 * $Id: quatern.c,v 1.6 2001/08/31 08:40:18 cocard Exp $
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
 *   AUTHORS: Cocard Marc  <cocard@geod.ethz.ch>
 *            Favey Etienne <favey@geod.ethz.ch>
 *
 *   HISTORY: 
 * $Log: quatern.c,v $
 * Revision 1.6  2001/08/31 08:40:18  cocard
 * *** empty log message ***
 *
 * Revision 1.5  2001/02/21 13:28:19  cocard
 * *** empty log message ***
 *
 * Revision 1.4  1999/09/08 09:12:18  favey
 * copyright info added
 *
 * Revision 1.3  1998/10/30 14:22:07  favey
 * changed config.h to global.h
 *
 * Revision 1.2  1998/09/28 15:34:11  cocard
 * warnings fixed
 *
 * Revision 1.1  1998/04/22 12:35:42  favey
 * Initial revision
 *
 * Revision 0.3  1997/12/24 12:12:03  favey
 * followup of V. 0.29
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
#include <string.h>
#include <math.h>
#include <ctype.h>

/*
 *  our special headers
 */
#include "global.h"
#include "mat3.h"
#include "quatern.h"

int q_add(double q1[4], double q2[4], double result[4])
{
  int i;
  for (i=0;i<4;i++) result[i] = q1[i]+q2[i];
  return 0;
} 

int q_mult(double q1[4], double q2[4], double result[4])
{
  result[0] =  q1[3]*q2[0] - q1[2]*q2[1] + q1[1]*q2[2] + q1[0]*q2[3]; 
  result[1] =  q1[2]*q2[0] + q1[3]*q2[1] - q1[0]*q2[2] + q1[1]*q2[3];
  result[2] = -q1[1]*q2[0] + q1[0]*q2[1] + q1[3]*q2[2] + q1[2]*q2[3];
  result[3] = -q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] + q1[3]*q2[3];
  return 0;  
}

int q_conjugate(double q[4], double result[4])
{
  int i;
  for (i=0;i<3;i++) result[i]=-q[i];
  result[3]=q[3];
  return 0;   
}

int q_scalarmult (double q[4], double factor, double result[4])
{
  int i;
  for (i=0;i<4;i++) result[i]=q[i]*factor; 
  return 0;
}

double q_norm(double q[4])
{
  int i;
  double result=0.0;
  for (i=0;i<4;i++) result+=q[i]*q[i];
  return result;
}

int q_inverse(double q[4], double result[4])
{
  int i;
  double factor=1/q_norm(q);
  for (i=0;i<3;i++) result[i] = -factor*q[i];
  result[3] = factor*q[3];
  return 0;
}

int rotvec_to_quat(double rotvec[3], double angle, double result[4])
{
  int i;
  double x=sin(angle/2.0);
  for (i=0;i<3;i++) result[i] = x*rotvec[i];
  result[3] = cos(angle/2.0);
  return 0;
}

int rotmat_to_quat(double rotmat[3][3], double result[4])
{
  int err;
  double rotvec[3];
  double angle;
  err = rotmat_to_vec(rotmat, rotvec, &angle);
  if (err != 0) return err;
  return rotvec_to_quat(rotvec, angle, result);
}

int quat_to_rotmat(double q[4], double result[3][3])
{
  double f = 2/q_norm(q);
  double x = f*q[0];
  double x2;

  result[1][1] = result[2][2] = 1 - x*q[0];
  result[0][1] = result[1][0] = x*q[1];
  result[0][2] = result[2][0] = x*q[2];
  result[2][1] = -(result[1][2] = x*q[3]);

  x = f*q[1];
  x2 = x*q[1];
  result[0][0] = 1 - x2;
  result[2][2] -= x2;
  x2 = x*q[2];
  result[1][2] += x2;
  result[2][1] += x2;
  x2 = x*q[3];
  result[2][0] += x2;
  result[0][2] -= x2;

  x = f*q[2];
  x2 = x*q[2];
  result[0][0] -= x2;
  result[1][1] -= x2;
  x2 = x*q[3];
  result[0][1] += x2;
  result[1][0] -= x2;

  return 0;
}

int q_rotate(double x[3], double q[4], double result[3])
{
  double r[3];
  double x1, x2;
  int i;

  /* this allows to use the same memspace 
     for x and result in the function call */
  memcpy(r, x, sizeof(double)*3);

  x1 = q[0]*q[0];
  x2 = q[1]*q[1];
  result[1] = -(result[0]=x1-x2);
  result[2] = -x1-x2;
  x1 = q[2]*q[2];
  x2 = q[3]*q[3];
  result[2] += (x1+x2);
  x1=x2-x1;
  result[0] += x1;
  result[1] += x1;
  for (i=0;i<3;i++) result[i] *= r[i];
    
  x1 = q[0]*q[1];
  x2 = q[2]*q[3];
  result[0] += 2*r[1]*(x1+x2);
  result[1] += 2*r[0]*(x1-x2);
  x1 = q[0]*q[2];
  x2 = q[1]*q[3];
  result[0] += 2*r[2]*(x1-x2);
  result[2] += 2*r[0]*(x1+x2);
  x1 = q[1]*q[2];
  x2 = q[0]*q[3];
  result[1] += 2*r[2]*(x1+x2);
  result[2] += 2*r[1]*(x1-x2);
 
  return 0; 
}
